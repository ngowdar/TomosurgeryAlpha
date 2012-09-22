using System;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Drawing;
using System.IO;

namespace TomosurgeryAlpha
{
   
    /// <summary>
    /// A collection class that contains the information about a group of RasterPath objects.
    /// Takes tumor object, compresses it into slices and creates a RasterPath object 
    /// for each slice.
    /// </summary>
    public class PathSet
    {
        public event ProgressChangedEventHandler PathsetWorkerProgressChanged;
        public event RunWorkerCompletedEventHandler PathsetWorkerCompleted;
        public event ProgressChangedEventHandler RasterPathWorkerProgress;
        public BackgroundWorker PSworker;
        public Boolean DoseModifiable;
        public static float[,] mask;
        public static int N = 40;
        public static int edgepadding = N/4;
        public static int tumorflag = 10;
        public static int CSflag = 2;
        public int NumSlices;
        public ArrayList RasterPaths;
        public int SliceThickness;
        public int RasterWidth;
        public int TolThickness;
        public float[][,] volume;
        public float[][,] dosespace;
        public int[] SlicePositions;
        public int[] boundaries; //6-term sequence (xstart, xend, ystart...etc.)
        public int X; public int Y; public int Z;

        public PathSet(float[][,] f, int sthick, int tolthick)
        {
            X = f[0].GetLength(0); Y = f[0].GetLength(1); Z = f.GetLength(0);
            boundaries = FindBoundaries(f);
            SliceThickness = sthick;
            TolThickness = tolthick;
            CalculateNumSlices();
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(f, SlicePositions[i], SliceThickness/2));
                RasterPaths.Add(rp);
            }
            //Slice thickness logic goes here:
                /* - Divide by slice thickness
                 * - Decide what to do with the extra bit
                 */
            
            //for (int i = 0; i < NumSlices; i++)
            //{
                /*For each slice, send the smaller float[][,] to the CompressSection()
                 *-An easy way to do this is to loop through x and y,
                 * then add up all the z for each pixel.
                 * -Then create a RasterPath object out of this slice.
                 * -Add this RasterPath object to some sort of collection
                 * that is stored as a variable in PathSet.
                 */
            //}
            volume = f;
            AttachHandlers();
        }

        private void AttachHandlers()
        {
            PSworker = new BackgroundWorker();            
            PSworker.WorkerReportsProgress = true;
            PSworker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(PSworker_RunWorkerCompleted);
            PSworker.ProgressChanged += new ProgressChangedEventHandler(PSworker_ProgressChanged);
            PSworker.DoWork += new DoWorkEventHandler(PSworker_DoWork);
        }

        void PSworker_DoWork(object sender, DoWorkEventArgs e)
        {
            double count = 0;
            for (int i = 0; i < NumSlices; i++)
            {                
                //((RasterPath)RasterPaths[i]).RPworker.RunWorkerAsync();
                ((RasterPath)RasterPaths[i]).OptimizeShotWeights();
                count++;
                PSworker.ReportProgress((int)(100*count/NumSlices));
                
            }
        }

        void PathSet_SliceWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (RasterPathWorkerProgress != null)
                RasterPathWorkerProgress.Invoke(null, e);
        }

        void PSworker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {              
            if (PathsetWorkerProgressChanged != null)
                PathsetWorkerProgressChanged.Invoke(null, e);
        }
        void PSworker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (PathsetWorkerCompleted != null)
                PathsetWorkerCompleted.Invoke(null, e);
        }
        
        public void RecalculateSlices(int sthick, int tolthick)
        {
            SliceThickness = sthick;
            TolThickness = tolthick;
            CalculateNumSlices();
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(volume, SlicePositions[i], SliceThickness / 2));
                RasterPaths.Add(rp);
            }
        }
        public void CalculatePlanDose(DoseKernel dk)
        {
            float[,] midplane = dk.Get2DSlice(N / 2);

        }
        private float[][,] GrabSlab(float[][,] f, int center)
        {
            int from = center - RasterWidth / 2;
            int to = center + RasterWidth / 2;
            int dist = to - from;
            float[][,] slab = new float[dist][,];
            for (int i = 0; i < dist; i++)
            {
                slab[i] = f[i + from];
            }
            return slab;
        }
        public void CalculateNumSlices()
        {
            int remainder_slices; int new_spacing;
            int padding;
            int[] zrasterpos;
            int zstart = boundaries[4];
            int zend = boundaries[5];
            //bool can_another_slice_fit = false;
            int StartPosition = zstart + TolThickness / 2;
            int initial_slice_estimate = ((zend + SliceThickness / 4) - (zstart - SliceThickness / 4)) / SliceThickness;
            int initial_remainder = ((zend + SliceThickness / 4) - (zstart - SliceThickness / 4)) % SliceThickness;

            /*If there is a significant remainder, add a slice and shift to equalize.
             * If there the remainder is small, just shift the slices */
            if (initial_remainder >= 0.7 * SliceThickness)
            {
                NumSlices = initial_slice_estimate + 1;
                zrasterpos = new int[NumSlices];
                padding = (SliceThickness - initial_remainder) / 2;
                zrasterpos[0] = (zstart + SliceThickness / 4) - padding;
                for (int i = 1; i < NumSlices; i++)
                {
                    zrasterpos[i] = zrasterpos[0] + i * SliceThickness;
                }

            }
            else
            {
                NumSlices = initial_slice_estimate;
                padding = initial_remainder / 2;
                zrasterpos = new int[NumSlices];
                zrasterpos[0] = (zstart + SliceThickness / 4) + padding;
                for (int i = 1; i < NumSlices; i++)
                {
                    zrasterpos[i] = zrasterpos[0] + i * SliceThickness;
                }
            }
            SlicePositions = zrasterpos;
        }
        public int[] FindBoundaries(float[][,] f)
        {
            int[] boundaries = new int[6];

            //Z boundaries
            Boolean z1 = false; Boolean z2 = false;
            for (int k = 0; k < Z; k++)
            {
                float[,] temp = f[k];
                float sum = Matrix.SumAll(temp);
                if (sum > 0)
                    if (z1 == false)
                    {
                        boundaries[4] = k;
                        z1 = true;
                    }
                    else
                        boundaries[5] = k;


                //X boundaries
                Boolean x1 = false; Boolean x2 = false;
                for (int i = 0; i < X; i++)
                    for (int j = 0; j < Y; j++)
                    {
                        if (!x1 && temp[i, j] >= 0)
                        {
                            x1 = true;
                            boundaries[0] = i;
                        }
                        if (!x2 && temp[X - i - 1, j] >= 0)
                        {
                            x2 = true;
                            boundaries[1] = X - i;
                        }
                        if (x1 && x2)
                            break;
                    }
                //Y boundaries
                Boolean y1 = false; Boolean y2 = false;
                
                    for (int j = 0; j < Y; j++)
                        for (int i = 0; i < X; i++)
                    {
                        if (!y1 && temp[i, j] >= 0)
                        {
                            y1 = true;
                            boundaries[2] = j;
                        }
                        if (!y2 && temp[i,Y - j - 1] > 0)
                        {
                            y2 = true;
                            boundaries[3] = Y - j;
                        }
                        if (y1 && y2)
                            break;
                    }
            }
            return boundaries;
        }
        public float[,] CompressSection(float[][,] f, int zt, int spt)
        {            
            float[,] squished = new float[X, Y];
            float[,] center = f[zt];
            int bz = zt - spt;
            int ez = zt + spt;            
                        
            for (int j = 0; j < Y; j++)                
                for (int i = 0; i < X; i++)
                    for (int k = 0; k < ez - bz; k++)
                    {
                        squished[i, j] += f[k][i, j];
                    }

            for (int i = 0; i < X; i++)
                for (int j = 0; j < Y; j++)
                    squished[i, j] -= center[i, j];

            squished = Matrix.Add(Matrix.ThresholdEq(squished, 2), center);
            double ddd = Matrix.SumAll(squished);
            for (int i = 0; i < squished.GetLength(0); i++)
                for (int j = 0; j < squished.GetLength(1); j++)
                    if (squished[i, j] > 0)
                        squished[i, j] = 1;
            ddd = Matrix.SumAll(squished);
            return squished;
        }

        #region Assembling Final Dose Computation

        /// <summary>
        /// The actual method that is called by the UI, which takes in a folder directory path, and 
        /// loads in each slicedose file and writes it to the global dosespace file.
        /// </summary>
        /// <param name="folderpath"></param>
        public void AssembleFinalDoseMatrix(string folderpath)
        {
            string path;
            int x = volume[0].GetLength(0); int y = volume[0].GetLength(1); int z = volume.GetLength(0);
            /* First, a for-loop to create the string name for each slice, which then loads each
             * slicedose file and in turn adds it to the result matrix.*/
            for (int s = 0; s < NumSlices; s++)
            {                
                path = string.Concat(folderpath, "slice_");
                path = string.Concat(path, s);
                float[] slicedose = ReadDoseFromFile(path);
                WriteSliceDoseToDoseSpace(slicedose, SlicePositions[s]);
            }

        }

        /// <summary>
        /// Takes in a slice 1D float matrix, and adds it to the global dosespace array for final calculation. Called by 
        /// AssembleFinalDoseMatrix().
        /// </summary>
        /// <param name="slicedose"></param>
        /// <param name="which_z"></param>
        private void WriteSliceDoseToDoseSpace(float[] slicedose, int which_z)
        {
            int TranslateZBy = SliceThickness / 2; // <- NEED TO CHANGE APPROPRIATELY
            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.
            for (int z = 0; z < SliceThickness; z++) 
            {
                int current_z = which_z - TranslateZBy + z;
                dosespace[current_z] = Matrix.Add(dosespace[current_z], GrabSlice(slicedose, z, dosespace[0].GetLength(0), dosespace[0].GetLength(1)));
            }
        }

        /// <summary>
        /// Takes in a 1D float matrix, grabs a 2D slice out based on which z position. Called by GrabSlice().
        /// </summary>
        /// <param name="input"></param>
        /// <param name="which"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private float[,] GrabSlice(float[] input, int which, int x, int y)
        {
            float[,] result = new float[x,y];
            for (int i = 0; i < x; i++)
                for (int j = 0; j < y; j++)
            {
                result[x, y] = input[which * x * y + j * x + i];
            }
            return result;
        }

        /// <summary>
        /// Given a path, loads in the associated slicedose 1D float matrix. Called by AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="loadpath"></param>
        /// <returns></returns>
        public float[] ReadDoseFromFile(string loadpath)
        {
            float[] d;
            using (FileStream fs = new FileStream(loadpath, FileMode.Open, FileAccess.Read))
            using (BinaryReader br = new BinaryReader(fs))
            {
                int x = br.ReadInt16(); int y = br.ReadInt16(); int z = br.ReadInt16();
                d = new float[x * y * z];
                for (int i = 0; i < d.GetLength(0); i++)
                    d[i] = br.ReadSingle();
            }
            return d;
        }
        #endregion


        #region Save SliceDoses to File
        /// <summary>
        /// Given a folder path and a DoseKernel, will create path names for each rasterpath in RasterPaths 
        /// and write each to the numbered path.
        /// </summary>
        /// <param name="dk"></param>
        /// <param name="folderpath"></param>
        public void CalculateSliceDosesAndWrite(DoseKernel dk, string folderpath)
        {
            string path;
            for (int s = 0; s < NumSlices; s++)
            {
                string date = System.DateTime.Now.ToShortTimeString();
                date = string.Concat(new string[]{folderpath,"\\", date, "\\slice_"});                
                path = string.Concat(date, s);
                RasterPath rp = (RasterPath)RasterPaths[s];
                rp.CalculateAndSaveSliceDose(dk, SliceThickness, path);
            }
            //Need to have some kind of confirmation event?
        }

        #endregion
    }

    

    /// <summary>
    /// Contains the shot paths for one slice, along with associated methods for
    /// finding and creating them.
    /// </summary>
    public class RasterPath
    {
        /*BACKGROUND LOGIC:
         *First find the x,y, and z boundaries.         * 
        * Use the x-width to find the number of lines, using LineSpacer().
        * Store the locations of these lines in a temp array, and set
        * the length of the ShotPoints array for this slice.
        * 
        * Now, loop through each line and find the boundary points along it
        * in the y-direction. Feed these points into ShotSpacer();
        * Use the output of both LineSpacer and ShotSpacer to fill in the
        * PointF[][] array for this slice.
        * 
        * Finally, have an adjustment method. 
        * AdjustLineSpacing(int width): which recalculates the line locations
        * AND the boundary positions on each line, and updates the PointF[][].
        * 
        * Try to do this via GPU for fast rendering.
        */
        
        #region Variables
        public static float[,] dosemidplane; //Dose midplane
        public static event ProgressChangedEventHandler SliceWorkerProgressChanged;
        public static event RunWorkerCompletedEventHandler SliceWorkerCompleted;
        public BackgroundWorker RPworker;
        public static int doseN = 161;     //Size of dose matrix
        public static int N;
        public static int X; public static int Y;
        public static bool DoseModifiable = false;
        public static float[,] mask;
        public int RasterWidth;
        public int StepSize;
        public int NumOfLines;
        public int NumOfShots;
        public double coverage;
        public int isovol;
        public int massvol;
        public int[] boundaries;
        float[,] slice;
        public float[,] ModdedSlice; //The "addition" layer for dose inhomogeneity.
        public float[] dosespace;
        public float[] weight;
        public Slice info;
        bool horizontal = false; //By default all lines will be vertical.
        bool optimized = false;
        PointF[][] shot_points; //A jagged array. One array of points for each line.
            //Note that this is meant to be direction independent (horizontal or vertical).
        PointF[] shots;
        #endregion

        #region Constructors
        public RasterPath(float[,] f)
        {
            SetParams(20, 20);
            slice = f;
            ModdedSlice = f;
            X = f.GetLength(0); Y = f.GetLength(1);
            FindAllShotPoints();
            CreateSliceInfo();
            InitWeightArray();
            AttachHandlers();
        }

        private void RePrioritizeDDS(float[,] mod)
        {
            //Add ModLayer onto binaryslice
            ModdedSlice = Matrix.Add(ModdedSlice, mod);
            Matrix.Normalize(ModdedSlice);
            
        }

        private void AttachHandlers()
        {
            RPworker = new BackgroundWorker();
            RPworker.ProgressChanged += new ProgressChangedEventHandler(RPworker_ProgressChanged);
            RPworker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(RPworker_RunWorkerCompleted);
            RPworker.WorkerReportsProgress = true;
            RPworker.DoWork += new DoWorkEventHandler(RPworker_DoWork);
        }
#endregion
        #region Background Worker Methods
        void RPworker_DoWork(object sender, DoWorkEventArgs e)
        {
            OptimizeShotWeights();
        }

        void RPworker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (SliceWorkerCompleted != null)
                SliceWorkerCompleted.Invoke(null, e);
        }

        void RPworker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (SliceWorkerProgressChanged != null)
                SliceWorkerProgressChanged.Invoke(null, e);
        }
        #endregion



        #region Optimization methods
        public void OptimizeShotWeights()
        {
            shots = ReturnSinglePoints();
            double Error = 1000; double coverage = 0.8;
            int index = 0;
            float[] w = new float[weight.GetLength(0)];
            float[] ds = new float[slice.GetLength(0)*slice.GetLength(1)];

            while (Error >= .0001 || index < 50)
            {
                //Add the weighted shots to the array
                Parallel.For(0, shots.GetLength(0), i =>
                {
                    if (weight[i] == 0)
                        weight[i] = 1.0f;
                    ds = IterateShot(ds, dosemidplane, i);
                });
                ds = Normalize(ds);
                Parallel.For(0, shots.GetLength(0), i =>
                {
                    float ratio = EvalShotWeightIteration(dosemidplane, ds, shots[i]);
                    w[i] = weight[i] * ratio;
                });
                w = Normalize(w);
                weight = Normalize(weight);                
                Error = FindError(weight, w);
                double cov = CalculateIterationCoverage(ds, 0.5f);
                if (index > 1 && cov < coverage)
                    break;
                else
                    coverage = cov;
                index++;
                weight = w;
                //REPORT PROGRESS HERE:
                string report = "Error: " + Error + "; Current iteration: " + index;
                SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, null));
                //RPworker.ReportProgress(2 * index);

            } //END of WHILE LOOP
            this.coverage = coverage;
            dosespace = ds;
            optimized = true;
        }
        
        private double FindError(float[] weight, float[] w)
        {
            double diff = 0;
            for (int i = 0; i < w.GetLength(0); i++)
            {
                diff += Math.Pow(Math.Abs((weight[i] - w[i])),2);
            }
            return Math.Sqrt(diff);
        }

        private float EvalShotWeightIteration(float[,] dosemidplane, float[] ds, PointF pf)
        {
            float[,] DStimesP = GetMultiplied_DS_Subset(ds, pf.X, pf.Y, dosemidplane);
            float[,] DDStimesP = GetMultiplied_DDS_Subset(slice, pf.X, pf.Y, dosemidplane);
            float dssum = Matrix.SumAll(DStimesP);
            float ddssum = Matrix.SumAll(DDStimesP);
            return (float)(ddssum / dssum);            
        }

        ///// <summary>
        ///// Calculates the full 3D matrix of dose, and converts it to a file.
        ///// </summary>
        ///// <param name="dk"></param>
        ///// <param name="slicethickness"></param>
        //private void Calculate_3D_Dose(DoseKernel dk, int slicethickness)
        //{
        //    int dosethickness = 3 * slicethickness + 1;
        //    int startz = (DoseKernel.N - dosethickness) / 2;
        //    int endz = startz + dosethickness;
        //    int actualx_start = 0;
        //    int actualy_start = 0;
        //    int actualz_start = 0;
        //    int ex = 0;
        //    int ey = 0;
        //    int ez = 0;
        //    float[][,] SliceSlab = new float[dosethickness][,];
        //    float[][,] DoseSlab = dk.GetDoseSlab(startz, endz);
        //    PointF center = new PointF(DoseSlab[0].GetLength(0) / 2, DoseSlab[0].GetLength(1) / 2);
        //    PointF[] startpts = new PointF[shots.GetLength(0)];
        //    for (int i = 0; i < shots.GetLength(0); i++)
        //    {
        //        //startpts[i] = GetStartingPoint(center, shots[i]);
        //        PointF temp = GetStartingPoint(center, shots[i]);
        //        PointF actualstart = temp;
        //        if (temp.X < 0)
        //            actualstart.X = 0 - (int)temp.X;
        //        else if (temp.X + DoseSlab[0].GetLength(0) > SliceSlab[0].GetLength(0))
        //            ex = SliceSlab[0].GetLength(0) - ((int)temp.X + DoseSlab[0].GetLength(0));
        //        if (temp.Y < 0)
        //            actualstart.Y = 0 - (int)temp.Y;
        //        else if (temp.Y + DoseSlab[0].GetLength(1) > SliceSlab[0].GetLength(1))
        //            ey = SliceSlab[0].GetLength(1) - ((int)temp.Y + DoseSlab[0].GetLength(1));
        //        //temp is the starting point within the TUMOR space (not the doseslab space)
        //        /*Need to check if the starting point is negative, then the starting location needs to change.                
        //         * Conditions:
        //         * 1) If starting point of dose slab is negative
        //         *   then need to start the doseslab loop below at the
        //         *   location that corresponds to the true 0 location within the tumor matrix.
        //         * 2) If ending point of dose slab is outside of matrix
        //         * then need to end the dose loop earlier (bc end is beyond matrix)
        //         * Do this by subtracting the endx/endy/endz variable from the doseslab conditions.                 
        //         * 
        //         * In plain english, if the start point in the tumor space is negative, make the
        //         * actual start point in the doseslab space equal to the negative shift amount
        //         * (i.e. move it forwards so that it aligns with 0,0)
        //         * If the start point (temp) plus the length of the dose slab is beyond the bounds
        //         * of the tumor matrix, then set the ending point within the doseslab loop to be the
        //         * coordinate that aligns with the last pixel in the tumorspace, which would just be 
        //         * the length of sliceslab - the total length, which will give the negative shift
        //         * amount to adjust the loop end by.
        //         */
        //        startpts[i] = actualstart;
        //    }

        //    //Assuming each shot has its own indexed weight, can just do one for loop to add all of the things.
        //    //Loop through each pixel in the dosekernel, and add it to the translated location relative to the shot center, 
        //    //weighted by the shotweight.
        //    for (int z = 0; z < DoseSlab.GetLength(0); z++)
        //        for (int doseY_index = 0; doseY_index < DoseSlab[0].GetLength(1)-ey; doseY_index++)
        //            for (int doseX_index = 0; doseX_index < DoseSlab[0].GetLength(0)-ex; doseX_index++)
        //        {
        //            Parallel.For(0, startpts.GetLength(0), (shotpoint) =>
        //            {
        //                //The 
        //                PointF topleft_dose = startpts[shotpoint];
        //                SliceSlab[z][(int)(topleft_dose.X + doseX_index), (int)(topleft_dose.Y + doseY_index)] = DoseSlab[z][doseX_index, doseY_index] * weight[shotpoint];
        //            });
        //        }
        //}

        public void CalculateAndSaveSliceDose(DoseKernel dk, int slicethickness, string savepath)
        {
            PointF[] startingpoints = GetStartingPoints(shots);
            int xsize = slice.GetLength(0); int ysize = slice.GetLength(1);
            int xmid = xsize / 2; int ymid = ysize / 2; int zmid = slicethickness / 2;
            float[] slicedose = new float[xsize * ysize * slicethickness];
            for (int k = 0; k < slicethickness; k++)
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                    {
                        for (int w = 0; w < shots.GetLength(0); w++)
                        {
                            PointF shot = shots[w];
                            PointF center = new PointF((N-1)/2,(N-1)/2);
                            //TODO: Need to find the bottom left pixel (starting position) of the dose in 
                            //terms of the slicedose coordinate system.
                            PointF firstdosepixel = FindFirstExistingDosePixel(shot);
                            PointF lastdosepixel = FindLastExistingDosePixel(shot,new PointF(xsize,ysize));
                            if (firstdosepixel.X > i || firstdosepixel.Y > j) //if the current dose pixel doesn't exist for the shot, continue
                                continue;
                            else                            
                                slicedose[k * xsize * ysize + ((int)shot.Y- (int)center.Y+j)* xsize + ((int)shot.X-(int)center.X+i)] = dk.ReturnSpecificDoseValue(i, j, k) * weight[w];
                        }
                    }
            WriteToFile(savepath, slicedose, xsize, ysize, slicethickness);
        }

        public PointF FindFirstExistingDosePixel(PointF shot)
        {
            PointF first = new PointF(0,0);
            PointF center = new Point((N - 1) / 2, (N - 1) / 2);
            float Xdiff = shot.X - center.X;
            float Ydiff = shot.Y - center.Y;
            
            if (Xdiff < 0) //i.e. first dosepixel is outside.
                first.X = center.X - shot.X;
            else 
                first.X = 0; //set first dose pixel to first corresponding real tumor pixel
            
            //Repeat above logic for Y coordinate.
            if (Ydiff >= 0)
                first.Y = 0;
            else if (Ydiff < 0)
                first.Y = center.Y - shot.Y;

            return first;
        }

        public PointF FindLastExistingDosePixel(PointF shot, PointF tumorsize)
        {
            PointF last = new PointF(0,0);
            PointF center = new PointF((N - 1) / 2, (N - 1) / 2);
            float XDistToTumorEdge = tumorsize.X - shot.X;
            float YDistToTumorEdge = tumorsize.Y - shot.Y;

            if (XDistToTumorEdge < center.X) //last dosepixel outside tumor
                last.X = center.X + XDistToTumorEdge;
            else  //last dosepixel inside tumor
                last.X = N - 1;

            //Repeat above logic for Y coordinate
            if (YDistToTumorEdge < center.Y)
                last.Y = center.Y + YDistToTumorEdge;
            else
                last.Y = N - 1;
            
            return last;


        }
        public PointF[] GetStartingPoints(PointF[] shots)
        {
            PointF[] output = new PointF[shots.GetLength(0)];
            for (int i = 0; i < shots.GetLength(0); i++)
            {
                output[i] = new PointF(shots[i].X - (N / 2), shots[i].Y - (N / 2));
            }
            return output;
        }

        public void WriteToFile(string savepath, float[] sd, int sizex, int sizey, int sizez)
        {
            using (FileStream fs = new FileStream(savepath, FileMode.OpenOrCreate, FileAccess.Write))
            using (BinaryWriter bw = new BinaryWriter(fs))
            {
                bw.Write(sizex); //first write the size
                bw.Write(sizey);
                bw.Write(sizez);
                foreach (float f in sd)
                    bw.Write(f);
            }
        }

        

        private PointF GetStartingPoint(PointF center, PointF coord)
        {
            PointF start = new PointF();
            /*The coord is the shot coordinate, which should correspond to the dead center of
             * the dose slab. If a shot is located at (0,0,0), then the actual starting pixel
             * is in the negative, because it is (coord - center) for find the starting point.
            */
            start.X = coord.X - center.X;
            start.Y = coord.Y - center.Y;
            return start;
        }

        public static float[,] GetMultiplied_DDS_Subset(float[,] slice, float px, float py, float[,] P)
        {
            //int startx = (int)px-(P.GetLength(0)-1)/2;
            //int starty = (int)py-(P.GetLength(1)-1)/2;
            float[,] output;

            // ADD DDS ADDITIONS HERE!!!
            if (DoseModifiable)
            {
                float[,] TempMod = Matrix.Add(mask,slice);                
                output = Matrix.MultiplySubset(TempMod, P, (int)px, (int)py);
            }
            else
                output = Matrix.MultiplySubset(slice, P, (int)px, (int)py);

            return output;
        }
        

        public static float[,] GetMultiplied_DS_Subset(float[] ds, float px, float py, float[,] P)
        {
            //int startx = (int)px - (P.GetLength(0) - 1) / 2;
            //int starty = (int)py - (P.GetLength(1) - 1) / 2;
            float[,] output = Matrix.MultiplySubset(ds, P, (int)px, (int)py, X, Y);
            return output;
        }

        public float[] IterateShot(float[] ds, float[,] dosemidplane, int whichshot)
        {
            int index;
            float ds_x; float ds_y;
            for (int i = 0; i < dosemidplane.GetLength(1); i++)
                for (int j = 0; j < dosemidplane.GetLength(0); j++)
                {
                    ds_x = ((PointF)shots[whichshot]).X - (doseN / 2) + j;
                    ds_y = ((PointF)shots[whichshot]).Y - (doseN / 2) + i;
                    index = (int)Math.Round(ds_x) + StructureSet.size * (int)Math.Round(ds_y);
                    ds[index] += dosemidplane[j, i] * weight[whichshot];
                }
            return ds;
        }

        private void NormalizeDose()
        {
            //Find max
            dosespace = Normalize(dosespace);
        }

        private float[] Normalize(float[] d)
        {
            float max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                if (d[i] > max)
                    max = d[i];
            for (int j = 0; j < d.GetLength(0); j++)
                d[j] = (float)(d[j] / max);
            return d;
        }

        public void Calculate2DCoverage(float iso)
        {
            float isovolume = 0; float both = 0; float tumor = 0;
            float dose; float div;
            for (int i = 0; i < slice.GetLength(0); i++)
                for (int j = 0; j < slice.GetLength(1); j++)
                {
                    dose = dosespace[i + (j * slice.GetLength(0))];
                    if (dose >= iso)
                    {
                        isovolume++;
                        if (slice[i, j] > 0)
                        {
                            tumor++;
                            both++;
                        }
                    }
                    else
                    {
                        if (slice[i, j] > 0)
                            tumor++;
                    }
                }
            div = (float)(both / tumor);
            coverage = div;
            isovol = (int)isovolume;
            massvol = (int)tumor;
        }

        public double CalculateIterationCoverage(float[] ds_temp, float iso)
        {
            double c = 0;
            float isovolume = 0; float both = 0; float tumor = 0; float dose;
            Parallel.For(0, shots.GetLength(0), i =>
                {
                    for (int j = 0; j < slice.GetLength(1); j++)
                    {
                        dose = ds_temp[i + j * slice.GetLength(0)];
                        if (dose >= iso)
                        {
                            isovolume++;
                            if (slice[i, j] > 0)
                            {
                                tumor++;
                                both++;
                            }
                        }
                        else
                        {
                            if (slice[i, j] > 0)
                                tumor++;
                        }
                    }
                });
            c = (float)(both / tumor);
            return c;
        }
        

        #endregion

        #region Preliminary Plan methods (paramters, shot locations, etc)
        public PointF[] ReturnSinglePoints()
        {
            PointF[] giant;
            int total = 0;
            for (int i = 0; i < shot_points.GetLength(0); i++)
                total += shot_points[i].GetLength(0);
            giant = new PointF[total];
            int counter = 0;
            for (int i = 0; i < shot_points.GetLength(0); i++)
                for (int j = 0; j < shot_points[i].GetLength(0); j++)
                {
                    giant[counter] = shot_points[i][j];
                    counter++;
                }
            return giant;
        }
        public int[] FindSliceBoundaries(float[,] temp)
        {
            int[] boundaries = new int[4];

            //X boundaries
            Boolean x1 = false; Boolean x2 = false;
            for (int i = 0; i < X; i++)
                for (int j = 0; j < Y; j++)
                {
                    if (!x1 && temp[i, j] > 0)
                    {
                        x1 = true;
                        boundaries[0] = i;
                    }
                    if (!x2 && temp[X - i - 1, j] > 0)
                    {
                        x2 = true;
                        boundaries[1] = X - i;
                    }
                    if (x1 && x2)
                        break;
                }
            //Y boundaries
            Boolean y1 = false; Boolean y2 = false;

            for (int j = 0; j < Y; j++)
                for (int i = 0; i < X; i++)
                {
                    if (!y1 && temp[i, j] > 0)
                    {
                        y1 = true;
                        boundaries[2] = j;
                    }
                    if (!y2 && temp[i, Y - j - 1] > 0)
                    {
                        y2 = true;
                        boundaries[3] = Y - j;
                    }
                    if (y1 && y2)
                        break;
                }

            return boundaries;
        }        
        public void ChangeParamsUpdatePoints(int stepsize, int rasterwidth)
        {
            SetParams(stepsize, rasterwidth);
            FindAllShotPoints();
        }
        private void SetParams(int stepsize, int rasterwidth)
        {
            StepSize = stepsize;
            RasterWidth = rasterwidth;
        }        
        /// <summary>
        /// For a given line position in x-coordinates, will find the y-boundaries.
        /// To be used by ShotSpacer function to get shot coordinates.
        /// </summary>
        /// <param name="linepos"></param>
        /// <returns></returns>
        private int[] LineBoundaries(int linepos)
        {
            int[] b = new int[2];
            float[] line = new float[slice.GetLength(1)];
            for (int i = 0; i < slice.GetLength(1); i++)
                line[i] = slice[linepos, i];
            Boolean y1=false;
            for (int y = 0; y < slice.GetLength(1); y++)
            {
                if (!y1 && slice[linepos, y] > 0)
                {
                    y1 = true;
                    b[0] = y;
                }
                else if (y1 && slice[linepos, y] == 0)
                {
                    b[1] = y;
                    break;
                }
            }
            return b;
        }        
        /// <summary>
        /// Umbrella method that determines the entire shot_points array
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public void FindAllShotPoints()
        {
            //Get slice boundaries first
            int[] boundaries = FindSliceBoundaries(slice);

            //Get line positions for the slice
            int[] lines = LineSpacer(boundaries[0], boundaries[1]);
            shot_points = new PointF[lines.GetLength(0)][];
            NumOfShots = 0;
            for (int i = 0; i < lines.GetLength(0); i++)
            {
                int[] ybounds = LineBoundaries(lines[i]);
                int[] lineshots = ShotSpacer(ybounds[0], ybounds[1]);
                PointF[] PF_shots = new PointF[lineshots.GetLength(0)];
                for (int m = 0; m < lineshots.GetLength(0); m++)
                {
                    PF_shots[m] = new PointF(lines[i], lineshots[m]);
                }
                shot_points[i] = PF_shots;
                NumOfShots += lineshots.GetLength(0);
            }
        }
        private void InitWeightArray()
        {
            weight = new float[NumOfShots];
            for (int i = 0; i < NumOfShots; i++)
                weight[i] = 1;
        }
        public int[] ShotSpacer(int ystart, int yend)
        {
            int[] shots;
            if ((yend - ystart) < (1.75 * StepSize))
            {
                shots = new int[1];
                shots[0] = (ystart + yend) / 2;
            }
            else
            {

                //Find shot spacing
                int edgepad = PathSet.edgepadding;
                int meat; int first; int last;
                int numshots; int newspacing;
                //Rest of Logic goes here...
                meat = yend - ystart - edgepad*2;
                first = ystart + edgepad;
                last = yend - edgepad;
                numshots = meat / StepSize;
                
                //If the meat section is cleanly divisible, -1 for the extra shot
                if (meat % StepSize == 0)
                    numshots = (meat / StepSize) - 1;
                //Find the new spacing
                newspacing = meat / (numshots + 1);
                //Add the first and last shots to the total number of shots
                numshots += 2;
                //Respace shots by newspacing distance.
                shots = new int[numshots];
                shots[0] = first;
                shots[numshots - 1] = last;
                for (int i = 1; i < numshots - 1; i++)
                    shots[i] = shots[i - 1] + newspacing;
            }
            //NumOfShots += shots.GetLength(0);
            return shots;

        }        
        /// <summary>
        /// Given a xstart and xend boundary of a slice, will calculate the number of lines.
        /// </summary>
        /// <param name="xstart"></param>
        /// <param name="xend"></param>
        /// <returns></returns>
        public int[] LineSpacer(int xstart, int xend)
        {            
            int[] lines;
            int edgepad = PathSet.edgepadding;
            int meat; int first; int last;
            int numlines; int newspacing;
            //Is there enough room for the two starting lines?
            if ((xend - xstart) < (1.75 * RasterWidth))
            {
                lines = new int[1];
                lines[0] = (xend + xstart) / 2; //Single line in dead center.
            }
            else
            {
                meat = xend - xstart - edgepad*2; //portion between the first and last lines.
                first = xstart + edgepad;
                last = xend - edgepad;
                numlines = meat / RasterWidth; //Get the maximum WHOLE lines that will fit
                if (meat % RasterWidth == 0)
                    numlines = (meat / RasterWidth) - 1;                
                newspacing = meat / (numlines + 1);
                numlines += 2; //Add the first and last lines
                //Respace by dividing the meat distance by the total number of spaces = (numlines +1)
                lines = new int[numlines];
                lines[0] = first;
                lines[numlines-1] = last;
                for (int i = 1; i < numlines - 1; i++)
                    lines[i] = lines[i - 1] + newspacing;
            }
            NumOfLines = lines.GetLength(0);
            return lines;
        }
        public void Calculate2DDoseSpace(float[,] dosemidplane)
        {
            /*Assuming that the shots are defined by their center coordinate,
             * the top left pixel (0,0) of the shot matrix can be found by
             * subtracting N/2 from the PointF X and Y found in the shots array.
             * 
             * Adding i or j to this will give the absolute position of the relative
             * i,j pixel values within the dose space array
             */
            
            shots = ReturnSinglePoints();
            if (!optimized)
                InitWeightArray();
            float ds_x; float ds_y;
            dosespace = new float[StructureSet.size * StructureSet.size];
            int index;
            //The first two arrays loop through the midplane
            for (int i = 0; i < dosemidplane.GetLength(1); i++)                            
                for (int j = 0; j < dosemidplane.GetLength(0); j++)
                {
                    //This loop calculates the same pixel [i,j] for each shot.
                    //Parallel.For(0, shots.GetLength(0), k =>
                    for (int k = 0; k < shots.GetLength(0); k++)
                    {
                        //Finds the coordinates relative to dosespace
                        ds_x = ((PointF)shots[k]).X - (doseN / 2) + j;
                        ds_y = ((PointF)shots[k]).Y - (doseN / 2) + i;

                        //Add the final result
                        index = (int)Math.Round(ds_x) + StructureSet.size * (int)Math.Round(ds_y);
                        dosespace[index] += dosemidplane[j, i] * weight[k];
                    
                    }//);
                }
            NormalizeDose();
            Calculate2DCoverage(0.5f);
        }                
        #endregion        
        
        
        /// <summary>
        /// Create a Slice struct unique to this slice that contains info
        /// to be used by display methods (like the listbox.
        /// </summary>
        internal void CreateSliceInfo()
        {
            info = new Slice();
            info.NumberOfLines = NumOfLines;
            info.NumberOfShots = NumOfShots;
            info.Coverage = coverage;            
        }
    }

    public struct Slice
    {
        public int SliceNumber { set; get; }
        public int NumberOfLines { set; get; }
        public int NumberOfShots { set; get; }
        public double Coverage { set; get; }
        public int SliceThickness { set; get; }
    }
}
