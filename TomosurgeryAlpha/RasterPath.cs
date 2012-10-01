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
using System.Diagnostics;
    

namespace TomosurgeryAlpha
{
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
        public PointF[][] shot_points; //A jagged array. One array of points for each line.
            //Note that this is meant to be direction independent (horizontal or vertical).
        public PointF[] shots;
        #endregion

        #region Constructors
        public RasterPath(float[,] f)
        {
            SetParams(20, 20);
            slice = f;
            ModdedSlice = f;
            X = f.GetLength(0); Y = f.GetLength(1);
            FindAllShotPoints();            
            InitWeightArray();
            AttachHandlers();
            CreateSliceInfo();
        }

        private void RePrioritizeDDS(float[,] mod)
        {
            //Add ModLayer onto binaryslice
            ModdedSlice = Matrix.Add(ModdedSlice, mod);
            Matrix.Normalize(ModdedSlice);
            
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
        private void AttachHandlers()
        {
            RPworker = new BackgroundWorker();
            RPworker.ProgressChanged += new ProgressChangedEventHandler(RPworker_ProgressChanged);
            RPworker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(RPworker_RunWorkerCompleted);
            RPworker.WorkerReportsProgress = true;
            RPworker.DoWork += new DoWorkEventHandler(RPworker_DoWork);
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
            int ds_x; int ds_y;
            dosespace = new float[StructureSet.size * StructureSet.size];
            int index;
            //The first two arrays loop through the midplane
            for (int j = 0; j < dosemidplane.GetLength(1); j++)                            
                for (int i = 0; i < dosemidplane.GetLength(0); i++)
                {
                    //This loop calculates the same pixel [i,j] for each shot.
                    //Parallel.For(0, shots.GetLength(0), k =>
                    for (int k = 0; k < shots.GetLength(0); k++)
                    {
                        //Finds the coordinates relative to dosespace
                        ds_x = (int)((PointF)shots[k]).X - ((doseN - 1) / 2) + i;
                        ds_y = (int)((PointF)shots[k]).Y - ((doseN - 1) / 2) + j;

                        //Add the final result
                        index = ds_x + StructureSet.size * ds_y;
                        dosespace[index] += dosemidplane[i, j] * weight[k];
                    
                    }//);
                }
            NormalizeDose();
            Calculate2DCoverage(0.5f);
        }                
        #endregion        
        
        #region Shot Optimization methods (Step 1)
          #region Helper methods for OptimizeShotWeight (Step 1)
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

        private float[] ReoptimizeShotWeights(float[] ds, float[] temp_weight)
        {            
            for (int shot = 0; shot < shots.GetLength(0); shot++)
            {
                PointF pf = shots[shot];
                float[,] DStimesP = GetMultiplied_DS_Subset(ds, pf.X, pf.Y, dosemidplane);
                float[,] DDStimesP = GetMultiplied_DDS_Subset(slice, pf.X, pf.Y, dosemidplane);
                float dssum = Matrix.SumAll(DStimesP);
                Debug.Assert(dssum > 0);
                float ddssum = Matrix.SumAll(DDStimesP);
                Debug.Assert(ddssum > 0);
                float ratio = (float)(ddssum / dssum);
                Debug.Assert(ratio > 0);
                temp_weight[shot] = weight[shot] * ratio;
                Debug.WriteLine("Old: " + weight[shot] + " * R: " + ratio + " = New: " + temp_weight[shot]);
            }
            temp_weight = Normalize(temp_weight);
            return temp_weight;
        }
        public void CalculateAndSaveSliceDose(DoseKernel dk, int dosecalcthickness, string savepath)
        {
            PointF[] startingpoints = GetStartingPoints(shots);
            N = dk.DKI.Size;
            int xsize = slice.GetLength(0); int ysize = slice.GetLength(1);
            int xmid = xsize / 2; int ymid = ysize / 2; int zmid = dosecalcthickness / 2;
            int StartingDoseSlice = ((N - 1) / 2) - zmid;
            float[] slicedose = new float[xsize * ysize * dosecalcthickness];
            for (int k = 0; k < dosecalcthickness; k++)
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                        for (int w = 0; w < shots.GetLength(0); w++)
                        {
                            PointF shot = shots[w];
                            PointF center = new PointF((N - 1) / 2, (N - 1) / 2);                            
                            PointF firstdosepixel = FindFirstExistingDosePixel(shot);
                            PointF lastdosepixel = FindLastExistingDosePixel(shot, new PointF(xsize, ysize));
                            if (i < firstdosepixel.X || j < firstdosepixel.Y) //if the current dose pixel doesn't exist for the shot, continue
                                continue;
                            else if (i > lastdosepixel.X || j > lastdosepixel.Y)
                                continue;
                            else
                            {
                                if (dk.ReturnSpecificDoseValue(i, j, k) * weight[w] > 0)
                                    slicedose[k * xsize * ysize + ((int)shot.Y - (int)center.Y + j) * xsize + ((int)shot.X - (int)center.X + i)] += dk.ReturnSpecificDoseValue(i, j, StartingDoseSlice+k) * weight[w];
                            }
                        }                    
            float f = dk.ReturnSpecificDoseValue(80, 80, 80);
            WriteToFile(savepath, slicedose, xsize, ysize, dosecalcthickness);
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
            using (FileStream fs = System.IO.File.Create(savepath))//new FileStream(savepath, FileMode.OpenOrCreate, FileAccess.Write))
                 using (StreamWriter bw = new StreamWriter(fs))
            {
                bw.WriteLine(sizex); //first write the size
                bw.WriteLine(sizey);
                bw.WriteLine(sizez);
                foreach (float f in sd)
                    bw.WriteLine(f);
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

        //public float[] IterateShot(float[] ds, float[,] dosemidplane, int whichshot)
        //{            
        //    //NOTE THIS IS REPLACED BY ITERATESHOTS()
        //    int ds_x; int ds_y;
        //    for (int j = 0; j < dosemidplane.GetLength(1); j++)
        //        for (int i = 0; i < dosemidplane.GetLength(0); i++)                
        //        {
        //            //Find index pixel relative to center shot coordinate
        //            ds_x = (int)((PointF)shots[whichshot]).X - ((N-1) / 2) + i;
        //            ds_y = (int)((PointF)shots[whichshot]).Y - ((N-1) / 2) + j;                    

        //            //Add the weighted dose pixel to the indexed location
        //            ds[(ds_y * StructureSet.size) + ds_x] += dosemidplane[i, j] * weight[whichshot];
        //        }            
            
        //    return ds;
        //}

        /// <summary>
        /// Adds a weighted dose midplane at each particular shot location to the dosespace matrix,
        /// based on the most recent weights. Called by OptimizeShotWeight().
        /// </summary>
        private float[] IterateShots(float[] ds)
        {
            for (int j = 0; j < dosemidplane.GetLength(1); j++)
                for (int i = 0; i < dosemidplane.GetLength(0); i++)
                    for (int whichshot = 0; whichshot < shots.GetLength(0); whichshot++)
                    {
                        //Find index pixel relative to center shot coordinate
                        int ds_x = (int)((PointF)shots[whichshot]).X - ((N - 1) / 2) + i;
                        int ds_y = (int)((PointF)shots[whichshot]).Y - ((N - 1) / 2) + j;

                        //Add the weighted dose pixel to the indexed location
                        ds[(ds_y * StructureSet.size) + ds_x] += dosemidplane[i, j] * weight[whichshot];
                    }
            ds = Normalize(ds);
            return ds;
        }

        private void NormalizeDose()
        {            
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
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)                
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

        public double CalculateIterationCoverage(float iso)
        {
            double c = 0;
            float isovolume = 0; double both = 0; double tumor = 0; float dose;
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
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
            c = both / tumor;
            return c;
        }
        #endregion
        
        /// <summary>
        /// Called from RPworker, which is called from within a slice loop in PathSet. 
        /// This method comprises the meat of Step 1 (shot-based) optimization, which
        /// finds shot weights within a single slice. Writes final result to the float[] ds matrix
        /// and calculates coverage too.
        /// </summary>        
        public void OptimizeShotWeights()
        {
            shots = ReturnSinglePoints();
            double Error = 1000; double coverage = 0.4;
            int index = 0;
            float[] temp_weight = new float[weight.GetLength(0)];
            for (int i = 0; i < shots.GetLength(0); i++)
            {
                temp_weight[i] = 1.0f;
                weight[i] = 1.0f;
            }            
            float[] ds = new float[slice.GetLength(0)*slice.GetLength(1)];

            while (Error >= .01 || index < 25)
            {
                ds = IterateShots(ds);
                Debug.Assert(ds.Max() > 0);                
                Debug.Assert(temp_weight.Max() > 0);
                temp_weight = ReoptimizeShotWeights(ds, temp_weight);                                
                Error = FindError(weight, temp_weight);
                ds = Normalize(ds);
                double temp_coverage = CalculateIterationCoverage(0.5f);

                if (temp_coverage < coverage) //Old version was index > 1 && cov < coverage
                    break;
                else
                {
                    coverage = Convert.ToDouble(temp_coverage);                    
                    weight = (float[])temp_weight.Clone();                    
                    string report = "Err: " + Error + "; Iter: " + index +"; Cov: " + coverage + ";";
                    Debug.WriteLine(report);
                    Debug.Write(weight);
                    index++;
                    SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, null));
                    //RPworker.ReportProgress(2 * index);
                }
            } //END of WHILE LOOP
            this.coverage = coverage;
            //dosespace = ds;
            optimized = true;
        }
        
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
        
        #endregion
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
