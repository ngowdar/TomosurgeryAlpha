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
        public static int RasterWidth;
        public static int StepSize;
        public static int ComparisonKernelSize = 20;
        public int NumOfLines;
        public int NumOfShots;
        public double coverage;
        public int isovol;
        public int massvol;
        public int WhichSlice;
        public int[] boundaries;
        float[,] slice;
        public float[,] DDS_slice;
        public float[,] ModdedSlice; //The "addition" layer for dose inhomogeneity.
        public float[] dosespace;
        public float[] weight;
        public Slice info;
        bool horizontal = false; //By default all lines will be vertical.
        bool optimized = false;
        public PointF[][] shot_points; //A jagged array. One array of points for each line.
            //Note that this is meant to be direction independent (horizontal or vertical).
        public PointF[] shots;
        public string[] ShotType;
        #endregion

        #region Constructors
        public RasterPath(float[,] f)
        {
            SetParams(StepSize, RasterWidth);
            slice = f;
            ModdedSlice = f;
            PrepareDDSFromSlice(f);
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

        private void PrepareDDSFromSlice(float[,] slice)
        {
            DDS_slice = new float[slice.GetLength(0), slice.GetLength(1)];
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
                {
                    if (slice[i, j] == 0)
                        DDS_slice[i, j] = PathSet.ToleranceDose;
                    if (slice[i, j] > 0)
                        DDS_slice[i, j] = PathSet.RxDose;
                }
            Debug.WriteLine("DDS_Slice sum: " + Matrix.SumAll(DDS_slice));
            
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
            ShotType = new string[total];
            int counter = 0;
            for (int i = 0; i < shot_points.GetLength(0); i++)
            {                
                for (int j = 0; j < shot_points[i].GetLength(0); j++)
                {
                    giant[counter] = shot_points[i][j];                   

                    if (i == 0 || i == (shot_points.GetLength(0) - 1))
                    {
                        if (j == 0 || j == (shot_points[i].GetLength(0) - 1))
                            ShotType[counter] = "corner";
                        else
                            ShotType[counter] = "edge";
                    }
                    else if (j == 0 || j == (shot_points[i].GetLength(0) - 1))
                    {
                        ShotType[counter] = "edge";
                    }
                    else
                        ShotType[counter] = "middle";
                    
                    counter++;
                }
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
            ComparisonKernelSize = (int)Math.Round(stepsize*1.4);
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
        private float[] InitWeightArray(float min, float max)
        {
            weight = new float[NumOfShots];
            
            for (int i = 0; i < NumOfShots; i++)
            {
                if (ShotType[i] == "edge" || ShotType[i] == "corner")
                    weight[i] = min;
                else
                    weight[i] = max;
            }
            return weight;
        }
        private float[] InitWeightArray()
        {
            weight = new float[NumOfShots];

            for (int i = 0; i < NumOfShots; i++)
            {
                weight[i] = 1;
            }
            return weight;
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
                int edgepad = PathSet.shot_edgepadding;
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
            int edgepad = PathSet.line_edgepadding;
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
                InitWeightArray(1, 0.7f);
            int ds_x; int ds_y;
            dosespace = new float[StructureSet.BIG_dim[0] * StructureSet.BIG_dim[1]];
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
                        if (ds_x < 0 || ds_y < 0)
                            continue;
                        if (ds_x >= StructureSet.BIG_dim[0] || ds_y >= StructureSet.BIG_dim[1])
                            continue;
                        //Add the final result
                        index = ds_x + StructureSet.BIG_dim[0] * ds_y;
                        dosespace[index] += dosemidplane[i, j] * weight[k];
                    
                    }//);
                }
                       //NormalizeDose();
            dosespace = Matrix.Normalize(dosespace);
            Calculate2DCoverage(0.5f);
        }                
        #endregion        
        
        #region Shot Optimization methods (Step 1)
          #region Helper methods for OptimizeShotWeight (Step 1)
       


        /// <summary>
        /// Called from RPworker, which is called from within a slice loop in PathSet. 
        /// This method comprises the meat of Step 1 (shot-based) optimization, which
        /// finds shot weights within a single slice. Writes final result to the float[] ds matrix
        /// and calculates coverage too.
        /// </summary>        
        public void OptimizeShotWeights()
        {
            bool IsCoverageBeingOptimized = false; //These simply mean is the coverage approaching 1.
            bool IsOverageBeingOptimized;
            shots = ReturnSinglePoints();
            double Error = 1000; double coverage = 1; double overage = 1;
            int index = 0;
            float[] temp_weight = new float[weight.GetLength(0)];
            weight = InitWeightArray(1, 0.7f);

            float[] ds = new float[slice.GetLength(0) * slice.GetLength(1)];
            for (int i = 0; i < ds.GetLength(0); i++)
                ds[i] = 0.0f;
            dosemidplane = Matrix.Normalize(dosemidplane);
            ds = PrepareDS(ds, weight, 1.0f);
            ds = Normalize(ds);

            double[] m = CalculateIterationCoverage(ds, 0.5f); //TODO: Take this out.            
            WriteFloatArray2BMP(ds, "starting_ds.bmp");
            WriteFloatArray2BMP(DDS_slice, "starting_DDS.bmp");
            double old_error;
            while (Error >= .01 && index < 40)
            {
                
                temp_weight = ReoptimizeShotWeights(ds);
                old_error = Convert.ToDouble(Error);
                Error = FindError(weight, temp_weight);

                Debug.Assert(ds.Max() > 0);
                Debug.Assert(temp_weight.Max() > 0);     

                //Re-prepare the DoseSpace matrix with the new weights
                ds = ReviseDS(ds, temp_weight);
                WriteFloatArray2BMP(ds, string.Concat(index, "_ds.bmp"));
                //ds = Normalize(ds);

                double[] measurements = CalculateIterationCoverage(ds, 0.5f);
                //1st value = tumor voxel total
                //2nd value = isovolume voxel total
                //3rd value = both tumor & >iso voxel total
                //4th value = pixels underdosed

                //Calculate Coverage = (tumor pixels covered by rx dose) / (tumor pixels)
                double temp_coverage = measurements[2] / measurements[0];
                if (index > 0 && (1 / temp_coverage) >= (1 / coverage))
                    IsCoverageBeingOptimized = true;
                

                //Calculate Overage = (total pixels covered by rx dose) / (tumor pixels)                
                double temp_overage = measurements[1] / measurements[0];

                //Need to make sure that with each iteration, coverage isn't going down, and overage isn't going up.

                if (index > 2 && (temp_coverage < coverage || temp_overage > overage)) //Old version was index > 1 && cov < coverage
                {
                    string r = "CAUTION: Bad trend in coverage!!! " + coverage + " --> " + temp_coverage;
                    Debug.WriteLine(r);
                    Debug.WriteLine("Overage: " + overage + "-->" + temp_overage);
                    break;
                }
                if (index > 10 && old_error < Error)
                {
                    string r = "WARNING: Error is increasing (" + old_error + " --> " + Error + "). Terminating...";
                    Debug.WriteLine(r);
                    //break;
                }
                else
                {
                    coverage = Convert.ToDouble(temp_coverage);
                    overage = Convert.ToDouble(temp_overage);
                    string status = "=";
                    if (IsCoverageBeingOptimized)
                        status = "Good optimization!";
                    //Debug.WriteLine("===============================");
                    //Debug.WriteLine("ITERATION: " + index); 
                    //WriteArrayAsList("Weight", temp_weight);                    
                    //Debug.WriteLine("Error: " + Math.Round(Error,2)); 
                    //Debug.WriteLine("Coverage: " + temp_coverage);
                    //Debug.WriteLine("Overage: " + temp_overage);
                    //Debug.WriteLine("========" + status + "=========");
                    
                    
                    weight = (float[])temp_weight.Clone();
                    index++;
                    SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, null));
                    //RPworker.ReportProgress(2 * index);
                }
            } //END of WHILE LOOP

            this.coverage = coverage;
            //dosespace = ds;
            optimized = true;
        }

        private double FindError(float[] weight, float[] w)
        {
            double diff = 0;
            for (int i = 0; i < w.GetLength(0); i++)
            {
                diff += Math.Pow(Math.Abs((weight[i] - w[i])), 2);
            }
            return Math.Sqrt(diff);
        }

        private double EvalShotWeightIteration(float[] ds, PointF pf)
        {
            float[,] DStimesP = Matrix.Subset(ds, DDS_slice.GetLength(0), DDS_slice.GetLength(1), (int)pf.X, (int)pf.Y, ComparisonKernelSize);
            float[,] DDStimesP = Matrix.Subset(DDS_slice, (int)pf.X, (int)pf.Y, ComparisonKernelSize);
            float max = Matrix.FindMax(DDStimesP);           

            double[] measurements = FindWindowCoverage(DStimesP, DDStimesP, PathSet.RxDose, PathSet.ToleranceDose);
            
            
            //ratio, RxVolume, TumorVol, LesionRx, Uncovered

            double ratio;
            double simplesum = measurements[0];
            double RxVolvsTumor = Convert.ToDouble(measurements[1]/measurements[2]); //Isovolume / Tumorvolume
            double BothvsTumor = Convert.ToDouble(measurements[3]/measurements[2]); //Covered tumor / Totaltumor
            double BothvsRxVol = Convert.ToDouble(measurements[3]/measurements[1]);
            double Importance_Factor = measurements[2] / (ComparisonKernelSize * ComparisonKernelSize);
                   
            
            //double ratio =  (double)(ddssum / dssum);
            //if (ratio < 0.2)
            //{
            //    string d = "mismatch_DS_" + pf.X + "_" + pf.Y + ".bmp";
            //    string dd = "mismatch_DDS_" + pf.X + "_" + pf.Y + ".bmp";
            //    WriteFloatArray2BMP(DStimesP, d);
            //    WriteFloatArray2BMP(DDStimesP, dd);
            //}

            //if (BothvsTumor < 0.98)
            //{
            //    double ratio1 = (1.0 / measurements[2]);
            //    double ratio2 = measurements[0];
            //    if (ratio1 > ratio2)
            //        ratio = ratio1;
            //    else
            //        ratio = ratio2;
            //}                
            //else
            //ratio = (1.0 / BothvsTumor);
            ratio = simplesum;
            //Debug.Assert(ratio > 0);
            return ratio;
        }

        private double[] FindWindowCoverage(float[,] ds_window, float[,] dds_window, double iso, double ToleranceDose)
        {
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            float dose; float dds_value; double Uncovered = 0;
            int x = ds_window.GetLength(0); int y = ds_window.GetLength(1); int z = ds_window.GetLength(0);            
                for (int j = 0; j < ds_window.GetLength(1); j++)
                    for (int i = 0; i < ds_window.GetLength(0); i++)
                    {
                        dose = ds_window[i, j];
                        dds_value = ds_window[i, j];
                        if (dds_value > ToleranceDose && dose < 0.5)
                            Uncovered++;
                        if (dose >= iso)
                        {
                            RxVolume++;
                            if (dds_value > ToleranceDose)
                            {
                                TumorVol++;
                                LesionRx++;
                            }
                        }
                        else if (dose < iso)
                        {
                            if (dds_value > ToleranceDose)
                                TumorVol++;
                        }
                    }
                double LesionRxoverTumorVol = LesionRx / TumorVol;
                double LesionRxoverRxVol = LesionRx / RxVolume;
                double Underdosed = (Uncovered / TumorVol) * 100;
                double sum_sd = Matrix.SumAll(ds_window);
                double sum_dds = Matrix.SumAll(dds_window);
                double ratio = sum_dds / sum_sd;
                return new double[5] { ratio, RxVolume, TumorVol, LesionRx, Uncovered };            
        }

        private float[] ReoptimizeShotWeights(float[] ds)
        {
            float[] tweight = (float[])weight.Clone();
            
            //float[] DS = PrepareDS(ds,tweight,1.0f);            
            //WriteFloatArray2BMP(Matrix.Normalize(DS), "wholeDS.bmp");
            for (int shot = 0; shot < shots.GetLength(0); shot++)
            {
                PointF pf = shots[shot];
                
                //double[] temp_r = CompareSlices(DStimesP, DDStimesP, false); //This was the elaborate rule-based algorithm                
                double ratio = EvalShotWeightIteration(ds, pf); //This comparison function just adds the DDS and DS and compares for a simple ratio.
                
                tweight[shot] = (float)(weight[shot] * ratio); // old weight multiplied by newest ratio.                
                //Debug.WriteLine("Old: " + weight[shot] + " * R: " + ratio + " = New: " + tweight[shot]);

            }
            tweight = Normalize(tweight);
            return tweight;
        }

        private void WriteArrayAsList(string prefix, float[] f)
        {
            string output = prefix + ": [" + Math.Round(f[0],2);
            for (int i = 0; i < f.GetLength(0); i++)
                output += ", " + Math.Round(f[i], 2);
            output += "]";
            Debug.WriteLine(output);
        }


        /// <summary>
        /// Prepares the DS matrix with the most recent shot weights. Optional normalization value
        /// to make the max value less than 1 if necessary.
        /// </summary>
        /// <param name="ds"></param>
        /// <param name="weight"></param>
        /// <param name="NormalizeValue"></param>
        /// <returns></returns>
        public float[] PrepareDS(float[] ds, float[] weight, float NormalizeValue)
        {
            shots = ReturnSinglePoints();
            int ds_x; int ds_y;
            //float[] ds = new float[StructureSet.size * StructureSet.size];
            //for (int i = 0; i < ds.GetLength(0); i++)
            //    ds[i] = 0f;
            int index = 0;
            //The first two arrays loop through the midplane
            for (int j = 0; j < doseN; j++)
                for (int i = 0; i < doseN; i++)
                {
                    //This loop calculates the same pixel [i,j] for each shot.
                    //Parallel.For(0, shots.GetLength(0), k =>
                    for (int k = 0; k < shots.GetLength(0); k++)
                    {
                        //Finds the coordinates relative to dosespace
                        ds_x = (int)((PointF)shots[k]).X - ((doseN - 1) / 2) + i;
                        ds_y = (int)((PointF)shots[k]).Y - ((doseN - 1) / 2) + j;

                        //Add the final result
                        index = ds_x + (StructureSet.BIG_dim[0] * ds_y);
                        ds[index] += dosemidplane[i, j] * weight[k];

                    }//);
                }
            //ds = Matrix.ScalarMultiply(Matrix.Normalize(ds), NormalizeValue); //Make the highest value equal to 0.6 to allow for more growth.
            return ds;
        }

        public float[] ReviseDS(float[] ds, float[] recent_weight)
        {
            shots = ReturnSinglePoints();
            int ds_x; int ds_y;
            
            //Recreate dosematrix
            ds = Matrix.Zero1DFloat(ds.GetLength(0));
            int index = 0;

            //The first two arrays loop through the midplane
            for (int j = 0; j < doseN; j++)
                for (int i = 0; i < doseN; i++)
                {
                    //This loop calculates the same pixel [i,j] for each shot.
                    //Parallel.For(0, shots.GetLength(0), k =>
                    for (int k = 0; k < shots.GetLength(0); k++)
                    {
                        //Finds the coordinates relative to dosespace
                        ds_x = (int)((PointF)shots[k]).X - ((doseN - 1) / 2) + i;
                        ds_y = (int)((PointF)shots[k]).Y - ((doseN - 1) / 2) + j;

                        //Finds appropriate index in ds
                        index = ds_x + (StructureSet.BIG_dim[0] * ds_y);

                        /*Each pixel may have dose contributions from multiple shots. Therefore, we cannot apply a global 
                         * weight change by setting the pixel = to the newest weight. We also cannot simply add the new weight,
                         * because this wouldn't alter the weight that is already applied. 
                         * 
                         * Therefore, the old weighted contribution must be removed (leaving any other dose contributions intact)
                         * and then the recent weight applied. In other words, only the incremental weight difference is applied (+ or -)*/
                        ds[index] += (dosemidplane[i, j] * (recent_weight[k]));
                        if (ds[index] < 0)
                            ds[index] = (dosemidplane[i, j] * recent_weight[k]);



                    }//);
                }
            //ds = Matrix.ScalarMultiply(Matrix.Normalize(ds), NormalizeValue); //Make the highest value equal to 0.6 to allow for more growth.
            return ds;
        }

        /// <summary>
        /// Method to compare the shot-weighting iteration DStimesP to the desired version, DDStimesP. 
        /// Outsourced to a separate method to try techniques like my "card-counting" strategy.
        /// </summary>
        /// <param name="DStimesP">First 2D float, or the DoseSpace</param>
        /// <param name="DDStimesP">Second 2D float, or the Desired DoseSpace</param>
        /// <param name="slice">Boolean value indicating if this should return the sum or the ratio</param>
        /// <returns></returns>
        public static double[] CompareSlices(float[,] DStimesP, float[,] DDStimesP, bool slice)
        {
            double LesionVolume_tally = 0;
            double nontumortally = 0;
            double Overdosed_Nonlesion_tally = 0;
            double Underdosed_Lesion_tally = 0;
            double TotalRX_tally = 0;
            double LesionRX_Covered_tally = 0;
            double uncovertally = 0;

            //Default threshholds
            float rx = PathSet.RxDose;
            float tolerance_dose = 0.95f;            
            float tol_background_dose = 0.35f;
            int tumorpixels = 0; int nontumorpixels = 0;
            

            for (int j = 0; j < DStimesP.GetLength(1); j++)
                for (int i = 0; i < DStimesP.GetLength(0); i++)
                {
                    float PixelDose = DStimesP[i,j];
                    float PixelType = DDStimesP[i, j]; //i.e. Tumor or non-tumor
                    /* Ground rules (I totally made these up):
                     * 1) If we are on a tumor pixel, then evaluate whether or not the DS pixel is above or below the RXdose value.
                     * If it is below, it is a +2. If it is above, it is a -1. If it's greater than the tolerance dose (default 0.8),
                     * then it's -2.
                     * 2) If we are on a non-tumor pixel, if dose is above 0.2, it's -0.5. If it's above 0.3, it's -1. 
                     * 
                     * The net result is that the final tally represents how badly the current location needs more dose. 
                     * If the net is negative, then the final (total / number of pixels) is a good fractional weight
                     
                     * 
                     * Alternate strategy: my own conformity index:
                     * (LesionVolume + underdosed count) / (Isovolume + overdosed count)
                     * 
                     
                     */                    

                    if (PixelType >= rx) //i.e. tumor pixel
                    {
                        tumorpixels++;
                        LesionVolume_tally++;
                        if (PixelDose >= rx) //if dose is correct for the tumor pixel
                        {
                            LesionRX_Covered_tally++;
                            TotalRX_tally++;                            
                        }
                        else if (PixelDose < rx)
                        {
                            Underdosed_Lesion_tally++;                                                        
                        }
                    }
                    else if (PixelType < rx) //i.e. a nontumor pixel
                    {
                        nontumorpixels++;
                        if (PixelDose >= rx)
                        {
                            Overdosed_Nonlesion_tally++;
                            TotalRX_tally++;
                        }
                    }
                    
                    
                    //if (DDStimesP[i, j] >= rx) //i.e. tumor pixel
                    //{
                    //    tumorpixels++;
                    //    if (value < rx)
                    //        tumortally += 2*(rx-value); //Optimum so far: 1
                    //    else if (value > tolerance_dose)
                    //        tumortally -= (value-tolerance_dose);

                    //}
                    //else //i.e. a nontumor pixel
                    //{
                    //    nontumorpixels++;
                    //    if (value >= rx)
                    //        nontumortally -= 2*(value-rx);
                    //    else if (value > tol_background_dose)
                    //        nontumortally -= (value-tol_background_dose);
                    //}
                }
            double totalsum = LesionVolume_tally + nontumortally;
            Debug.Assert(Underdosed_Lesion_tally == (LesionVolume_tally - LesionRX_Covered_tally));
            if (slice)
            {
                return new double[5] { LesionVolume_tally, TotalRX_tally, LesionRX_Covered_tally, Underdosed_Lesion_tally, Overdosed_Nonlesion_tally };
                //return totalsum;
            }
            else
            {
                double ratio = 1;
                //if (totalsum <= 0)
                //{
                //    ratio = ((totalsum * (-1)) / (double)(tumorpixels + nontumorpixels)); //i.e. ratio should be between 0 and 1.
                //    //ratio = covertally / dosetally;
                //    //and if the total sum is negative and high, the ratio should approach zero (i.e. don't need no more dose.
                //    if (ratio < 0)
                //    {
                //        Debug.Assert(ratio > 0);
                //        ratio = 0;
                //    }

                //}
                //else if (totalsum > 0)
                //{
                //    //ratio = (double)(totalsum + tumorpixels + nontumorpixels) / (double)(tumorpixels + nontumorpixels); //ratio is over 1
                //    ratio = (LesionVolume_tally + Underdosed_Lesion_tally) / (TotalRX_tally + Overdosed_Nonlesion_tally);
                //}
                if (TotalRX_tally <= 0)
                    TotalRX_tally = 1.0;
                double lesioncov = LesionRX_Covered_tally / LesionVolume_tally;
                if (lesioncov >= 0.9)
                    ratio = lesioncov * (LesionRX_Covered_tally / TotalRX_tally); //VanReits
                else
                    ratio = lesioncov;
                //Using VanReits conformity index for ratio:
                //ratio = (LesionRX_Covered_tally / LesionVolume_tally) * (LesionRX_Covered_tally / TotalRX_tally);
                return new double[6] { ratio, LesionVolume_tally, TotalRX_tally, LesionRX_Covered_tally, Underdosed_Lesion_tally, Overdosed_Nonlesion_tally };
            }
        }


        private void WriteFloatArray2BMP(float[,] temp, string p)
        {
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            float[,] temp2 = (float[,])Matrix.Normalize(temp).Clone(); 
            int color = 0;
            Bitmap b = new Bitmap(temp.GetLength(0), temp.GetLength(1));
            for (int j = 0; j < temp.GetLength(1); j++)
                for (int i = 0; i < temp.GetLength(0); i++)
                {
                    color = (int)(temp2[i, j] * 255);
                    b.SetPixel(i, j, Color.FromArgb(color, color, color));
                }
            b.Save(path);
        }
        private void WriteFloatArray2BMP(float[] temp, string p)
        {
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            float[] temp2 = (float[])Matrix.Normalize(temp).Clone(); 
            int color = 0;
            //int size = (int)Math.Sqrt(temp.GetLength(0));
            Bitmap b = new Bitmap(X, Y);
            for (int j = 0; j < Y; j++)
                for (int i = 0; i < X; i++)
                {
                    color = (int)(temp2[i+(j*X)] * 255);
                    b.SetPixel(i, j, Color.FromArgb(color, color, color));
                }
            b.Save(path);
        }

        
//TODO: 
        public void CalculateAndSaveSliceDose(DoseKernel dk, int dosecalcthickness, string savepath)
        {
            Stopwatch s = new Stopwatch();
            s.Start();
            //PointF[] startingpoints = GetStartingPoints(shots);
            N = dk.DKI.Size;
            int xsize = DDS_slice.GetLength(0); int ysize = DDS_slice.GetLength(1);
            int xmid = xsize / 2; int ymid = ysize / 2; int zmid = dosecalcthickness / 2;
            int StartingDoseSlice = ((N - 1) / 2) - zmid;
            float[] slicedose = new float[xsize * ysize * dosecalcthickness];
            for (int k = 0; k < dosecalcthickness; k++)
                {
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            Parallel.For(0, shots.GetLength(0), (w) =>
                            {
                                PointF shot = shots[w];
                                PointF center = new PointF((N - 1) / 2, (N - 1) / 2);
                                //PointF FDP = FindFirstExistingDosePixel(shot, new PointF(xsize, ysize));
                                //PointF LDP = FindLastExistingDosePixel(shot, new PointF(xsize, ysize));
                                PointF FDP = new PointF(shot.X - center.X, shot.Y - center.Y);
                                float dose = dk.ReturnSpecificDoseValue(i, j, StartingDoseSlice + k) * weight[w];
                                int index = (k * xsize * ysize) + (((int)FDP.Y + j) * xsize) + ((int)FDP.X + i);
                                slicedose[index] += dose;
                                //
                            });
                }
            //s.Stop(); Debug.WriteLine("Calculate and save doses takes: " + s.ElapsedMilliseconds);
            //s.Reset(); s.Start();
            //float f = dk.ReturnSpecificDoseValue(80, 80, 80);
            WriteToFile(savepath, slicedose, xsize, ysize, dosecalcthickness);
            //s.Stop(); Debug.WriteLine("WriteToFile() takes: " + s.ElapsedMilliseconds);
        }

        public PointF FindFirstExistingDosePixel(PointF shot, PointF dims)
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
            //float XDistToTumorEdge = tumorsize.X - shot.X;
            //float YDistToTumorEdge = tumorsize.Y - shot.Y;

            float XDist = shot.X + center.X;
            float YDist = shot.Y + center.Y;

            if (XDist >= tumorsize.X) //last dosepixel outside tumor
                last.X = shot.X + center.X - XDist;
            else  //last dosepixel inside tumor
                last.X = N - 1;

            //Repeat above logic for Y coordinate
            if (XDist >= tumorsize.Y)
                last.Y = shot.Y + center.Y - YDist;
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

        public static float[,] GetDDS_Subset(float[,] DDS_slice, float px, float py, float[,] P)
        {
            //int startx = (int)px-(P.GetLength(0)-1)/2;
            //int starty = (int)py-(P.GetLength(1)-1)/2;
            float[,] output;

            // ADD DDS ADDITIONS HERE!!!
            if (DoseModifiable)
            {
                float[,] TempMod = Matrix.Add(mask, DDS_slice);
                output = Matrix.Subset(DDS_slice, (int)px, (int)py, ComparisonKernelSize);
                //output = Matrix.MultiplySubset(TempMod, P, (int)px, (int)py);
            }
            else
                output = Matrix.Subset(DDS_slice, (int)px, (int)py, ComparisonKernelSize);

            return output;
        }
        

        public static float[,] GetMultiplied_DS_Subset(float[] ds, float px, float py, float[,] P)
        {
            //int startx = (int)px - (P.GetLength(0) - 1) / 2;
            //int starty = (int)py - (P.GetLength(1) - 1) / 2;
            float[,] output = Matrix.MultiplySubset(ds, P, (int)px, (int)py, ComparisonKernelSize,ComparisonKernelSize);
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
        private float[] RepopulateSliceDS(float[] ds)
        {
            float[] norm_weight = Normalize(weight);
            ds = new float[ds.GetLength(0)];
            for (int k = 0; k < ds.GetLength(0); k++)
                ds[k] = 0.0f;
            for (int j = 0; j < dosemidplane.GetLength(1); j++)
                for (int i = 0; i < dosemidplane.GetLength(0); i++)
                    for (int whichshot = 0; whichshot < shots.GetLength(0); whichshot++)
                    {
                    //TODO: Make this loop parallel
                        //Find index pixel relative to center shot coordinate
                        int ds_x = (int)((PointF)shots[whichshot]).X - ((N - 1) / 2) + i;
                        int ds_y = (int)((PointF)shots[whichshot]).Y - ((N - 1) / 2) + j;

                        //Add the weighted dose pixel to the indexed location
                        ds[(ds_y * StructureSet.BIG_dim[0]) + ds_x] = dosemidplane[i, j] * norm_weight[whichshot];
                    }
            //ds = Normalize(ds);
            
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

        public double[] CalculateIterationCoverage(float[] DS, float iso)
        {
            //float[] ds = Normalize(DS);            
            double c = 0;
            float isovolume = 0; double both = 0; double tumor = 0; float dose; double uncovered = 0;
            float dds_value;
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
                {
                    dose = DS[i + (j * slice.GetLength(0))];
                    dds_value = slice[i, j];
                    if (dds_value > 0 && dose < 0.5)
                        uncovered++;
                    if (dose >= iso)
                    {
                        isovolume++;
                        if (dds_value > 0)
                        {
                            tumor++;
                            both++;
                        }
                    }
                    else if (dose < iso)
                    {
                        if (dds_value > 0)
                            tumor++;
                    }

                }
            c = (isovolume / tumor);
            return new double[4]{tumor, isovolume, both, uncovered};
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
