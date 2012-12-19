﻿using System;
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
        public static event ProgressChangedEventHandler ShotWeightsProgressHandler;
        public static event RunWorkerCompletedEventHandler SliceWorkerCompleted;
        public BackgroundWorker RPworker;
        public static int doseN = 161;     //Size of dose matrix
        public static int N;
        public static int X; public static int Y;
        public static bool DoseModifiable = false;
        public static float[,] mask;
        public static int RasterWidth;
        public static int StepSize;
        public static int ComparisonKernelSize;
        public int LineSidePadding;
        public int LineEdgePadding;
        public int[] Lines;
        public int NumOfLines;
        public int NumOfShots;
        public double coverage;
        public int isovol;
        public int massvol;
        public int WhichSlice;
        public int[] boundaries;
        float[,] slice;
        public ArrayList OverdosePoints;
        public float[,] DDS_slice;
        public float[,] ModdedSlice; //The "addition" layer for dose inhomogeneity.
        public float[] dosespace;
        public double[] ShotWeights;
        public Slice info;
        bool horizontal = false; //By default all lines will be vertical.
        bool optimized = false;
        public PointF[][] shot_points; //A jagged array. One array of points for each line.
            //Note that this is meant to be direction independent (horizontal or vertical).
        public PointF[] shots;
        public string[] ShotType;
        #endregion

        #region Constructors
        public RasterPath(float[,] tumor, float[,] combined, int stepsize, int rasterwidth, int edgepad, int sidepad)
        {
            SetParams(stepsize, rasterwidth, edgepad, sidepad);
            //ComparisonKernelSize = 120;
            //LineEdgePadding = (int)Math.Round(0.5 * StepSize);
            //LineSidePadding = (int)Math.Round(0.5 * StepSize);
            slice = tumor;
            ModdedSlice = tumor;
            //slice = PrepareDDSFromSlice(f);
            DDS_slice = PrepareDDSFromSlice(DilateSlice(DilateSlice(tumor)));
            WriteFloatArray2BMP(DDS_slice, "dilated_slice.bmp");
            WriteFloatArray2BMP(tumor, "non_dilatedslice.bmp");
            //DDS_slice = PrepareDDSFromSlice(DilateSlice(tumor));
            //DDS_slice = PrepareDDSFromSlice(f);
            
            //PrepareDDSFromSlice(f);
            X = tumor.GetLength(0); Y = tumor.GetLength(1);
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

        

        private float[,] PrepareDDSFromSlice(float[,] slice)
        {
            float[,] dds_slice = new float[slice.GetLength(0), slice.GetLength(1)];
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
                {
                    if (slice[i, j] == 0)
                        dds_slice[i, j] = PathSet.ToleranceDose;
                    else if (slice[i, j] == 1)
                        dds_slice[i, j] = PathSet.RxDose;
                    else if (slice[i, j] > 1)
                        dds_slice[i, j] = 0.01f;
                }
            return dds_slice;
            
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
        public void ChangeParamsUpdatePoints(int stepsize, int rasterwidth, int edgepad, int sidepad)
        {
            SetParams(stepsize, rasterwidth, edgepad, sidepad);
            FindAllShotPoints();
        }
        private void SetParams(int stepsize, int rasterwidth, int edgepad, int sidepad)
        {            
            StepSize = stepsize;            
            //ComparisonKernelSize = (int)Math.Round((decimal)stepsize);
            //ComparisonKernelSize = 80;
            //if (ComparisonKernelSize < 10)
            //    ComparisonKernelSize = 10;
            RasterWidth = rasterwidth;
            LineEdgePadding = edgepad;
            LineSidePadding = sidepad;            
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
            boundaries = FindSliceBoundaries(slice);
            //TODO: jhg
            //Get line positions for the slice
            if (Lines == null)
                Lines = LineSpacer(boundaries[0], boundaries[1], LineSidePadding);
            
                WriteArrayAsList("Line locations: ", (int[])Lines.Clone());
                shot_points = new PointF[Lines.GetLength(0)][];
                NumOfShots = 0;
                for (int i = 0; i < Lines.GetLength(0); i++)
                {
                    int[] ybounds = LineBoundaries(Lines[i]);
                    int[] lineshots = ShotSpacer(ybounds[0], ybounds[1]);
                    PointF[] PF_shots = new PointF[lineshots.GetLength(0)];
                    for (int m = 0; m < lineshots.GetLength(0); m++)
                    {
                        PF_shots[m] = new PointF(Lines[i], lineshots[m]);
                    }
                    shot_points[i] = PF_shots;
                    NumOfShots += lineshots.GetLength(0);
                }
            
        }
        private double[] InitWeightArray(float min, float max)
        {
            ShotWeights = new double[NumOfShots];
            
            for (int i = 0; i < NumOfShots; i++)
            {
                if (ShotType[i] == "edge" || ShotType[i] == "corner")
                    ShotWeights[i] = min;
                else
                    ShotWeights[i] = max;
            }
            return ShotWeights;
        }
        private double[] InitWeightArray()
        {
            ShotWeights = new double[NumOfShots];

            for (int i = 0; i < NumOfShots; i++)
            {
                ShotWeights[i] = 1;
            }
            return ShotWeights;
        }
        public int[] ShotSpacer(int ystart, int yend)
        {
            int[] shots;
            //if ((yend - ystart) < (1.75 * StepSize))
            //{
            //    shots = new int[1];
            //    shots[0] = (ystart + yend) / 2;
            //}
            //else
            //{

                //Find shot spacing
                int edgepad = LineEdgePadding;
                int meat; int first; int last;
                int numshots; int newspacing;
                //Rest of Logic goes here...
                meat = yend - ystart - edgepad*2;
                first = ystart + edgepad;
                last = yend - edgepad;
                if (last - first < edgepad)
                {
                    edgepad = (last - first);
                    first = ystart + edgepad;
                    last = yend - edgepad;
                    meat = last - first;
                }
                numshots = meat / StepSize;
                
                //If the meat section is cleanly divisible, -1 for the extra shot
                if (meat % StepSize == 0)
                    numshots = (meat / StepSize) - 1;
                //Find the new spacing
                newspacing = meat / (numshots + 1);
                if (newspacing < 2)
                    newspacing = 2;
                //Add the first and last shots to the total number of shots
                numshots += 2;
                //Respace shots by newspacing distance.
                shots = new int[numshots];
                shots[0] = first;
                shots[numshots - 1] = last;
                for (int i = 1; i < numshots - 1; i++)
                {
                    if (shots[i - 1] + newspacing > last)
                        break;
                    else
                        shots[i] = shots[i - 1] + newspacing;
                }
            //}
            //NumOfShots += shots.GetLength(0);
            return shots;

        }        
        /// <summary>
        /// Given a xstart and xend boundary of a slice, will calculate the number of lines.
        /// </summary>
        /// <param name="xstart"></param>
        /// <param name="xend"></param>
        /// <returns></returns>
        public int[] LineSpacer(int xstart, int xend, int edgepad)
        {            
            int[] lines;
            //int edgepad = LineSidePadding;
            int meat; int first; int last;
            int numlines; int newspacing;
            //Is there enough room for the two starting lines?
            if ((xend - xstart) < (1.25 * RasterWidth))
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

        public void Calculate2DDoseSpace(float[,] dosemidplane, double[] shotweights)
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
                InitWeightArray(1, 1.0f);
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
                        if (index < 0 || index > dosespace.GetLength(0))
                            continue;
                        else
                        {
                            if (shotweights == null)
                                dosespace[index] += (float)(dosemidplane[i, j] * ShotWeights[k]);
                            else
                                dosespace[index] += (float)(dosemidplane[i, j] * shotweights[k]);
                        }
                    
                    }//);
                }


              
            
            //Matrix.Normalize(ref dosespace);
            /* Instead of Normalizing dosespace (which reduces all shot weights)
             * search for local maxima, and reduce the involved shot weights 
             * proportionally to make the hotspot equal to 1.0. 
             * 
             * Two approaches: Global weight reduction or "top contributor" reduction.
             * 
             * Choice A: Global Weight Reduction
             * For instance, if two shots are 0.5 and 0.75 weighted, respectively,
             * and the maximum dose point is 1.5:
             * ==> 0.5*valueA + 0.75*valueB = 1.5
             * =====> k * (1.5) = 1.0 ==> k = 1 / 1.5 = 0.66667
             * =====> 0.6667 should be multiplied by shot A and shot B's weight.
             * 
             * Choice B: Top Contributor reduction
             * 1) Find the location of the hotspot pixel (i.e. the max dose)
             * 2) Find which shots contribute to the max dose
             * 3) Since the calculation is 161x161 pixels, chances are many of the shots will 
             * contribute, so sort the list and pick the biggest 3 contributors.
             * 4) Subtract the remaining contributing weights from 1.0 to find the total contribution
             * of the 3 top shots.
             * 5) I.E. if the dose was 1.5, and top 3 contribute 0.9, then need to reduce
             * the weighted sum of the top 3 shots to 0.4
             */

            Calculate2DCoverage(0.5f);
        }                
        #endregion        
        
        #region Shot Optimization methods (Step 1)
          
       


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
            double[] temp_weight = new double[ShotWeights.GetLength(0)];
            ShotWeights = InitWeightArray(1, 1);
            
            float[,] temp_slice = PrepareDDSFromSlice(slice);
            float[] ds = new float[slice.GetLength(0) * slice.GetLength(1)];
            for (int i = 0; i < ds.GetLength(0); i++)
                ds[i] = 0.0f;
            dosemidplane = Matrix.Normalize(dosemidplane);

            ds = PrepareDS(ShotWeights, 1.0f);

            double max = ds.Max();
            //ds = Normalize(ds);
            //double[] restricted_weights = new double[shots.GetLength(0)];

            //while (max > 1.0)
            //{
            //    double[] dummy_priority = new double[2] { 0, 0 };
            //    restricted_weights = ShotWeightPostProcess(ds, dummy_priority);
            //    ds = PrepareDS(ds, restricted_weights, 1.0f);
            //    weight = (double[])restricted_weights.Clone();
            //    max = ds.Max();
            //}

            //weight = (double[])restricted_weights.Clone();

            double[] SliceCoverageValues = CalculateIterationCoverage(ds, temp_slice, 0.5f); //TODO: Take this out.            
            //WriteFloatArray2BMP(ds, "starting_ds.bmp");
            //WriteFloatArray2BMP(DDS_slice, "starting_DDS.bmp");
            int MultiplierChoice = 2;
            if (index == 0)
                Error = 1000;

            bool BruteForce = false;

            double old_error = 1000;
            while (Error >= .015 && index < 40)
            {
                if (BruteForce == true)
                {
                    temp_weight = (double[])ShotWeights.Clone();
                    int max_index = Array.IndexOf(temp_weight, temp_weight.Max());
                    while (max > 1.0)
                    {
                        max_index = Array.IndexOf(temp_weight, temp_weight.Max());
                        temp_weight[max_index] = temp_weight[max_index]*0.9;
                        ds = PrepareDS(temp_weight, 1.0f);
                        ShotWeights = (double[])temp_weight.Clone();
                        max = ds.Max();
                    }
                    //for (int i = 0; i < temp_weight.GetLength(0); i++)
                    //    temp_weight[i] = temp_weight[i] * 0.92;
                }
                else
                    temp_weight = ReoptimizeShotWeights(ds);
                
                old_error = Convert.ToDouble(Error);

                //max = ds.Max();
                
                

                Debug.Assert(ds.Max() > 0);
                Debug.Assert(temp_weight.Max() > 0);     

                //Re-prepare the DoseSpace matrix with the new weights
                max = ds.Max();
                //ds = ReviseDS(ds, temp_weight);
                double[] measurements_before = CalculateIterationCoverage(ds, DDS_slice, 0.5f);
                ds = PrepareDS(temp_weight, 1.0f);
                //WriteFloatArray2BMP(ds, string.Concat(index, "_ds.bmp"));
                //ds = Normalize(ds);
                double max2 = ds.Max();
                Debug.WriteLine("Change in DS max dose value: " + Math.Round(max,2) + " ==> " + Math.Round(max2, 2));
                double[] measurements = CalculateIterationCoverage(ds, DDS_slice, 0.5f);
                //1st value = tumor voxel total
                //2nd value = isovolume voxel total
                //3rd value = both tumor & >iso voxel total
                //4th value = pixels underdosed

                //Calculate Coverage = (tumor pixels covered by rx dose) / (tumor pixels)
                double IterationCoverage = measurements[2] / measurements[0];
                double PercentUnderdosed = (measurements[3] / measurements[0]) * 100;
                if (index > 0 && (1 / IterationCoverage) >= (1 / coverage))
                    IsCoverageBeingOptimized = true;

                //Calculate Overage = (total pixels covered by rx dose) / (tumor pixels)                
                double temp_overage = measurements[1] / measurements[0];
                Error = FindError(ShotWeights, temp_weight);
                ShotWeights = (double[])temp_weight.Clone();
                if (index > 3 && IterationCoverage >= 1.0)
                {
                    //if the newest iteration max is > 1.0 and the max is increasing, switch to BruteForce methods
                    if (max2 > 1.0 && max2 > max)
                    {
                        BruteForce = true;
                        index++;
                        Error = 1;
                        continue;
                    }

                    //if BruteForce is on and max drops below 1, switch back to regular
                    else if (BruteForce == true && max2 <= 1.0)
                    {
                        BruteForce = false;
                        coverage = Convert.ToDouble(IterationCoverage);
                        ShotWeights = (double[])temp_weight.Clone();
                        overage = Convert.ToDouble(temp_overage);
                        index++;
                        SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, ShotWeights));
                        continue;
                    }
                    //if max is <= 1 and error is < 0.05, stop
                    else if (max2 <= 1.0 && Error < 0.05)
                        break;                    

                    else if (IterationCoverage < coverage && max2 <= 1.0) //coverage reversing/oscillating
                    {
                        Debug.WriteLine("Stopped bc coverage > 98%, coverage starting to decrease");
                        Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " + temp_overage);
                        index++;
                        break;
                    }

                    else if ((Math.Abs(old_error - Error) < .005) && max2 <= 1.0) //error isn't changing that much
                    {
                        Debug.WriteLine("Stopped bc error difference negligible");
                        Debug.WriteLine("Error: " + old_error + " --> " + Error);
                        Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " + temp_overage);
                        index++;
                        break;
                    }
                    else
                    {
                        coverage = Convert.ToDouble(IterationCoverage);
                        ShotWeights = (double[])temp_weight.Clone();
                        overage = Convert.ToDouble(temp_overage);
                        index++;
                        SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, ShotWeights));                    
                        continue;
                    }                    
                }
                else
                {
                    coverage = Convert.ToDouble(IterationCoverage);
                    ShotWeights = (double[])temp_weight.Clone();
                    overage = Convert.ToDouble(temp_overage);
                    string status = "=";
                    if (IsCoverageBeingOptimized)
                        status = "Good optimization!";
                    index++;
                    SliceWorkerProgressChanged.Invoke(null, new ProgressChangedEventArgs(2 * index, ShotWeights));                    
                        continue;
                }
            } //END of WHILE LOOP

            this.coverage = coverage;
            //dosespace = ds;
            optimized = true;
        }

        #region Helper methods for OptimizeShotWeight (Step 1)
        
        /// <summary>
        /// Given a global index location, this returns an array of the particular local
        /// index of each shot that contributes to the global location.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public PointF[] FindShotsRelatedToPoint(int index)
        {            
            PointF[] locations = new PointF[shots.GetLength(0)];
            for (int i = 0; i < shots.GetLength(0); i++)
            {
                PointF pf = GlobalToLocalShotIndex(shots[i], index);
                if (pf.IsEmpty)
                    continue;
                else
                    locations[i] = pf;                
            }
            return locations;
        }

        /// <summary>
        /// Given a global index location, this returns an array of the particular local
        /// index of each shot that contributes to the global location.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public PointF[] FindShotsRelatedToPoint(int x, int y)
        {
            int index = (y * X) + x;
            return FindShotsRelatedToPoint(index);
        }

        /// <summary>
        /// Called by FindShotsRelatedToPoint(), finds the relative location
        /// of a point in the coordinate system centered around a particular shot location.
        /// If relative location is outside the boundaries of the shot, then return
        /// a blank PointF.
        /// </summary>
        /// <param name="pf"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        private PointF GlobalToLocalShotIndex(PointF pf, int index)
        {
            int midlength = dosemidplane.GetLength(0) / 2;
            
                
            int startpoint = ((int)(pf.Y - midlength) * X) + ((int)pf.X - midlength);
            //Find start and stop boundaries of shot in global coordinates
            int startX = (int)pf.X - midlength;
            int startY = (int)pf.Y - midlength;
            if (startX < 0)
                startX = 0;
            if (startY < 0)
                startY = 0;
            int endX = (int)pf.X + midlength;
            int endY = (int)pf.Y + midlength;
            if (endX >= X)
                endX = X - 1;
            if (endY >= Y)
                endY = Y - 1;
            //Convert index to global x,y
            int shotx = (index % X);
            int shoty = (index - shotx) / X;
            int localx=0; int localy=0;
            //See if point is within shot boundaries
            if (shotx >= startX && shotx < endX)
                if (shoty >= startY && shoty < endY)
                {
                    localx = shotx - startX;
                    localy = shoty - startY;
                    return new PointF(localx, localy);
                }
                else
                    return new PointF();
            else
                return new PointF();
        }

        

        public double[] ShotWeightPostProcess(float[] ds, double[] iteration_weight)
        {
            float max = 0;
            int maxindex = 0;
            for (int j = 0; j < Y; j++)
                for (int i = 0; i < X; i++)
                {
                    float value = ds[j * X + i];
                    if (value > max)
                    {
                        max = value;
                        maxindex = (j * X) + i;
                    }
                }
            PointF[] RelativeIndices = FindShotsRelatedToPoint(maxindex);
            double[] ContributingWeights = new double[shots.GetLength(0)];
            double[] Values = new double[shots.GetLength(0)];
            int biggestcontributor = 0; double temp = 0;
            for (int i = 0; i < shots.GetLength(0); i++)
            {
                PointF pf = RelativeIndices[i];
                if (pf.IsEmpty)
                {
                    Values[i] = 0;
                    ContributingWeights[i] = 0;
                    continue;
                }
                else
                {
                    if (iteration_weight[0] > 0)
                        ContributingWeights[i] = iteration_weight[i];
                    else
                        ContributingWeights[i] = ShotWeights[i];

                    Values[i] = dosemidplane[(int)pf.X, (int)pf.Y];
                    if (Values[i] > temp)
                    {
                        temp = Values[i];
                        biggestcontributor = i;
                    }                    
                }
            }

            double sum = 0;
            for (int i = 0; i < Values.GetLength(0); i++)
                sum += Values[i];
            double multiplier = 1.0 / sum;
            double[] weights = new double[shots.GetLength(0)];
            double[] priorities = null;
            int ShotThatNeedsWeightReduction = biggestcontributor;
            if (iteration_weight[0] != 0)
            {
                priorities = new double[shots.GetLength(0)];
                double distance = 1000;
                for (int i = 0; i < shots.GetLength(0); i++)
                {
                     if (ContributingWeights[i] != 0)
                     {
                    PointF center = shots[i];
                    PointF pf = RelativeIndices[i];
                    double dist = Math.Sqrt(Math.Pow(pf.X - center.X, 2) + Math.Pow(pf.Y - center.Y, 2));
                    if (dist < distance)
                    { 
                        distance = dist;
                        priorities[i] = dist;
                    }
                     }                        
                    else
                        priorities[i] = 1000;
                }
                ShotThatNeedsWeightReduction = Array.IndexOf(priorities, priorities.Min());
            }

            for (int s = 0; s < shots.GetLength(0); s++)
            {
                if (ContributingWeights[s] != 0) //<- Check to see if this slice is important
                {
                    //If its the important slice, check if the adjacent slices are involved.
                    if (s == ShotThatNeedsWeightReduction)
                    {
                        double threshold = 0.3 * Values[s];
                        //Only play with weights if adjacent slice priorities are < 30% of the important slice
                        //if ((s + 1) < shots.GetLength(0) && Values[s + 1] <= threshold && Values[s + 1] != 0)
                        //    ContributingWeights[s + 1] = ContributingWeights[s + 1] * 0.9;
                        
                        //if ((s - 1) >= 0 && Values[s - 1] <= threshold && Values[s - 1] != 0)
                        //    weights[s - 1] = ContributingWeights[s - 1] * multiplier * 0.9;                        
                        
                        //weights[s] = ContributingWeights[s] * multiplier * 0.9;
                        weights[s] = ContributingWeights[s] * multiplier;
                    }
                    else
                        weights[s] = ContributingWeights[s] * multiplier;
                }
                else
                    weights[s] = ShotWeights[s];                
            }
            return weights;
        }

        private double FindError(double[] weight, double[] w)
        {
            double diff = 0;
            for (int i = 0; i < w.GetLength(0); i++)
            {
                diff += Math.Pow(Math.Abs((weight[i] - w[i])), 2);
            }
            return Math.Sqrt(diff);
        }

        private void CreateRedWindowBitmap(float[,] fp, PointF pf, string name)
        {
            float[,] FramePic = (float[,])fp.Clone();
            int halfwindow = ComparisonKernelSize / 2;
            int xstart = (int)pf.X - halfwindow;
            int xend = (int)pf.X + halfwindow;
            int ystart = (int)pf.Y - halfwindow;
            int yend = (int)pf.Y + halfwindow;
            if (xstart < 0)
                xstart = 0;
            else if (xend > DDS_slice.GetLength(0))
                xend = DDS_slice.GetLength(0) - 1;
            if (ystart < 0)
                ystart = 0;
            else if (yend > DDS_slice.GetLength(1))
                yend = DDS_slice.GetLength(1) - 1;

            for (int i = xstart; i <= xend; i++)
            {
                FramePic[i, ystart] = 100;
                FramePic[i, yend] = 100;
            }
            for (int j = ystart; j <= yend; j++)
            {
                FramePic[xstart, j] = 100;
                FramePic[xend, j] = 100;
            }
            string p = "OMG_" + name + "Currentwindow_at_" + pf.X + "_" + pf.Y + ".bmp";
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            int color = 0;
            Bitmap b = new Bitmap(FramePic.GetLength(0), FramePic.GetLength(1));
            for (int j = 0; j < FramePic.GetLength(1); j++)
                for (int i = 0; i < FramePic.GetLength(0); i++)
                {
                    if (FramePic[i, j] == 100)
                        b.SetPixel(i, j, Color.FromArgb(255, 0, 0));
                    else
                    {
                        color = (int)(FramePic[i, j] * 255);
                        if (color > 255)
                            color = 255;
                        b.SetPixel(i, j, Color.FromArgb(color, color, color));
                    }
                }
            b.Save(path);
        }

        private void CreateRedWindowBitmap(float[] ds, PointF pf, string name)
        {          

            float[,] FramePic = new float[DDS_slice.GetLength(0), DDS_slice.GetLength(1)];
            for (int j = 0; j < FramePic.GetLength(1); j++)
                for (int i = 0; i < FramePic.GetLength(0); i++)
                {
                    FramePic[i, j] = ds[(j * DDS_slice.GetLength(0)) + i];                                        
                }
            CreateRedWindowBitmap(FramePic, pf, name);            
        }

        private double EvalShotWeightIteration(float[] ds, PointF pf, int MultiplierChoice)
        {
            float[,] DStimesP = Matrix.ThresholdEq(Matrix.Subset(ds, DDS_slice.GetLength(0), DDS_slice.GetLength(1), (int)pf.X, (int)pf.Y, ComparisonKernelSize),0.5f);
            float[,] DDStimesP = Matrix.Subset(DDS_slice, (int)pf.X, (int)pf.Y, ComparisonKernelSize);
                  
            float max = 0; PointF maxpoint;

            double[] measurements = FindWindowCoverage(DStimesP, DDStimesP, PathSet.RxDose, PathSet.ToleranceDose);
            //ratio, RxVolume, TumorVol, LesionRx, Uncovered
            double DDStoDS_Sum_ratio = measurements[0];
            double RxVolume = measurements[1];
            double TumorVol = measurements[2];
            double BothVol = measurements[3];
            double Uncovered = measurements[4];


            //WriteFloatArray2BMP(DStimesP, "OMG_orig_DStimesP.bmp");
            //WriteFloatArray2BMP(DDStimesP, "OMG_orig_DDS.bmp");
            //CreateRedWindowBitmap(ds, pf, "DS");
            //WriteArrayAsList("Original measurements: ", measurements);
            double simplesum = measurements[0];
            float maximum = Matrix.FindMax(DStimesP);
            double RxVolvsTumor = Convert.ToDouble(measurements[1] / measurements[2]); //Isovolume / Tumorvolume
            double Coverage = Convert.ToDouble(measurements[3] / measurements[2]); //Covered tumor / Totaltumor
            double BothvsRxVol = Convert.ToDouble(measurements[3] / measurements[1]);
            double BothvsTumor = Convert.ToDouble(measurements[3] / measurements[2]);
            double Importance_Factor = measurements[2] / (ComparisonKernelSize * ComparisonKernelSize);
            double[] Conform_Indices = new double[4] { simplesum, RxVolvsTumor, BothvsTumor, BothvsRxVol };

            //Debug.WriteLine("Measurements are: ratio, RxVolume, TumorVol, LesionRx, Uncovered");

            //Trying to see what happens to measurements when (1.0 / BothVsTumor) is applied
            double[] measurements1 = FindWindowCoverage(Matrix.ScalarMultiply(DStimesP, (float)(1.0 / BothvsTumor)), DDStimesP, PathSet.RxDose, PathSet.ToleranceDose);            

            //When (1.0 / RxVolvsTumor) is applied...
            double[] measurements2 = FindWindowCoverage(Matrix.ScalarMultiply(DStimesP, (float)(1.0 / RxVolvsTumor)), DDStimesP, PathSet.RxDose, PathSet.ToleranceDose);
            


            for (int j = 0; j < DStimesP.GetLength(1); j++)
                for (int i = 0; i < DStimesP.GetLength(0); i++)
                    if (DStimesP[i, j] > max)
                    {
                        max = DStimesP[i, j];
                        maxpoint = new PointF(i, j);
                    }            
            
            double ratio = 1.0;            
            //if (BothvsTumor == 1.0 && RxVolvsTumor > 1.0)
            //    ratio = 1.0 / RxVolvsTumor;
            //else
            //    ratio = 1.0 / BothvsTumor;
            int which = CompareImprovements(measurements1, measurements2);
            switch (which)
            {
                case 1:
                    ratio = 1.0 / BothvsTumor;
                    break;
                case 2:
                    ratio = 1.0 / RxVolvsTumor;
                    break;
                case 3:
                    ratio = simplesum;
                    break;
            }           
            #region old
            //double ratio = simplesum;
            //if (max > 1.0)
            //    ratio = ratio;// *(1 / max);


            //if (Coverage < 0.9)
            //{
            //    double BiggestIndex = Conform_Indices.Max();
            //    if (Conform_Indices[MultiplierChoice] < BiggestIndex)
            //    {
            //        MultiplierChoice = Array.LastIndexOf(Conform_Indices, BiggestIndex);

            //        switch (MultiplierChoice)
            //        {
            //            case (0):
            //                ratio = simplesum;
            //                break;
            //            case (1):
            //                ratio = 1.0 / RxVolvsTumor;
            //                break;
            //            case (2):
            //                ratio = 1.0 / BothvsTumor;
            //                break;
            //            case (3):
            //                ratio = 1.0 / BothvsRxVol;
            //                break;
            //            default:
            //                ratio = 1.0 / BothvsTumor;
            //                break;
            //        }
            //    }
            //}
            //ratio = simplesum;
            //ratio = Convert.ToDouble(measurements[2] / measurements[3]); // 1 / (both / tumor)     
            #endregion
            return ratio;
        }

        private int CompareImprovements(double[] m1, double[] m2)
        {
            double sumratio1 = m1[0];
            double RxVolume1 = m1[1];
            double TumorVol1 = m1[2];
            double BothVol1= m1[3];
            double Uncovered1 = m1[4];
            double cov1 = m1[3] / m1[2];

            double sumratio2 = m2[0];
            double RxVolume2 = m2[1];
            double TumorVol2 = m2[2];
            double BothVol2 = m2[3];
            double Uncovered2 = m2[4];
            double cov2 = m2[3] / m2[2];

            if (cov1 > cov2)
                return 1;
            else if (cov2 > cov1)
                return 2;
            else
            {                
                if (sumratio1 <= sumratio2)
                    return 3;
                else
                    return 2;
            }
        }

        private double[] FindWindowCoverage(float[,] ds_window, float[,] dds_window, double iso, double ToleranceDose)
        {
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            double max = Matrix.FindMax(ds_window); double alt_iso = iso*max;
            float dose; float dds_value; double Uncovered = 0;
            int x = ds_window.GetLength(0); int y = ds_window.GetLength(1); int z = ds_window.GetLength(0);            
                for (int j = 0; j < ds_window.GetLength(1); j++)
                    for (int i = 0; i < ds_window.GetLength(0); i++)
                    {
                        dose = ds_window[i, j];
                        dds_value = dds_window[i, j];
                        if (dds_value > ToleranceDose)
                        {
                            TumorVol++;
                            if (dose >= alt_iso)
                            {
                                RxVolume++;
                                LesionRx++;
                            }
                            else if (dose < alt_iso)
                            {
                                Uncovered++;
                            }
                        }
                        else if (dds_value <= ToleranceDose)
                        {
                            if (dose >= alt_iso)
                                RxVolume++;
                        }
                    }



                    //    //Dose is > Rxdose?
                    //    if (dose >= iso)
                    //    {
                    //        RxVolume++; // <-- +1 to the RxVolume
                    //        //Check if pixel is a tumor pixel...
                    //        if (dds_value > ToleranceDose)
                    //        {
                    //            TumorVol++; // <-- if so, +1 to tumor volume and lesionrx
                    //            LesionRx++;
                    //        }
                    //    }
                    //    //If dose < iso, then doesn't get added to rxvolume
                    //    else if (dose < iso)
                    //    {
                    //        if (dds_value > ToleranceDose) // <-- still need to check tumor status.
                    //        {
                    //            TumorVol++; // <-- +1 to tumor volume
                    //            Uncovered++; // <-- +1 to "underdosed" voxel count
                    //        }
                    //    }
                    //}
                double BothvsTumorVol = LesionRx / TumorVol;
                double BothvsRxVol = LesionRx / RxVolume;
                double Underdosed = (Uncovered / TumorVol) * 100;
                double sum_sd = Matrix.SumAll(ds_window) / 2; //<- ds_window is thresholded to 0s and 1s. Halfing it will make them all 0.5.
                double sum_dds = Matrix.SumAll(dds_window);
                double min_ds = Matrix.FindMin(ds_window);
                double min_dds = Matrix.FindMin(dds_window);
                double ratio = sum_dds / sum_sd;

                
                //if (ratio < (1 / LesionRxoverTumorVol))
                //    ratio = (1.0 / LesionRxoverTumorVol);
                //else if (ratio < (1.0 / LesionRxoverRxVol))
                //    ratio = (1.0 / LesionRxoverRxVol);
                return new double[5] { ratio, RxVolume, TumorVol, LesionRx, Uncovered };            
        }

        private double[] ReoptimizeShotWeights(float[] ds)
        {
            double[] tweight = (double[])ShotWeights.Clone();
            int MultiplierChoice = 2;
            float max = ds.Max();
            //float[] DS = PrepareDS(ds,tweight,1.0f);            
            //WriteFloatArray2BMP(Matrix.Normalize(DS), "wholeDS.bmp");
            Debug.WriteLine("==========SLICE " + WhichSlice + "=========");
            for (int shot = 0; shot < shots.GetLength(0); shot++)
            {
                PointF pf = shots[shot];
                CreateRedWindowBitmap(DDS_slice, pf, "DDS");
                //double[] temp_r = CompareSlices(DStimesP, DDStimesP, false); //This was the elaborate rule-based algorithm                
                double ratio = EvalShotWeightIteration(ds, pf, MultiplierChoice); //This comparison function just adds the DDS and DS and compares for a simple ratio.

                //TODO: Limiting shotweight to 1.0
                //if (weight[shot] * ratio >= 1.0)
                //    tweight[shot] = 1.0f;
                //else
                if (ratio > 0)
                {
                    if (ratio > 1.0 && (ShotWeights[shot] * ratio) <= 1.0)
                        tweight[shot] = (ShotWeights[shot] * ratio);
                    else if (ratio > 1.0 && (ShotWeights[shot] * ratio) > 1.0)
                        tweight[shot] = 1.0;
                    else
                        tweight[shot] = ShotWeights[shot] * ratio;                    
                }
                 
                //else
                //    tweight[shot] = ShotWeights[shot]*ratio;
                // old weight multiplied by newest ratio.                
                Debug.WriteLine("Shot " + shot + ": " + Math.Round(ShotWeights[shot],2) + " ==> " + Math.Round(tweight[shot],2));
                
            }
            //tweight = Normalize(tweight);
            return tweight;
        }

        private void WriteArrayAsList(string prefix, decimal[] f)
        {
            string output = prefix + ": => \n 0: " + Math.Round((decimal)f[0],2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += "\n " + i + ": " + Math.Round((decimal)f[i], 2);
            //output += "]";
            Debug.WriteLine(output);
        }

        private void WriteArrayAsList(string prefix, double[] f)
        {
            string output = prefix + ": => \n 0: " + Math.Round((double)f[0], 2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += "\n " + i + ": " + Math.Round((double)f[i], 2);
            //output += "]";
            Debug.WriteLine(output);
        }

        private void WriteArrayAsList(string prefix, int[] f)
        {
            string output = prefix + ": [" + Math.Round((double)f[0],2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += ", " + Math.Round((double)f[i], 2);
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
        public float[] PrepareDS(double[] weight, float NormalizeValue)
        {
            //WriteFloatArray2BMP(ds, "before_ds.bmp");
            shots = ReturnSinglePoints();
            int ds_x; int ds_y;
            //ds = new float[StructureSet.size * StructureSet.size];
            float[] ds = new float[StructureSet.BIG_dim[0] * StructureSet.BIG_dim[1]];
            for (int i = 0; i < ds.GetLength(0); i++)
                ds[i] = 0.0f;
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
                        index = ds_x + (StructureSet.BIG_dim[0] * ds_y);
                        if (index < 0 || index >= ds.GetLength(0))
                            continue;
                        else
                        {
                            //Add the final result
                            float currentval = ds[index];
                            float newval = currentval + (float)(dosemidplane[i,j] * weight[k]);
                            ds[index] = newval;
                        }

                    }//);
                }
            //ds = Matrix.ScalarMultiply(Matrix.Normalize(ds), NormalizeValue); //Make the highest value equal to 0.6 to allow for more growth.
            //WriteFloatArray2BMP(ds, "ds.bmp");
            return ds;
        }

        public float[] ReviseDS(float[] ds, double[] recent_weight)
        {
            //shots = ReturnSinglePoints();
            int ds_x; int ds_y;
            
            //Recreate dosematrix
            float[] DS = new float[ds.GetLength(0)];
            DS = Matrix.Zero1DFloat(ds.GetLength(0));
            int index = 0;
            OverdosePoints = new ArrayList();
            //The first two arrays loop through the midplane
            for (int j = 0; j < doseN; j++)
                for (int i = 0; i < doseN; i++)
                {
                    //This loop calculates the same pixel [i,j] for each shot.
                    //Parallel.For(0, shots.GetLength(0), k =>
                    for (int k = 0; k < shots.GetLength(0); k++)
                    {
                        //Finds the coordinates relative to dosespace
                        ds_x = (int)(((PointF)shots[k]).X - ((doseN - 1) / 2) + i);
                        ds_y = (int)(((PointF)shots[k]).Y - ((doseN - 1) / 2) + j);

                        //Finds appropriate index in ds
                        index = ds_x + (StructureSet.BIG_dim[0] * ds_y);
                        if (index < 0 || index >= ds.GetLength(0))
                            continue;
                        else
                        {
                            /*Each pixel may have dose contributions from multiple shots. Therefore, we cannot apply a global 
                             * weight change by setting the pixel = to the newest weight. We also cannot simply add the new weight,
                             * because this wouldn't alter the weight that is already applied. 
                             * 
                             * Therefore, the old weighted contribution must be removed (leaving any other dose contributions intact)
                             * and then the recent weight applied. In other words, only the incremental weight difference is applied (+ or -)*/
                            float add = (float)(dosemidplane[i, j] * (recent_weight[k]));
                            float current_value = (float)ds[index];
                            float new_value = (float)(current_value + add);
                            if (new_value < 0)
                                DS[index] = (float)(current_value + add);
                            else if ((new_value) > 1.0)
                                OverdosePoints.Add(new double[]{index, k, ds[index], i, j, dosemidplane[i,j], recent_weight[k], ds[index]+add});
                            
                            DS[index] = new_value;
                        }

                    }  //);
                }
            //ds = Matrix.ScalarMultiply(Matrix.Normalize(ds), NormalizeValue); //Make the highest value equal to 0.6 to allow for more growth.
            /*OVERDOSE POINTS ARRAY LEGEND:
             * 0: Index
             * 1: Shot number
             * 2: Existing dose before addition
             * 3: i-coordinate in dose matrix
             * 4: j-coordinate in dose matrix
             * 5: Dose value
             * 6: Weight
             * 
             * 
             */
            OverdosePoints.TrimToSize();
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
            float[] temp2 = (float[])temp.Clone();
            Matrix.Normalize(ref temp2);
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
                            for (int w = 0; w < shots.GetLength(0); w++)
                            {
                                PointF shot = shots[w];
                                PointF center = new PointF((N - 1) / 2, (N - 1) / 2);
                                //PointF FDP = FindFirstExistingDosePixel(shot, new PointF(xsize, ysize));
                                //PointF LDP = FindLastExistingDosePixel(shot, new PointF(xsize, ysize));
                                PointF FDP = new PointF(shot.X - center.X, shot.Y - center.Y);
                                float dose;
                                int index = (k * xsize * ysize) + (((int)FDP.Y + j) * xsize) + ((int)FDP.X + i);
                                if (index < 0 || index >= slicedose.GetLength(0))
                                    continue;
                                else
                                {
                                    dose = (float)(dk.ReturnSpecificDoseValue(i, j, StartingDoseSlice + k) * ShotWeights[w]);
                                    slicedose[index] += dose;
                                }
                                    
                                //
                            }
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
            double[] norm_weight = Normalize(ShotWeights);
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
                        ds[(ds_y * StructureSet.BIG_dim[0]) + ds_x] = (float)(dosemidplane[i, j] * norm_weight[whichshot]);
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
            float[] dd = new float[d.GetLength(0)];
            float max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                if (d[i] > max)
                    max = d[i];
            for (int j = 0; j < d.GetLength(0); j++)
                dd[j] = (float)(d[j] / max);
            return dd;
        }
        private double[] Normalize(double[] d)
        {
            double max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                if (d[i] > max)
                    max = d[i];
            for (int j = 0; j < d.GetLength(0); j++)
                d[j] = (double)(d[j] / max);
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

        public float[,] DilateSlice(float[,] p)
        {
            float[,] output = (float[,])p.Clone();
            float top;
            float right;
            float left; float bottom;
            for (int j = 1; j < p.GetLength(1) - 1; j++)
                for (int i = 1; i < p.GetLength(0) - 1; i++)
                {
                    if (p[i, j] > 0)
                        continue;
                    else
                    {
                        top = p[i, j - 1];
                        left = p[i - 1, j];
                        right = p[i + 1, j];
                        bottom = p[i, j + 1];
                        
                        
                        float sum = top + left + right + bottom;
                        if (sum > 0)
                            output[i, j] = 1;
                    }
                }
            return output;
        }

        public double[] CalculateIterationCoverage(float[] ds, float[,] s, float iso)
        {
            /*Return double[] output:
             *0: tumor
             *1: isovolume
             *2: both
             *3: uncovered 
             */
            float max = ds.Max();
            //float[] DS = (float[])Normalize(ds).Clone();     //Just added 12/4       
            float[] DS = (float[])ds.Clone();
            float min = Matrix.FindMin(slice);
            max = DS.Max();
            double c = 0;
            
            float isovolume = 0; double both = 0; double tumor = 0; float dose; double uncovered = 0;
            float dds_value;
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
                {
                    dose = DS[i + (j * slice.GetLength(0))];
                    dds_value = slice[i, j];
                    if (dds_value > 0 && dose < iso)
                        uncovered++;
                    if (dose >= (iso*max))
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
