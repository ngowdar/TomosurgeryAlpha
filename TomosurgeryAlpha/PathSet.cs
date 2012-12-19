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
using OpenCLTemplate;
using Cloo;

namespace TomosurgeryAlpha
{
    

    /// <summary>
    /// A collection class that contains the information about a group of RasterPath objects.
    /// Takes tumor object, compresses it into slices and creates a RasterPath object 
    /// for each slice.
    /// </summary>
    public class PathSet
    {
        #region Variables        
        public Boolean DoseModifiable;
        public Boolean Dilate = true;
        public bool ShotsWeighted = false;
        public string folderpath = ActiveDirectory; //The active directory where slice data will be written.
        public DoseKernel DK; //'DoseKernel', the object representing the dose matrix.
        public StructureSet SS; //'SS' represents the binary tumor volume from the DICOM.
        public static float[,] mask;
        public static int N;
        public static int LineEdgePadding;
        public static int LineSidePadding;
        public static int tumorflag = 10;
        public static int CSflag = 2;
        public static float ToleranceDose = 0.1f;
        public static float RxDose = 0.5f;
        public static float CSdose = 0.05f;
        public static int NumSlices; //Set in step 1
        public static int TumorVolCount;
        public static int CSvolCount;
        public ArrayList RasterPaths; //Collection of objects representing each slice
        public static string ActiveDirectory;
        
        public static int DCT; //static version of dosecalcthickness
        public static int StepSize;
        public static int RasterWidth;
        public int SliceThickness;        
        public int TolThickness;
        public float[][,] volume; //Presumably the structure set binary volume?
        public float[][,] DDS;
        public float[][,] dds_nondilated;
        public float[][,] DoseSpace; //The actual final dosespace matrix containing the real dose.
        public float max;
        public int[] ReferencePoint;
        public double[] ReferenceWeights;
        public double[] ReferenceValues;
        public int[] SlicePositions; //The centered z-indexes of each slice location
        public double[] SliceWeights; //Weights, from 0 to 1, of each slice. Calculated from Step 2.
        public double[][] OldSliceCoverage;
        public double[][] SliceCoverage;
        public double[] OldCoverage;
        public double[] TotalCoverage;
        public double[] NonnormalizedCoverage;
        public int[] boundaries; //6-term sequence (xstart, xend, ystart...etc.)
        public int X; public int Y; public int Z;
        public event ProgressChangedEventHandler PathsetWorkerProgressChanged;
        public event RunWorkerCompletedEventHandler PathsetWorkerCompleted;
        public event ProgressChangedEventHandler RasterPathWorkerProgress;
        public event RunWorkerCompletedEventHandler OptimizationWorkerCompleted;
        public event ProgressChangedEventHandler OptimizationWorkerProgress;
        public event RunWorkerCompletedEventHandler SliceweightWorkerCompleted;
        public event ProgressChangedEventHandler SliceweightWorkerProgress;
        public BackgroundWorker PS_1_ShotOptimize_worker; // <- called when "Optimize" button is clicked 
        public BackgroundWorker PS_3_SliceWeightOpt_worker; 
        public BackgroundWorker PS_2_CalcDose_worker; //<- called when the Calc/Save/Dose button is clicked.
        #endregion




        #region Constructors
        public PathSet(float[][,] f, int sthick, int tolthick, int padding, DoseKernel dk, StructureSet ss)
        {
            
            X = f[0].GetLength(0); Y = f[0].GetLength(1); Z = f.GetLength(0);
            DK = dk; SS = ss;
            N = dk.dose.GetLength(0);
            boundaries = FindBoundaries(f);
            
            SliceThickness = sthick;
            DCT = SliceThickness * 2;
            TolThickness = tolthick;
            CalculateNumSlices(padding);
            
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(f, SlicePositions[i], SliceThickness / 2), CompressSection(ss.fj_Combined, SlicePositions[i], SliceThickness / 2), StepSize, RasterWidth, LineSidePadding, LineEdgePadding);
                rp.WhichSlice = i;
                RasterPaths.Add(rp);
            }            
            volume = f;
            AttachHandlers();
        }
        public PathSet(string dosespace_path)
        {
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1); int z = DoseSpace.GetLength(0);
            DoseSpace = GPU.BackTo3D(ReadDoseSpaceFromFile(dosespace_path),x,y,z);            
        }

        public void WritePathSetInfoToReport()
        {
            Analysis.AddLineToReport("==============PLAN SUMMARY================");
            Analysis.AddLineToReport("Slice Thickness: " + SliceThickness);
            Analysis.AddLineToReport("Dose Calculation Thickness: " + DCT);
            Analysis.AddLineToReport("Number of Slices: " + NumSlices);
            Analysis.AddLineToReport(WriteArrayAsList("Tumor boundaries: ", boundaries));

            Analysis.AddLineToReport("==========================================");
        }

        private void AttachHandlers()
        {

            PS_1_ShotOptimize_worker = new BackgroundWorker();
            PS_1_ShotOptimize_worker.WorkerReportsProgress = true;
            PS_1_ShotOptimize_worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(PS_1_ShotOptimize_RunWorkerCompleted);
            PS_1_ShotOptimize_worker.ProgressChanged += new ProgressChangedEventHandler(PS_1_ShotOptimize_ProgressChanged);
            PS_1_ShotOptimize_worker.DoWork += new DoWorkEventHandler(PS_1_ShotOptimize_DoWork);
            PS_3_SliceWeightOpt_worker = new BackgroundWorker();
            PS_3_SliceWeightOpt_worker.WorkerReportsProgress = true;
            PS_3_SliceWeightOpt_worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(PS_3_SliceOptimize_worker_RunWorkerCompleted);
            PS_3_SliceWeightOpt_worker.ProgressChanged += PS_3_SliceOptimize_worker_ProgressChanged;
            PS_3_SliceWeightOpt_worker.DoWork += PS_3_SliceOptimize_worker_DoWork;
            PS_2_CalcDose_worker = new BackgroundWorker();
            PS_2_CalcDose_worker.WorkerReportsProgress = true;
            PS_2_CalcDose_worker.RunWorkerCompleted += PS_2_CalcDose_worker_RunWorkerCompleted;
            PS_2_CalcDose_worker.ProgressChanged += PS_2_CalcDose_worker_ProgressChanged;
            PS_2_CalcDose_worker.DoWork += PS_2_CalcDose_worker_DoWork;
           
        }

        
        #endregion


        #region All Background Worker methods
        #region PS_1 Shot Optimization Background worker methods
        void PS_1_ShotOptimize_DoWork(object sender, DoWorkEventArgs e)
        {
            double count = 0;
            for (int i = 0; i < NumSlices; i++)
            {
                //((RasterPath)RasterPaths[i]).RPworker.RunWorkerAsync();
                Debug.WriteLine("SLICE " + Convert.ToString(i) + ":");
                ((RasterPath)RasterPaths[i]).OptimizeShotWeights();
                count++;
                PS_1_ShotOptimize_worker.ReportProgress((int)(100 * count / NumSlices));
            }
        }
        void PS_1_ShotOptimize_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (PathsetWorkerProgressChanged != null)
                PathsetWorkerProgressChanged.Invoke(null, e);
        }
        void PS_1_ShotOptimize_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            ShotsWeighted = true;
            if (PathsetWorkerCompleted != null)
                PathsetWorkerCompleted.Invoke(null, e);
            
        }
        #endregion

        void PS_2_CalcDose_worker_DoWork(object sender, DoWorkEventArgs e)
        {            
                
            CalculateAndWriteSliceDoses(DK, folderpath);            
            AssembleDoseSpaceFromFiles();
            //PS_2_CalcDose_worker.ReportProgress(50);
        }
        void PS_2_CalcDose_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (OptimizationWorkerProgress != null)
                OptimizationWorkerProgress.Invoke(null, e);
        }
        void PS_2_CalcDose_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            Debug.WriteLine("Shot weighting optimization complete. Review slice coverage before continuing...");
            if (OptimizationWorkerCompleted != null)
                OptimizationWorkerCompleted.Invoke(null, e);

            PS_3_SliceWeightOpt_worker.RunWorkerAsync();

        }
        void PS_3_SliceOptimize_worker_DoWork(object sender, DoWorkEventArgs e)
        {
            //PathsetWorkerProgressChanged += Opt_PathSet_PathsetWorkerProgressChanged;
            OptimizeSliceWeights();
        }

        
        void PS_3_SliceOptimize_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (SliceweightWorkerProgress != null)
                SliceweightWorkerProgress.Invoke(null, e);
        }
        void PS_3_SliceOptimize_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (SliceweightWorkerCompleted != null)
                SliceweightWorkerCompleted.Invoke(null, e);
        }
        

        void PathSet_SliceWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (RasterPathWorkerProgress != null)
                RasterPathWorkerProgress.Invoke(null, e);
        }
        
        #endregion

        #region STEP 1: Shot Optimization (per slice)
        /* TODO: Step-by-step walkthrough
         * 1) CreateDoseMatrix() is called
         * 2) 
         */

        /// <summary>
        /// Given a folder path and a DoseKernel, will create path names for each rasterpath in RasterPaths 
        /// and write each resulting slicedose to path specified.
        /// </summary>
        /// <param name="dk"></param>
        /// <param name="folderpath"></param>
        public void CalculateAndWriteSliceDoses(DoseKernel dk, string folderpath)
        {
            DK = dk; this.folderpath = folderpath;
            string path;
            string subfolder = System.IO.Path.Combine(folderpath, DateTime.Now.ToString("yyyyMMddHHmmssfff"));
            System.IO.Directory.CreateDirectory(subfolder);
            ActiveDirectory = subfolder;
            this.folderpath = subfolder;
            int progress = 0;
            for (int s = 0; s < NumSlices; s++)
            {
                string filename = string.Concat("slice_", s);
                path = System.IO.Path.Combine(subfolder, filename);
                RasterPath rp = (RasterPath)RasterPaths[s];
                rp.CalculateAndSaveSliceDose(dk, DCT, path);
                progress = (int)Math.Round((double)(50 / NumSlices));
                PS_2_CalcDose_worker.ReportProgress(progress);
            }
            //TODO: Need to have some kind of confirmation event?
        }
        /// <summary>
        /// Frontend method that calls PS_InitialDose_worker backgroundworker. 
        /// Background worker runs CalculateSliceDosesAndWrite() and AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="dk"></param>
        /// <param name="path"></param>
        public void CreateDoseMatrix(DoseKernel dk, string path)
        {
            this.folderpath = path;
            this.DK = dk;
            PS_2_CalcDose_worker.RunWorkerAsync();

        }
       
        /// <summary>
        /// Takes in a 1D float matrix, grabs a 2D slice out based on which z position. Called by GrabSlice().
        /// </summary>
        /// <param name="input"></param>
        /// <param name="which"></param>
        /// <param name="xsize"></param>
        /// <param name="ysize"></param>
        /// <returns></returns>
        public static float[,] GrabSlice(float[] input, int which, int xsize, int ysize)
        {
            float[,] result = new float[xsize, ysize];
            for (int j = 0; j < ysize; j++)
                for (int i = 0; i < xsize; i++)
                {
                    int index = (which * xsize * ysize) + (j * xsize) + i;
                    result[i, j] = input[index];
                }
            return result;
        }
       
        
        public void RecalculateSlices(int sthick, int tolthick)
        {
            SliceThickness = sthick;
            TolThickness = tolthick;
            CalculateNumSlices((int)Math.Round(SliceThickness*0.4));
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(SS.fj_Tumor, SlicePositions[i], SliceThickness / 2), CompressSection(SS.fj_Combined, SlicePositions[i], SliceThickness / 2), StepSize, RasterWidth, LineEdgePadding, LineSidePadding);
                RasterPaths.Add(rp);
            }
        }      
        public void CalculateNumSlices(int padding)
        {
            int remainder_slices; int new_spacing;            
            int[] zrasterpos;
            int zstart = boundaries[4];
            int zend = boundaries[5];
            int meat; int first; int last; int numslices; int newspacing;
            //Is there enough for 2 slices?
            if ((zend - zstart) < (1.75 * SliceThickness))
                zrasterpos = new int[1] { (zend + zstart) / 2 };
            else
            {
                meat = zend - zstart - (padding * 2);
                first = zstart + padding;
                last = zend - padding;
                numslices = meat / SliceThickness;
                if (meat % SliceThickness == 0)
                    numslices = (meat / SliceThickness) - 1;
                newspacing = meat / (numslices + 1);
                numslices += 2;

                zrasterpos = new int[numslices];
                zrasterpos[0] = first;
                zrasterpos[numslices - 1] = last;
                for (int i = 1; i < numslices - 1; i++)
                    zrasterpos[i] = zrasterpos[i - 1] + newspacing;
                this.SliceThickness = newspacing;
            }
            NumSlices = zrasterpos.GetLength(0);
            SlicePositions = zrasterpos;
            DCT = SliceThickness * 2;
            
            
            //int StartPosition = zstart + SliceThickness / 2;
            //int initial_slice_estimate = ((zend + SliceThickness / 4) - (zstart - SliceThickness / 4)) / SliceThickness;
            ////int initial_slice_estimate = (int)Math.Round((double)(zend - zstart) / SliceThickness);
            //int initial_remainder = ((zend + SliceThickness / 4) - (zstart - SliceThickness / 4)) % SliceThickness;

            ///*If there is a significant remainder, add a slice and shift to equalize.
            // * If there the remainder is small, just shift the slices */
            //if (initial_remainder >= 0.7 * SliceThickness)
            //{
            //    NumSlices = initial_slice_estimate + 1;
            //    zrasterpos = new int[NumSlices];
            //    padding = (SliceThickness - initial_remainder) / 2;
            //    //zrasterpos[0] = (zstart + (SliceThickness / 2)) - padding;
            //    zrasterpos[0] = (zstart + (SliceThickness / 4)) - padding;
            //    for (int i = 1; i < NumSlices; i++)
            //    {
            //        zrasterpos[i] = zrasterpos[0] + i * SliceThickness;
            //    }                

            //}
            //else
            //{
            //    NumSlices = initial_slice_estimate;
            //    padding = initial_remainder / 2;
            //    zrasterpos = new int[NumSlices];
            //    //zrasterpos[0] = (zstart + (SliceThickness / 2)) - padding;
            //    zrasterpos[0] = (zstart + (SliceThickness / 4)) - padding;
            //    for (int i = 1; i < NumSlices; i++)
            //    {
            //        zrasterpos[i] = zrasterpos[0] + i * SliceThickness;
            //    }                
            //}
            //SlicePositions = zrasterpos;
        }
        public int[] FindBoundaries(float[][,] f)
        {
            int[] boundaries = new int[6]; //xstart, xend, ystart, yend, zstart, zend
            double dsum = (double)Matrix.SumAll(f);
            //Z boundaries
            Boolean z1 = false; Boolean z2 = false;
            for (int k = 0; k < Z; k++)
            {
                float[,] temp = f[k];
                float sum = Matrix.SumAll(temp);
                if (sum > 0)
                {
                    if (z1 == false)
                    {
                        boundaries[4] = k;
                        z1 = true;
                    }
                    else
                        boundaries[5] = k;
                }


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
        
        
        #endregion

        #region STEP 2: Slice Optimization
        /* TODO: Step-by-step walkthrough
         * 1)
         * 2)
         * 3) 
         */

        /// <summary>
        /// Main method of Step 2. Using a error-minimization strategy, the weights are adjusted
        /// via DDS/DS ratios, just as in the shot-optimization of step 1.
        /// </summary>
        public void OptimizeSliceWeights()
        {
            //folderpath = ActiveDirectory;
            Debug.Assert(folderpath != null);
            double Error = 1000; int index = 0; double coverage = 0.8; double overage = 1;

            //Initialize the sliceweight matrices (temp_weights and SliceWeights)
            SliceWeights = new double[NumSlices];
            double[] temp_weights = new double[NumSlices];
            for (int i = 0; i < SliceWeights.GetLength(0); i++)
            {                
                SliceWeights[i] = 1.0;
                temp_weights[i] = 1.0;                
            }
            
            dds_nondilated = (float[][,])PrepareDDS(SS.fj_Tumor).Clone();
            DDS = (float[][,])PrepareEdgeEnhanced_DDS(SS.fj_Tumor).Clone();
            //DDS = (float[][,])PrepareDDS(SS.fj_Tumor).Clone();

            DoseSpace = PrepareWeighted_DS(SliceWeights, folderpath);
            float max = Matrix.FindMax(DoseSpace);
            double[] restricted_weights = new double[NumSlices];
            
            //while (max > 1.0)
            //{
            //    double[] dummy_priority = new double[2]{0,0};
            //    restricted_weights = SliceWeightPostProcess(DoseSpace, folderpath, dummy_priority);
            //    DoseSpace = PrepareWeighted_DS(restricted_weights, folderpath);
            //    SliceWeights = (double[])restricted_weights.Clone();
            //    max = Matrix.FindMax(DoseSpace);
            //}
            
            //SliceWeights = (double[])restricted_weights.Clone();
            //DoseSpace = PrepareWeighted_DS_GPU(SliceWeights, folderpath, DDS);
            //WriteFloatArray2BMP(DoseSpace[DoseSpace.GetLength(0) / 2], "PS_487_ds_midplane_" + index + ".bmp");
            
            TotalCoverage = new double[5];
            NonnormalizedCoverage = new double[5];
            NonnormalizedCoverage = FindTotalCoverage(0.5, DoseSpace, dds_nondilated, SliceWeights);
            
            
            //{ TumorVol, RxVolume, LesionRx, Uncovered };

            int MultiplierChoice = 1;
            if (index == 0)
                Error = 1000;
            double old_error = 1000;

            while (Error >= .01 && index < 20)
            {
                temp_weights = ReOptimizeSliceWeights(dds_nondilated, MultiplierChoice);
                //for (int i = 0; i < temp_weights.GetLength(0); i++)
                //    if (temp_weights[i] > restricted_weights[i])
                //        temp_weights[i] = restricted_weights[i];

                old_error = Convert.ToDouble(Error);

                //Re-prepare dosespace with latest iteration of sliceweight
                ClearDosespace();
                DoseSpace = PrepareWeighted_DS(temp_weights, folderpath); //TODO: This method is time-consuming, make this GPU?
                
                //Restrict Slice Weights to avoid overdosing
                max = Matrix.FindMax(DoseSpace);


                while (max > 1.0)
                {
                    restricted_weights = SliceWeightPostProcess(DoseSpace, folderpath, temp_weights);
                    DoseSpace = PrepareWeighted_DS(restricted_weights, folderpath);
                    temp_weights = (double[])restricted_weights.Clone();
                    max = Matrix.FindMax(DoseSpace);
                }

                //if (index > 0)
                    Error = FindError(SliceWeights, temp_weights);
                
                
                
                //DoseSpace = PrepareWeighted_DS_GPU(temp_weights, folderpath, DDS);
                //WriteFloatArray2BMP(DoseSpace[DoseSpace.GetLength(0) / 2], "PS_508_ds_midplane_" + index + ".bmp");
                //DoseSpace = Matrix.Normalize(DoseSpace);            
                TotalCoverage = FindTotalCoverage(RxDose, DoseSpace, dds_nondilated, temp_weights);
                //FakeCoverage = FindTotalCoverage(RxDose, DoseSpace, DDS, temp_weights);
                /*
                 * TumorVol, RxVolume, LesionRx, Uncovered, Overdosed
                 * 0: Tumor Volume
                 * 1: Rx Volume
                 * 2: Tumor Vol Covered by Rx
                 * 3: Underdosed voxels
                 * 4: Overdosed pixels
                 */
                PrintIterationSummary(OldCoverage, TotalCoverage, index, Error);
                double RTOG = TotalCoverage[1] / TotalCoverage[0];
                double BothvsTumor = TotalCoverage[2] / TotalCoverage[0];
                double BothvsRx = TotalCoverage[2] / TotalCoverage[1];
                double[] Conform_Indices = new double[3] { (1/RTOG), (1/BothvsTumor), (1/BothvsRx) };
                double IterationCoverage = TotalCoverage[2] / TotalCoverage[0];
                
                double PercentUnderdosed = (TotalCoverage[3] / TotalCoverage[0]) * 100;
                double TempOverage = TotalCoverage[1]/TotalCoverage[0];
                //index++;
                if (index > 1)
                {
                    if (IterationCoverage >= 1.0 && IterationCoverage < coverage) //coverage reversing/oscillating
                    {
                        Debug.WriteLine("Stopped bc coverage > 98%, coverage starting to decrease");
                        Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " + TempOverage);
                        break;
                    }
                    
                    else if ((Math.Abs(old_error - Error) < .005)) //error isn't changing that much
                    {
                        if (IterationCoverage < 0.9)
                        {
                            double BiggestIndex = Conform_Indices.Max();
                            if (Conform_Indices[MultiplierChoice] < BiggestIndex)
                                MultiplierChoice = Array.LastIndexOf(Conform_Indices, BiggestIndex);
                            Debug.WriteLine("Reached error asymptote with current index, coverage not great.");
                            Debug.WriteLine("SWITCHING INDEX! MultiplierChoice: " + MultiplierChoice);
                            continue;
                        }
                        else
                        {
                            Debug.WriteLine("Stopped bc error difference negligible");
                            Debug.WriteLine("Error: " + old_error + " --> " + Error);
                            Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " + TempOverage);
                            break;
                        }
                    }
                    //else if (IterationCoverage < coverage)
                    //{
                    //    if (IterationCoverage < 0.9)
                    //    {
                    //        double BiggestIndex = Conform_Indices.Max();
                    //        if (Conform_Indices[MultiplierChoice] < BiggestIndex)
                    //        {
                    //            MultiplierChoice = Array.LastIndexOf(Conform_Indices, BiggestIndex);
                    //            Debug.WriteLine("Reached error asymptote with current index, coverage not great.");
                    //            Debug.WriteLine("SWITCHING INDEX! MultiplierChoice: " + MultiplierChoice);
                    //            coverage = Convert.ToDouble(IterationCoverage);
                    //            continue;
                    //        }
                    //        else
                    //        {
                    //            Debug.WriteLine("Best possible. Still sucks. Try improving shot coverage first.");
                    //            break;
                    //        }
                    //    }
                    //    else
                    //    {
                    //        Debug.WriteLine("Stopped bc coverage starting to reverse");
                    //        Debug.WriteLine("Error: " + old_error + " --> " + Error);
                    //        Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " + TempOverage);
                    //        break;
                    //    }
                    //}
                    else
                    {
                        coverage = Convert.ToDouble(IterationCoverage);
                        SliceWeights = (double[])temp_weights.Clone();
                        index++;
                        PS_3_SliceWeightOpt_worker.ReportProgress(index * 5);
                        continue;
                    }
                }
                else
                {
                    coverage = Convert.ToDouble(IterationCoverage);
                    SliceWeights = (double[])temp_weights.Clone();
                    index++;
                    PS_3_SliceWeightOpt_worker.ReportProgress(index * 5);
                    continue;
                }
            } // <- end of while loop
            Debug.WriteLine("While Loop ended. Error: " + Error + "; Index: " + index);
            PS_3_SliceWeightOpt_worker.ReportProgress(100);
            FindMaxDose();
        }

        private void PrintIterationSummary(double[][] OldCoverage, double[][] SliceCoverage)
        {
            for (int i = 0; i < NumSlices; i++)
            {
                Debug.WriteLine("====================================================");
                Debug.WriteLine("====================SLICE " + (i + 1) + "===================");
                Debug.WriteLine("Sum Ratio: " + Math.Round(OldCoverage[i][0], 2) + " -> " + Math.Round(SliceCoverage[i][0], 2));
                Debug.WriteLine("Isovolume: " + Math.Round(OldCoverage[i][1], 2) + " -> " + Math.Round(SliceCoverage[i][1], 2));
                Debug.WriteLine("Tumor volume coverage: " + Math.Round(100 * (OldCoverage[i][2] / OldCoverage[i][3]), 2) + "% -> " + Math.Round(100 * (SliceCoverage[i][2] / SliceCoverage[i][3]), 2) + "%");
                Debug.WriteLine("Underdosed %: " + Math.Round(100 * (OldCoverage[i][4] / OldCoverage[i][2]), 2) + "% -> " + Math.Round(100 * (SliceCoverage[i][4] / SliceCoverage[i][2]), 2) + "%");
                Debug.WriteLine("====================================================+");
            }
        }
        private void PrintIterationSummary(double[] OldCoverage, double[] TotalCoverage, int iteration, double error)
        {
            
                Debug.WriteLine("====================================================");
                Debug.WriteLine("====================Iteration " + iteration + "===================");
                //Debug.WriteLine("Sum Ratio: " + Math.Round(OldCoverage[0], 2) + " -> " + Math.Round(TotalCoverage[0], 2));
                Debug.WriteLine("Isovolume: " + Math.Round(OldCoverage[1], 2) + " -> " + Math.Round(TotalCoverage[1], 2));
                Debug.WriteLine("Tumor volume coverage: " + Math.Round(100 * (OldCoverage[2] / OldCoverage[0]), 2) + "% -> " + Math.Round(100 * (TotalCoverage[2] / TotalCoverage[0]), 2) + "%");
                Debug.WriteLine("Underdosed %: " + Math.Round(100 * (OldCoverage[3] / OldCoverage[0]), 2) + "% -> " + Math.Round(100 * (TotalCoverage[3] / TotalCoverage[0]), 2) + "%");
                Debug.WriteLine("Overdosed voxels: " + Math.Round(OldCoverage[4], 2) + " -> " + Math.Round(TotalCoverage[4], 2));
                Debug.WriteLine("ERROR: " + Math.Round(error, 2));
                Debug.WriteLine("====================================================+");
            
        }

        private string WriteArrayAsList(string prefix, double[] f)
        {
            string output = prefix + ": [" + Math.Round(f[0], 2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += ", " + Math.Round(f[i], 2);
            output += "]";
            Debug.WriteLine(output);
            return output;
        }

        private string WriteArrayAsList(string prefix, int[] f)
        {
            string output = prefix + ": [" + f[0];
            for (int i = 1; i < f.GetLength(0); i++)
                output += ", " + f[i];
            output += "]";
            Debug.WriteLine(output);
            return output;
        }

        private void ReviseWeighted_DS(double[] recent_weights, double[] old_weights, string folderpath)
        {
           if (GPU.GPUenabled)
            {
                Stopwatch gputimer = new Stopwatch();
                gputimer.Start();
                int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1); int z = DoseSpace.GetLength(0);
                float[] weighted_slicedose = new float[NumSlices * DCT * DoseSpace[0].GetLength(0) * DoseSpace[0].GetLength(1)];
                
                GPU.PrepareDoseSpace(weighted_slicedose, SlicePositions, recent_weights, new int[3] { x, y, z }, DCT, folderpath);
                gputimer.Stop();
                Debug.WriteLine("ReviseWeighted_DS takes " + gputimer.Elapsed);
                //TODO: Uncoment below==============
            }
           else
           {
               Parallel.For(0, NumSlices, (s) =>
               {
                   float[] slicedose = LoadSliceDose(s, folderpath);
                   //Debug.WriteLine(slicedose.Sum());
                   slicedose = Matrix.ScalarMultiply(slicedose, (float)recent_weights[s]);
                   ClearDosespace();
                   DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);
               });
           }
                //==========================================
        }

        private float[][,] PrepareWeighted_DS_Dynamic(double[] SliceWeights, float[][,] ods)
        {
            if (GPU.GPUenabled)
            {
                for (int z = 0; z < NumSlices; z++)
                {
                    int startz = FindStartZ(z);
                    int endz = FindEndZ(startz);
                    for (int k = 0; k < DCT; k++)
                    {
                        ods[k + startz] = Matrix.Add(ods[k + startz], GPU.ScalarMultiply(ods[k + startz], (float)SliceWeights[z]));
                    }
                }
            }
            else
            {
                Parallel.For(0, NumSlices, (z) =>
                    {
                        int startz = FindStartZ(z);
                        int endz = FindEndZ(startz);
                        for (int k = 0; k < DCT; k++)
                        {
                            ods[k + startz] = Matrix.Add(ods[k + startz], Matrix.ScalarMultiply(ods[k + startz], (float)SliceWeights[z]));
                        }

                    });
            }
            return ods;
        }

        private double[] ReOptimizeSliceWeights(float[][,] dds, int MultiplierChoice)
        {
           double[] tweight = (double[])SliceWeights.Clone();
            if (SliceCoverage != null)
                OldSliceCoverage = (double[][])SliceCoverage.Clone();
            else
                SliceCoverage = new double[NumSlices][];
           int count = 0;

           for (int s = 0; s < NumSlices; s++)
           {
               float[][,] ds_slab = GrabSlab(DoseSpace,PathSet.DCT, SlicePositions[s]);
               //WriteFloatArray2BMP(ds_slab[ds_slab.GetLength(0)/2],String.Concat(s,"_ds_slab.bmp"));
               float[][,] dds_slab = GrabSlab(dds, PathSet.DCT, SlicePositions[s]);
               //WriteFloatArray2BMP(dds_slab[dds_slab.GetLength(0)/2], String.Concat(s, "_dds_slab.bmp"));
               double[] measurements = FindSliceDoseCoverage(Matrix.ThresholdEq(ds_slab, 0.5f),s,0.5,dds_slab);

               /* 0: ratio
                * 1: RxVol
                * 2: TumorVol
                * 3: Both (LesionRx)
                * 4: Underdosed
                * 5: Overdosed
                */ 
              
               double simplesum_ratio = measurements[0];
               double RxVolvsTumor = measurements[1]/measurements[2];
               double BothvsTumor = measurements[3]/measurements[2];
               double BothvsRxVol = measurements[3]/measurements[1];
               double Overdosed = measurements[5];
               float bvt = (float)(1.0 / BothvsTumor);
               float simplesum = (float)Math.Round(simplesum_ratio, 3);
               float rvt = (float)(1.0 / RxVolvsTumor);

               double[] m1 = FindSliceDoseCoverage(Matrix.ThresholdEq(Matrix.ScalarMultiply(ds_slab, bvt), 0.5f), s, 0.5, dds_slab);
               double[] m2 = FindSliceDoseCoverage(Matrix.ThresholdEq(Matrix.ScalarMultiply(ds_slab, rvt), 0.5f), s, 0.5, dds_slab);

               MultiplierChoice = CompareImprovements(m1, m2);

               double ratio;
               switch (MultiplierChoice)
               {   
                   case(1):
                       ratio = bvt;
                       break;
                   case(2):
                       ratio = rvt;
                       break;
                   case(3):
                       ratio = simplesum_ratio;
                       break;
                   default:
                       ratio = 1.0 / BothvsTumor;
                       break;
               }
               //double ratio = 1.0 / RxVolvsTumor;
               //double ratio = simplesum_ratio;               
               Debug.WriteLine("Slice: " + s + " Tumor coverage: " + BothvsTumor);
               Debug.WriteLine("(Isovol/Tumor)x(Isotumor/Tumor): " + ratio + "; 1//BothVsTumor: " + (1.0 / BothvsTumor) + "; 1//RxVolVsTumor: " + (1.0 / RxVolvsTumor));
               Debug.WriteLine("Currently using BothVsTumor");
               
               //if (SliceWeights[s] * ratio > 1.0)
               //{
               //    tweight[s] = 1.0;
               //    continue;
               //    Debug.WriteLine("Ratio greater than 1: " + ratio);
               //}
               //else
               if (ratio > 0)
               {
                   if (ratio > 1.0 && (SliceWeights[s] * ratio) <= 1.0)
                       tweight[s] = SliceWeights[s] * ratio;
                   else if (ratio > 1.0 && (SliceWeights[s] * ratio) > 1.0)
                       tweight[s] = SliceWeights[s] * ratio;
                   else
                       tweight[s] = SliceWeights[s] * ratio;
               }               
               
               SliceCoverage[s] = measurements;
               count++;
           }
           //tweight = Normalize(tweight);
           return tweight;
        }

        private int CompareImprovements(double[] m1, double[] m2)
        {
            double sumratio1 = m1[0];
            double RxVolume1 = m1[1];
            double TumorVol1 = m1[2];
            double BothVol1 = m1[3];
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

        private double[] Normalize(double[] d)
        {
           double max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                if (d[i] > max)
                    max = d[i];
            for (int j = 0; j < d.GetLength(0); j++)
                d[j] = (float)(d[j] / max);
            return d;
        }

        
        public static float[][,] PrepareDDS(float[][,] dds_slice)
        {
            float[][,] pDDS = new float[dds_slice.GetLength(0)][,];
            TumorVolCount = 0; CSvolCount = 0;
            int _debug_tumorcount = 0;
            for (int k = 0; k < dds_slice.GetLength(0); k++)
            {
                pDDS[k] = new float[dds_slice[0].GetLength(0), dds_slice[0].GetLength(1)];
                for (int j = 0; j < dds_slice[0].GetLength(1); j++)
                    for (int i = 0; i < dds_slice[0].GetLength(0); i++)
                    {
                        Debug.Assert(dds_slice[k][i, j] < 2);
                        if (dds_slice[k][i, j] <= (ToleranceDose*2))
                            pDDS[k][i, j] = ToleranceDose;                        
                        else if (dds_slice[k][i,j] > (ToleranceDose*2))
                        {
                            TumorVolCount++;
                            pDDS[k][i, j] = RxDose;
                        }
                    }
            }
            //Debug.WriteLine("DDS counts >= 2: " + _debug_tumorcount);
            //Debug.WriteLine("DDS tumor counts: " + TumorVolCount);
            return pDDS;
        }

        public float[][,] PrepareEdgeEnhanced_DDS(float[][,] dds)
        {
            float[][,] EE_DDS;
            double sum = Matrix.SumAll(dds);

            if (Dilate == true)
            {
                EE_DDS = PrepareDDS(DilateDDS(dds));
                EE_DDS = PrepareDDS(DilateDDS(EE_DDS));
                EE_DDS = PrepareDDS(DilateDDS(EE_DDS));                
            }
            else
            {
                EE_DDS = PrepareDDS(dds);                
            }
            return EE_DDS;
        }

        private float[][,] DilateDDS(float[][,] p)
        {
            float[][,] output = new float[p.GetLength(0)][,];
            /*
             *  0          A (i,j-1)          0
             *  B (i-1,j)     C(i,j)          D (i+1,j)
             *  0           E(i, j+1)        0
             * 
             * prev slice = (k-1,i,j), next slice (k+1,i,j)
             */

            float top; 
            float right; 
            float left; float bottom; float back; float front;
            //Starting at 1 instead of 0-index, to save having to worry about edges. Tumor wont go that far anyway.
            output[0] = (float[,])p[0].Clone();
            output[p.GetLength(0) - 1] = (float[,])p[0].Clone();

            for (int k = 1; k < p.GetLength(0)-1; k++)
            {
                output[k] = (float[,])p[k].Clone();
                for (int j = 1; j < p[k].GetLength(1) - 1; j++)
                    for (int i = 1; i < p[k].GetLength(0) - 1; i++)
                    {
                        if (p[k][i, j] > 0)
                            continue;
                        else
                        {
                            top = p[k][i, j - 1];
                            left = p[k][i - 1, j];
                            right = p[k][i + 1, j];
                            bottom = p[k][i, j + 1];
                            back = p[k - 1][i, j];
                            front = p[k + 1][i, j];
                            float sum = top + left + right + bottom + back + front;
                            if (sum > 0)
                                output[k][i, j] = 1;
                        }
                    }
            }
            return output;
        }

        public float[][,] RetrieveFinalizedDoseSpace()
        {                        
            DoseSpace = PrepareWeighted_DS(SliceWeights, folderpath);
            FindMaxDoseVoxel();
            FindReferenceContributions(ReferencePoint, folderpath);            
            //DoseSpace = PrepareWeighted_DS_GPU(SliceWeights, folderpath, DDS);
           // WriteFloatArray2BMP(DoseSpace[DoseSpace.GetLength(0) / 2], "PS_778_ds_midplane.bmp");
            //Debug.WriteLine("Before Normalization: ");
            //FindTotalCoverage(0.5, DoseSpace, DDS, SliceWeights);
            //DoseSpace = Matrix.Normalize(DoseSpace);
            //Debug.WriteLine("After Normalization: ");
            //FindTotalCoverage(0.5, DoseSpace, DDS, SliceWeights);
            //MakeMidplaneXPic(DoseSpace, "X_view_Dosespace.bmp");
            //MakeMidplaneYpic(DoseSpace, "Y_view_Dosespace.bmp");
            //MakeMidplaneZpic(DoseSpace, "Z_view_Dosespace.bmp");

            return DoseSpace;
        }

        

        

        private void MakeMidplaneXPic(float[][,] DoseSpace, string p)
        {
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            float[,] X = new float[y, z];
            float[,] Xdds = new float[y,z];
            for (int k = 0; k < z; k++)
                for (int j = 0; j < y; j++)
                {
                    X[j,k] = DoseSpace[k][x / 2, j];
                    Xdds[j,k] = dds_nondilated[k][x / 2, j];
                }
            WriteFloatArray2BMP(X, p);
            WriteFloatArray2BMP(Xdds, "X_midplane_DS.bmp");
            MakeComparisonPic(X, Xdds, "X_midplane_comparison.bmp");
        }

        private void MakeMidplaneYpic(float[][,] DoseSpace, string p)
        {
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            float[,] Y = new float[DoseSpace[0].GetLength(0), DoseSpace.GetLength(0)];
            float[,] Ydds = new float[DoseSpace[0].GetLength(0), DoseSpace.GetLength(0)];
            for (int k = 0; k < z; k++)
                for (int i = 0; i < x; i++)
                {
                    Y[i,k] = DoseSpace[k][i, y / 2];
                    Ydds[i,k] = dds_nondilated[k][i, y / 2];
                }
            WriteFloatArray2BMP(Y, p);
            WriteFloatArray2BMP(Ydds, "Y_midplane_DS.bmp");
            MakeComparisonPic(Y, Ydds, "Y_midplane_comparison.bmp");
        }

        private void MakeMidplaneZpic(float[][,] DoseSpace, string p)
        {
            float[,] Z = DoseSpace[DoseSpace.GetLength(0) / 2];
            float[,] Zdds = dds_nondilated[DDS.GetLength(0) / 2];
            WriteFloatArray2BMP(Z, p);
            WriteFloatArray2BMP(Zdds, "Z_midplane_DS.bmp");
            MakeComparisonPic(Z, Zdds, "Z_midplane_comparison.bmp");
        }



        private float[][,] PrepareWeighted_DS_GPU(double[] weights, string subfolder, float[][,] DDS)
        {
            float[][,] weighted_slicedoses = new float[DoseSpace.GetLength(0)][,];
            bool restrictweight = false;
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < weighted_slicedoses.GetLength(0); k++)
                weighted_slicedoses[k] = Matrix.Zeroes(x, y);
            PathSet.DCT = DCT;


            //wSD = GPU.PrepareDoseSpace(wSD, SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
            float[] wSD = GPU.WeightOriginalDS(SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
            WriteFloatArray2BMP(GrabSlice(wSD, DoseSpace.GetLength(0)/2, x, y), "wSD.bmp");
            
            float[][,] GPUsd = GPU.BackTo3D(wSD, x, y, z);
            return GPUsd;


        }
        /// <summary>
        /// Initializes the DS matrix using a preliminary weight set and the original slicedoses.
        /// </summary>
        /// <param name="weights"></param>
        /// <param name="subfolder"></param>
        private float[][,] PrepareWeighted_DS(double[] weights, string subfolder)
        {
            float[][,] weighted_slicedoses = new float[DoseSpace.GetLength(0)][,];
            
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < weighted_slicedoses.GetLength(0); k++)
                weighted_slicedoses[k] = Matrix.Zeroes(x,y);
            //PathSet.DCT = DCT;
            #region GPU
            if (GPU.GPUenabled)
            {
                float[][,] GPUsd = new float[DoseSpace.GetLength(0)][,];
                Stopwatch gputime = new Stopwatch();
                gputime.Start();
                float[] wSD = new float[x * y * z];
                for (int i = 0; i < wSD.GetLength(0); i++)
                    wSD[i] = 0.0f;
                
                //wSD = GPU.PrepareDoseSpace(wSD, SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
                wSD = GPU.WeightOriginalDS(SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
                gputime.Stop();
                Debug.WriteLine("GPU PrepareWeighted_DS time: " + gputime.Elapsed);
                GPUsd = GPU.BackTo3D(wSD,x,y,z);
                return GPUsd;
             
            }
#endregion
            else
            {
                //Stopwatch cputimer = new Stopwatch();
                //cputimer.Start();
                Debug.Write("Preparing weighted DS.");
                ArrayList OverdosePoints = new ArrayList();
                for (int s = 0; s < NumSlices; s++)
                {
                    int whichz = SlicePositions[s] - (DCT / 2); // <- NEED TO CHANGE APPROPRIATELY
                    int StartAt = 0;
                    int EndAt = DCT;
                    if (whichz < 0)
                    {
                        StartAt += (-1) * (whichz);
                        whichz = 0;
                    }
                    if ((whichz + DCT) >= weighted_slicedoses.GetLength(0))
                        EndAt = (weighted_slicedoses.GetLength(0) - whichz);    

                    //int whichz = SlicePositions[s] - (DCT / 2);
                    float[] slicedose = LoadSliceDose(s, subfolder);
                    
                    //Parallel.For(0, DCT, (k) =>
                    //    {
                    //        for (int j = 0; j < y; j++)
                    //            for (int i = 0; i < x; i++)
                    //            {
                    //                weighted_slicedoses[whichz + k][i, j] += (slicedose[(k * x * y) + (j * x) + i] * (float)weights[s]);
                    //            }
                    //    });
                    for (int k = StartAt; k < EndAt; k++)
                        {
                            Parallel.For(0, y, (j) =>
                            {
                                for (int i = 0; i < x; i++)
                                {
                                    float currentdose = Convert.ToSingle(weighted_slicedoses[whichz + k][i, j]);
                                    float originaldose = slicedose[(k * x * y) + (j * x) + i];
                                    float sliceweight = (float)weights[s];
                                    float updated_dose = Convert.ToSingle(currentdose) + (originaldose * sliceweight);
                                    //if (currentdose > 0.0 && updated_dose > 1.0)
                                    //{
                                    //    //Add coordinate index, value, and sliceweight to the overdose arraylist for post-processing
                                    //    float global_index = ((whichz + k) * x * y) + (j * x) + i;
                                    //    float slice_index = (k * x * y) + (j * x) + i;
                                    //    float[] overdose = new float[7] { global_index, updated_dose, slice_index, originaldose, s, sliceweight, currentdose };
                                    //    //OverdosePoints.Add(overdose);
                                    //    //OverdosePoints.TrimToSize();
                                    //}
                                    weighted_slicedoses[whichz + k][i, j] = updated_dose;
                                }
                            });                        
                                
                            OverdosePoints.TrimToSize();
                        }
                    //
                    
                    
                    int progress = Convert.ToInt16((double)(s * 100.0) / (double)NumSlices);
                    //PS_3_SliceWeightOpt_worker.ReportProgress(progress);
                    Debug.Write(".");
                }
                Debug.Write("done.");
                Debug.WriteLine("Post-processing weights to ensure no overdosing...");
                
                return weighted_slicedoses;
            }
        }

        private PointF[] GetAllShotsArray()
        {
            //Get total shotcount
            int shotcount = 0;
            foreach (RasterPath rp in RasterPaths)
            {
                shotcount += rp.shots.GetLength(0);
            }
            
            //Populate array with all shots
            PointF[] AllShots = new PointF[shotcount];
            int counter = 0;
            foreach (RasterPath rp in RasterPaths)
            {
                for (int i = 0; i < rp.shots.GetLength(0); i++)
                {
                    AllShots[counter] = rp.shots[i];
                    counter++;
                }
            }
            return AllShots;
        }

        public void AddSliceDoseToDoseSpace(float[][,] DS, int which_slice)
        {
            PointF[] shots = ((RasterPath)RasterPaths[which_slice]).shots;
            double[] weights = ((RasterPath)RasterPaths[which_slice]).ShotWeights;
            int ds_x = 0;
            int ds_y = 0;
            int ds_z = 0;
            int index = 0;
            //Loop through dosekernel and fill in each location
            for (int k = 0; k < N; k++)
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                    {
                        for (int shot = 0; shot < shots.GetLength(0); shot++)
                        {
                            ds_x = (int)(shots[shot].X - ((N-1)/2) + i);
                            ds_y = (int)(shots[shot].Y - ((N-1)/2) + j);
                            ds_z = (int)(SlicePositions[which_slice] - ((N-1)/2) + k);
                            index = ds_z * (StructureSet.BIG_dim[0] * StructureSet.BIG_dim[1]) + (ds_y * StructureSet.BIG_dim[0]) + ds_x;
                            if (ds_x < 0 || ds_y < 0 || ds_z < 0)
                                continue;
                            else if (ds_x >= StructureSet.BIG_dim[0] || ds_y >= StructureSet.BIG_dim[1] || ds_z >= StructureSet.BIG_dim[2])
                                continue;
                            else if (index < 0 || index > (StructureSet.BIG_dim[0] * StructureSet.BIG_dim[1] * StructureSet.BIG_dim[2]))
                                continue;
                            else
                            {
                                DS[ds_z][ds_x, ds_y] += (float)(DK.dose[k][(j * N) + i] * weights[shot]);
                            }
                        }
                    }
        }

        public void Calculate3DDoseSpace(DoseKernel dk)
        {
            //PointF[] AllShots = GetAllShotsArray();
            float[][,] DoseSpace = new float[StructureSet.BIG_dim[2]][,];
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
            {
                //Create a new 2D slice.
                float[,] layer = Matrix.Zeroes(StructureSet.BIG_dim[0], StructureSet.BIG_dim[1]);
                DoseSpace[k] = (float[,])layer.Clone();                
            }
            for (int i = 0; i < SlicePositions.GetLength(0); i++)
                AddSliceDoseToDoseSpace(DoseSpace, i);            
        }

        public double[] FindMaxPointsInEachSlice(float[][,] ds)
        {
            double[] MaxValues = new double[NumSlices];
            int StartAt = 0; int EndAt = 0;
            for (int s = 0; s < NumSlices; s++)
            {
                StartAt = SlicePositions[s] - (DCT / 2);
                EndAt = StartAt + DCT;
                double max = 0;
                if (StartAt < 0)
                    StartAt = 0;
                if (EndAt > ds.GetLength(0))
                    EndAt = ds.GetLength(0);
                for (int k = StartAt; k < EndAt; k++)
                        for (int j = 0; j < ds[0].GetLength(1); j++)
                            for (int i = 0; i < ds[0].GetLength(0); i++)
                            {
                                float value = ds[k][i, j];
                                if (value > max)
                                    max = value;
                            }
                MaxValues[s] = max;
            }
            return MaxValues;
        }

        public double[] GetSliceWeightRestrictions(float[][,] ds, double[] weight, string subfolder)
        {
            double[] restrictweight = new double[weight.GetLength(0)];
            double[] maxvalues = FindMaxPointsInEachSlice(ds);
            double max = maxvalues.Max();

            while (max > 1.0)
            {
                
                for (int i = 0; i < NumSlices; i++)
                {
                    if (maxvalues[i] == max)
                        restrictweight[i] = weight[i] / max;
                    else
                        restrictweight[i] = weight[i];
                }

                ds = PrepareWeighted_DS(restrictweight, subfolder);

                maxvalues = FindMaxPointsInEachSlice(ds);
                max = maxvalues.Max();

            }
            return restrictweight;
        }

        public double[] SliceWeightPostProcess(float[][,] ds, string subfolder, double[] iteration_weight)
        {
            float max = 0; int maxindex = 0;
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < ds.GetLength(0); k++)
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        float value = ds[k][i, j];
                        if (value > max)
                        {
                            max = value;
                            maxindex = (k * x * y) + (j * x) + i;
                        }
                    }
            //First, make a list of the index ranges of each slice.
            int[] SliceIndexStartPositions = new int[NumSlices];
            //SliceIndexStartPositions[0] = 0;
            for (int s = 0; s < NumSlices; s++)
                SliceIndexStartPositions[s] = (SlicePositions[s] - (DCT / 2)) * x * y;


            int[] RelativeIndices = new int[NumSlices];
            for (int s = 0; s < NumSlices; s++)
            {
                RelativeIndices[s] = (maxindex - SliceIndexStartPositions[s]);
                if (RelativeIndices[s] < 0)
                    RelativeIndices[s] = 0;
                else if (RelativeIndices[s] > (DCT * x * y))
                    RelativeIndices[s] = 0;
            }

            
            //Third, find the weight contributions of each slice to maxvalue
            double[] ContributingWeights = new double[NumSlices];
            double[] Values = new double[NumSlices];
            for (int s = 0; s < NumSlices; s++)
            {
                if (RelativeIndices[s] > 0)
                {
                    if (iteration_weight[0] > 0)
                        ContributingWeights[s] = iteration_weight[s];
                    else
                        ContributingWeights[s] = SliceWeights[s];
                    float[] slicedose = LoadSliceDose(s, subfolder);
                    Values[s] = slicedose[RelativeIndices[s]];
                }
                else
                {
                    ContributingWeights[s] = 0;
                    Values[s] = 0;
                }

            }      

            //ContributingWeights = Matrix.Normalize(ContributingWeights);

            /*To find a multiplier without disturbing the relative proportions, 
             *  k(uv1 + uv2 + uv3) = 1.0
             *  k = 1.0 / (uv1 + uv2 + uv3)
             *  k = "multiplier"
             */
            double sum = 0;
            for (int i = 0; i < Values.GetLength(0); i++)
                sum += (Values[i] * ContributingWeights[i]);
            double multiplier = 1.0 / sum;
            double[] weights = new double[NumSlices];

            //double[] priorities = null;
            //int SliceThatNeedsWeightBadly = 0;
            //if (iteration_weight[0] != 0) //if the iteration_weight array isn't a dummy variable
            //{
            //    priorities = new double[NumSlices];
            //    for (int i = 0; i < NumSlices; i++)
            //    {
            //        //If priorities is POSITIVE, that means the weight needs to go up for adequate coverage.
            //        //If NEGATIVE, thats fine, weight can go down.
            //        if (ContributingWeights[i] != 0)
            //            priorities[i] = iteration_weight[i] - ContributingWeights[i];
            //        else
            //            priorities[i] = 0;
            //    }
            //    SliceThatNeedsWeightBadly = Array.IndexOf(priorities, priorities.Max());
            //}

            //Change only the contributing weights. In effect, this "restricts them" to a certain weight.
            for (int s = 0; s < NumSlices; s++)
            {
                //if (ContributingWeights[s] != 0) //<- Check to see if this slice is important
                //{
                //    //If its the important slice, check if the adjacent slices are involved.
                //    if (priorities != null && priorities[s] != 0 && s == SliceThatNeedsWeightBadly)
                //    {
                //        double threshold = 0.3 * priorities[s];
                //        //Only play with weights if adjacent slice priorities are < 30% of the important slice
                //        if ((s+1) < NumSlices && priorities[s + 1] <= threshold && priorities[s+1] != 0)
                //        {
                //            weights[s] = ContributingWeights[s] * multiplier * 1.1; // Add 10% more weight
                //            ContributingWeights[s + 1] = ContributingWeights[s + 1] * 0.9;
                //        }
                //        if ((s - 1) >= 0 && priorities[s - 1] <= threshold && priorities[s - 1] != 0)
                //        {
                //            weights[s] = ContributingWeights[s] * multiplier * 1.1; // Add 10% more weight
                //            weights[s - 1] = ContributingWeights[s - 1] * multiplier * 0.9;
                //        }
                //    }
                //    else
                //    weights[s] = ContributingWeights[s] * multiplier;
                //    Debug.WriteLine("Restricting slice " + s + " weight: " + SliceWeights[s] + " --> " + weights[s]);
                //}
                //else
                //    weights[s] = SliceWeights[s];

                weights[s] = SliceWeights[s] * multiplier;
            }
            return weights;
        }



        public void FindReferenceContributions(int[] refpoint, string subfolder)
        {
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            int refindex = (refpoint[2]*x*y)+(refpoint[1]*x)+refpoint[0];
            ReferenceWeights = new double[NumSlices];
            ReferenceValues = new double[NumSlices];
            for (int s = 0; s < NumSlices; s++)
            {
                ReferenceWeights[s] = 0.0;
                ReferenceValues[s] = 0.0;
                int whichz = SlicePositions[s] - (DCT / 2);
                float[] slicedose = LoadSliceDose(s, subfolder);                
                Parallel.For(0, DCT, (k) =>
                {                    
                        for (int j = 0; j < y; j++)                        
                            for (int i = 0; i < x; i++)
                            {
                                int tempindex = ((whichz+k) * x * y) + (j * x) + i;
                                if (tempindex == refindex)
                                {
                                    ReferenceWeights[s] = SliceWeights[s];
                                    ReferenceValues[s] = slicedose[(k * x * y) + (j * x) + i];
                                }
                                else
                                    continue;                                
                            }
                });
            }
            Debug.WriteLine(WriteArrayAsList("Reference Point: ", ReferencePoint));
            Debug.WriteLine("Slice contributions to reference point: ");
            for (int s = 0; s < NumSlices; s++)
            {
                Debug.WriteLine("Slice " + s + ": " + ReferenceValues[s] + " x (Weight: " + ReferenceWeights[s] + ") = " + Math.Round(ReferenceWeights[s] * ReferenceValues[s], 3));
            }
        }

        public void FindDoseContributionsToReferencePoint(double[] normweights, double doserate, double max)
        {
            double[] DoseContributions = new double[NumSlices];
            double sum = 0;
            for (int s = 0; s < NumSlices; s++)
            {
                DoseContributions[s] = ReferenceValues[s] * normweights[s];
                sum += DoseContributions[s];
                Debug.WriteLine("Slice " + s + ": " + Math.Round(ReferenceValues[s],3) + " x (Weight: " + Math.Round(normweights[s],3) + ") = " + Math.Round(DoseContributions[s], 3));                
            }
            Debug.WriteLine("Doserate: " + doserate);
            Debug.WriteLine("Max Dose: " + max);
            Debug.WriteLine(WriteArrayAsList("Contributions: ", DoseContributions));
            Debug.WriteLine("Total = " + Math.Round(sum, 3));
        }

        public double[] FindSliceWeightsForMaxDose(double max)
        {
            double sum = 0;
            double weight_multiplier = 0;
            for (int s = 0; s < NumSlices; s++)
                sum += ReferenceValues[s] * ReferenceWeights[s];
            weight_multiplier = (max / sum);
            double[] NormalizedWeights = new double[NumSlices];
            for (int s = 0; s < NumSlices; s++)
                NormalizedWeights[s] = SliceWeights[s] * weight_multiplier;
            return NormalizedWeights;
        }

        private float[] LoadSliceDose(int which_slice, string subfolder)
        {
            string filename = String.Concat("slice_", which_slice);
            string path = System.IO.Path.Combine(subfolder, filename);
            float[] d = ReadSliceDoseFromFile(path);
            Matrix.Normalize(ref d);
            return d;
        }

        private int FindStartZ(int s)
        {
            int start = SlicePositions[s] - (DCT/2);
            if (start < 0)
                start = 0;
            return start;
        }
        private int FindEndZ(int startz)
        {
            int end = startz + DCT;
            if (end >= DoseSpace.GetLength(0))
                end = DoseSpace.GetLength(0)-1;
            return end;
        }




        /// <summary>
        /// Finds fraction of tumor that is covered by the desired dose. Called from the OptimizeSliceWeights()
        /// function of Step 2.
        /// </summary>
        /// <param name="iso"></param>
        /// <param name="Tumor"></param>
        /// <returns></returns>
        public double[] FindTotalCoverage(double iso, float[][,] ds, float[][,] Tumor, double[] recent_weight)
        {
            if (TotalCoverage != null)
                OldCoverage = (double[])TotalCoverage.Clone();
            float[][,] normds = Matrix.Normalize(ds);
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            float dose; float dds_value; double Uncovered = 0; 
            for (int k = 0; k < ds.GetLength(0); k++)                            
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        dose = normds[k][i, j];
                        dds_value = Tumor[k][i, j];                        
                        if (dose >= iso)
                        {
                            RxVolume++;
                            if (dds_value > ToleranceDose*2)
                            {
                                TumorVol++;
                                LesionRx++;
                            }
                            else if (dds_value < ToleranceDose * 2)
                            {
                                double overage = 2;
                            }
                        }
                        else if (dose < iso)
                        {
                            if (dds_value > ToleranceDose*2)
                            {
                                TumorVol++;
                                Uncovered++;
                            }
                        }
                    }            

            // (isovolume / tumor);
            double LesionRxoverTumorVol = LesionRx / TumorVol;
            double LesionRxoverRxVol = LesionRx / RxVolume;
            double Underdosed = (Uncovered / TumorVol) * 100;
            double Overdosed = TumorVol / RxVolume;
            string weights = "SliceWeights: " + recent_weight[0];
            for (int i = 1; i < recent_weight.GetLength(0); i++)
                weights += ", " + recent_weight[i];                      
            return new double[5] { TumorVol, RxVolume, LesionRx, Uncovered, Overdosed };
        }

        private double[] FindSliceDoseCoverage(float[][,] ds, int which_slice, double iso, float[][,] DDS)
        {
            //float[][,] sd = Matrix.Normalize(ds);
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            float dose; float dds_value; double Uncovered = 0;
            int x = DDS[0].GetLength(0);int y = DDS[0].GetLength(1); int z = DDS.GetLength(0);
            for (int k = 0; k < DDS.GetLength(0); k++)
                for (int j = 0; j < DDS[0].GetLength(1); j++)
                    for (int i = 0; i < DDS[0].GetLength(0); i++)
                    {
                        dose = ds[k][i,j];
                        dds_value = DDS[k][i, j];
                        if (dose < iso)
                        {
                            if (dds_value > ToleranceDose)
                            {
                                TumorVol++;
                                Uncovered++;
                            }
                        }
                        else if (dose >= iso)
                        {
                            RxVolume++;
                            if (dds_value > ToleranceDose)
                            {
                                TumorVol++;
                                LesionRx++;
                            }
                        }
                       
                    }           
            double LesionRxoverTumorVol = LesionRx / TumorVol;
            double LesionRxoverRxVol = LesionRx / RxVolume;
            double Underdosed = (Uncovered / TumorVol);
            double Overdosed = RxVolume / TumorVol;
            //Debug.WriteLine("SLICE " + which_slice);
            //Debug.WriteLine("==============================");
            //Debug.WriteLine("LesionRx/TumorVol: " + LesionRxoverTumorVol);
            //Debug.WriteLine("LesionRx/RxVol: " + LesionRxoverRxVol);
            //Debug.WriteLine("Percentage underdosed: " + Underdosed);
            double sum_sd = Matrix.SumAll(ds);
            double sum_dds = Matrix.SumAll(DDS);            
            double ratio = sum_dds / sum_sd;
            //return new double[5] { ratio, RxVolume / TumorVol, LesionRxoverTumorVol, LesionRxoverRxVol, Uncovered };
            return new double[6] { ratio, RxVolume, TumorVol, LesionRx, Underdosed, Overdosed };
        }


        private void ClearDosespace()
        {
            //DoseSpace = new float[DoseSpace.GetLength(0)][,];
            for (int i = 0; i < DoseSpace.GetLength(0); i++)
                DoseSpace[i] = Matrix.Zeroes(DoseSpace[0].GetLength(0), DoseSpace[0].GetLength(1));            

        }
        
        /// <summary>
        /// Called from EvaluateIterationForThisSlice(). Gets the DDS slab from the SS tumor at the specified location (slicecenter_z)
        /// </summary>
        /// <param name="SS"></param>
        /// <param name="slicecenter_z"></param>
        /// <param name="doseslab"></param>
        /// <returns></returns>
        //public double GetSumOfMultiplied_DDS_Subset(StructureSet SS, int slicecenter_z, float[][,] doseslab)
        //{
        //    int dstart = (N - DCT) / 2; int dend = dstart + DCT;

        //    float[][,] DDSslab = GrabSlab(SS.fj_Tumor, DCT, slicecenter_z);

        //    //Multiply the elements and return.
        //    return Matrix.SumAll(Matrix.MultiplyElements(doseslab, DDSslab));
        //}
        /// <summary>
        /// Called from EvaluateIterationForThisSlice(). Gets the DS slab from the SS tumor at the specified location (slicecenter_z)
        /// </summary>
        /// <param name="dosespace"></param>
        /// <param name="z"></param>
        /// <param name="doseslab"></param>
        /// <returns></returns>
        //public double GetSumOfMultiplied_DS_Subset(float[][,] dosespace, int z, float[][,] doseslab)
        //{
        //    //Get the dose slab
        //    int dstart = (N - DCT) / 2; int dend = dstart + DCT;

        //    //Convert the 1D float ds matrix into the jagged DS matrix, and then get the slab.            
        //    float[][,] DSslab = GrabSlab(dosespace, DCT, z);

        //    return Matrix.SumAll(Matrix.MultiplyElements(doseslab, DSslab));
        //}
        private double FindError(double[] weight, double[] w)
        {
            double diff = 0;
            for (int i = 0; i < w.GetLength(0); i++)
            {
                diff += Math.Pow(Math.Abs((weight[i] - w[i])), 2);
            }
            return Math.Sqrt(diff);
        }
        /// <summary>
        /// Called by many Step 2 methods above. Takes a float[][,] matrix, a centered z-index, and desired thickness and returns
        /// an appropriate section, depending on boundaries.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="thickness"></param>
        /// <param name="z_center"></param>
        /// <returns></returns>
       public static float[][,] GrabSlab(float[][,] f, int thickness, int z_center)
        {
            //Find the start and end points based on the thickness
            int temp_start = z_center - (thickness / 2);
            int temp_end = z_center + (thickness / 2);

            //Check if the endpoints exist, change them if necessary 
            if (temp_start < 0)
                temp_start = 0;
            if (temp_end >= f.GetLength(0))
                temp_end = f.GetLength(0) - 1;

            //Initialize new matrix, fill, return
            float[][,] slab = new float[temp_end - temp_start][,];
            for (int i = 0; i < (temp_end - temp_start); i++)
            {
                slab[i] = f[i + temp_start];
            }
            return slab;
        }

        public static float[][,] GrabSlab(float[][,] f, int temp_start, int temp_end, bool startend)
        {
            //Find the start and end points based on the thickness
            
            //Initialize new matrix, fill, return
            float[][,] slab = new float[temp_end - temp_start][,];
            for (int i = 0; i < (temp_end - temp_start); i++)
            {
                slab[i] = f[i + temp_start];
            }
            return slab;
        }


        private void WriteFloatArray2BMP(float[,] temp, string p)
        {
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            float max = Matrix.FindMax(temp);
            int color = 0;
            Bitmap b = new Bitmap(temp.GetLength(1), temp.GetLength(0));
            for (int j = 0; j < temp.GetLength(0); j++)
                for (int i = 0; i < temp.GetLength(1); i++)
                {
                    color = (int)((temp[j, i]/max) * 255);
                    b.SetPixel(i, j, Color.FromArgb(color, color, color));
                }
            b.Save(path);
        }

        private void MakeComparisonPic(float[,] temp, float[,] dds, string p)
        {
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            float[,] temp2 = (float[,])Matrix.Normalize(temp).Clone();
            int color = 0;
            Bitmap b = new Bitmap(temp.GetLength(0), temp.GetLength(1));
            for (int j = 0; j < temp.GetLength(1); j++)
                for (int i = 0; i < temp.GetLength(0); i++)
                {
                    float dose = temp[i, j]; float pixel = dds[i, j];
                    if (pixel > ToleranceDose)
                    {
                        if (dose >= RxDose)
                            b.SetPixel(i, j, Color.Green);
                        else if (dose < RxDose)
                            b.SetPixel(i, j, Color.Blue);
                    }
                    else if (pixel <= ToleranceDose)
                    {
                        if (dose >= RxDose)
                            b.SetPixel(i, j, Color.Red);
                    }
                }
            b.Save(path);
        }

        public static void RestrictWeight(double[] d)
        {

        }

        
        #endregion

        #region Helper methods that Read/Write dose to files

        public void FindMaxDose()
        {
            max = Matrix.FindMax(DoseSpace);
            ReferencePoint = new int[3];
            ReferencePoint = FindMaxDoseVoxel();
        }

        public int[] FindMaxDoseVoxel()
        {
            int[] point = new int[3]; float max = 0;
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
                for (int j = 0; j < DoseSpace[0].GetLength(1); j++)
                    for (int i = 0; i < DoseSpace[0].GetLength(0); i++)
                        {
                            float dose = DoseSpace[k][i,j];
                            if (dose > max)
                            {
                                max = dose;
                                point = new int[3]{i, j, k};
                            }
                        }
            return point;
        }



        /// <summary>
        /// Takes in a slice 1D float matrix, and adds it to the global dosespace array for final calculation. Called by 
        /// AssembleFinalDoseMatrix(). If dosespace is null, it will create the slice, else it will add
        /// it to the existing one.
        /// </summary>
        /// <param name="slicedose"></param>
        /// <param name="which_z_slice"></param>
        public float[][,] WriteSliceDoseToDoseSpace(float[] slicedose, float[][,] output_ds, int which_slice)
        {            
            int TranslateZBy = SlicePositions[which_slice] - (DCT/2); // <- NEED TO CHANGE APPROPRIATELY
            int StartAt = 0;
            int EndAt = DCT;
            if (TranslateZBy < 0)
            {
                StartAt += (-1) * (TranslateZBy);
                TranslateZBy = 0;
            }
            if ((TranslateZBy + DCT) >= output_ds.GetLength(0))
                EndAt = (output_ds.GetLength(0) - TranslateZBy);            
                
            //int startz = SlicePositions[which_slice] - DCT/2;

            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.

            if (GPU.GPUenabled)
                return output_ds;
               //return GPU.AddToDoseSpace(slicedose, output_ds, TranslateZBy, 1.0f);
            //else
                
            
            //Parallel.For(0, DCT, (z) =>
                //Parallel.For(StartAt, EndAt, (z) =>
            for (int z = StartAt; z < EndAt; z++)
                    {
                        int current_z = TranslateZBy + z;                        
                        if (output_ds[current_z] == null)
                            output_ds[current_z] = Matrix.Zeroes(output_ds[0].GetLength(0), output_ds[0].GetLength(1));                        
                        output_ds[current_z] = Matrix.Add(output_ds[current_z], GrabSlice(slicedose, z, volume[0].GetLength(0), volume[0].GetLength(1)));
                        int dssum = Matrix.SumAll(output_ds[current_z]);
                        //WriteFloatArray2BMP(output_ds[current_z], "ds_z_" + z + "_" + dssum + ".bmp");
                        //WriteFloatArray2BMP(GrabSlice(slicedose, z, output_ds[0].GetLength(0), output_ds[0].GetLength(1)), "sd_z_" + z + ".bmp");
                    }//);
                return output_ds;
            

            
        }
        /// <summary>
        /// Given a path, loads in the associated slicedose 1D float matrix. Called by AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="loadpath"></param>
        /// <returns></returns>
        public static float[] ReadSliceDoseFromFile(string loadpath)
        {
            float[] d;
            using (FileStream fs = new FileStream(loadpath, FileMode.Open, FileAccess.Read))
            using (StreamReader br = new StreamReader(fs))
            {
                int x = Convert.ToInt16(br.ReadLine()); int y = Convert.ToInt16(br.ReadLine()); int z = Convert.ToInt16(br.ReadLine());
                d = new float[x * y * z];
                for (int i = 0; i < d.GetLength(0); i++)
                    d[i] = Convert.ToSingle(br.ReadLine());
            }
            float sum = d.Sum();
            return d;
        }
        /// <summary>
        /// End of STEP 1. The actual method that is called by the UI, which takes in a folder directory path, and 
        /// loads in each slicedose file and writes it to the global dosespace file.
        /// </summary>
        /// <param name="folderpath"></param>
        public void AssembleDoseSpaceFromFiles()
        {
            string subfolder = ActiveDirectory;
            int progress = 0;
            int x = volume[0].GetLength(0); int y = volume[0].GetLength(1); int z = volume.GetLength(0);
            /* First, a for-loop to create the string name for each slice, which then loads each
             * slicedose file and in turn adds it to the result matrix.*/
            DoseSpace = new float[z][,];
            for (int i = 0; i < z; i++)
                DoseSpace[i] = Matrix.Zeroes(x, y);
            for (int s = 0; s < NumSlices; s++)
            {                
                float[] slicedose = LoadSliceDose(s, subfolder);
                DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);
                progress = (int)((double)s * 50.0 / NumSlices);
                PS_2_CalcDose_worker.ReportProgress(progress);
            }           

            //Write dosespace to file.
            WriteDoseSpaceToFile("OriginalDS.txt");
            PS_2_CalcDose_worker.ReportProgress(100);
        }
        private void NormalizeDosespace()
        {
            float max = 0;
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
                if (Matrix.FindMax(DoseSpace[k]) > max)
                    max = Matrix.FindMax(DoseSpace[k]);
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
                for (int j = 0; j < DoseSpace[0].GetLength(1); j++)
                    for (int i = 0; i < DoseSpace[0].GetLength(0); i++)
                        DoseSpace[k][i, j] = DoseSpace[k][i, j] / max;

        }
        /// <summary>
        /// Called from AssembleDoseSpaceFromFiles() to write the output to a file. Note the header that 
        /// is written first containing the list of dimension order.
        /// </summary>
        public void WriteDoseSpaceToFile(string filename)
        {
            string subfolder = ActiveDirectory;            
            string path = System.IO.Path.Combine(subfolder, filename);
            using (FileStream fs = new FileStream(path, FileMode.Create, FileAccess.Write))
            using (StreamWriter bw = new StreamWriter(fs))
            {
                //Write dimensions and plan info
                bw.WriteLine("Dimensions (x,y,z):");
                bw.WriteLine(DoseSpace[0].GetLength(0));
                bw.WriteLine(DoseSpace[0].GetLength(1));
                bw.WriteLine(DoseSpace.GetLength(0));
                bw.WriteLine("Number of Slices:");
                bw.WriteLine(NumSlices);
                bw.WriteLine("Dose data:");
                for (int k = 0; k < DoseSpace.GetLength(0); k++)
                    for (int j = 0; j < DoseSpace[0].GetLength(1); j++)
                        for (int i = 0; i < DoseSpace[0].GetLength(0); i++)
                            bw.WriteLine(DoseSpace[k][i, j]);
            }
            //System.Windows.MessageBox.Show(string.Concat("Dosespace written to: ", path));
        }
        /// <summary>
        /// Counterpart to WriteDoseSpaceToFile(). Loads in a dosespace and sets it to the main variable
        /// when given a path.
        /// </summary>
        /// <param name="path"></param>
        public static float[] ReadDoseSpaceFromFile(string filename)
        {
            float[] ds;
            string subfolder = ActiveDirectory;
            string path = System.IO.Path.Combine(subfolder, filename);
            using (FileStream fs = new FileStream(path, FileMode.Open, FileAccess.Read))
            using (StreamReader br = new StreamReader(fs))
            {
                br.ReadLine();
                int x = Convert.ToInt16(br.ReadLine()); int y = Convert.ToInt16(br.ReadLine()); int z = Convert.ToInt16(br.ReadLine());
                br.ReadLine();
                NumSlices = Convert.ToInt16(br.ReadLine());
                br.ReadLine();
                ds = new float[z*x*y];
                for (int k = 0; k < ds.GetLength(0); k++)
                {   
                    ds[k] = Convert.ToSingle(br.ReadLine());                    
                }
            }
            return ds;
        }

        #endregion
    }
}