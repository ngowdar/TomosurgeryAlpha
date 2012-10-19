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
        public string folderpath = ActiveDirectory; //The active directory where slice data will be written.
        public DoseKernel DK; //'DoseKernel', the object representing the dose matrix.
        public StructureSet SS; //'SS' represents the binary tumor volume from the DICOM.
        public static float[,] mask;
        public static int N;
        public static int shot_edgepadding;
        public static int line_edgepadding;
        public static int tumorflag = 10;
        public static int CSflag = 2;
        public static float ToleranceDose = 0.1f;
        public static float RxDose = 0.5f;
        public static float CSdose = 0.05f;
        public static int NumSlices; //Set in step 1
        public int TumorVolCount;
        public int CSvolCount;
        public ArrayList RasterPaths; //Collection of objects representing each slice
        public static string ActiveDirectory;
        public int DoseCalculationThickness; //How much of the dose kernel to use in calculation
        public static int DCT; //static version of dosecalcthickness
        public static int StepSize;
        public static int RasterWidth;        
        public int SliceThickness;        
        public int TolThickness;
        public float[][,] volume; //Presumably the structure set binary volume?
        public float[][,] DDS;
        public float[][,] DoseSpace; //The actual final dosespace matrix containing the real dose.
        public int[] SlicePositions; //The centered z-indexes of each slice location
        public double[] SliceWeights; //Weights, from 0 to 1, of each slice. Calculated from Step 2.
        public double[] SliceCoverage;
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
        public PathSet(float[][,] f, int sthick, int tolthick, DoseKernel dk, StructureSet ss)
        {
            
            X = f[0].GetLength(0); Y = f[0].GetLength(1); Z = f.GetLength(0);
            DK = dk; SS = ss;
            N = dk.dose.GetLength(0);
            boundaries = FindBoundaries(f);
            
            SliceThickness = sthick;
            DCT = SliceThickness * 2;
            TolThickness = tolthick;
            CalculateNumSlices();
            
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(f, SlicePositions[i], SliceThickness / 2));
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
            if (PathsetWorkerCompleted != null)
                PathsetWorkerCompleted.Invoke(null, e);
            //Calculate the slicedoses so far...
            CreateDoseMatrix(DK, folderpath); //<- This calls second background worker "PS_Two"
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
            Debug.WriteLine("Shot weighting optimization complete. Now commencing slice weighting...");
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
            CalculateNumSlices();
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(volume, SlicePositions[i], SliceThickness / 2));
                RasterPaths.Add(rp);
            }
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
                //Set an initial weight discrepancy if the slice is the first or last
                if (i == 0 || i == (SliceWeights.GetLength(0) - 1))
                    SliceWeights[i] = 1;
                else
                    SliceWeights[i] = 0.7;
                //temp_weights[i] = 1.0;                
            }
            

            //Start a while loop, and first create an updated DS matrix with the newest weights

            //UNDONE: Prepare DDS Tumor matrix 
            DDS = PrepareDDS(SS.fj_Tumor);
            //float[][,] OriginalDS = ReadDoseSpaceFromFile("OriginalDS.txt");


            //Reset the dosespace
            //ClearDosespace();
            DoseSpace = PrepareWeighted_DS(SliceWeights, folderpath, DDS);
            DoseSpace = Matrix.Normalize(DoseSpace);
            double[] initial_cov = FindTotalCoverage(0.5, DoseSpace, DDS, SliceWeights);
            //{ TumorVol, RxVolume, LesionRx, Uncovered };
            


            double old_error;
            while (Error >= .05 && index < 5)
            {
                
                temp_weights = ReOptimizeSliceWeights(DDS);
                old_error = Convert.ToDouble(Error);
                Error = FindError(SliceWeights, temp_weights);
                
                              
                //Re-prepare dosespace with latest iteration of sliceweight
                ClearDosespace();
                DoseSpace = PrepareWeighted_DS(temp_weights, folderpath, DDS); //TODO: This method is time-consuming, make this GPU?
                //DoseSpace = PrepareWeighted_DS_Dynamic(SliceWeights, OriginalDS);
                //DoseSpace = Matrix.Normalize(DoseSpace);            
                double[] result = FindTotalCoverage(RxDose, DoseSpace, DDS, temp_weights);
                double IterationCoverage = result[2] / result[0];
                //index++;

                if (index > 2 && IterationCoverage < coverage)
                {
                    string r = "CAUTION: Coverage is Decreasing!!! " + coverage + " --> " + IterationCoverage;
                    Debug.WriteLine(r);
                    break;
                }
                if (index > 2 && old_error < Error)
                {
                    string r = "WARNING: Error is increasing (" + old_error + " --> " + Error + "). Terminating...";
                    Debug.WriteLine(r);
                    //break;
                }
                else
                {
                    coverage = Convert.ToDouble(IterationCoverage);

                    Debug.WriteLine("Iteration: " + index + " Error: " + Error);
                    Debug.WriteLine("===============================");
                    Debug.WriteLine("ITERATION: " + index);
                    WriteArrayAsList("Weight", temp_weights);
                    Debug.WriteLine("Error: " + Math.Round(Error, 2));
                    Debug.WriteLine("Coverage: " + IterationCoverage);                    
                    Debug.WriteLine("===============================");


                    SliceWeights = (double[])temp_weights.Clone();
                    index++;
                }
                PS_3_SliceWeightOpt_worker.ReportProgress(index * 10);
            } // <- end of while loop
            
            PS_3_SliceWeightOpt_worker.ReportProgress(100);
        }

        private void WriteArrayAsList(string prefix, double[] f)
        {
            string output = prefix + ": [" + Math.Round(f[0], 2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += ", " + Math.Round(f[i], 2);
            output += "]";
            Debug.WriteLine(output);
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

        /// <summary>
        /// Helper method for ReviseWeighted_DS(), which takes both old and new sliceweights, removes the old contribuiton and adds the new one.
        /// </summary>
        /// <param name="slicedose"></param>
        /// <param name="old_weights"></param>
        /// <param name="recent_weights"></param>
        /// <param name="which_slice"></param>
        /// <returns></returns>
        private float[][,] ReviseSliceWeightContribution(float[] slicedose, double old_weight, double recent_weight, int which_slice)
        {
            int TranslateZBy = SlicePositions[which_slice] - (DCT / 2); // <- NEED TO CHANGE APPROPRIATELY
            //Check if this is correct ^

            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.

            //if (GPU.GPUenabled)
            //{
            //    float incremental_weight = (float)(recent_weight - old_weight);
                
            //    return GPU.PrepareDoseSpace(,slicedose, DoseSpace, TranslateZBy, 1.0f);
            //}
            //else
            //{
                Parallel.For(0, DCT, (z) =>
                {
                    float[,] incremental_slice = GrabSlice(slicedose, z, volume[0].GetLength(0), volume[0].GetLength(1));
                    int current_z = TranslateZBy + z;
                    for (int j = 0; j < incremental_slice.GetLength(1); j++)
                        for (int i = 0; i < incremental_slice.GetLength(0); i++)
                        {
                            float current = DoseSpace[current_z][i, j];
                            float subtract = (float)(incremental_slice[i, j] * old_weight);
                            float add = (float)(incremental_slice[i, j] * recent_weight);
                            current = current - subtract + add;
                            DoseSpace[current_z][i, j] = current;
                            Debug.Assert(DoseSpace[current_z][i, j] >= 0);
                        }
                });
                return DoseSpace;
           // }
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

        private double[] ReOptimizeSliceWeights(float[][,] dds)
        {
           double[] tweight = (double[])SliceWeights.Clone();
           int count = 0;

           for (int s = 0; s < NumSlices; s++)
           {
               float[][,] ds_slab = GrabSlab(DoseSpace,SliceThickness, SlicePositions[s]);
               float[][,] dds_slab = GrabSlab(dds, SliceThickness, SlicePositions[s]);
               double[] measurements = FindSliceDoseCoverage(ds_slab,s,0.5,dds_slab);
               
               //double DDS_slicesum = Matrix.SumAll(dds_slab);
               //double DS_slicesum = Matrix.SumAll(ds_slab);
               //int startz = FindStartZ(s); //Just takes starting z - DCT / 2
               //int endz = FindEndZ(startz);
               //double ratio = CompareSlicedoses(DoseSpace, dds, startz, startz + DCT);
               //double ratio = CompareSliceSums(DoseSpace, dds, startz, endz);
               //for (int z = 0; z < DCT; z++)
               //{
               //    /* The DDS matrix is element-multiplied by the newest dosespace
               //     * and the sum is added to DDS_slicesum. The DS matrix is also
               //     * multiplied by the latest iteration (squaring it?)*/
               //    float[,] DDS_slice = DDS[startz + z]; //Debug.WriteLine("Sum of DDS: " + Matrix.SumAll(DDS_slice));
               //    float[,] DS_slice = DoseSpace[startz + z]; //Debug.WriteLine("Sum of DS: " + Matrix.SumAll(DS_slice));
               //    //DDS_slicesum += Matrix.SumAll(Matrix.MultiplyElements(DDS_slice, OriginalDS[startz + z]));
               //    //DS_slicesum += Matrix.SumAll(Matrix.MultiplyElements(DS_slice, OriginalDS[startz + z]));
               //    DDS_slicesum += Matrix.SumAll(DDS_slice);
               //    DS_slicesum += Matrix.SumAll(DS_slice);
               //}
               //Debug.Assert(DDS_slicesum != 0);
               //Debug.Assert(DS_slicesum != 0);               
               double ratio = measurements[0];
               double RxVolvsTumor = measurements[1];
               double BothvsTumor = measurements[2];
               double BothvsRxVol = measurements[3];
               tweight[s] = SliceWeights[s] / RxVolvsTumor;
               //Debug.WriteLine("Ratio=" + ratio + "; Weight: " + SliceWeights[s] + "-->" + tweight[s]);
               //Debug.WriteLine("DDS / DS sum: " + DDS_slicesum + " / " + DS_slicesum);
               count++;
           }
           //tweight = Normalize(tweight);
           return tweight;
        }

        private double CompareSliceSums(float[][,] DoseSpace, float[][,] dds, int startz, int endz)
        {
            double ds_sum = 0;
            double dds_sum = 0;

            Parallel.For(0, endz - startz, (i) =>
                {
                    float[,] ds_slice = DoseSpace[i + startz];
                    float[,] dds_slice = dds[i + startz];
                    ds_sum += Matrix.SumAll(ds_slice);
                    dds_sum += Matrix.SumAll(dds_slice);
                });
            return (double)(dds_sum / ds_sum);
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

        /// <summary>
        /// "Card counting" comparison method that uses the same method from the
        /// shot weighting algorithm, CompareSlices(), but generalizes it to 3
        /// dimensions
        /// </summary>
        /// <param name="ds"></param>
        /// <param name="dds"></param>
        /// <param name="startz"></param>
        /// <param name="endz"></param>
        /// <returns></returns>
        private double CompareSlicedoses(float[][,] ds, float[][,] dds, int startz, int endz)
        {
            double total_tally = 0;           
            double ratio = 0;
            double LesionVol = 0;
            double TotalRx = 0;
            double LesionRx = 0;
            double underdosed = 0;
            double overdosed = 0;
            Parallel.For(0, DCT, (k) =>
                {
                    //total_tally += RasterPath.CompareSlices(ds[k + startz], dds[k + startz], true);
                    double[] d = RasterPath.CompareSlices(ds[k + startz], dds[k + startz], true);
                    LesionVol += d[0];
                    TotalRx += d[1];
                    LesionRx += d[2];
                    underdosed += d[3];
                    overdosed += d[4];
                });            

            //if (total_tally <= 0)
            //{
            //    ratio = ((total_tally * (-1)) / (double)(ds[0].GetLength(0) * ds[0].GetLength(1) * DCT));
            //}
            //else if (total_tally > 0)
            //    ratio = (double)(1 + (total_tally / (ds[0].GetLength(0) * ds[0].GetLength(1) * DCT)));
            ratio = (LesionVol + underdosed) / (TotalRx + overdosed);

            return ratio;
        }
        public float[][,] PrepareDDS(float[][,] dds_slice)
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
                        if (dds_slice[k][i, j] == 0)
                            pDDS[k][i, j] = ToleranceDose;
                        //else if (dds_slice[k][i, j] == CSflag)
                        //{
                        //    CSvolCount++;
                        //    pDDS[k][i, j] = CSdose;
                        //}
                        else if (dds_slice[k][i,j] > 0)
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

        public float[][,] RetrieveFinalizedDoseSpace()
        {                        
            DoseSpace = PrepareWeighted_DS(SliceWeights, folderpath, DDS);
            Debug.WriteLine("Before Normalization: ");
            FindTotalCoverage(0.5, DoseSpace, DDS, SliceWeights);
            DoseSpace = Matrix.Normalize(DoseSpace);
            Debug.WriteLine("After Normalization: ");
            FindTotalCoverage(0.5, DoseSpace, DDS, SliceWeights);
            return DoseSpace;
        }

        /// <summary>
        /// Initializes the DS matrix using a preliminary weight set and the original slicedoses.
        /// </summary>
        /// <param name="weights"></param>
        /// <param name="subfolder"></param>
        private float[][,] PrepareWeighted_DS(double[] weights, string subfolder, float[][,] DDS)
        {
            float[][,] weighted_slicedoses = new float[DoseSpace.GetLength(0)][,];
            
            int x = DoseSpace[0].GetLength(0); int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < weighted_slicedoses.GetLength(0); k++)
                weighted_slicedoses[k] = Matrix.Zeroes(x,y);
            PathSet.DCT = DCT;
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
            else
            {
                //Stopwatch cputimer = new Stopwatch();
                //cputimer.Start();
                for (int s = 0; s < NumSlices; s++)
                {                    
                    int whichz = SlicePositions[s] - (DCT / 2);
                    float[] slicedose = LoadSliceDose(s, subfolder);
                    slicedose = Matrix.Normalize(slicedose);
                    //FindSliceDoseCoverage(Matrix.ScalarMultiply(slicedose,(float)SliceWeights[s]), s, 0.5, DDS);
                    //Debug.WriteLineIf(slicedose.GetLength(0) != weighted_slicedoses.GetLength(0)*x*y, "slicedose not equal to weighted_slicedose in PrepareWeighted_DS()!!!!!");
                    //Debug.WriteLine(slicedose.Sum());
                    Parallel.For(0, DCT, (k) =>
                        {
                            for (int j = 0; j < y; j++)
                                for (int i = 0; i < x; i++)
                                {
                                    weighted_slicedoses[whichz + k][i, j] += (slicedose[(k * x * y) + (j * x) + i] * (float)weights[s]);
                                }
                        });
                }
                //Just took these out 10/15/2012, put back in if above doesn't work.
                //slicedose = Matrix.ScalarMultiply(slicedose, (float)weights[s]);                        
                //DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);
                //cputimer.Stop(); Debug.WriteLine("CPU time for PrepareWeightedDS: " + cputimer.Elapsed);
                //Debug.WriteLine("CPU_sum: " + Matrix.SumAll(weighted_slicedoses));
                //Debug.WriteLine("GPU_sum: " + Matrix.SumAll(GPUsd));
                return weighted_slicedoses;
            }
        }

        private float[] LoadSliceDose(int which_slice, string subfolder)
        {
            string filename = String.Concat("slice_", which_slice);
            string path = System.IO.Path.Combine(subfolder, filename);
            return Matrix.Normalize(ReadSliceDoseFromFile(path));
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
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            float dose; float dds_value; double Uncovered = 0;
            for (int k = 0; k < ds.GetLength(0); k++)                            
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        dose = ds[k][i, j];
                        dds_value = Tumor[k][i, j];
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

            // (isovolume / tumor);
            double LesionRxoverTumorVol = LesionRx / TumorVol;
            double LesionRxoverRxVol = LesionRx / RxVolume;
            double Underdosed = (Uncovered / TumorVol) * 100;
            string weights = "SliceWeights: " + recent_weight[0];
            for (int i = 1; i < recent_weight.GetLength(0); i++)
                weights += ", " + recent_weight[i];           
            Debug.WriteLine("**********Combined Coverage**************");
            Debug.WriteLine(weights);
            Debug.WriteLine("LesionRx/TumorVol: " + LesionRxoverTumorVol);
            Debug.WriteLine("LesionRx/RxVol: " + LesionRxoverRxVol);
            Debug.WriteLine("Percentage underdosed: " + Underdosed);
            WriteArrayAsList("Coverage: ", new double[4] { TumorVol, RxVolume, LesionRx, Uncovered });
            Debug.WriteLine("========================================");
            return new double[4] { TumorVol, RxVolume, LesionRx, Uncovered };
        }

        private double[] FindSliceDoseCoverage(float[][,] sd, int which_slice, double iso, float[][,] DDS)
        {
            double Coverage = 0; double TumorVol = 0; double LesionRx = 0; double RxVolume = 0;
            float dose; float dds_value; double Uncovered = 0;
            int x = DDS[0].GetLength(0);int y = DDS[0].GetLength(1); int z = DDS.GetLength(0);
            for (int k = 0; k < DDS.GetLength(0); k++)
                for (int j = 0; j < DDS[0].GetLength(1); j++)
                    for (int i = 0; i < DDS[0].GetLength(0); i++)
                    {
                        dose = sd[k][i,j];
                        dds_value = DDS[k][i, j];
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
            //Debug.WriteLine("SLICE " + which_slice);
            //Debug.WriteLine("==============================");
            //Debug.WriteLine("LesionRx/TumorVol: " + LesionRxoverTumorVol);
            //Debug.WriteLine("LesionRx/RxVol: " + LesionRxoverRxVol);
            //Debug.WriteLine("Percentage underdosed: " + Underdosed);
            double sum_sd = Matrix.SumAll(sd);
            double sum_dds = Matrix.SumAll(DDS);
            double ratio = sum_dds / sum_sd;
            return new double[5] { ratio, RxVolume / TumorVol, LesionRxoverTumorVol, LesionRxoverRxVol, Uncovered };
            
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
        private static float[][,] GrabSlab(float[][,] f, int thickness, int z_center)
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



        
        #endregion

        #region Helper methods that Read/Write dose to files

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
            //int startz = SlicePositions[which_slice] - DCT/2;

            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.

            if (GPU.GPUenabled)
                return output_ds;
               //return GPU.AddToDoseSpace(slicedose, output_ds, TranslateZBy, 1.0f);
            //else
                
            {
                Parallel.For(0, DCT, (z) =>
                    {
                        int current_z = TranslateZBy + z;
                        if (output_ds[current_z] == null)
                            output_ds[current_z] = Matrix.Zeroes(output_ds[0].GetLength(0), output_ds[0].GetLength(1));                        
                        output_ds[current_z] = Matrix.Add(output_ds[current_z], GrabSlice(slicedose, z, volume[0].GetLength(0), volume[0].GetLength(1)));
                        int dssum = Matrix.SumAll(output_ds[current_z]);
                        //WriteFloatArray2BMP(output_ds[current_z], "ds_z_" + z + "_" + dssum + ".bmp");
                        //WriteFloatArray2BMP(GrabSlice(slicedose, z, output_ds[0].GetLength(0), output_ds[0].GetLength(1)), "sd_z_" + z + ".bmp");
                    });
                return output_ds;
            }

            
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
            }
            //Normalize Dosespace
            //NormalizeDosespace();

            //Write dosespace to file.
            WriteDoseSpaceToFile("OriginalDS.txt");
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