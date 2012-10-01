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
        public static int edgepadding = 10;
        public static int tumorflag = 10;
        public static int CSflag = 2;
        public static float ToleranceDose = 0.2f;
        public static float RxDose = 0.8f;
        public static float CSdose = 0.1f;
        public int NumSlices; //Set in step 1
        public int TumorVolCount;
        public int CSvolCount;
        public ArrayList RasterPaths; //Collection of objects representing each slice
        public static string ActiveDirectory;
        public int DoseCalculationThickness; //How much of the dose kernel to use in calculation
        public static int DCT; //static version of dosecalcthickness
        public int SliceThickness;
        public int RasterWidth;
        public int TolThickness;
        public float[][,] volume; //Presumably the structure set binary volume?
        public float[][,] DoseSpace; //The actual final dosespace matrix containing the real dose.
        public int[] SlicePositions; //The centered z-indexes of each slice location
        public double[] SliceWeights; //Weights, from 0 to 1, of each slice. Calculated from Step 2.
        public int[] boundaries; //6-term sequence (xstart, xend, ystart...etc.)
        public int X; public int Y; public int Z;
        public event ProgressChangedEventHandler PathsetWorkerProgressChanged;
        public event RunWorkerCompletedEventHandler PathsetWorkerCompleted;
        public event ProgressChangedEventHandler RasterPathWorkerProgress;
        public event RunWorkerCompletedEventHandler OptimizationWorkerCompleted;
        public event ProgressChangedEventHandler OptimizationWorkerProgress;
        public event RunWorkerCompletedEventHandler SliceweightWorkerCompleted;
        public BackgroundWorker PS_ShotOptimize_worker; // <- called when "Optimize" button is clicked 
        public BackgroundWorker PS_SliceOptimize_worker; 
        public BackgroundWorker PS_CalcDose_worker; //<- called when the Calc/Save/Dose button is clicked.
        #endregion




        #region Constructors
        public PathSet(float[][,] f, int sthick, int tolthick, DoseKernel dk, StructureSet ss)
        {
            X = f[0].GetLength(0); Y = f[0].GetLength(1); Z = f.GetLength(0);
            DK = dk; SS = ss;
            N = dk.dose.GetLength(0);
            boundaries = FindBoundaries(f);
            SliceThickness = sthick;
            DoseCalculationThickness = SliceThickness * 2;
            TolThickness = tolthick;
            CalculateNumSlices();
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                RasterPath rp = new RasterPath(CompressSection(f, SlicePositions[i], SliceThickness / 2));
                RasterPaths.Add(rp);
            }            
            volume = f;
            AttachHandlers();
        }
        public PathSet(string dosespace_path)
        {
            DoseSpace = ReadDoseSpaceFromFile(dosespace_path);
        }
        private void AttachHandlers()
        {

            PS_ShotOptimize_worker = new BackgroundWorker();
            PS_ShotOptimize_worker.WorkerReportsProgress = true;
            PS_ShotOptimize_worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(PS_ShotOptimize_RunWorkerCompleted);
            PS_ShotOptimize_worker.ProgressChanged += new ProgressChangedEventHandler(PS_ShotOptimize_ProgressChanged);
            PS_ShotOptimize_worker.DoWork += new DoWorkEventHandler(PS_ShotOptimize_DoWork);
            PS_SliceOptimize_worker = new BackgroundWorker();
            PS_SliceOptimize_worker.WorkerReportsProgress = true;
            PS_SliceOptimize_worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(PS_SliceOptimize_worker_RunWorkerCompleted);
            PS_SliceOptimize_worker.ProgressChanged += PS_SliceOptimize_worker_ProgressChanged;
            PS_SliceOptimize_worker.DoWork += PS_SliceOptimize_worker_DoWork;
            PS_CalcDose_worker = new BackgroundWorker();
            PS_CalcDose_worker.WorkerReportsProgress = true;
            PS_CalcDose_worker.RunWorkerCompleted += PS_CalcDose_worker_RunWorkerCompleted;
            PS_CalcDose_worker.ProgressChanged += PS_CalcDose_worker_ProgressChanged;
            PS_CalcDose_worker.DoWork += PS_CalcDose_worker_DoWork;
           
        }

        
        #endregion


        #region All Background Worker methods

        void PS_ShotOptimize_DoWork(object sender, DoWorkEventArgs e)
        {
            double count = 0;
            for (int i = 0; i < NumSlices; i++)
            {
                //((RasterPath)RasterPaths[i]).RPworker.RunWorkerAsync();
                Debug.WriteLine("SLICE " + Convert.ToString(i) + ":");
                ((RasterPath)RasterPaths[i]).OptimizeShotWeights();
                count++;
                PS_ShotOptimize_worker.ReportProgress((int)(100 * count / NumSlices));
            }
        }
        void PS_ShotOptimize_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (PathsetWorkerProgressChanged != null)
                PathsetWorkerProgressChanged.Invoke(null, e);
        }
        void PS_ShotOptimize_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (PathsetWorkerCompleted != null)
                PathsetWorkerCompleted.Invoke(null, e);

            //Calculate the slicedoses so far...
            CreateDoseMatrix(DK, folderpath);

            
        }

        void PS_CalcDose_worker_DoWork(object sender, DoWorkEventArgs e)
        {
            CalculateAndWriteSliceDoses(DK, folderpath);
            AssembleDoseSpaceFromFiles();
        }
        void PS_CalcDose_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            throw new NotImplementedException();
        }
        void PS_CalcDose_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            Debug.WriteLine("Shot weighting optimization complete. Now commencing slice weighting...");
            if (OptimizationWorkerCompleted != null)
                OptimizationWorkerCompleted.Invoke(null, e);

            PS_SliceOptimize_worker.RunWorkerAsync();

        }
        void PS_SliceOptimize_worker_DoWork(object sender, DoWorkEventArgs e)
        {
            OptimizeSliceWeights();
        }
        void PS_SliceOptimize_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            throw new NotImplementedException();
        }
        void PS_SliceOptimize_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
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
            for (int s = 0; s < NumSlices; s++)
            {
                string filename = string.Concat("slice_", s);
                path = System.IO.Path.Combine(subfolder, filename);
                RasterPath rp = (RasterPath)RasterPaths[s];
                rp.CalculateAndSaveSliceDose(dk, DoseCalculationThickness, path);
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
            PS_CalcDose_worker.RunWorkerAsync();

        }
       
        /// <summary>
        /// Takes in a 1D float matrix, grabs a 2D slice out based on which z position. Called by GrabSlice().
        /// </summary>
        /// <param name="input"></param>
        /// <param name="which"></param>
        /// <param name="xsize"></param>
        /// <param name="ysize"></param>
        /// <returns></returns>
        private float[,] GrabSlice(float[] input, int which, int xsize, int ysize)
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
            double Error = 1000; int index = 0; double coverage = 0.8;
            SliceWeights = new double[NumSlices];
            double[] oldweights = new double[NumSlices];
            for (int i = 0; i < SliceWeights.GetLength(0); i++)
            {
                oldweights[i] = 1.0;
                SliceWeights[i] = 1.0;
            }
            //Start a while loop, and first create an updated DS matrix with the newest weights

            //UNDONE: Prepare DDS Tumor matrix 
            float[][,] DDS = PrepareDDS(SS.fj_Tumor);
            float[][,] OriginalDS = ReadDoseSpaceFromFile("OriginalDS.txt");
            

            while (Error >= .0001 && index <= 10)
            {
                //Reset the dosespace
                ClearDosespace();
                
                //Re-prepare dosespace with latest iteration of sliceweight
                PrepareWeighted_DS(SliceWeights, folderpath); //TODO: This method is time-consuming, make this GPU?
                
                //Evaluate each slice against the DDS slice
                for (int s = 0; s < NumSlices; s++)
                {
                    double DDS_slicesum = 0;
                    double DS_slicesum = 0;
                    int startz = FindStartZ(s);
                    int endz = FindEndZ(startz);
                    for (int z = 0; z < DoseCalculationThickness; z++)
                    {
                        /* The DDS matrix is element-multiplied by the newest dosespace
                         * and the sum is added to DDS_slicesum. The DS matrix is also
                         * multiplied by the latest iteration (squaring it?)
                         */
                        float[,] DDS_slice = DDS[startz + z]; Debug.WriteLine("Sum of DDS: " + Matrix.SumAll(DDS_slice));
                        float[,] DS_slice = DoseSpace[startz+z]; Debug.WriteLine("Sum of DS: " + Matrix.SumAll(DS_slice));
                        DDS_slicesum += Matrix.SumAll(Matrix.MultiplyElements(DDS_slice, OriginalDS[startz+z]));
                        DS_slicesum += Matrix.SumAll(Matrix.MultiplyElements(DS_slice, OriginalDS[startz+z]));
                    }
                    Debug.Assert(DDS_slicesum != 0);
                    Debug.Assert(DS_slicesum != 0);
                    double ratio = DDS_slicesum / DS_slicesum;
                    SliceWeights[s] = oldweights[s] * ratio;
                    Debug.WriteLine("Ratio=" + ratio + "; Weight: " + oldweights[s] + "-->" + SliceWeights[s]);
                }
                Error = FindError(SliceWeights, oldweights);

                Debug.WriteLine("Iteration: " + index + " Error: " + Error);


                Debug.WriteLine("Error:" + Error);
                double IterationCoverage = FindCoverage(RxDose/2, SS.fj_Tumor);
                index++;


                //
                if (index > 2 && IterationCoverage < coverage)
                {
                    SliceWeights = (double[])oldweights.Clone();
                    NormalizeDosespace();
                    break;
                }
                else
                {
                    coverage = IterationCoverage;
                    for (int s = 0; s < NumSlices; s++)
                    {
                        //Set the old weight = to the current weight
                        oldweights[s] = Convert.ToDouble(SliceWeights[s]);
                    }
                }
                //PS_ CalcDose_worker.ReportProgress(index);
            } // <- end of while loop
        }
        private float[][,] PrepareDDS(float[][,] dds_slice)
        {
            float[][,] pDDS = new float[dds_slice.GetLength(0)][,];
            TumorVolCount = 0; CSvolCount = 0;
            int _debug_tumorcount = 0;
            for (int k = 0; k < dds_slice.GetLength(0); k++)
            {
                pDDS[k] = new float[dds_slice[0].GetLength(0), dds_slice[0].GetLength(0)];
                for (int j = 0; j < dds_slice[0].GetLength(1); j++)
                    for (int i = 0; i < dds_slice[0].GetLength(0); i++)
                    {
                        Debug.Assert(dds_slice[k][i, j] < 2);
                        if (dds_slice[k][i, j] == 0)
                            pDDS[k][i, j] = ToleranceDose;
                        else if (dds_slice[k][i, j] == CSflag)
                        {
                            CSvolCount++;
                            pDDS[k][i, j] = CSdose;
                        }
                        else
                        {
                            TumorVolCount++;
                            if (dds_slice[k][i, j] >= 2)
                            {
                                Debug.WriteLine(dds_slice[k][i, j]);
                                _debug_tumorcount++;
                            }
                            pDDS[k][i, j] = RxDose;
                        }
                    }
            }
            //Debug.WriteLine("DDS counts >= 2: " + _debug_tumorcount);
            //Debug.WriteLine("DDS tumor counts: " + TumorVolCount);
            return pDDS;
        }
        
        
        /// <summary>
        /// Called from Step 2, uses the most recent round of optimized weights to read in each initial slice
        /// dose file (found in 'subfolder' path), multiply it by its respective weight, then add it to slicedose.
        /// </summary>
        /// <param name="weights"></param>
        /// <param name="subfolder"></param>
        private void PrepareWeighted_DS(double[] weights, string subfolder)
        {
            for (int s = 0; s < NumSlices; s++)
            {
                string filename = String.Concat("slice_", s);
                string path = System.IO.Path.Combine(subfolder, filename);
                float[] slicedose = ReadSliceDoseFromFile(path);
                //Debug.WriteLine(slicedose.Sum());
                slicedose = Matrix.ScalarMultiply(slicedose, (float)weights[s]);                
                DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);                
            }
        }

        

        private int FindStartZ(int s)
        {
            int start = SlicePositions[s] - (DoseCalculationThickness/2);
            if (start < 0)
                start = 0;
            return start;
        }
        private int FindEndZ(int startz)
        {
            int end = startz + (DoseCalculationThickness/2);
            if (end >= DoseSpace.GetLength(0))
                end = DoseSpace.GetLength(0)-1;
            return end;
        }




        /// <summary>
        /// Finds fraction of tumor that is covered by the desired dose. Called from the OptimizeSliceWeights()
        /// function of Step 2.
        /// </summary>
        /// <param name="desired_dose"></param>
        /// <param name="Tumor"></param>
        /// <returns></returns>
        private double FindCoverage(double desired_dose, float[][,] Tumor)
        {
            double Coverage = 0; double tumor = 0;
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
                for (int j = 0; j < DoseSpace[0].GetLength(1); j++)
                    for (int i = 0; i < DoseSpace[0].GetLength(0); i++)
                        if (Tumor[k][i, j] > 0) //if inside a tumor voxel
                        {
                            tumor++;
                            if (DoseSpace[k][i, j] >= desired_dose)
                                Coverage++;
                        }
                        else
                        {
                            continue;
                        }
            return (Coverage / tumor);
        }
       

        private void ClearDosespace()
        {
            DoseSpace = Matrix.ScalarMultiply(DoseSpace, 0.0f);

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
        //    int dstart = (N - DoseCalculationThickness) / 2; int dend = dstart + DoseCalculationThickness;

        //    float[][,] DDSslab = GrabSlab(SS.fj_Tumor, DoseCalculationThickness, slicecenter_z);

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
        //    int dstart = (N - DoseCalculationThickness) / 2; int dend = dstart + DoseCalculationThickness;

        //    //Convert the 1D float ds matrix into the jagged DS matrix, and then get the slab.            
        //    float[][,] DSslab = GrabSlab(dosespace, DoseCalculationThickness, z);

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
            int temp_start = z_center - thickness / 2;
            int temp_end = z_center + thickness / 2;

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
            int TranslateZBy = SlicePositions[which_slice] - (DoseCalculationThickness/2); // <- NEED TO CHANGE APPROPRIATELY
            int startz = SlicePositions[which_slice] - DoseCalculationThickness/2;

            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.
            for (int z = 0; z < DoseCalculationThickness; z++)
            {
                int current_z = TranslateZBy + z;
                if (output_ds[current_z] == null)
                    output_ds[current_z] = Matrix.Zeroes(volume[0].GetLength(0),volume[0].GetLength(1));               
                output_ds[current_z] = Matrix.Add(DoseSpace[current_z], GrabSlice(slicedose, z, volume[0].GetLength(0), volume[0].GetLength(1)));                
            }
            return output_ds;
        }
        /// <summary>
        /// Given a path, loads in the associated slicedose 1D float matrix. Called by AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="loadpath"></param>
        /// <returns></returns>
        public float[] ReadSliceDoseFromFile(string loadpath)
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
                string filename = String.Concat("slice_", s);
                string path = System.IO.Path.Combine(subfolder, filename);
                float[] slicedose = ReadSliceDoseFromFile(path);
                DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);
            }
            //Normalize Dosespace
            NormalizeDosespace();

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
        public float[][,] ReadDoseSpaceFromFile(string filename)
        {
            float[][,] ds;
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
                ds = new float[z][,];
                for (int k = 0; k < z; k++)
                {
                    float[,] slice = new float[x, y];
                    for (int j = 0; j < y; j++)
                        for (int i = 0; i < x; i++)
                            slice[i, j] = Convert.ToSingle(br.ReadLine());
                    ds[k] = slice;
                }
            }
            return ds;
        }

        #endregion
    }
}