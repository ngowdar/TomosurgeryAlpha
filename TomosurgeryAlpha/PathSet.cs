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
    ///   A collection class that contains the information about a group of RasterPath objects.
    ///   Takes tumor object, compresses it into slices and creates a RasterPath object 
    ///   for each slice.
    /// </summary>
    public class PathSet
    {
        #region Variables        

        public static float[,] mask;
        private static int N;
        public static int LineEdgePadding;
        public static int LineSidePadding;
        public static int tumorflag = 10;
        public static int CSflag = 2;
        public const float ToleranceDose = 0.2f;
        public const float RxDose = 0.5f;
        public const float CSdose = 0.05f;
        public static int NumSlices; //Set in step 1
        private static int TumorVolCount = 0;
        private static int CSVolCount = 0;
        
        public static string ActiveDirectory;

        public static int DCT; //static version of dosecalcthickness
        public static int StepSize;
        public static int RasterWidth;
        private float[][,] DDS;
        public static DoseKernel DK; //'DoseKernel', the object representing the dose matrix.
        private Boolean Dilate = true;
        public Boolean DoseModifiable;
        public float[][,] DoseSpace; //The actual final dosespace matrix containing the real dose.
        private double[] NonnormalizedCoverage;
        private double[] OldCoverage;
        private double[][] OldSliceCoverage;
        public BackgroundWorker Ps1ShotOptimizeWorker; // <- called when "Optimize" button is clicked 
        private BackgroundWorker Ps2CalcDoseWorker; //<- called when the Calc/Save/Dose button is clicked.
        private BackgroundWorker Ps3SliceWeightOptWorker;
        public ArrayList RasterPaths; //Collection of objects representing each slice
        private int[] ReferencePoint;
        private double[] ReferenceValues;
        private double[] ReferenceWeights;
        public StructureSet SS; //'SS' represents the binary tumor volume from the DICOM.
        public bool ShotsWeighted = false;
        private double[][] SliceCoverage;
        public int[] SlicePositions; //The centered z-indexes of each slice location
        private int SliceThickness;
        public static int SuggestedThickness;
        public static SuggThick ST;
        public double[] SliceWeights; //Weights, from 0 to 1, of each slice. Calculated from Step 2.
        private int TolThickness;
        private double[] TotalCoverage;
        public readonly int X;
        public readonly int Y;
        public readonly int Z;
        private readonly int[] boundaries; //6-term sequence (xstart, xend, ystart...etc.)
        private float[][,] dds_nondilated;
        public string folderpath = ActiveDirectory; //The active directory where slice data will be written.
        public float Max;
        public float[][,] Volume; //Presumably the structure set binary volume?
        public int VolumeZSize;
        public event ProgressChangedEventHandler PathsetWorkerProgressChanged;
        public event ProgressChangedEventHandler RasterPathWorkerProgress;
        public event RunWorkerCompletedEventHandler PathsetWorkerCompleted;
        public event RunWorkerCompletedEventHandler OptimizationWorkerCompleted;
        public event ProgressChangedEventHandler OptimizationWorkerProgress;
        public event RunWorkerCompletedEventHandler SliceweightWorkerCompleted;
        public event ProgressChangedEventHandler SliceweightWorkerProgress;

        #endregion

        #region Constructors

        public PathSet(float[][,] f, int sthick, int tolthick, DoseKernel dk, StructureSet ss)
        {
            X = f[0].GetLength(0);
            Y = f[0].GetLength(1);
            Z = f.GetLength(0);
            DK = dk;
            SS = ss;
            N = dk.dose.GetLength(0);
            boundaries = FindBoundaries(f);
            Volume = f;
            VolumeZSize = f.GetLength(0);

/*
            SliceThickness = sthick;
            DCT = SliceThickness*2;
            TolThickness = tolthick;*/
            //CalculateNumSlices();
            if (dk.DKI.Name == "8mmKernel")
                CalculateOptimumSliceThickness(32, 44);
            else if (dk.DKI.Name == "4mmKernel")
                CalculateOptimumSliceThickness(16, 24);
            ST = new SuggThick(SuggestedThickness);

        }

        public void ProcessVolume(int sthick, int tolthick)
        {
            SliceThickness = sthick;
            DCT = SliceThickness * 2;
            TolThickness = tolthick;
            ApplySliceThickness(sthick);
            RasterPaths = new ArrayList();
            for (var i = 0; i < NumSlices; i++)
            {
                /*var rp = new RasterPath(CompressSection(f, SlicePositions[i], SliceThickness/2),
                                        CompressSection(ss.fj_Combined, SlicePositions[i], SliceThickness/2), StepSize,
                                        RasterWidth, LineSidePadding, LineEdgePadding) {WhichSlice = i};*/
                var rp = new RasterPath(Volume, SlicePositions[i], sthick / 2, StepSize, RasterWidth, LineEdgePadding,
                                        LineSidePadding);
                rp.WhichSlice = i;
                RasterPaths.Add(rp);
            }
            
            AttachHandlers();
        }

        public PathSet(string dosespace_path)
        {
            var x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            DoseSpace = GPU.BackTo3D(ReadDoseSpaceFromFile(dosespace_path), x, y, z);
        }

        #endregion

        #region STEP 1: Shot Optimization (per slice)

        /// <summary>
        ///   Given a folder path and a DoseKernel, will create path names for each rasterpath in RasterPaths 
        ///   and write each resulting slicedose to path specified.
        /// </summary>
        /// <param name="dk"> </param>
        /// <param name="temppath"> </param>
        private void CalculateAndWriteSliceDoses(DoseKernel dk, string temppath)
        {
            DK = dk;
            folderpath = temppath;
            string subfolder = Path.Combine(temppath, DateTime.Now.ToString("yyyyMMddHHmmssfff"));
            Directory.CreateDirectory(subfolder);
            ActiveDirectory = subfolder;
            folderpath = subfolder;
            for (int s = 0; s < NumSlices; s++)
            {
                string filename = string.Concat("slice_", s);
                string path = Path.Combine(subfolder, filename);
                var rp = (RasterPath) RasterPaths[s];
                rp.CalculateAndSaveSliceDose(dk, DCT, path);
                var progress = (int) Math.Round((50.0/NumSlices));
                Ps2CalcDoseWorker.ReportProgress(progress);
            }
        }

        /// <summary>
        ///   Frontend method that calls PS_InitialDose_worker backgroundworker. 
        ///   Background worker runs CalculateSliceDosesAndWrite() and AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="dk"> </param>
        /// <param name="path"> </param>
        public void CreateDoseMatrix(DoseKernel dk, string path)
        {
            folderpath = path;
            DK = dk;
            Ps2CalcDoseWorker.RunWorkerAsync();
        }

        public void RecalculateSlices(int sthick, int tolthick)
        {
            SliceThickness = sthick;
            TolThickness = tolthick;
            CalculateNumSlices();
            RasterPaths = new ArrayList();
            for (int i = 0; i < NumSlices; i++)
            {
                /*var rp = new RasterPath(CompressSection(SS.fj_Tumor, SlicePositions[i], SliceThickness/2),
                                        CompressSection(SS.fj_Combined, SlicePositions[i], SliceThickness/2), StepSize,
                                        RasterWidth, LineEdgePadding, LineSidePadding);*/
                var rp = new RasterPath(SS.fj_Tumor, SlicePositions[i], SliceThickness / 2, StepSize, RasterWidth, LineEdgePadding,
                                        LineSidePadding);
                RasterPaths.Add(rp);
            }
        }

        private void CalculateOptimumSliceThickness(int low, int high)
        {
            int zstart = boundaries[4];
            int zend = boundaries[5];
            int meat = zend - zstart;
            int best_thick = low;
            int best_remainder = high;
            int best_numslices = 0;
            
            int BestShiftAmount = 100;
            Debug.WriteLine("Testing thicknesses of " + low + " to " + high + ":");
            Debug.WriteLine("Total z-length: " + meat);
            for (int i = low; i <= high; i++)
            {
                var padding = (int) Math.Round(0.5*i);
                int first = zstart + padding;
                int last = zend - padding;
                int numslices = meat / i;
                int remainder = ((zend - zstart) % i);
                int ShiftAmount = 0;
                if (remainder <= Math.Round(0.25 * i))
                    ShiftAmount = (-1) * remainder;
                else
                {
                    ShiftAmount = (i - remainder) / 2;
                    numslices = numslices + 1;
                }
                Debug.WriteLine("Testing " + i + ": ");
                Debug.WriteLine("Slices: " + numslices + ", Remainder: " + remainder);
                Debug.WriteLine("Shift Amount: " + ShiftAmount);
                
                if (remainder == 0)
                {
                    best_thick = i;
                    best_numslices = numslices;
                    best_remainder = remainder;
                }
                else if (remainder < best_remainder)
                {
                    if (ShiftAmount > 0 && ShiftAmount < BestShiftAmount)
                    {
                        BestShiftAmount = ShiftAmount;
                        best_thick = i;
                        best_numslices = numslices;
                        best_remainder = remainder;
                    }
                }
            }
            
            Debug.WriteLine("Best thickness is: " + best_thick);
           
            SuggestedThickness = best_thick;

        }

        private void ApplySliceThickness(int i)
        {
            int[] zrasterpos;
            int zstart = boundaries[4];
            int zend = boundaries[5];
            int meat = zend - zstart;

            var padding = (int)Math.Round(0.5 * i);
            int first = zstart + padding;
            int last = zend - padding;
            int numslices = meat / i;
            int remainder = ((zend - zstart) % i);
            int ShiftAmount = 0;

            if (remainder <= Math.Round(0.25 * i))
                ShiftAmount = (-1) * remainder;
            else
            {
                ShiftAmount = (i - remainder) / 2;
                numslices = numslices + 1;
            }
            int pad = (int) Math.Round(0.5*i);
            int ffirst = zstart + pad;
            zrasterpos = new int[numslices];
            if (remainder != 0)
                zrasterpos[0] = ffirst - ShiftAmount;
            else
                zrasterpos[0] = ffirst;
            for (int j = 1; j < numslices; j++)
                zrasterpos[j] = zrasterpos[0] + (j * i);

            NumSlices = zrasterpos.GetLength(0);
            SlicePositions = (int[])zrasterpos.Clone();
            SliceThickness = i;
            DCT = SliceThickness * 2;

        }

        private void CalculateNumSlices()
        {
            var padding = (int)Math.Round(0.5 * SliceThickness);
            int[] zrasterpos;
            int zstart = boundaries[4];
            int zend = boundaries[5];
            //Is there enough for 2 slices?
            if ((zend - zstart) < (1.75*SliceThickness))
                zrasterpos = new[] {(zend + zstart)/2};
            else
            {
               /* int meat = zend - zstart - (padding*2);
                int first = zstart + padding;
                int last = zend - padding;
                int numslices = meat/SliceThickness;
                if (meat % SliceThickness == 0)
                    numslices = (meat/SliceThickness) - 1;
                int newspacing = meat/(numslices + 1);
                numslices += 2;*/

                int meat = zend - zstart;
                int first = zstart + padding;
                int last = zend - padding;
                int numslices = meat/SliceThickness;
                int remainder = ((zend - zstart) % SliceThickness);
                int ShiftAmount = 0;
                if (remainder <= Math.Round(0.25 * SliceThickness))
                    ShiftAmount = (-1)*remainder;
                else
                {
                    ShiftAmount = (SliceThickness - remainder) / 2;
                    numslices = numslices + 1;    
                }
                zrasterpos = new int[numslices];
                zrasterpos[0] = first - ShiftAmount;
                
                for (int i = 1; i < numslices; i++)
                    zrasterpos[i] = zrasterpos[0] + (i*SliceThickness);
                    /*zrasterpos[i] = zrasterpos[i - 1] + newspacing;
                this.SliceThickness = newspacing;*/
            }
            NumSlices = zrasterpos.GetLength(0);
            SlicePositions = (int[])zrasterpos.Clone();
            DCT = SliceThickness*2;
        }

        public int[] FindBoundaries(float[][,] f)
        {
            var boundaries = new int[6]; //xstart, xend, ystart, yend, zstart, zend
            double dsum = Matrix.SumAll(f);
            //Z boundaries
            Boolean z1 = false;
            Boolean z2 = false;
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
                Boolean x1 = false;
                Boolean x2 = false;
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
                Boolean y1 = false;
                Boolean y2 = false;

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

        #endregion

        #region STEP 2: Slice Optimization

        /// <summary>
        ///   Main method of Step 2. Using a error-minimization strategy, the weights are adjusted
        ///   via DDS/DS ratios, just as in the shot-optimization of step 1.
        /// </summary>
        public void OptimizeSliceWeights()
        {
            //path = ActiveDirectory;
            Debug.Assert(folderpath != null);
            double Error = 1000;
            int index = 0;
            double coverage = 0.8;

            //Initialize the sliceweight matrices (temp_weights and SliceWeights)
            SliceWeights = new double[NumSlices];
            var temp_weights = new double[NumSlices];
            for (int i = 0; i < SliceWeights.GetLength(0); i++)
            {
                SliceWeights[i] = 1.0;
                temp_weights[i] = 1.0;
            }

            dds_nondilated = (float[][,]) PrepareDDS(SS.fj_Tumor).Clone();
            DDS = (float[][,]) PrepareEdgeEnhanced_DDS(SS.fj_Tumor).Clone();
            //DDS = (float[][,])PrepareDDS(SS.fj_Tumor).Clone();


            //Use initial sliceweights and find restrictions, test it out first.
            DoseSpace = PrepareWeighted_DS(SliceWeights, folderpath);
            double max = Matrix.FindMax(DoseSpace);
            Debug.WriteLine("Max value before sliceweight restriction: " + max);
            double[] restricted_weights = GetSliceWeightRestrictions(DoseSpace, SliceWeights, folderpath);
            DoseSpace = PrepareWeighted_DS(restricted_weights, folderpath);
            max = Matrix.FindMax(DoseSpace);
            Debug.WriteLine("Max value AFTER restriction: " + max);
            TotalCoverage = FindTotalCoverage(0.5, DoseSpace, dds_nondilated);
            double RxVolvsTumor = Math.Round(TotalCoverage[1]/TotalCoverage[0],3);
            double BVT = Math.Round(TotalCoverage[2]/TotalCoverage[0],3);
            Debug.WriteLine("Coverage with RestrictedWeights: ");
            Debug.WriteLine("Rx Volume / Tumor Volume (RTOG): " + RxVolvsTumor);
            Debug.WriteLine("Isodose Coverage: " + BVT);

            //while (max > 1.0)
            //{
            //    double[] dummy_priority = new double[2]{0,0};
            //    restricted_weights = SliceWeightPostProcess(DoseSpace, path, dummy_priority);
            //    DoseSpace = PrepareWeighted_DS(restricted_weights, path);
            //    SliceWeights = (double[])restricted_weights.Clone();
            //    max = Matrix.FindMax(DoseSpace);
            //}

            //SliceWeights = (double[])restricted_weights.Clone();
            //DoseSpace = PrepareWeighted_DS_GPU(SliceWeights, path, DDS);
            //WriteFloatArray2BMP(DoseSpace[DoseSpace.GetLength(0) / 2], "PS_487_ds_midplane_" + index + ".bmp");

            
            //NonnormalizedCoverage = new double[5];
            //NonnormalizedCoverage = FindTotalCoverage(0.5, DoseSpace, dds_nondilated);


            //TumorVol, RxVolume, LesionRx, Uncovered, Overdosed

            int MultiplierChoice = 1;
            if (index == 0)
                Error = 1000;
            double old_error = 1000;


            //TODO: Need to fix this whole loop. Currently goes to infinity.
            while (Error >= .01 && index < 5)
            {
                temp_weights = ReOptimizeSliceWeights(dds_nondilated, MultiplierChoice);
                //for (int i = 0; i < temp_weights.GetLength(0); i++)
                //    if (temp_weights[i] > restricted_weights[i])
                //        temp_weights[i] = restricted_weights[i];

                old_error = Convert.ToDouble(Error);

                //Re-prepare dosespace with latest iteration of sliceweight
                ClearDosespace();
                DoseSpace = PrepareWeighted_DS(temp_weights, folderpath);
                    //TODO: This method is time-consuming, make this GPU?

                //Restrict Slice Weights to avoid overdosing
                max = Matrix.FindMax(DoseSpace);


                if (max > 1.0)
                {
                    for (int i = 0; i < temp_weights.GetLength(0); i++)
                    {
                        if (temp_weights[i] > restricted_weights[i])
                        {
                            temp_weights[i] = restricted_weights[i];
                        }
                    }
                    DoseSpace = PrepareWeighted_DS(temp_weights, folderpath);
                    max = Matrix.FindMax(DoseSpace);
                }

                //while (max > 1.0)
                //{
                    //restricted_weights = SliceWeightPostProcess(DoseSpace, folderpath, temp_weights);
                    //DoseSpace = PrepareWeighted_DS(restricted_weights, folderpath);
                    //temp_weights = (double[]) restricted_weights.Clone();
                    //max = Matrix.FindMax(DoseSpace);
                //}

                //if (index > 0)
                Error = FindError(SliceWeights, temp_weights);


                //DoseSpace = PrepareWeighted_DS_GPU(temp_weights, path, DDS);
                //WriteFloatArray2BMP(DoseSpace[DoseSpace.GetLength(0) / 2], "PS_508_ds_midplane_" + index + ".bmp");
                //DoseSpace = Matrix.Normalize(DoseSpace);            
                TotalCoverage = FindTotalCoverage(RxDose, DoseSpace, dds_nondilated);
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
                double RTOG = TotalCoverage[1]/TotalCoverage[0];
                double BothvsTumor = TotalCoverage[2]/TotalCoverage[0];
                double BothvsRx = TotalCoverage[2]/TotalCoverage[1];
                var Conform_Indices = new double[3] {(1/RTOG), (1/BothvsTumor), (1/BothvsRx)};
                double IterationCoverage = TotalCoverage[2]/TotalCoverage[0];

                double PercentUnderdosed = (TotalCoverage[3]/TotalCoverage[0])*100;
                double TempOverage = TotalCoverage[1]/TotalCoverage[0];
                //index++;
                if (index > 1)
                {
                    if (IterationCoverage >= 1.0 && IterationCoverage < coverage) //coverage reversing/oscillating
                    {
                        Debug.WriteLine("Stopped bc coverage > 98%, coverage starting to decrease");
                        Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " +
                                        TempOverage);
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
                            Debug.WriteLine("Index: " + index + "; Coverage: " + IterationCoverage + "; Overage: " +
                                            TempOverage);
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
                        SliceWeights = (double[]) temp_weights.Clone();
                        index++;
                        Ps3SliceWeightOptWorker.ReportProgress(index*5);
                        continue;
                    }
                }
                else
                {
                    coverage = Convert.ToDouble(IterationCoverage);
                    SliceWeights = (double[]) temp_weights.Clone();
                    index++;
                    Ps3SliceWeightOptWorker.ReportProgress(index*5);
                    continue;
                }
            } // <- end of while loop
            Debug.WriteLine("While Loop ended. Error: " + Error + "; Index: " + index);
            Ps3SliceWeightOptWorker.ReportProgress(100);
            FindMaxDose();
        }

        private double[] ReOptimizeSliceWeights(float[][,] dds, int MultiplierChoice)
        {
            var tweight = (double[]) SliceWeights.Clone();
            if (SliceCoverage != null)
                OldSliceCoverage = (double[][]) SliceCoverage.Clone();
            else
                SliceCoverage = new double[NumSlices][];
            int count = 0;

            for (int s = 0; s < NumSlices; s++)
            {
                float[][,] ds_slab = GrabSlab(DoseSpace, DCT, SlicePositions[s]);
                //WriteFloatArray2BMP(ds_slab[ds_slab.GetLength(0)/2],String.Concat(s,"_ds_slab.bmp"));
                float[][,] dds_slab = GrabSlab(dds, DCT, SlicePositions[s]);
                //WriteFloatArray2BMP(dds_slab[dds_slab.GetLength(0)/2], String.Concat(s, "_dds_slab.bmp"));
                double[] measurements = FindSliceDoseCoverage(Matrix.ThresholdEq(ds_slab, 0.5f), s, 0.5, dds_slab);

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
                var bvt = (float) (1.0/BothvsTumor);
                var simplesum = (float) Math.Round(simplesum_ratio, 3);
                var rvt = (float) (1.0/RxVolvsTumor);

                double[] m1 = FindSliceDoseCoverage(Matrix.ThresholdEq(Matrix.ScalarMultiply(ds_slab, bvt), 0.5f), s,
                                                    0.5, dds_slab);
                double[] m2 = FindSliceDoseCoverage(Matrix.ThresholdEq(Matrix.ScalarMultiply(ds_slab, rvt), 0.5f), s,
                                                    0.5, dds_slab);

                MultiplierChoice = CompareImprovements(m1, m2);

                double ratio;
                switch (MultiplierChoice)
                {
                    case (1):
                        ratio = bvt;
                        break;
                    case (2):
                        ratio = rvt;
                        break;
                    case (3):
                        ratio = simplesum_ratio;
                        break;
                    default:
                        ratio = 1.0/BothvsTumor;
                        break;
                }
                //double ratio = 1.0 / RxVolvsTumor;
                //double ratio = simplesum_ratio;               
                Debug.WriteLine("Slice: " + s + " Tumor coverage: " + BothvsTumor);
                Debug.WriteLine("(Isovol/Tumor)x(Isotumor/Tumor): " + ratio + "; 1//BothVsTumor: " + (1.0/BothvsTumor) +
                                "; 1//RxVolVsTumor: " + (1.0/RxVolvsTumor));
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
                    if (ratio > 1.0 && (SliceWeights[s]*ratio) <= 1.0)
                        tweight[s] = SliceWeights[s]*ratio;
                    else if (ratio > 1.0 && (SliceWeights[s]*ratio) > 1.0)
                        tweight[s] = SliceWeights[s]*ratio;
                    else
                        tweight[s] = SliceWeights[s]*ratio;
                }

                SliceCoverage[s] = measurements;
                count++;
            }
            //tweight = Normalize(tweight);
            return tweight;
        }

        #region Prepare DDS / DS Volumes

        /// <summary>
        ///   Initializes the DS matrix using a preliminary weight set and the original slicedoses.
        /// </summary>
        /// <param name="weights"> </param>
        /// <param name="subfolder"> </param>
        private float[][,] PrepareWeighted_DS(double[] weights, string subfolder)
        {
            var weighted_slicedoses = new float[DoseSpace.GetLength(0)][,];

            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < weighted_slicedoses.GetLength(0); k++)
                weighted_slicedoses[k] = Matrix.Zeroes(x, y);
            //PathSet.DCT = DCT;

            #region GPU

            if (GPU.GPUenabled)
            {
                var GPUsd = new float[DoseSpace.GetLength(0)][,];
                var gputime = new Stopwatch();
                gputime.Start();
                var wSD = new float[x*y*z];
                for (int i = 0; i < wSD.GetLength(0); i++)
                    wSD[i] = 0.0f;

                //wSD = GPU.PrepareDoseSpace(wSD, SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
                wSD = GPU.WeightOriginalDS(SlicePositions, weights, new int[3] {x, y, z}, DCT, subfolder);
                gputime.Stop();
                Debug.WriteLine("GPU PrepareWeighted_DS time: " + gputime.Elapsed);
                GPUsd = GPU.BackTo3D(wSD, x, y, z);
                return GPUsd;
            }
                #endregion

            else
            {
                //Stopwatch cputimer = new Stopwatch();
                //cputimer.Start();
                Debug.Write("Preparing weighted DS.");
                var OverdosePoints = new ArrayList();
                for (int s = 0; s < NumSlices; s++)
                {
                    int GlobalStartZ = SlicePositions[s] - (DCT / 2); // <- Global start z
                    int GlobalEndZ = GlobalStartZ + DCT;
                    int LocalSliceStart = 0;
                    int LocalSliceEnd = DCT;
                    if (GlobalStartZ < 0)
                    {
                        LocalSliceStart = Math.Abs(GlobalStartZ);
                        GlobalStartZ = 0;
                        LocalSliceEnd = DCT;
                    }
                    else if (GlobalEndZ > DoseSpace.GetLength(0))
                    {
                        LocalSliceEnd = DoseSpace.GetLength(0) - GlobalStartZ;
                        GlobalEndZ = DoseSpace.GetLength(0);
                        
                    }
                    
                    
                    float[] slicedose = LoadSliceDose(s, subfolder);

                    //Parallel.For(0, DCT, (k) =>
                    //    {
                    //        for (int j = 0; j < y; j++)
                    //            for (int i = 0; i < x; i++)
                    //            {
                    //                weighted_slicedoses[whichz + k][i, j] += (slicedose[(k * x * y) + (j * x) + i] * (float)weights[s]);
                    //            }
                    //    });
                    
                    //for (int k = LocalSliceStart; k < LocalSliceEnd; k++)
                    Parallel.For(LocalSliceStart, LocalSliceEnd, (k) =>
                                                                     {
                                                                         Parallel.For(0, y, (j) =>
                                                                                                {
                                                                                                    for (int i = 0;
                                                                                                         i < x;
                                                                                                         i++)
                                                                                                    {
                                                                                                        float
                                                                                                            currentdose
                                                                                                                =
                                                                                                                Convert.
                                                                                                                    ToSingle
                                                                                                                    (weighted_slicedoses
                                                                                                                         [
                                                                                                                             GlobalStartZ +
                                                                                                                             (k -
                                                                                                                              LocalSliceStart)
                                                                                                                         ]
                                                                                                                         [
                                                                                                                             i,
                                                                                                                             j
                                                                                                                         ]);
                                                                                                        float
                                                                                                            originaldose
                                                                                                                =
                                                                                                                slicedose
                                                                                                                    [
                                                                                                                        (k*
                                                                                                                         x*
                                                                                                                         y) +
                                                                                                                        (j*
                                                                                                                         x) +
                                                                                                                        i
                                                                                                                    ];
                                                                                                        var sliceweight
                                                                                                            =
                                                                                                            (float)
                                                                                                            weights[s];
                                                                                                        float
                                                                                                            updated_dose
                                                                                                                =
                                                                                                                Convert.
                                                                                                                    ToSingle
                                                                                                                    (currentdose) +
                                                                                                                (originaldose*
                                                                                                                 sliceweight);
                                                                                                        weighted_slicedoses
                                                                                                            [
                                                                                                                GlobalStartZ +
                                                                                                                (k -
                                                                                                                 LocalSliceStart)
                                                                                                            ][i, j] =
                                                                                                            updated_dose;
                                                                                                    }
                                                                                                });

                                                                         OverdosePoints.TrimToSize();
                                                                     });
                    //


                    int progress = Convert.ToInt16((s*100.0)/NumSlices);
                    //PS_3_SliceWeightOpt_worker.ReportProgress(progress);
                    Debug.Write(".");
                }
                Debug.Write("done.");
                Debug.WriteLine("Post-processing weights to ensure no overdosing...");

                return weighted_slicedoses;
            }
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
                        ods[k + startz] = Matrix.Add(ods[k + startz],
                                                     GPU.ScalarMultiply(ods[k + startz], (float) SliceWeights[z]));
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
                                                       ods[k + startz] = Matrix.Add(ods[k + startz],
                                                                                    Matrix.ScalarMultiply(
                                                                                        ods[k + startz],
                                                                                        (float) SliceWeights[z]));
                                                   }
                                               });
            }
            return ods;
        }

        public static float[][,] PrepareDDS(float[][,] dds_slice)
        {
            var pDDS = new float[dds_slice.GetLength(0)][,];
            int _debug_tumorcount = 0;
            for (int k = 0; k < dds_slice.GetLength(0); k++)
            {
                pDDS[k] = new float[dds_slice[0].GetLength(0),dds_slice[0].GetLength(1)];
                for (int j = 0; j < dds_slice[0].GetLength(1); j++)
                    for (int i = 0; i < dds_slice[0].GetLength(0); i++)
                    {
                        Debug.Assert(dds_slice[k][i, j] < 2);
                        float val = dds_slice[k][i, j];
                        /*if (dds_slice[k][i, j] <= (ToleranceDose*2))
                            pDDS[k][i, j] = ToleranceDose;
                        else if (dds_slice[k][i, j] > (ToleranceDose*2))
                        {
                            TumorVolCount++;
                            pDDS[k][i, j] = RxDose;
                        }*/
                        if (val < 1)
                            pDDS[k][i, j] = ToleranceDose;
                        else if (val == 1)
                        {
                            TumorVolCount++;
                            pDDS[k][i, j] = RxDose;
                        }
                        else if (val > 1)
                        {
                            CSVolCount++;
                            pDDS[k][i, j] = CSdose;
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

            if (Dilate)
            {
                EE_DDS = PrepareDDS(DilateDDS(DilateDDS(dds)));
                //EE_DDS = PrepareDDS(DilateDDS(EE_DDS));
                //EE_DDS = PrepareDDS(DilateDDS(EE_DDS));
            }
            else
            {
                EE_DDS = PrepareDDS(dds);
            }
            return EE_DDS;
        }

        private float[][,] DilateDDS(float[][,] p)
        {
            var output = new float[p.GetLength(0)][,];
            /*
             *  0          A (i,j-1)          0
             *  B (i-1,j)     C(i,j)          D (i+1,j)
             *  0           E(i, j+1)        0
             * 
             * prev slice = (k-1,i,j), next slice (k+1,i,j)
             */

            float top;
            float right;
            float left;
            float bottom;
            float back;
            float front;
            //Starting at 1 instead of 0-index, to save having to worry about edges. Tumor wont go that far anyway.
            output[0] = (float[,]) p[0].Clone();
            output[p.GetLength(0) - 1] = (float[,]) p[0].Clone();

            for (int k = 1; k < p.GetLength(0) - 1; k++)
            {
                output[k] = (float[,]) p[k].Clone();
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
            //DoseSpace = PrepareWeighted_DS_GPU(SliceWeights, path, DDS);
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

        private float[][,] PrepareWeighted_DS_GPU(double[] weights, string subfolder, float[][,] DDS)
        {
            var weighted_slicedoses = new float[DoseSpace.GetLength(0)][,];
            bool restrictweight = false;
            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < weighted_slicedoses.GetLength(0); k++)
                weighted_slicedoses[k] = Matrix.Zeroes(x, y);
            DCT = DCT;


            //wSD = GPU.PrepareDoseSpace(wSD, SlicePositions, weights, new int[3] { x, y, z }, DCT, subfolder);
            float[] wSD = GPU.WeightOriginalDS(SlicePositions, weights, new int[3] {x, y, z}, DCT, subfolder);
            WriteFloatArray2BMP(GrabSlice(wSD, DoseSpace.GetLength(0)/2, x, y), "wSD.bmp");

            float[][,] GPUsd = GPU.BackTo3D(wSD, x, y, z);
            return GPUsd;
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
            var AllShots = new PointF[shotcount];
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

        #endregion

        #region Sliceweight Post-Optimization Processing

        private double[,] FindMaxPointsInEachSlice(float[][,] temp_ds)
        {
            double[,] MaxValues = new double[NumSlices, 4];
            for (int s = 0; s < NumSlices; s++)
            {
                int GlobalStartZ = SlicePositions[s] - (DCT / 2); // <- Global start z
                int GlobalEndZ = GlobalStartZ + DCT;
                int LocalSliceStart = 0;
                int LocalSliceEnd = DCT;
                if (GlobalStartZ < 0)
                {
                    LocalSliceStart = Math.Abs(GlobalStartZ);
                    GlobalStartZ = 0;
                    LocalSliceEnd = DCT;
                }
                else if (GlobalEndZ > DoseSpace.GetLength(0))
                {
                    LocalSliceEnd = DoseSpace.GetLength(0) - GlobalStartZ;
                    GlobalEndZ = DoseSpace.GetLength(0);

                }

               
                double max = 0;
                int max_x = -1;
                int max_y = -1;
                int max_z = -1;
                for (int k = LocalSliceStart; k < LocalSliceEnd; k++)
                    for (int j = 0; j < temp_ds[0].GetLength(1); j++)
                        for (int i = 0; i < temp_ds[0].GetLength(0); i++)
                        {
                            float value = temp_ds[GlobalStartZ + (k - LocalSliceStart)][i, j];
                            if (value > max)
                            {
                                max = value;
                                max_x = i;
                                max_y = j;
                                max_z = k;
                            }
                        }
                MaxValues[s,0] = max;
                MaxValues[s, 1] = max_x;
                MaxValues[s, 2] = max_y;
                MaxValues[s, 3] = max_z;
            }
            return MaxValues;
        }

        private double[] GetSliceWeightRestrictions(float[][,] ds, double[] weight, string subfolder)
        {
            var restrictweight = (double[])weight.Clone();
            double[,] slicemaxvalues = FindMaxPointsInEachSlice(ds);
            double maxval = 0;

            //Find the biggest max value to treat first.
            for (int i = 0; i < slicemaxvalues.GetLength(0); i++)
                if (slicemaxvalues[i, 0] > maxval)
                    maxval = slicemaxvalues[i, 0];

            double maxdose = Matrix.FindMax(ds);
            int counter = 0;
            bool stop = false;

            #region WhileLoop
            /*while (maxval > 1.05)
            {
                int[] WhichSlices = new int[NumSlices];
                int bestcovslice = 0;
                int bestrtogslice = 0;
                double bestcoverage = 0;
                double mostoverdosed = 0;
                double bestcovslice_val = 0;
                double bestrtogslice_val = 0;
                double[] impt_d = new double[5];
                double[] cov_d = new double[5];
                double[] rtog_d = new double[5];
                

                float best_unweighted_value = 0;
                int mostimportantslice = 0;

                for (int i = 0; i < NumSlices; i++)
                {
                    double slicevalue = slicemaxvalues[i, 0];
                    int z = (int)slicemaxvalues[i, 3];
                    int x = (int)slicemaxvalues[i, 1];
                    int y = (int)slicemaxvalues[i, 2];
                    
                    //Check if this slice contains the max value.
                    const double EPSILON = 0.001;
                    if (Math.Abs(slicevalue - maxdose) < EPSILON)
                    {
                        //Find the slicecoverage of this slice, store it
                        //Find both coverage and Rx/Tumor values
                        //{tumor, rx, both};
                        //float[][,] slicedose = Matrix.ConvertTo3D(LoadSliceDose(i, subfolder), DoseSpace[0].GetLength(0), DoseSpace[0].GetLength(1),
                          //                   DoseSpace.GetLength(0));
                        //float value = slicedose[z][x, y];
                        int index = z * DoseSpace[0].GetLength(0) * DoseSpace[0].GetLength(1) + y * DoseSpace[0].GetLength(0) + x;
                        float[] sd = LoadSliceDose(i, subfolder);
                        float unweighted_slicevalue = sd[index];
                        double[] d = FindSliceCoverage(ds, i, 0.5, DDS);
                        double cov = d[2] / d[0];
                        double rtog = d[1] / d[0];
                        double currentweight = restrictweight[i];
                        if (unweighted_slicevalue > best_unweighted_value)
                        {
                            best_unweighted_value = unweighted_slicevalue;
                            mostimportantslice = i;
                            impt_d = new double[] { i, unweighted_slicevalue, currentweight, cov, rtog };
                        }                        
                        if (cov > bestcoverage)
                        {   
                            bestcoverage = cov;
                            bestcovslice = i;
                            bestcovslice_val = unweighted_slicevalue;
                            cov_d = new double[] { i, unweighted_slicevalue, currentweight, cov, rtog };
                        }
                        if (rtog > mostoverdosed)
                        {
                            mostoverdosed = rtog;
                            bestrtogslice = i;
                            bestrtogslice_val = unweighted_slicevalue;
                            rtog_d = new double[] { i, unweighted_slicevalue, currentweight, cov, rtog };
                        }
                    }
                }
                double current_weight = restrictweight[mostimportantslice];
                double current_val = 0;
                double chosenslice = 0;
                if (impt_d[0] != rtog_d[0])
                {
                    if ((rtog_d[1] * rtog_d[2]) > (impt_d[1] * impt_d[2]))
                    {
                        chosenslice = rtog_d[0];
                        current_weight = rtog_d[2];
                        current_val = rtog_d[1];
                    }
                    else
                    {
                        chosenslice = impt_d[0];
                        current_weight = impt_d[2];
                        current_val = impt_d[1];
                    }
                }
                else if (impt_d[0] != cov_d[0])
                {
                    if ((cov_d[1] * cov_d[2]) > (impt_d[1] * impt_d[2]))
                    {
                        chosenslice = cov_d[0];
                        current_weight = cov_d[2];
                        current_val = cov_d[1];
                    }
                    else
                    {
                        chosenslice = impt_d[0];
                        current_weight = impt_d[2];
                        current_val = impt_d[1];
                    }
                }
                else
                {
                    chosenslice = impt_d[0];
                    current_weight = impt_d[2];
                    current_val = impt_d[1];
                }
                


                //Now that you have the slice that is closest to the culprit point, find out how much to reduce it by.                

                double leftover = (current_val * current_weight) - (maxval - 1.0);
                double contribution = current_val * current_weight;
                double percent_contribution = contribution / maxval;
                double desired_multiplier = ((contribution) - (maxval - 1.0)) / (contribution);
                double desired_weight = current_weight * desired_multiplier;
                if (stop == false)
                    desired_weight = ((current_weight - desired_weight) * 0.3) + desired_weight;
                //desired_weight = 1.0 / (current_val + (maxval - current_val));
                
                double desired_contribution = desired_weight * current_val;
                //desired_weight = (1.0 - leftover) * current_weight;
                //double newval = desired_weight * best_unweighted_value;


                /* Idea here is to reduce the sliceweight of slices containing the max value. However,
                 * rather than reduce the weight of all slices at once, its better to reduce one and then
                 * see how it affects the other slices and max values. To choose which slice to begin with,
                 * check to see which slice has the "most to lose" (i.e. had great coverage to start with, or 
                 * had an excess of RTOG compared to tumor vol (RTOG/tumorvol).
                 * 
                 #1#
                Debug.WriteLine("Lowering Slice " + mostimportantslice + " from " + restrictweight[mostimportantslice] + " to " + desired_weight);
                //restrictweight[bestrtogslice] = restrictweight[bestrtogslice]/maxval;
                restrictweight[(int)chosenslice] = desired_weight;
                counter++;
                if (restrictweight[bestrtogslice] < 0.005)
                {
                    restrictweight[bestrtogslice] = 0;
                    break;
                }
                ds = PrepareWeighted_DS(restrictweight, subfolder);
                maxdose = Matrix.FindMax(ds);
                Debug.WriteLine("Current max dose: " + maxdose);
                slicemaxvalues = FindMaxPointsInEachSlice(ds);
                maxval = 0;
                for (int i = 0; i < slicemaxvalues.GetLength(0); i++)
                    if (slicemaxvalues[i, 0] > maxval)
                        maxval = slicemaxvalues[i, 0];
                if (counter >= 5)
                    stop = true;
            }*/
            #endregion

            #region Another Method
            while (maxval > 1.05)
            {
                int[] whichslices = new int[NumSlices];
                double[] contributing_values = new double[NumSlices];
                for (int i = 0; i < NumSlices; i++)
                {
                    double slicevalue = slicemaxvalues[i, 0];
                    int z = (int)slicemaxvalues[i, 3];
                    int x = (int)slicemaxvalues[i, 1];
                    int y = (int)slicemaxvalues[i, 2];

                    const double EPSILON = 0.001;
                    if (Math.Abs(slicevalue - maxdose) < EPSILON)
                        whichslices[i] = 1;
                    else
                        whichslices[i] = 0;
                }

                double biggest_contributing_value = 0;
                int major_contributing_slice = -1;
                double desired_weight_multiplier = 1.0;

                for (int i = 0; i < NumSlices; i++)
                {
                    if (whichslices[i] > 0)
                    {
                        double slicevalue = slicemaxvalues[i, 0];
                        int z = (int)slicemaxvalues[i, 3];
                        int x = (int)slicemaxvalues[i, 1];
                        int y = (int)slicemaxvalues[i, 2];

                        int index = z * DoseSpace[0].GetLength(0) * DoseSpace[0].GetLength(1) + y * DoseSpace[0].GetLength(0) + x;
                        float[] sd = LoadSliceDose(i, subfolder);
                        float unweighted_slicevalue = sd[index];
                        double weighted_contribution = unweighted_slicevalue * restrictweight[i];
                        if (weighted_contribution > biggest_contributing_value)
                        {
                            biggest_contributing_value = weighted_contribution;
                            major_contributing_slice = i;
                            //double desired_contribution = weighted_contribution - ((maxval - 1.0) / maxval);
                            double desired_contribution = weighted_contribution / maxval;
                            desired_weight_multiplier = desired_contribution / weighted_contribution;
                        }
                        contributing_values[i] = weighted_contribution;
                    }
                    else
                        continue;
                }
                if (major_contributing_slice >= 0)
                    restrictweight[major_contributing_slice] = restrictweight[major_contributing_slice] * desired_weight_multiplier;

                ds = PrepareWeighted_DS(restrictweight, subfolder);
                maxdose = Matrix.FindMax(ds);
                Debug.WriteLine("Current max dose: " + maxdose);
                slicemaxvalues = FindMaxPointsInEachSlice(ds);
                maxval = 0;
                for (int i = 0; i < slicemaxvalues.GetLength(0); i++)
                    if (slicemaxvalues[i, 0] > maxval)
                        maxval = slicemaxvalues[i, 0];
            }
            #endregion
            WriteArrayAsList("RestrictedWeights: ", restrictweight);
            return restrictweight;
        }

        private double[] SliceWeightPostProcess(float[][,] ds, string subfolder, double[] iteration_weight)
        {
            float max = 0;
            int maxindex = 0;
            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            for (int k = 0; k < ds.GetLength(0); k++)
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        float value = ds[k][i, j];
                        if (value > max)
                        {
                            max = value;
                            maxindex = (k*x*y) + (j*x) + i;
                        }
                    }
            //First, make a list of the index ranges of each slice.
            var SliceIndexStartPositions = new int[NumSlices];
            //SliceIndexStartPositions[0] = 0;
            for (int s = 0; s < NumSlices; s++)
                SliceIndexStartPositions[s] = (SlicePositions[s] - (DCT/2))*x*y;


            var RelativeIndices = new int[NumSlices];
            for (int s = 0; s < NumSlices; s++)
            {
                RelativeIndices[s] = (maxindex - SliceIndexStartPositions[s]);
                if (RelativeIndices[s] < 0)
                    RelativeIndices[s] = 0;
                else if (RelativeIndices[s] > (DCT*x*y))
                    RelativeIndices[s] = 0;
            }


            //Third, find the weight contributions of each slice to maxvalue
            var ContributingWeights = new double[NumSlices];
            var Values = new double[NumSlices];
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
                sum += (Values[i]*ContributingWeights[i]);
            double multiplier = 1.0/sum;
            var weights = new double[NumSlices];

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

                weights[s] = SliceWeights[s]*multiplier;
            }
            return weights;
        }

        public void FindReferenceContributions(int[] refpoint, string subfolder)
        {
            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            int refindex = (refpoint[2]*x*y) + (refpoint[1]*x) + refpoint[0];
            ReferenceWeights = new double[NumSlices];
            ReferenceValues = new double[NumSlices];
            for (int s = 0; s < NumSlices; s++)
            {
                ReferenceWeights[s] = 0.0;
                ReferenceValues[s] = 0.0;
                int whichz = SlicePositions[s] - (DCT/2);
                float[] slicedose = LoadSliceDose(s, subfolder);
                Parallel.For(0, DCT, (k) =>
                                         {
                                             for (int j = 0; j < y; j++)
                                                 for (int i = 0; i < x; i++)
                                                 {
                                                     int tempindex = ((whichz + k)*x*y) + (j*x) + i;
                                                     if (tempindex == refindex)
                                                     {
                                                         ReferenceWeights[s] = SliceWeights[s];
                                                         ReferenceValues[s] = slicedose[(k*x*y) + (j*x) + i];
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
                Debug.WriteLine("Slice " + s + ": " + ReferenceValues[s] + " x (Weight: " + ReferenceWeights[s] + ") = " +
                                Math.Round(ReferenceWeights[s]*ReferenceValues[s], 3));
            }
        }

        public void FindDoseContributionsToReferencePoint(double[] normweights, double doserate, double max)
        {
            var DoseContributions = new double[NumSlices];
            double sum = 0;
            for (int s = 0; s < NumSlices; s++)
            {
                DoseContributions[s] = ReferenceValues[s]*normweights[s];
                sum += DoseContributions[s];
                Debug.WriteLine("Slice " + s + ": " + Math.Round(ReferenceValues[s], 3) + " x (Weight: " +
                                Math.Round(normweights[s], 3) + ") = " + Math.Round(DoseContributions[s], 3));
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
                sum += ReferenceValues[s]*ReferenceWeights[s];
            weight_multiplier = (max/sum);
            var NormalizedWeights = new double[NumSlices];
            for (int s = 0; s < NumSlices; s++)
                NormalizedWeights[s] = SliceWeights[s]*weight_multiplier;
            return NormalizedWeights;
        }

        #endregion

        #region Coverage and Evaluation Functions

        /// <summary>
        ///   Finds fraction of tumor that is covered by the desired dose. Called from the OptimizeSliceWeights()
        ///   function of Step 2.
        /// </summary>
        /// <param name="iso"> </param>
        /// <param name="ds"> </param>
        /// <param name="Tumor"> </param>
        /// <returns> </returns>
        private double[] FindTotalCoverage(double iso, float[][,] ds, float[][,] Tumor)
        {
            if (TotalCoverage != null)
                OldCoverage = (double[]) TotalCoverage.Clone();
            //float[][,] normds = Matrix.Normalize(ds);
            double Coverage = 0;
            double TumorVol = 0;
            double LesionRx = 0;
            double RxVolume = 0;
            float dose;
            float dds_value;
            double Uncovered = 0;
            for (int k = 0; k < ds.GetLength(0); k++)
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        dose = ds[k][i, j];
                        dds_value = Tumor[k][i, j];
                        if (dose >= iso)
                        {
                            RxVolume++;
                            if (dds_value > ToleranceDose*2)
                            {
                                TumorVol++;
                                LesionRx++;
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
            double Overdosed = TumorVol/RxVolume;
            return new double[5] {TumorVol, RxVolume, LesionRx, Uncovered, Overdosed};
        }

        private double[] FindSliceCoverage(float[][,] DS, int which_slice, double iso, float[][,] dds)
        {
            int position = SlicePositions[which_slice];
            int halfdct = DCT / 2;
            int startpoint = position - (DCT / 2);
            int endpoint = startpoint + DCT;
            if (startpoint < 0)
                startpoint = 0;
            if (endpoint > DS.GetLength(0))
                endpoint = DS.GetLength(0) - 1;
            int tumorsum = 0;
            int rxsum = 0;
            int bothsum = 0;
            float dose = 0;
            float ddsval = 0;
            for (int k = startpoint; k <= endpoint; k++)
                    for (int j = 0; j < DS[0].GetLength(1); j++)
                                for (int i = 0; i < DS[0].GetLength(0); i++)
                                {
                                    dose = DS[k][i, j];
                                    ddsval = dds[k][i, j];
                                    if (dose < iso)
                                    {
                                        if (ddsval > ToleranceDose)
                                            tumorsum++;
                                    }
                                    else if (dose >= iso)
                                    {
                                        rxsum++;
                                        if (ddsval > ToleranceDose)
                                        {
                                            tumorsum++;
                                            bothsum++;
                                        }

                                    }

                                }
            return new double[]{tumorsum, rxsum, bothsum};
        }

        private double[] FindSliceDoseCoverage(float[][,] sd, int which_slice, double iso, float[][,] DDS)
        {
            //float[][,] sd = Matrix.Normalize(ds);
            double Coverage = 0;
            double TumorVol = 0;
            double LesionRx = 0;
            double RxVolume = 0;
            float dose;
            float dds_value;
            double Uncovered = 0;
            int x = DDS[0].GetLength(0);
            int y = DDS[0].GetLength(1);
            int z = DDS.GetLength(0);
            for (int k = 0; k < DDS.GetLength(0); k++)
                for (int j = 0; j < DDS[0].GetLength(1); j++)
                    for (int i = 0; i < DDS[0].GetLength(0); i++)
                    {
                        dose = sd[k][i, j];
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
            double LesionRxoverTumorVol = LesionRx/TumorVol;
            double LesionRxoverRxVol = LesionRx/RxVolume;
            double Underdosed = (Uncovered/TumorVol);
            double Overdosed = RxVolume/TumorVol;
            //Debug.WriteLine("SLICE " + which_slice);
            //Debug.WriteLine("==============================");
            //Debug.WriteLine("LesionRx/TumorVol: " + LesionRxoverTumorVol);
            //Debug.WriteLine("LesionRx/RxVol: " + LesionRxoverRxVol);
            //Debug.WriteLine("Percentage underdosed: " + Underdosed);
            double sum_sd = Matrix.SumAll(sd);
            double sum_dds = Matrix.SumAll(DDS);
            double ratio = sum_dds/sum_sd;
            //return new double[5] { ratio, RxVolume / TumorVol, LesionRxoverTumorVol, LesionRxoverRxVol, Uncovered };
            return new double[6] {ratio, RxVolume, TumorVol, LesionRx, Underdosed, Overdosed};
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

        private int CompareImprovements(double[] m1, double[] m2)
        {
            double sumratio1 = m1[0];
            double RxVolume1 = m1[1];
            double TumorVol1 = m1[2];
            double BothVol1 = m1[3];
            double Uncovered1 = m1[4];
            double cov1 = m1[3]/m1[2];

            double sumratio2 = m2[0];
            double RxVolume2 = m2[1];
            double TumorVol2 = m2[2];
            double BothVol2 = m2[3];
            double Uncovered2 = m2[4];
            double cov2 = m2[3]/m2[2];

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

        #endregion

        #endregion

        #region Helper Methods (Matrix handling, I/O, Background Workers)

        #region Matrix Functions (Normalize, Max, etc.)

        public void FindMaxDose()
        {
            Max = Matrix.FindMax(DoseSpace);
            ReferencePoint = new int[3];
            ReferencePoint = FindMaxDoseVoxel();
        }

        public int[] FindMaxDoseVoxel()
        {
            var point = new int[3];
            float max = 0;
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
                for (int j = 0; j < DoseSpace[0].GetLength(1); j++)
                    for (int i = 0; i < DoseSpace[0].GetLength(0); i++)
                    {
                        float dose = DoseSpace[k][i, j];
                        if (dose > max)
                        {
                            max = dose;
                            point = new int[3] {i, j, k};
                        }
                    }
            return point;
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
                        DoseSpace[k][i, j] = DoseSpace[k][i, j]/max;
        }

        private double[] Normalize(double[] d)
        {
            double max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                if (d[i] > max)
                    max = d[i];
            for (int j = 0; j < d.GetLength(0); j++)
                d[j] = (float) (d[j]/max);
            return d;
        }

        public static float[][,] GrabSlab(float[][,] f, int thickness, int z_center)
        {
            //Find the start and end points based on the thickness
            int temp_start = z_center - (thickness/2);
            int temp_end = z_center + (thickness/2);

            //Check if the endpoints exist, change them if necessary 
            if (temp_start < 0)
                temp_start = 0;
            if (temp_end >= f.GetLength(0))
                temp_end = f.GetLength(0) - 1;

            //Initialize new matrix, fill, return
            var slab = new float[temp_end - temp_start][,];
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
            var slab = new float[temp_end - temp_start][,];
            for (int i = 0; i < (temp_end - temp_start); i++)
            {
                slab[i] = f[i + temp_start];
            }
            return slab;
        }

        public static float[,] GrabSlice(float[] input, int which, int xsize, int ysize)
        {
            var result = new float[xsize,ysize];
            for (int j = 0; j < ysize; j++)
                for (int i = 0; i < xsize; i++)
                {
                    int index = (which*xsize*ysize) + (j*xsize) + i;
                    result[i, j] = input[index];
                }
            return result;
        }

        private void ClearDosespace()
        {
            //DoseSpace = new float[DoseSpace.GetLength(0)][,];
            for (int i = 0; i < DoseSpace.GetLength(0); i++)
                DoseSpace[i] = Matrix.Zeroes(DoseSpace[0].GetLength(0), DoseSpace[0].GetLength(1));
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
                end = DoseSpace.GetLength(0) - 1;
            return end;
        }

        public float[,] CompressSection(float[][,] f, int zt, int spt)
        {
            var squished = new float[X,Y];
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

        public void AddSliceDoseToDoseSpace(float[][,] DS, int which_slice)
        {
            PointF[] shots = ((RasterPath) RasterPaths[which_slice]).shots;
            double[] weights = ((RasterPath) RasterPaths[which_slice]).ShotWeights;
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
                            ds_x = (int) (shots[shot].X - ((N - 1)/2) + i);
                            ds_y = (int) (shots[shot].Y - ((N - 1)/2) + j);
                            ds_z = (SlicePositions[which_slice] - ((N - 1)/2) + k);
                            index = ds_z*(StructureSet.BIG_dim[0]*StructureSet.BIG_dim[1]) +
                                    (ds_y*StructureSet.BIG_dim[0]) + ds_x;
                            if (ds_x < 0 || ds_y < 0 || ds_z < 0)
                                continue;
                            else if (ds_x >= StructureSet.BIG_dim[0] || ds_y >= StructureSet.BIG_dim[1] ||
                                     ds_z >= StructureSet.BIG_dim[2])
                                continue;
                            else if (index < 0 ||
                                     index >
                                     (StructureSet.BIG_dim[0]*StructureSet.BIG_dim[1]*StructureSet.BIG_dim[2]))
                                continue;
                            else
                            {
                                DS[ds_z][ds_x, ds_y] += (float) (DK.dose[k][(j*N) + i]*weights[shot]);
                            }
                        }
                    }
        }

        public void Calculate3DDoseSpace(DoseKernel dk)
        {
            //PointF[] AllShots = GetAllShotsArray();
            var DoseSpace = new float[StructureSet.BIG_dim[2]][,];
            for (int k = 0; k < DoseSpace.GetLength(0); k++)
            {
                //Create a new 2D slice.
                float[,] layer = Matrix.Zeroes(StructureSet.BIG_dim[0], StructureSet.BIG_dim[1]);
                DoseSpace[k] = (float[,]) layer.Clone();
            }
            for (int i = 0; i < SlicePositions.GetLength(0); i++)
                AddSliceDoseToDoseSpace(DoseSpace, i);
        }

        #endregion

        #region Reading/Writing to Files and BMPs

        public void WritePathSetInfoToReport()
        {
            Analysis.AddLineToReport("==============PLAN SUMMARY================");
            Analysis.AddLineToReport("Slice Thickness: " + SliceThickness);
            Analysis.AddLineToReport("Dose Calculation Thickness: " + DCT);
            Analysis.AddLineToReport("Number of Slices: " + NumSlices);
            Analysis.AddLineToReport(WriteArrayAsList("Tumor boundaries: ", boundaries));

            Analysis.AddLineToReport("==========================================");
        }

        private void PrintIterationSummary(double[][] OldCoverage, double[][] SliceCoverage)
        {
            for (int i = 0; i < NumSlices; i++)
            {
                Debug.WriteLine("====================================================");
                Debug.WriteLine("====================SLICE " + (i + 1) + "===================");
                Debug.WriteLine("Sum Ratio: " + Math.Round(OldCoverage[i][0], 2) + " -> " +
                                Math.Round(SliceCoverage[i][0], 2));
                Debug.WriteLine("Isovolume: " + Math.Round(OldCoverage[i][1], 2) + " -> " +
                                Math.Round(SliceCoverage[i][1], 2));
                Debug.WriteLine("Tumor volume coverage: " + Math.Round(100*(OldCoverage[i][2]/OldCoverage[i][3]), 2) +
                                "% -> " + Math.Round(100*(SliceCoverage[i][2]/SliceCoverage[i][3]), 2) + "%");
                Debug.WriteLine("Underdosed %: " + Math.Round(100*(OldCoverage[i][4]/OldCoverage[i][2]), 2) + "% -> " +
                                Math.Round(100*(SliceCoverage[i][4]/SliceCoverage[i][2]), 2) + "%");
                Debug.WriteLine("====================================================+");
            }
        }

        private void PrintIterationSummary(double[] OldCoverage, double[] TotalCoverage, int iteration, double error)
        {
            Debug.WriteLine("====================================================");
            Debug.WriteLine("====================Iteration " + iteration + "===================");
            //Debug.WriteLine("Sum Ratio: " + Math.Round(OldCoverage[0], 2) + " -> " + Math.Round(TotalCoverage[0], 2));
            Debug.WriteLine("Isovolume: " + Math.Round(OldCoverage[1], 2) + " -> " + Math.Round(TotalCoverage[1], 2));
            Debug.WriteLine("Tumor volume coverage: " + Math.Round(100*(OldCoverage[2]/OldCoverage[0]), 2) + "% -> " +
                            Math.Round(100*(TotalCoverage[2]/TotalCoverage[0]), 2) + "%");
            Debug.WriteLine("Underdosed %: " + Math.Round(100*(OldCoverage[3]/OldCoverage[0]), 2) + "% -> " +
                            Math.Round(100*(TotalCoverage[3]/TotalCoverage[0]), 2) + "%");
            Debug.WriteLine("Overdosed voxels: " + Math.Round(OldCoverage[4], 2) + " -> " +
                            Math.Round(TotalCoverage[4], 2));
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

        private void MakeMidplaneXPic(float[][,] DoseSpace, string p)
        {
            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            var X = new float[y,z];
            var Xdds = new float[y,z];
            for (int k = 0; k < z; k++)
                for (int j = 0; j < y; j++)
                {
                    X[j, k] = DoseSpace[k][x/2, j];
                    Xdds[j, k] = dds_nondilated[k][x/2, j];
                }
            WriteFloatArray2BMP(X, p);
            WriteFloatArray2BMP(Xdds, "X_midplane_DS.bmp");
            MakeComparisonPic(X, Xdds, "X_midplane_comparison.bmp");
        }

        private void MakeMidplaneYpic(float[][,] DoseSpace, string p)
        {
            int x = DoseSpace[0].GetLength(0);
            int y = DoseSpace[0].GetLength(1);
            int z = DoseSpace.GetLength(0);
            var Y = new float[DoseSpace[0].GetLength(0),DoseSpace.GetLength(0)];
            var Ydds = new float[DoseSpace[0].GetLength(0),DoseSpace.GetLength(0)];
            for (int k = 0; k < z; k++)
                for (int i = 0; i < x; i++)
                {
                    Y[i, k] = DoseSpace[k][i, y/2];
                    Ydds[i, k] = dds_nondilated[k][i, y/2];
                }
            WriteFloatArray2BMP(Y, p);
            WriteFloatArray2BMP(Ydds, "Y_midplane_DS.bmp");
            MakeComparisonPic(Y, Ydds, "Y_midplane_comparison.bmp");
        }

        private void MakeMidplaneZpic(float[][,] DoseSpace, string p)
        {
            float[,] Z = DoseSpace[DoseSpace.GetLength(0)/2];
            float[,] Zdds = dds_nondilated[DDS.GetLength(0)/2];
            WriteFloatArray2BMP(Z, p);
            WriteFloatArray2BMP(Zdds, "Z_midplane_DS.bmp");
            MakeComparisonPic(Z, Zdds, "Z_midplane_comparison.bmp");
        }

        private void WriteFloatArray2BMP(float[,] temp, string p)
        {
            string path = Path.Combine(ActiveDirectory, p);
            float max = Matrix.FindMax(temp);
            int color = 0;
            var b = new Bitmap(temp.GetLength(1), temp.GetLength(0));
            for (int j = 0; j < temp.GetLength(0); j++)
                for (int i = 0; i < temp.GetLength(1); i++)
                {
                    color = (int) ((temp[j, i]/max)*255);
                    b.SetPixel(i, j, Color.FromArgb(color, color, color));
                }
            b.Save(path);
        }

        private void MakeComparisonPic(float[,] temp, float[,] dds, string p)
        {
            string path = Path.Combine(ActiveDirectory, p);
            var temp2 = (float[,]) Matrix.Normalize(temp).Clone();
            int color = 0;
            var b = new Bitmap(temp.GetLength(0), temp.GetLength(1));
            for (int j = 0; j < temp.GetLength(1); j++)
                for (int i = 0; i < temp.GetLength(0); i++)
                {
                    float dose = temp[i, j];
                    float pixel = dds[i, j];
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

        private float[] LoadSliceDose(int which_slice, string subfolder)
        {
            string filename = String.Concat("slice_", which_slice);
            string path = Path.Combine(subfolder, filename);
            float[] d = ReadSliceDoseFromFile(path);
            Matrix.Normalize(ref d);
            return d;
        }
        

        /// <summary>
        ///   Takes in a slice 1D float matrix, and adds it to the global dosespace array for final calculation. Called by 
        ///   AssembleFinalDoseMatrix(). If dosespace is null, it will create the slice, else it will add
        ///   it to the existing one.
        /// </summary>
        /// <param name="slicedose"> </param>
        /// <param name="which_z_slice"> </param>
        public float[][,] WriteSliceDoseToDoseSpace(float[] slicedose, float[][,] output_ds, int which_slice)
        {
            int TranslateZBy = SlicePositions[which_slice] - (DCT/2); // <- NEED TO CHANGE APPROPRIATELY
            int StartAt = 0;
            int EndAt = DCT;
            if (TranslateZBy < 0)
            {
                StartAt += (-1)*(TranslateZBy);
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
                output_ds[current_z] = Matrix.Add(output_ds[current_z],
                                                  GrabSlice(slicedose, z, Volume[0].GetLength(0), Volume[0].GetLength(1)));
                int dssum = Matrix.SumAll(output_ds[current_z]);
                //WriteFloatArray2BMP(output_ds[current_z], "ds_z_" + z + "_" + dssum + ".bmp");
                //WriteFloatArray2BMP(GrabSlice(slicedose, z, output_ds[0].GetLength(0), output_ds[0].GetLength(1)), "sd_z_" + z + ".bmp");
            } //);
            return output_ds;
        }

        /// <summary>
        ///   Given a path, loads in the associated slicedose 1D float matrix. Called by AssembleFinalDoseMatrix()
        /// </summary>
        /// <param name="loadpath"> </param>
        /// <returns> </returns>
        public static float[] ReadSliceDoseFromFile(string loadpath)
        {
            float[] d;
            using (var fs = new FileStream(loadpath, FileMode.Open, FileAccess.Read))
            using (var br = new StreamReader(fs))
            {
                int x = Convert.ToInt16(br.ReadLine());
                int y = Convert.ToInt16(br.ReadLine());
                int z = Convert.ToInt16(br.ReadLine());
                d = new float[x*y*z];
                for (int i = 0; i < d.GetLength(0); i++)
                    d[i] = Convert.ToSingle(br.ReadLine());
            }
            float sum = d.Sum();
            return d;
        }

        /// <summary>
        ///   End of STEP 1. The actual method that is called by the UI, which takes in a folder directory path, and 
        ///   loads in each slicedose file and writes it to the global dosespace file.
        /// </summary>
        /// <param name="folderpath"> </param>
        public void AssembleDoseSpaceFromFiles()
        {
            string subfolder = ActiveDirectory;
            int progress = 0;
            int x = Volume[0].GetLength(0);
            int y = Volume[0].GetLength(1);
            int z = Volume.GetLength(0);
            /* First, a for-loop to create the string name for each slice, which then loads each
             * slicedose file and in turn adds it to the result matrix.*/
            DoseSpace = new float[z][,];
            for (int i = 0; i < z; i++)
                DoseSpace[i] = Matrix.Zeroes(x, y);
            for (int s = 0; s < NumSlices; s++)
            {
                float[] slicedose = LoadSliceDose(s, subfolder);
                DoseSpace = WriteSliceDoseToDoseSpace(slicedose, DoseSpace, s);
                progress = (int) (s*50.0/NumSlices);
                Ps2CalcDoseWorker.ReportProgress(progress);
            }

            //Write dosespace to file.
            WriteDoseSpaceToFile("OriginalDS.txt");
            Ps2CalcDoseWorker.ReportProgress(100);
        }

        /// <summary>
        ///   Called from AssembleDoseSpaceFromFiles() to write the output to a file. Note the header that 
        ///   is written first containing the list of dimension order.
        /// </summary>
        public void WriteDoseSpaceToFile(string filename)
        {
            string subfolder = ActiveDirectory;
            string path = Path.Combine(subfolder, filename);
            using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
            using (var bw = new StreamWriter(fs))
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
        ///   Counterpart to WriteDoseSpaceToFile(). Loads in a dosespace and sets it to the main variable
        ///   when given a path.
        /// </summary>
        /// <param name="path"> </param>
        public static float[] ReadDoseSpaceFromFile(string filename)
        {
            float[] ds;
            string subfolder = ActiveDirectory;
            string path = Path.Combine(subfolder, filename);
            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read))
            using (var br = new StreamReader(fs))
            {
                br.ReadLine();
                int x = Convert.ToInt16(br.ReadLine());
                int y = Convert.ToInt16(br.ReadLine());
                int z = Convert.ToInt16(br.ReadLine());
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

        #region All Background Worker methods

        private void AttachHandlers()
        {
            Ps1ShotOptimizeWorker = new BackgroundWorker();
            Ps1ShotOptimizeWorker.WorkerReportsProgress = true;
            Ps1ShotOptimizeWorker.RunWorkerCompleted += PS_1_ShotOptimize_RunWorkerCompleted;
            Ps1ShotOptimizeWorker.ProgressChanged += PS_1_ShotOptimize_ProgressChanged;
            Ps1ShotOptimizeWorker.DoWork += PS_1_ShotOptimize_DoWork;
            Ps3SliceWeightOptWorker = new BackgroundWorker();
            Ps3SliceWeightOptWorker.WorkerReportsProgress = true;
            Ps3SliceWeightOptWorker.RunWorkerCompleted += PS_3_SliceOptimize_worker_RunWorkerCompleted;
            Ps3SliceWeightOptWorker.ProgressChanged += PS_3_SliceOptimize_worker_ProgressChanged;
            Ps3SliceWeightOptWorker.DoWork += PS_3_SliceOptimize_worker_DoWork;
            Ps2CalcDoseWorker = new BackgroundWorker();
            Ps2CalcDoseWorker.WorkerReportsProgress = true;
            Ps2CalcDoseWorker.RunWorkerCompleted += PS_2_CalcDose_worker_RunWorkerCompleted;
            Ps2CalcDoseWorker.ProgressChanged += PS_2_CalcDose_worker_ProgressChanged;
            Ps2CalcDoseWorker.DoWork += PS_2_CalcDose_worker_DoWork;
        }

        #region PS_1 Shot Optimization Background worker methods

        private void PS_1_ShotOptimize_DoWork(object sender, DoWorkEventArgs e)
        {
            double count = 0;
            for (int i = 0; i < NumSlices; i++)
            {
                //((RasterPath)RasterPaths[i]).RPworker.RunWorkerAsync();
                Debug.WriteLine("SLICE " + Convert.ToString(i) + ":");
                ((RasterPath) RasterPaths[i]).OptimizeShotWeights();
                count++;
                Ps1ShotOptimizeWorker.ReportProgress((int) (100*count/NumSlices));
            }
        }

        private void PS_1_ShotOptimize_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (PathsetWorkerProgressChanged != null)
                PathsetWorkerProgressChanged.Invoke(null, e);
        }

        private void PS_1_ShotOptimize_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            ShotsWeighted = true;
            if (PathsetWorkerCompleted != null)
                PathsetWorkerCompleted.Invoke(null, e);
        }

        #endregion

        #region PS_2 CalcDose BackgroundWorker methods

        private void PS_2_CalcDose_worker_DoWork(object sender, DoWorkEventArgs e)
        {
            CalculateAndWriteSliceDoses(DK, folderpath);
            AssembleDoseSpaceFromFiles();
            //PS_2_CalcDose_worker.ReportProgress(50);
        }

        private void PS_2_CalcDose_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (OptimizationWorkerProgress != null)
                OptimizationWorkerProgress.Invoke(null, e);
        }

        private void PS_2_CalcDose_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            Debug.WriteLine("Shot weighting optimization complete. Review slice coverage before continuing...");
            if (OptimizationWorkerCompleted != null)
                OptimizationWorkerCompleted.Invoke(null, e);

            Ps3SliceWeightOptWorker.RunWorkerAsync();
        }

        #endregion

        #region PS_3 SliceOptimize BackgroundWorker methods

        private void PS_3_SliceOptimize_worker_DoWork(object sender, DoWorkEventArgs e)
        {
            //PathsetWorkerProgressChanged += Opt_PathSet_PathsetWorkerProgressChanged;
            OptimizeSliceWeights();
        }


        private void PS_3_SliceOptimize_worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (SliceweightWorkerProgress != null)
                SliceweightWorkerProgress.Invoke(null, e);
        }

        private void PS_3_SliceOptimize_worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (SliceweightWorkerCompleted != null)
                SliceweightWorkerCompleted.Invoke(null, e);
        }


        private void PathSet_SliceWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (RasterPathWorkerProgress != null)
                RasterPathWorkerProgress.Invoke(null, e);
        }

        #endregion

        #endregion

        #endregion
    }

    public class SuggThick
    {
        public int thickness { get; set; }
        public string Name { get; set; }
        public SuggThick(int thick)
        {
            thickness = thick;
            Name = "Suggested Thickness: " + thickness;
        }
    }
}