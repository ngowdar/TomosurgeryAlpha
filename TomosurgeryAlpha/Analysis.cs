﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.IO;

namespace TomosurgeryAlpha
{
    public static class Analysis
    {
        public static List<AnalysisInfo> AIList = new List<AnalysisInfo>();
        public static float[] DVH;
        public static string SessionName;
        private static double RX;
        private static double TolDose;
        private static double tumorvolume;
        private static double isovol;
        private static double isotumorvol;
        private static double csvol;
        private static double csrxvol;
        private static double maskvolume;
        private static double maskcoveragebyrx;
        private static double rtog_cov;
        public static double maxdose;
        public static double doserate;
        public static double beamontime;
        public static DICOMDoseFile ddf;

        public static void CreateReportFile(string name)
        {
            string datetimeString = string.Format("{0:yyyy-MM-dd_hh-mm-ss-tt}.txt", DateTime.Now);
            string path = string.Concat(SessionName, "_", name);
            path = Path.Combine(PathSet.ActiveDirectory, path);
            var fs = new FileStream(path, FileMode.Create);
            var tw = new StreamWriter(fs);
            tw.WriteLine("====================" + SessionName + " Report" + "====================");
            tw.Close();
            fs.Close();
        }

        public static string CreateFinalReportFile()
        {
            string datetimeString = string.Format("{0:yyyy-MM-dd_hh-mm-ss-tt}.txt", DateTime.Now);
            string path = string.Concat(SessionName, "_SummaryReport.txt");
            path = Path.Combine(PathSet.ActiveDirectory, path);
            var fs = new FileStream(path, FileMode.Create);
            var tw = new StreamWriter(fs);
            tw.WriteLine("====================" + SessionName + " Summary Report" + "====================");

            tw.WriteLine();

            tw.WriteLine("================PARAMETERS===================");
            tw.WriteLine("Raster Width: " + RasterPath.RasterWidth);
            tw.WriteLine("Shot Step Size: " + RasterPath.StepSize);
            tw.WriteLine("Slice Thickness: " + PathSet.ST);
            tw.WriteLine("Dose Calculation Thickness: " + PathSet.DCT);
            tw.WriteLine("Prescription dose: " + PathSet.RxDose);
            tw.WriteLine("Tolerance dose: " + PathSet.ToleranceDose);

            tw.WriteLine("================ANALYSIS===================");
            tw.WriteLine("Lesion Size: " + tumorvolume);
            tw.WriteLine("Rx Volume: " + isovol);
            tw.WriteLine("Tumor Volume covered by Rx: " + isotumorvol);
            tw.WriteLine("Shell Volume: " + maskvolume);
            tw.WriteLine("Percent volume covered by shell: " + (maskcoveragebyrx / maskvolume));
            tw.WriteLine(" ");
            tw.WriteLine("----------------Indexes------------------");
            tw.WriteLine("COVERAGE: " + Math.Round((isotumorvol/tumorvolume),2));
            tw.WriteLine("RTOG-CONFORMITY: " + Math.Round((isovol/tumorvolume),2));
            tw.WriteLine("LOMAX-SCHEIB: " + Math.Round((isotumorvol/isovol),2));
            tw.WriteLine("VANT REIT: " + Math.Round(((isotumorvol * isotumorvol) / (tumorvolume * isovol)),2));
            tw.WriteLine("RTOG-COVERAGE: " + Math.Round(rtog_cov,2));
            tw.WriteLine(" ");
            tw.WriteLine("----------------Time Cost--------------------");
            tw.WriteLine("Assuming central doserate of " + doserate + ":");
            tw.WriteLine("" + PathSet.TotalNumShots + " shots at a max of " + maxdose + "Gy...");
            tw.WriteLine(" ==> approximately " + Math.Round(beamontime,2) + "mins");
            tw.WriteLine("===========================================");
            
            tw.Close();
            fs.Close();

            return path;
        }


        public static void AddLineToReport(string filename, string text)
        {
            string name = string.Concat(SessionName, "_", filename);
            string path = Path.Combine(PathSet.ActiveDirectory, name);
            if (File.Exists(path))
            {
                var fs = new FileStream(path, FileMode.Append);
                var tw = new StreamWriter(fs);
                tw.WriteLine(text);
                tw.Close();
                fs.Close();
            }
            else
            {
                CreateReportFile(filename);
                AddLineToReport(filename, text);
            }
        }

        public static double FindRTOGCoverage(float[][,] tumor, float[][,] ndose, float rxdose)
        {
            //float[][,] ndose = Matrix.Normalize(dose);
            double max = FindMaxDose(ndose);
            var rx = (float) (max*rxdose);
            float min = 10000;
            for (int k = 0; k < tumor.GetLength(0); k++)
                for (int j = 0; j < tumor[0].GetLength(1); j++)
                    for (int i = 0; i < tumor[0].GetLength(0); i++)
                    {
                        float tumorstate = tumor[k][i, j];
                        float d = ndose[k][i, j];
                        if (tumorstate > 0)
                        {
                            if (d < min)
                                min = d;
                        }
                    }
            Debug.WriteLine("Max dose: " + max + "; Min dose: " + min + "; Min/Rx = " + (min/rx));
            return min/rx;
        }

        public static void AnalyzeDICOMdose(DICOMDoseFile ddf, StructureSet SS)
        {
            float[][,] tumor;
            float[][,] cs;
            SS.GetCSOnly(out cs, StructureSet.originalTumor);
            SS.GetTumorOnly(out tumor, StructureSet.originalTumor);
            float[][,] dose = ddf.OriginalDose;
            int tumor_xlength = tumor[0].GetLength(0);
            int dose_xlength = dose[0].GetLength(0);
            double maxdose = FindMaxDose(dose);
            var rxdose = (float)(0.5 * maxdose);
            var toldose = (float)(0.3 * maxdose);
            dose = Matrix.Normalize(dose);
            double isovol = 0;
            var ai = new AnalysisInfo();
            if (tumor_xlength != dose_xlength)
            {
                Debug.WriteLine("Tumor file doesn't match! Analzying dose alone...");
                for (int k = 0; k < dose.GetLength(0); k++)
                    for (int j = 0; j < dose[0].GetLength(1); j++)
                        for (int i = 0; i < dose[0].GetLength(0); i++)
                        {
                            double d = dose[k][i, j];
                            if (d >= rxdose)
                                isovol++;
                        }
                ai.Rx_Volume = isovol;
            }

            else
            {
                int zmid1 = (dose.GetLength(0)/2);
                int zmid2 = tumor.GetLength(0)/2;

                //Matrix.WriteFloatArray2BMP(Matrix.ThresholdEq(dose[zmid1], 0.5f), "DS_halfslice.bmp");
                //Matrix.WriteFloatArray2BMP(tumor[zmid2], "Tumor_halfslice.bmp");

                double isotumorvol = 0;
                double toldosevol = 0;
                double tolcsvol = 0;
                double tumorvol = Matrix.SumAll(Matrix.Normalize(tumor));
                double csvol = Matrix.SumAll(Matrix.Normalize(cs));
                for (int k = 0; k < tumor.GetLength(0); k++)
                    for (int j = 0; j < tumor[0].GetLength(1); j++)
                        for (int i = 0; i < tumor[0].GetLength(0); i++)
                        {
                            double t = tumor[k][i, j];
                            double d = dose[k][i, j];
                            double c = cs[k][i, j];

                            if (d >= rxdose)
                            {
                                isovol++;
                                if (t > 0)
                                    isotumorvol++;
                            }
                            else if (d >= toldose)
                            {
                                toldosevol++;
                                if (c > 0)
                                    tolcsvol++;
                            }
                        }
                GetDVH(dose, tumor);
                
                Debug.WriteLine("DVH Volume covered by 50% dose: " + DVH[50] + "//" + DVH[1] + " = " + DVH[50]/DVH[1]);
                ai.CS_Volume = csvol;
                ai.CS_overdose = tolcsvol;
                ai.LesionVolume = tumorvol;
                ai.Rx_Volume = isovol;
                ai.RxLesion_Volume = isotumorvol;
                ai.Coverage = isotumorvol / tumorvol;
                ai.RTOG = (isovol / tumorvol);
                ai.LomaxScheib = isotumorvol / isovol;
                ai.VantReits = (isotumorvol * isotumorvol) / (tumorvol * isovol);
                ai.TestName = String.Concat("original_", DateTime.Now.ToShortTimeString());
                rtog_cov = FindRTOGCoverage(tumor, dose, rxdose);
                ai.RTOGCoverage = Math.Round(rtog_cov,2);
            }
            AIList.Add(ai);
        }

        private static void GetDVH(float[][,] origdose, float[][,] tumor)
        {
            float[][,] dose = Matrix.Normalize(origdose);
            DVH = new float[100];
            DVH = Matrix.Zero1DFloat(100);
            for (int k = 0; k < dose.GetLength(0); k++)
                for (int j = 0; j < dose[0].GetLength(1); j++)
                    for (int i = 0; i < dose[0].GetLength(0); i++)
                    {
                        float tumor_state = tumor[k][i, j];
                        float dose_value = dose[k][i, j];
                        if (tumor_state > 0)
                        {
                            var d = (int) (Math.Round(100.0*dose_value));
                            for (int x = 0; x < d; x++)
                            {
                                DVH[x] += 1;
                            }
                        }
                    }
        }

        private static double FindMaxDose(float[][,] dose)
        {
            double max = 0;
            for (int k = 0; k < dose.GetLength(0); k++)
                for (int j = 0; j < dose[0].GetLength(1); j++)
                    for (int i = 0; i < dose[0].GetLength(0); i++)
                    {
                        if (dose[k][i, j] > max)
                            max = dose[k][i, j];
                    }
            return max;
        }

        public static void RunAnalysis(PathSet PS, StructureSet SS, double rxlevel)
        {
            //PS.DoseSpace = Matrix.Normalize(PS.DoseSpace);
            var ai = new AnalysisInfo();
            double maxdose = FindMaxDose(PS.DoseSpace);
            var rxdose = (float) (rxlevel*maxdose);
            var toldose = (float) (0.3*maxdose);
            RX = rxdose;
            TolDose = PathSet.ToleranceDose*maxdose;
            int startingz = PS.SlicePositions[0] - (PathSet.DCT/2);
            //Debug
            //float[][,] DDS = PS.PrepareDDS(SS.fj_Tumor);
            //AnalyzeLesionCoverage(DDS,SS.fj_Tumor,startingz);

            AnalyzeLesionCoverage(PS.DoseSpace, SS.fj_Tumor, SS.fj_CS, startingz);
            AnalyzeMaskCoverage(PS.DoseSpace, Matrix.PrepareDDSMask(SS.fj_Tumor, 3), rxlevel);
            ai.RxLevel = RX;
            ai.LesionVolume = tumorvolume;
            ai.Rx_Volume = isovol;
            ai.RxLesion_Volume = isotumorvol;
            ai.Coverage = isotumorvol/tumorvolume;
            ai.RTOG = isovol/tumorvolume;
            ai.LomaxScheib = isotumorvol/isovol;
            rtog_cov = FindRTOGCoverage(SS.fj_Tumor, PS.DoseSpace, (float)rxlevel);
            ai.RTOGCoverage = rtog_cov;
            ai.VantReits = (isotumorvol*isotumorvol)/(tumorvolume*isotumorvol);
            ai.MaskCoverage = maskcoveragebyrx;
            ai.MaskVolume = maskvolume;
            ai.CS_Volume = csvol;
            ai.CS_overdose = csrxvol;
            ai.TestName = DateTime.Now.ToShortTimeString();
            ai.CreateInfoString();
            AIList.Add(ai);
            //if (ddf != null)
            //{
            //    AnalysisInfo DICOM_dose = new AnalysisInfo();
            //    DICOM_dose.LesionVolume = lesionvolume;
            //    DICOM_dose.RxLevel = RX;
            //    //AnalyzeLesionCoverage(ddf.GetJaggedDoseArray(), SS.fj_Tumor, startingz);                                
            //    DICOM_dose.Rx_Volume = isovol;
            //    DICOM_dose.RxLesion_Volume = isotumorvol;
            //    DICOM_dose.RTOG = isovol / lesionvolume;
            //    DICOM_dose.LomaxScheib = isotumorvol / isovol;
            //    DICOM_dose.VantReits = (isotumorvol * isotumorvol) / (lesionvolume * isovol);
            //    DICOM_dose.TestName = "DICOMdose";
            //    AIList.Add(DICOM_dose);
            //}
        }


        public static void RunAnalysis(DICOMDoseFile ddf, float[][,] tumor, float[][,] cs, double rxlevel)
        {
            //TODO: Figure out what size slab of tumor to take such that it lines up with the beginning of the dose matrix.

            /*Problem:
             * 
             * The structure tumor is by definition going to be a smaller matrix 
             * than the dose matrix, because the dose has to cover outside the tumor 
             * boundaries.
             *                   62                52        
             * So, dose start --> |================|--------|&&&&&&&&&&&... and so on
             *                    ^doseoffset      ^SS      ^where tumor starts
             * So, (dose-offset - structureoffset) gives how far the 
             * SS matrix is in DICOM coordinates, and zstart gives how far before tumor edge starts.
             *              * 
             *  So, (doseoffset - SSoffset)*4-3 gives where the SS matrix starts
             *  Then, that + the zstart gives where the tumor starts. 
             */
            int[] zends = StructureSet.FindZBoundaries(tumor);
            //The dose DICOM offset - the structureset's offset = where the structure set starts in the dosespace
            //Add the starting boundary for the z-slice to get the absolute start location for the tumor within the 
            //dosespace matrix.
            int startingz = (int) Math.Round(DICOMDoseFile.doseoffset[2] - StructureSet.f_SSoffset[2]) + zends[0];
            float[][,] t = PathSet.GrabSlab(tumor, zends[0], zends[1], true);
            t = PathSet.PrepareDDS(t, 1.0f);


            var DICOM_dose = new AnalysisInfo();
            AnalyzeLesionCoverage(ddf.EnlargeAndInterpolate(ddf.OriginalDose), t, cs, (int) ddf.ZStart);
            AnalyzeMaskCoverage(ddf.EnlargeAndInterpolate(ddf.OriginalDose), Matrix.PrepareDDSMask(t, 3),
                                rxlevel);
            DICOM_dose.RxLevel = RX;
            DICOM_dose.LesionVolume = tumorvolume;
            DICOM_dose.Rx_Volume = isovol;
            DICOM_dose.RxLesion_Volume = isotumorvol;
            DICOM_dose.RTOG = isovol/tumorvolume;
            DICOM_dose.Coverage = isotumorvol/tumorvolume;
            DICOM_dose.LomaxScheib = isotumorvol/isovol;
            DICOM_dose.VantReits = (isotumorvol*isotumorvol)/(tumorvolume*isovol);
            DICOM_dose.RTOGCoverage = FindRTOGCoverage(tumor, ddf.EnlargeAndInterpolate(ddf.OriginalDose),
                                                       (float) rxlevel);
            DICOM_dose.TestName = "DICOMdose";
            AIList.Add(DICOM_dose);
        }


        //l


        private static void AnalyzeLesionCoverage(float[][,] ds, float[][,] tumor, float[][,] cs, int startingz)
        {
            //ds = Matrix.Normalize(ds);

            double maxdose = Matrix.FindMax(ds);
            double rxdose = RX;
            //tumorvolume = Matrix.SumAll(Matrix.Normalize(tumor));
            isovol = 0;
            isotumorvol = 0;
            csvol = 0;
            double toldosevol = 0;
            int zmid1 = (ds.GetLength(0)/2);
            int zmid2 = tumor.GetLength(0)/2;
            Matrix.WriteFloatArray2BMP(Matrix.ThresholdEq(ds[zmid1], 0.5f), "DS_halfslice.bmp");
            Matrix.WriteFloatArray2BMP(tumor[zmid2], "Tumor_halfslice.bmp");
            for (int k = 0; k < tumor.GetLength(0); k++)
                for (int j = 0; j < tumor[0].GetLength(1); j++)
                    for (int i = 0; i < tumor[0].GetLength(0); i++)
                    {
                        double t = tumor[k][i, j];
                        double d = ds[k][i, j];
                        double c = cs[k][i, j];
                        if (c > 0)
                            csvol++;

                        if (d >= rxdose)
                        {
                            isovol++;
                            if (t > TolDose)
                            {
                                isotumorvol++;
                                tumorvolume++;
                            }
                            if (c > 0)
                                csrxvol++;
                        }
                        else if (d < rxdose && t > TolDose)
                        {
                            toldosevol++;
                            tumorvolume++;
                            if (c > 0)
                                csrxvol++;
                        }
                    }
            GetDVH(ds, tumor);
            Debug.WriteLine("DVH Volume covered by 50% dose: " + DVH[50] + "//" + DVH[1] + " = " + DVH[50]/DVH[1]);
        }

        public static void AnalyzeMaskCoverage(float[][,] ds, float[][,] mask, double rxlevel)
        {
            //ds = Matrix.Normalize(ds);
            maskvolume = 0;
            maskcoveragebyrx = 0;
            int zmid2 = mask.GetLength(0)/2;
            Matrix.WriteFloatArray2BMP(mask[zmid2], "Tumor_halfslice.bmp");
            for (int k = 0; k < ds.GetLength(0); k++)
            {
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        if (mask[k][i, j] > 0)
                        {
                            maskvolume++;
                            if (ds[k][i, j] >= RX)
                                maskcoveragebyrx++;
                        }
                    }
            }
        }


        internal static void AddAnalysisReport(string filename)
        {
            AddLineToReport(filename, "================ANALYSIS===================");
            AddLineToReport(filename, "Lesion Size: " + tumorvolume);
            AddLineToReport(filename, "Rx Volume: " + isovol);
            AddLineToReport(filename, "Tumor Volume covered by Rx: " + isotumorvol);
            AddLineToReport(filename, "Shell Volume: " + maskvolume);
            AddLineToReport(filename, "Percent volume covered by shell: " + (maskcoveragebyrx/maskvolume));
            AddLineToReport(filename, "===========================================");
        }
    }

    public struct AnalysisInfo
    {
        public string TestName { get; set; }
        public double RxLevel { get; set; }
        public double Coverage { get; set; }
        public double RTOGCoverage { get; set; }
        public double RTOG { get; set; }
        public double LomaxScheib { get; set; }
        public double VantReits { get; set; }
        public double LesionVolume { get; set; }
        public double Rx_Volume { get; set; }
        public double RxLesion_Volume { get; set; }
        public double MaskVolume { get; set; }
        public double MaskCoverage { get; set; }
        public double CS_Volume { get; set; }
        public double CS_overdose { get; set; }

        internal void CreateInfoString()
        {
            string s = TestName + ": ";
            s += "\n Lesion Size: " + LesionVolume;
            s += "\n Isovolume Size: " + Rx_Volume;
            s += "\n Lesion Covered by Isodose: " + RxLesion_Volume;
            if (CS_Volume > 0)
                s += "\n Percent CS covered by isodose: " + (CS_overdose/CS_Volume);
        }
    }
}