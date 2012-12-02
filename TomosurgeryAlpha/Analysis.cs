﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace TomosurgeryAlpha
{
    public static class Analysis
    {
        public static List<AnalysisInfo> AIList = new List<AnalysisInfo>(); 
        public static float[] DVH;
        static double RX;
        static double TolDose;
        static double lesionvolume;
        static double totalvolcoveredbyrx;
        static double lesioncoveragebyrx;
        static double maskvolume;
        static double maskcoveragebyrx;
        public static DICOMDoseFile ddf;

        public static void CreateReportFile(string name)
        {
            string path = "report.txt";
            path = System.IO.Path.Combine(PathSet.ActiveDirectory, path);
            FileStream fs = new FileStream(path, FileMode.Create);
            StreamWriter tw = new StreamWriter(fs);
            tw.WriteLine("====================" + name + " Report" + "====================");
            tw.Close();
            fs.Close();
        }

        public static void AddLineToReport(string s)
        {
            string path = "report.txt";
            path = System.IO.Path.Combine(PathSet.ActiveDirectory, path);
            bool exists = System.IO.File.Exists(path);
            if (System.IO.File.Exists(path))
            {
                FileStream fs = new FileStream(path, FileMode.Append);
                StreamWriter tw = new StreamWriter(fs);
                tw.WriteLine(s);
                tw.Close();
                fs.Close();
            }
            else
            {
                CreateReportFile("Unnamed");
                AddLineToReport(s);
            }
        }

        public static void AnalyzeDICOMdose(DICOMDoseFile ddf, StructureSet SS)
        {
            
            float[][,] tumor;            
            float[][,] cs;
            SS.GetCSOnly(out cs, StructureSet.originalTumor);
            SS.GetTumorOnly(out tumor, StructureSet.originalTumor);
            float[][,] dose = Matrix.TransposeMatrix(DICOMDoseFile.OriginalDose);
            double maxdose = FindMaxDose(dose);
            float rxdose = (float)(0.5 * maxdose); float toldose = (float)(0.3 * maxdose);
            //dose = Matrix.Normalize(dose);            

            double isovol = 0; double isotumorvol = 0;
            double toldosevol = 0; double tolcsvol = 0;
            double tumorvol = Matrix.SumAll(Matrix.Normalize(tumor));
            double csvol = Matrix.SumAll(Matrix.Normalize(cs));
            for (int k = 0; k < tumor.GetLength(0); k++)
                for (int j = 0; j < tumor[0].GetLength(1); j++)
                    for (int i = 0; i < tumor[0].GetLength(0); i++)
                    {
                        double t = tumor[k][i,j];
                        double d = dose[k][i, j]*maxdose;
                        double c = cs[k][i,j];

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
            System.Diagnostics.Debug.WriteLine("DVH Volume covered by 50% dose: " + DVH[50] + "//" + DVH[1] + " = " + DVH[50]/DVH[1]);


            AnalysisInfo ai = new AnalysisInfo();
            ai.CS_Volume = csvol;
            ai.CS_overdose = tolcsvol;
            ai.LesionVolume = tumorvol;
            ai.Rx_Volume = isovol;
            ai.RxLesion_Volume = isotumorvol;
            ai.RTOG = (double)(isovol / tumorvol);
            ai.LomaxScheib = isotumorvol / isovol;
            ai.VantReits = (isotumorvol * isotumorvol) / (tumorvol * isovol);
            ai.TestName = String.Concat("original_", System.DateTime.Now.ToShortTimeString());
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
                            int d = (int)(Math.Round(100.0 * dose_value));
                            for (int x = 0; x < d; x++)
                            {
                                DVH[x] += 1;
                            }
                        }
                    }
        }

        private static float GetVolumeCoveredByIso(float iso)
        {
            float totalvol = 0; float isovol = 0;
            int d = (int)Math.Round(100*iso,0);
            for (int i = 0; i < 100; i++)
            {
                if (DVH[i] < d)
                {
                    totalvol++;
                    isovol++;
                }
                else
                    totalvol++;
            }
            return (isovol / totalvol);

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
            AnalysisInfo ai = new AnalysisInfo();
            RX = rxlevel;
            TolDose = PathSet.ToleranceDose;
            int startingz = PS.SlicePositions[0] - (PathSet.DCT / 2);
            //Debug
            //float[][,] DDS = PS.PrepareDDS(SS.fj_Tumor);
            //AnalyzeLesionCoverage(DDS,SS.fj_Tumor,startingz);

            AnalyzeLesionCoverage(PS.DoseSpace, SS.fj_Tumor, startingz);
            AnalyzeMaskCoverage(PS.DoseSpace, Matrix.PrepareDDSMask(SS.fj_Tumor,3), rxlevel);
            ai.RxLevel = RX;
            ai.LesionVolume = lesionvolume;
            ai.Rx_Volume = totalvolcoveredbyrx;
            ai.RxLesion_Volume = lesioncoveragebyrx;
            ai.RTOG = totalvolcoveredbyrx / lesionvolume;
            ai.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
            ai.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
            ai.TestName = System.DateTime.Now.ToShortTimeString();
            ai.CreateInfoString();
            AIList.Add(ai);
            //if (ddf != null)
            //{
            //    AnalysisInfo DICOM_dose = new AnalysisInfo();
            //    DICOM_dose.LesionVolume = lesionvolume;
            //    DICOM_dose.RxLevel = RX;
            //    //AnalyzeLesionCoverage(ddf.GetJaggedDoseArray(), SS.fj_Tumor, startingz);                                
            //    DICOM_dose.Rx_Volume = totalvolcoveredbyrx;
            //    DICOM_dose.RxLesion_Volume = lesioncoveragebyrx;
            //    DICOM_dose.RTOG = totalvolcoveredbyrx / lesionvolume;
            //    DICOM_dose.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
            //    DICOM_dose.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
            //    DICOM_dose.TestName = "DICOMdose";
            //    AIList.Add(DICOM_dose);
            //}
        }


        public static void RunAnalysis(float[][,] dosespace, float[][,] tumor, double rxlevel)
        {
            double rxdoselevel;
            double rtogindex;
            double lomaxscheib;
            double vantreits;
            int startingz = 20;
            //dosespace = Matrix.Normalize(dosespace);
            AnalysisInfo ai = new AnalysisInfo();
            RX = rxlevel;            
            AnalyzeLesionCoverage(dosespace, tumor, startingz);
            AnalyzeMaskCoverage(dosespace, Matrix.PrepareDDSMask(tumor, 3), rxlevel);
            ai.RxLevel = RX;
            ai.LesionVolume = lesionvolume;
            ai.Rx_Volume = totalvolcoveredbyrx;
            ai.RxLesion_Volume = lesioncoveragebyrx;
            ai.RTOG = totalvolcoveredbyrx / lesionvolume;
            ai.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
            ai.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
            ai.TestName = System.DateTime.Now.ToShortTimeString();
            AIList.Add(ai);
            if (ddf != null)
            {
                AnalysisInfo DICOM_dose = new AnalysisInfo();
                AnalyzeLesionCoverage(ddf.EnlargeAndInterpolate(DICOMDoseFile.OriginalDose), tumor, startingz);
                DICOM_dose.RxLevel = RX;
                DICOM_dose.LesionVolume = lesionvolume;
                DICOM_dose.Rx_Volume = totalvolcoveredbyrx;
                DICOM_dose.RxLesion_Volume = lesioncoveragebyrx;
                DICOM_dose.RTOG = totalvolcoveredbyrx / lesionvolume;
                DICOM_dose.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
                DICOM_dose.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
                DICOM_dose.TestName = "DICOMdose";
                AIList.Add(DICOM_dose);
            }
        }

        

        public static void RunAnalysis(DICOMDoseFile ddf, float[][,] tumor, double rxlevel)
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
            int startingz = (int)Math.Round(DICOMDoseFile.doseoffset[2] - StructureSet.f_SSoffset[2]) + zends[0];
            float[][,] t = PathSet.GrabSlab(tumor, zends[0], zends[1], true);            
            t = PathSet.PrepareDDS(t);

            
            AnalysisInfo DICOM_dose = new AnalysisInfo();
            AnalyzeLesionCoverage(ddf.EnlargeAndInterpolate(DICOMDoseFile.OriginalDose), t, (int)ddf.ZStart);
            AnalyzeMaskCoverage(ddf.EnlargeAndInterpolate(DICOMDoseFile.OriginalDose), Matrix.PrepareDDSMask(t, 3), rxlevel);
            DICOM_dose.RxLevel = RX;
            DICOM_dose.LesionVolume = lesionvolume;
            DICOM_dose.Rx_Volume = totalvolcoveredbyrx;
            DICOM_dose.RxLesion_Volume = lesioncoveragebyrx;
            DICOM_dose.RTOG = totalvolcoveredbyrx / lesionvolume;
            DICOM_dose.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
            DICOM_dose.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
            DICOM_dose.TestName = "DICOMdose";
            AIList.Add(DICOM_dose);
        }
        

        //l

        private static void FindTumorVolume(PathSet PS)
        {
            //UNDONE: A menu-based option to show the volume of the current structure set.
            //Make sure it uses a binary matrix, then multiply it by the standard reference dose
            //using element-multiply, and THEN count it. This will the most realistic.

        }
        

        private static void AnalyzeLesionCoverage(float[][,] ds, float[][,] tumor, int startingz)
        {
            //ds = Matrix.Normalize(ds);
            lesionvolume = 0; totalvolcoveredbyrx = 0; lesioncoveragebyrx = 0;
            int zmid1 = (ds.GetLength(0) / 2); int zmid2 = tumor.GetLength(0) / 2;
            Matrix.WriteFloatArray2BMP(ds[zmid1], "DS_halfslice.bmp");
            Matrix.WriteFloatArray2BMP(tumor[zmid2], "Tumor_halfslice.bmp");
            for (int k = 0; k < ds.GetLength(0); k++)
            {
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        if (ds[k][i, j] >= RX)
                        {
                            totalvolcoveredbyrx++;
                            if (tumor[k][i, j] > TolDose)
                            {
                                lesioncoveragebyrx++;
                                lesionvolume++;
                            }
                        }
                        else if (ds[k][i, j] < RX && tumor[k][i, j] > TolDose)
                            lesionvolume++;                        
                    }
            }
            
        }
        
        public static void AnalyzeMaskCoverage(float[][,] ds, float[][,] mask, double rxlevel)
        {
            //ds = Matrix.Normalize(ds);
            maskvolume = 0; maskcoveragebyrx = 0;
            int zmid2 = mask.GetLength(0) / 2;            
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




        internal static void AddAnalysisReport()
        {            
            AddLineToReport("================ANALYSIS===================");
            AddLineToReport("Lesion Size: " + lesionvolume);
            AddLineToReport("Rx Volume: " + totalvolcoveredbyrx);
            AddLineToReport("Tumor Volume covered by Rx: " + lesioncoveragebyrx);
            AddLineToReport("Shell Volume: " + maskvolume);
            AddLineToReport("Percent volume covered by shell: " + (maskcoveragebyrx / maskvolume));
            AddLineToReport("===========================================");
        }
    }

    public struct AnalysisInfo
    {
        public string TestName { get; set; }
        public double RxLevel { get; set; }
        public double LesionVolume { get; set; }
        public double Rx_Volume { get; set; }
        public double RxLesion_Volume { get; set; }
        public double RTOG { get; set; }
        public double LomaxScheib { get; set; }
        public double VantReits { get; set; }
        public double CS_Volume { get; set; }
        public double CS_overdose { get; set; }

        internal void CreateInfoString()
        {
            string s = TestName + ": ";
            s += "\n Lesion Size: " + LesionVolume;
            s += "\n Isovolume Size: " + Rx_Volume;
            s += "\n Lesion Covered by Isodose: " + RxLesion_Volume;
            if (CS_Volume > 0)
                s += "\n Percent CS covered by isodose: " + (CS_overdose / CS_Volume);

        }
    }
}
