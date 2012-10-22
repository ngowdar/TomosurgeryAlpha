using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TomosurgeryAlpha
{
    public static class Analysis
    {
        public static List<AnalysisInfo> AIList = new List<AnalysisInfo>();        
        static double RX;
        static double TolDose;
        static double lesionvolume;
        static double totalvolcoveredbyrx;
        static double lesioncoveragebyrx;
        public static DICOMDoseFile ddf;

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
                DICOM_dose.LesionVolume = lesionvolume;
                DICOM_dose.RxLevel = RX;
                AnalyzeLesionCoverage(ddf.GetJaggedDoseArray(), SS.fj_Tumor, startingz);                                
                DICOM_dose.Rx_Volume = totalvolcoveredbyrx;
                DICOM_dose.RxLesion_Volume = lesioncoveragebyrx;
                DICOM_dose.RTOG = totalvolcoveredbyrx / lesionvolume;
                DICOM_dose.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
                DICOM_dose.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
                DICOM_dose.TestName = "DICOMdose";
                AIList.Add(DICOM_dose);
            }
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
                AnalyzeLesionCoverage(ddf.GetJaggedDoseArray(), tumor, startingz);
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
            double rxdoselevel;
            double rtogindex;
            double lomaxscheib;
            double vantreits;
            AnalysisInfo DICOM_dose = new AnalysisInfo();
            AnalyzeLesionCoverage(ddf.GetJaggedDoseArray(), tumor, (int)ddf.ZStart);
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
            ds = Matrix.Normalize(ds);
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
        
    }
}
