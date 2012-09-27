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
        static double lesionvolume;
        static double totalvolcoveredbyrx;
        static double lesioncoveragebyrx;

        public static void RunAnalysis(PathSet PS, StructureSet SS, double rxlevel)
        {
            double rxdoselevel;            
            double rtogindex;
            double lomaxscheib;
            double vantreits;
            PS.DoseSpace = Matrix.Normalize(PS.DoseSpace);
            AnalysisInfo ai = new AnalysisInfo();
            RX = rxlevel;
            int startingz = PS.SlicePositions[0] - (PS.DoseCalculationThickness / 2);
            AnalyzeLesionCoverage(PS.DoseSpace,SS.fj_Tumor,startingz);
            ai.RxLevel = RX;
            ai.LesionVolume = lesionvolume;
            ai.Rx_Volume = totalvolcoveredbyrx;
            ai.RxLesion_Volume = lesioncoveragebyrx;
            ai.RTOG = totalvolcoveredbyrx / lesionvolume;
            ai.LomaxScheib = lesioncoveragebyrx / totalvolcoveredbyrx;
            ai.VantReits = (lesioncoveragebyrx * lesioncoveragebyrx) / (lesionvolume * totalvolcoveredbyrx);
            ai.TestName = System.DateTime.Now.ToShortTimeString();
            AIList.Add(ai);
        }

        

        

        public static void RunAnalysis(float[][,] dosespace, float[][,] tumor, double rxlevel)
        {
            double rxdoselevel;
            double rtogindex;
            double lomaxscheib;
            double vantreits;
            int startingz = 20;
            dosespace = Matrix.Normalize(dosespace);
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
        }

        private static void AnalyzeLesionCoverage(float[][,] ds, float[][,] tumor, int startingz)
        {
            
            lesionvolume = 0; totalvolcoveredbyrx = 0; lesioncoveragebyrx = 0;
            
            for (int k = 0; k < ds.GetLength(0); k++)
            {
                for (int j = 0; j < ds[0].GetLength(1); j++)
                    for (int i = 0; i < ds[0].GetLength(0); i++)
                    {
                        if (ds[k][i, j] >= RX)
                               totalvolcoveredbyrx++;
                        if (tumor[k+startingz][i, j] > 0)                                                            
                                lesionvolume++;
                        if (tumor[k + startingz][i, j] > 0 && ds[k][i, j] >= RX)
                            lesioncoveragebyrx++;
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
