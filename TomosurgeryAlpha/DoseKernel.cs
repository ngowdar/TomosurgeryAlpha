﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Collections;
using System.Runtime.Serialization.Formatters.Binary;

namespace TomosurgeryAlpha
{
    public class DoseKernel
    {
        float[][] dose;
        public static int N;
        int[] sectorsizes;
        double resolution;
        double[] location;
        public float[,] midplane;
        public DoseKernelInfo DKI;

        public DoseKernel(string path)
        {
            string hpath = GetHeaderFileName(path);
            FileInfo fi = new FileInfo(path);
            
            DKI = new DoseKernelInfo();
            if (hpath != null && fi.Exists)
            {
                LoadHeader(hpath);
                LoadDose(path);
                FindMidplane();
                DKI.Name = fi.Name;
                DKI.GetListboxInfo();
                SetDosemidplaneForOpt();
            }            
        }

        private void FindMidplane()
        {
            float[] temp = dose[N / 2];
            midplane = new float[N,N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    midplane[i, j] = temp[j * N + i];            
        }

        private string GetHeaderFileName(string path)
        {
            if (path.EndsWith("txt") == true)
            {
                path = Path.ChangeExtension(path, ".h");
                //path.Remove(count - 4,3);
                //String.Concat(path, "h");
                //path.Replace("txt", "h");
                //path.TrimEnd(("txt").ToCharArray());
                return path;
            }
            else
            {
                System.Windows.MessageBox.Show("You didn't choose a text file!");
                path = null;
            }
            return path;

        }

        private static ArrayList BinaryDeserialize(string path)
        {
            ArrayList output = null;
            using (FileStream str = File.OpenRead(path))
            {
                BinaryFormatter bf = new BinaryFormatter();
                output = (ArrayList)bf.Deserialize(str);
            }
            return output;
        }
        
        private void LoadDose(string path)
        {
            StreamReader t = new StreamReader(path);
            dose = new float[N][];
            float[] temp;
            for (int j = 0; j < N; j++)
            {
                temp = new float[N * N];
                for (int i = 0; i < N * N; i++)
                {
                    temp[i] = (float)Convert.ToDecimal(t.ReadLine());
                    //dose[j][i] = (float)Convert.ToDecimal(t.ReadLine());
                }
                dose[j] = temp;
            }
        }

        private void LoadHeader(string hpath)
        {
            StreamReader head = new StreamReader(hpath);
            head.ReadLine();
            string sectorconfig = head.ReadLine();
            DKI.SectorConfig = Convert.ToInt32(sectorconfig);
            SetSectors(sectorconfig);
            head.ReadLine();
            resolution = Convert.ToDouble(head.ReadLine());
            DKI.Resolution = resolution;
            head.ReadLine();
            //resolution = Convert.ToDouble(res.Remove(0, 12));
            //head.ReadLine();
            location = new double[3];
            location[0] = Convert.ToDouble(head.ReadLine());
            location[1] = Convert.ToDouble(head.ReadLine());
            location[2] = Convert.ToDouble(head.ReadLine());
            DKI.ShotLocation = "" + location[0] + ", " + location[1] + ", " + location[2];
            N = Convert.ToInt16(head.ReadLine());
            DKI.Size = N;
        }

        private void SetSectors(string sectorconfig)
        {
            sectorsizes = new int[8];
            for (int i = 0; i < 8; i++)
            {
                sectorsizes[i] = Convert.ToInt16(sectorconfig[i]);
            }
        }

        public float[] GetSlice(int z)
        {
            return dose[z];
        }

        public float[,] Get2DSlice(int z)
        {
            float[,] r = new float[N, N];
            float[] temp = dose[z];

            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                    r[i, j] = temp[i * j];
            return r;
        }

        private void SetDosemidplaneForOpt()
        {
            float[,] dmp = new float[N, N];
            for (int j = 0; j < N; j++)
                for (int i = 0; i < N; i++)
                {
                    float z = 0;
                    for (int k = 0; k < N; k++)
                        z += dose[k][j * N + i];
                    dmp[i, j] = z;
                }
            RasterPath.dosemidplane = dmp;
        }
    }

    public struct DoseKernelInfo
    {
        public string Name { get; set; }
        public int Size { get; set; }
        public double Resolution { get; set; }
        public int SectorConfig { get; set; }
        public string ShotLocation { get; set; }
        public string Info { get; set; }

        public void GetListboxInfo()
        {
            Info = "DOSE: " + Name + ", " + Size + "x" + Size + ", Res: " + Resolution + "mm/px, Config: " + SectorConfig;            
        }
    }
}