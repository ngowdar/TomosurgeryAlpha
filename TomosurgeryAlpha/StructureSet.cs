using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using openDicom;
using System.IO;
using System.Collections.ObjectModel;
using System.Collections;
using System.Windows.Media.Imaging;
using System.ComponentModel;

namespace TomosurgeryAlpha
{
    public class StructureSet
    {
        public event ProgressChangedEventHandler StructureProgress;
        public event RunWorkerCompletedEventHandler StructureDone;
        public BackgroundWorker Structure_bw;
        public static string s_dictionarypath;
        static openDicom.Registry.DataElementDictionary dd;
        public openDicom.File.DicomFile difile;
        
        public static float f_global_xoffset;
        public static float f_global_yoffset;
        public static float f_global_zoffset;

        public float[] f_offset;
        public float[][] f_structurearray;
        public static int size;
        public ArrayList pointdata;
        public ArrayList pointF_data;
        public ArrayList shiftedpointdata;
        public ArrayList structure_names;
        public StructureInfo SI;
        public float[][,] fj_Tumor;
        public float[][,] fj_CS;

        public StructureSet(string header, string tumor)
        {
            float[][,] BinaryVolume = ReadArrayFromFile(tumor, header);
            BinaryVolume = EnlargeTumor(BinaryVolume, 160);
            f_structurearray = CreateArray(BinaryVolume);
            fj_Tumor = GetTumorOnly(BinaryVolume);
            fj_CS = GetCSOnly(BinaryVolume);
            size = f_structurearray.GetLength(0);
        }

        public StructureSet(float[][] tumorobj)
        {            
            f_structurearray = tumorobj;
            //fj_Tumor = tumorobj;
            SI = new StructureInfo();

        }

        private float[][] CreateArray(float[][,] BinaryVolume)
        {
            float[][] temp = new float[BinaryVolume.GetLength(0)][];
            LinMatrix bv = new LinMatrix(BinaryVolume);
            for (int i = 0; i < BinaryVolume.GetLength(0); i++)
                temp[i] = bv.GrabLinearSlice(i);
            return temp;
        }

        private float[][,] GetCSOnly(float[][,] d)
        {
            float[][,] cs = d;
            for (int k = 0; k < d.GetLength(0); k++)
                for (int i = 0; i < d[0].GetLength(0); i++)
                    for (int j = 0; j < d[0].GetLength(1); j++)
                    {
                        if (d[k][i, j] > 0 && d[k][i,j] < 5)
                            cs[k][i, j] = 1;
                        else
                            cs[k][i, j] = 0;
                    }
            return cs;
        }

        private float[][,] GetTumorOnly(float[][,] d)
        {
            float[][,] tumor = d;
            for (int k = 0; k < d.GetLength(0); k++)
                for (int i = 0; i < d[0].GetLength(0); i++)
                    for (int j = 0; j < d[0].GetLength(1); j++)
                    {
                        if (d[k][i, j] > 10)
                            tumor[k][i, j] = 1;
                        else
                            tumor[k][i, j] = 0;
                    }
            return tumor;
        }

        private float[][,] EnlargeTumor(float[][,] dd, int p)
        {   
            float[][,] d = Matrix.EnlargeAndCenter(dd, p, p / 2, p / 2, p / 2);
            return d;
        }

        
        public float[][,] ReadArrayFromFile(string fpath, string hpath)
        {
            //Header file is stored as a *.h file. Tumor is just text.
            float[][,] d; int x; int y; int z;
            using (System.IO.StreamReader header = new System.IO.StreamReader(hpath))
            {
                //Read in the size of the tumor
                y = Convert.ToInt16(header.ReadLine());
                x = Convert.ToInt16(header.ReadLine());
                z = Convert.ToInt16(header.ReadLine());

                //Read in the global img_offset vector (this is the top-left corner, and sets the DICOM coordinate frame)
                StructureSet.f_global_xoffset = (float)Convert.ToDecimal(header.ReadLine());
                StructureSet.f_global_yoffset = (float)Convert.ToDecimal(header.ReadLine());
                StructureSet.f_global_zoffset = (float)Convert.ToDecimal(header.ReadLine());

                //Read in the dose_offset vector (this is the top-left corner of the cropped window showing the tumor)
                f_offset = new float[3];
                f_offset[0] = (float)Convert.ToDouble(header.ReadLine());
                f_offset[1] = (float)Convert.ToDouble(header.ReadLine());
                f_offset[2] = (float)Convert.ToDouble(header.ReadLine());
            }

            using (System.IO.StreamReader file = new System.IO.StreamReader(fpath))
            {

                d = new float[z][,]; float[,] temp;
                for (int k = 0; k < d.GetLength(0); k++)
                {
                    temp = new float[x, y];
                    for (int j = 0; j < y; j++)
                        for (int i = 0; i < x; i++)
                            temp[i,j] = (float)Convert.ToDecimal(file.ReadLine());
                    d[k] = temp;
                }
            }
            return d;
        }

        
    }
    
    public struct StructureInfo
    {
        public string Name { get; set; }
        public string Info { get; set; }
        public int Size { get; set; }

        public void CreateInfo()
        {
            Info = "STRUCTURE: " + Name + ", " + Size + "x" + Size + "x" + Size;
        }

        public void CreateTestInfo(int radius)
        {
            Name = "Test Sphere r=" + radius;
            CreateInfo();
        }
    }
}
