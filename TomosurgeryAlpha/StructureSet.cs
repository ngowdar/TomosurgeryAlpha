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
using System.Diagnostics;

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
        public static int padsize;
        public string headerpath;
        public string tumorpath;
        public float[] f_offset;
        public float[][] f_structurearray;
        public static int size;
        public static int[] SS_dim;
        public static int[] BIG_dim;
        public ArrayList pointdata;
        public ArrayList pointF_data;
        public ArrayList shiftedpointdata;
        public ArrayList structure_names;
        public StructureInfo SI;
        public float[][,] fj_Tumor;
        public float[][,] fj_CS;


        public StructureSet(string header, string tumor)
        {
            if (DoseKernel.N != null)
                padsize = DoseKernel.N/2; // <- need to make this (N-1)/2 once N has been established.
            float[] LinearVolume = Read1DArrayFromFile(tumor, header);
            float[][,] BinaryVolume = EnlargeTumor(LinearVolume, padsize);
            f_structurearray = CreateArray(BinaryVolume, padsize);
            fj_Tumor = GetTumorOnly(BinaryVolume);
            fj_CS = GetCSOnly(BinaryVolume);
            size = f_structurearray.GetLength(0);
            SI = new StructureInfo();
            headerpath = header;
            tumorpath = tumor;

            
            ////TODO: Remove this, just a test
            //float[] d = new float[5] { 1, 2, 3, 4, 5 };
            //Matrix.WriteArrayAsList("d: ", d);
            //float[] dd = new float[5] { 6, 7, 8, 9, 10 };            
            //float[][] ff = new float[2][] { Matrix.InterpolateSingleVector(d), Matrix.InterpolateSingleVector(dd) };
            
            //float[,] awesome = Matrix.InterpolateVectors(ff);
            //for (int i = 0; i < awesome.GetLength(0); i++)
            //    Debug.Write(awesome[i, 1]);
            //float[,] awesome = Matrix.LinearlyInterpolateSlices(new float[,]{{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
            //for (int i = 0; i < awesome.GetLength(1); i++)
            //{
            //    Debug.Write(" " + awesome[0, i]);
            //    Debug.WriteLine(" ");
            //    Debug.Write(" " + awesome[1, i]);
            //    Debug.WriteLine(" ");
            //    Debug.Write(" " + awesome[2, i]);
            //    Debug.WriteLine(" ");
            //    Debug.Write(" " + awesome[3, i]);
            //    Debug.WriteLine(" ");
            //    Debug.Write(" " + awesome[4, i]);
            //}

        }

        public StructureSet(float[][] tumorobj)
        {            
            f_structurearray = tumorobj;
            //fj_Tumor = tumorobj;
            SI = new StructureInfo();

        }

        private float[][] CreateArray(float[][,] BV, int padsize)
        {
            float[][] temp = new float[BV.GetLength(0)][];
            
            //Fill in with zeroes first.
            for (int i = 0; i < temp.GetLength(0); i++)
                temp[i] = Matrix.Zero1DFloat(BV[0].GetLength(0) * BV[0].GetLength(1));

            //Fill in tumor portion.
            for (int i = 0; i < SS_dim[2]; i++)
                temp[padsize + i] = Convert2D_to1D(BV[padsize + i]);
            return temp;            
        }

        private float[] Convert2D_to1D(float[,] A)
        {
            float[] B = new float[A.GetLength(0)*A.GetLength(1)];
            for (int j = 0; j < A.GetLength(1); j++)
                for (int i = 0; i < A.GetLength(0); i++)
                    B[(j * A.GetLength(0)) + i] = A[i, j];
            return B;

        }

        private float[][,] GetCSOnly(float[][,] d)
        {
            float[][,] cs = new float[d.GetLength(0)][,];
            for (int i = 0; i < d.GetLength(0); i++)
                cs[i] = Matrix.Zeroes(d[0].GetLength(0), d[0].GetLength(1));
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
            float[][,] tumor = new float[d.GetLength(0)][,];
            for (int i = 0; i < d.GetLength(0); i++)
                tumor[i] = Matrix.Zeroes(d[0].GetLength(0), d[0].GetLength(1));
            for (int k = 0; k < d.GetLength(0); k++)
                for (int i = 0; i < d[0].GetLength(0); i++)
                    for (int j = 0; j < d[0].GetLength(1); j++)
                    {
                        float value = d[k][i, j];
                        if (value > 0)
                            tumor[k][i, j] = 1.0f;
                        else
                            tumor[k][i, j] = 0.0f;
                    }
            return tumor;
        }

        private float[][,] EnlargeTumor(float[] dd, int p)
        {   
            float[][,] d = Matrix.TransposeMatrix(Matrix.EnlargeAndCenter(dd, p, SS_dim[0], SS_dim[1], SS_dim[2]));

            Enlarge_SS_Dim(d);
            return d;
        }

        private void Enlarge_SS_Dim(float[][,] d)
        {
            BIG_dim = new int[3];
            BIG_dim[0] = d[0].GetLength(0);
            BIG_dim[1] = d[0].GetLength(1);
            BIG_dim[2] = d.GetLength(0);
        }

        public float[] Read1DArrayFromFile(string fpath, string hpath)
        {            
            SS_dim = new int[3];
            using (System.IO.StreamReader header = new System.IO.StreamReader(hpath))
            {
                SS_dim[0] = Convert.ToInt16(header.ReadLine());
                SS_dim[1] = Convert.ToInt16(header.ReadLine());
                SS_dim[2] = Convert.ToInt16(header.ReadLine());

                //Read in the global img_offset vector (this is the top-left corner, and sets the DICOM coordinate frame)
                StructureSet.f_global_xoffset = (float)Convert.ToDecimal(header.ReadLine());
                StructureSet.f_global_yoffset = (float)Convert.ToDecimal(header.ReadLine());
                StructureSet.f_global_zoffset = (float)Convert.ToDecimal(header.ReadLine());

            }

            float[] d = new float[SS_dim[0] * SS_dim[1] * SS_dim[2]];

            using (System.IO.StreamReader file = new System.IO.StreamReader(fpath))
            {                
                for (int k = 0; k < d.GetLength(0); k++)
                {         
                    d[k] = (float)Convert.ToDecimal(file.ReadLine());
                }
            }
            headerpath = hpath;
            tumorpath = fpath;
            return d;
        }
        
        public float[][,] ReadArrayFromFile(string fpath, string hpath, bool IsDoseFile)
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
                if (IsDoseFile)
                {
                    f_offset = new float[3];
                    f_offset[0] = (float)Convert.ToDouble(header.ReadLine());
                    f_offset[1] = (float)Convert.ToDouble(header.ReadLine());
                    f_offset[2] = (float)Convert.ToDouble(header.ReadLine());
                }
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
            tumorpath = fpath;
            headerpath = hpath;
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
