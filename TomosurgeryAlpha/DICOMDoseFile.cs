﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Collections;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using openDicom;

namespace TomosurgeryAlpha
{
    public class DICOMDoseFile
    {
        openDicom.File.DicomFile DF;
        public StructureSet SS;
        public static string dictionarypath;
        public string path;        
        static openDicom.Registry.DataElementDictionary dd;
        byte[] pixeldata;
        public static float[][,] OriginalDose;
        //public static float[] dose;
        public static float[][,] Dose;
        public float ZStart;
        public static ushort rows;
        public static ushort columns;
        public ushort bitsalloc;
        static public long numframes;
        public static UInt16 maxdose;
        public ushort scaling;
        public static decimal[] doseoffset = new decimal[3];
        public static int[] structurebounds;
        public string dosesavepath;

        public DICOMDoseFile(string path, bool IsDICOM)
        {
            if (IsDICOM)
            {
                this.path = path;
                SetDictionaryPath();
                DF = new openDicom.File.DicomFile(path);
                OriginalDose = Matrix.Normalize(ExtractDoseData());                
                Dose = EnlargeAndInterpolate(OriginalDose);

                dosesavepath = System.IO.Path.Combine(PathSet.ActiveDirectory, "ExtractedDoseFile.bin");

                //WriteDoseToFile(dosesavepath);
                DEBUG_WriteFileSummary();
            }
            else
            {
                ReadDoseFromFile(path);
                Dose = EnlargeAndInterpolate(OriginalDose);
            }

        }

        private void SetDictionaryPath()
        {
            if (dictionarypath == null)
            {
                dictionarypath = CreateDictionaryFile(TomosurgeryAlpha.Properties.Resources.dicomdictionary);
                dd = new openDicom.Registry.DataElementDictionary(dictionarypath, openDicom.Registry.DictionaryFileFormat.BinaryFile);
            }
        }

        public decimal[] AssociateStructureSet(StructureSet ss)
        {
            SS = ss;
            SS.AssociateDose(this);
            decimal[] offset = new decimal[3];
            for (int i = 0; i < 3; i++)
                offset[i] = (decimal)SS.f_offset[i] - (decimal)doseoffset[i];
            structurebounds = new int[6];
            structurebounds = SS.FindAllAxisBoundaries();

            //TODO: Find relative structure bounds
            decimal[] zbound = new decimal[2] { structurebounds[4], structurebounds[5] };
            zbound[0] += offset[2];
            zbound[1] += offset[2];

            return zbound;

            
        }

        public void ReadDoseFromFile(string path)
        {
            using (FileStream fs = new FileStream(path, FileMode.Open, FileAccess.Read))
            using (StreamReader br = new StreamReader(fs))
            {
                rows = Convert.ToUInt16(br.ReadLine());
                columns = Convert.ToUInt16(br.ReadLine());
                numframes = Convert.ToInt32(br.ReadLine());
                doseoffset[0] = Convert.ToDecimal(br.ReadLine());
                doseoffset[1] = Convert.ToDecimal(br.ReadLine());
                doseoffset[2] = Convert.ToDecimal(br.ReadLine());
                scaling = Convert.ToUInt16(br.ReadLine());
                ZStart = (float)Convert.ToDecimal(br.ReadLine());
                int readlength = Convert.ToInt32(br.ReadLine());
                OriginalDose = new float[numframes][,];
                for (int k = 0; k < (int)numframes; k++)
                    for (int j = 0; j < rows; j++)
                        for (int i = 0; i < columns; i++)
                        {
                            OriginalDose[k][i, j] = (float)Convert.ToDecimal(br.ReadLine());
                        }
            }
            OriginalDose = Matrix.Normalize(OriginalDose);
            Debug.WriteLine("Successfully loaded dose file.");
            DEBUG_WriteFileSummary();

        }

        public string CreateDictionaryFile(byte[] b)
        {
            //Create path
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, "tempdict.bin");

            //Writing byte array to a file
            FileStream fs = new FileStream(path, FileMode.OpenOrCreate);

            BinaryWriter bw = new BinaryWriter(fs);
            bw.Write(b);
            bw.Close();

            DICOMImageFile.s_dictionarypath = path;
            DICOMImageSet.s_dictionarypath = path;
            dictionarypath = path;
            return path;
            //DICOMRT.dictionarypath = path;
            //DICOMdose.dictionarypath = path;
        }  

        private void DEBUG_WriteFileSummary()
        {
            Debug.WriteLine("Rows x Columns x NumFrames: " + rows + " x " + columns + " x " + numframes);
            Debug.WriteLine("Offset Vector: < " + doseoffset[0] + ", " + doseoffset[1] + ", " + doseoffset[2] + " > ");
            Debug.WriteLine("Scaling: " + scaling);
            Debug.WriteLine("ZStart: " + ZStart);
            Debug.WriteLine("Saved to: " + dosesavepath);
            Debug.WriteLine("Maximum value: " + Matrix.FindMax(OriginalDose));
            Debug.WriteLine("------------------------------");
        }

        public void WriteDoseToFile(string path)
        {
            //rows
            //columns
            //numframes
            //offsets:
            //x
            //y
            //z
            //scaling
            //zstart
            //total length
            //10 elements in***
            //start dose data here.

            FileStream ff;
            if (System.IO.File.Exists(path))
                ff = new FileStream(path, FileMode.Truncate, FileAccess.Write);
            else
                ff = new FileStream(path, FileMode.Create, FileAccess.Write);
            using (ff)
            using (StreamWriter ds = new StreamWriter(ff))
            {
                ds.WriteLine(rows);
                ds.WriteLine(columns);
                ds.WriteLine(numframes);
                ds.WriteLine(doseoffset[0]);
                ds.WriteLine(doseoffset[1]);
                ds.WriteLine(doseoffset[2]);
                ds.WriteLine(scaling);
                ds.WriteLine(ZStart);
                ds.WriteLine(OriginalDose.GetLength(0));
                for (int k = 0; k < OriginalDose.GetLength(0); k++)
                    for (int j = 0; j < OriginalDose[0].GetLength(1); j++)
                        for (int i = 0; i < OriginalDose[0].GetLength(0); i++)
                            ds.WriteLine(OriginalDose[k][i,j]);
            }
        }

        public float[][,] EnlargeAndInterpolate(float[][,] dose)
        {
            float[][,] r = new float[numframes][,];            
            for (int k = 0; k < numframes; k++)
            {                
                r[k] = Matrix.LinearlyInterpolateSlices(Matrix.LinearlyInterpolateSlices(dose[k]));
            }
            r = InterpolateInbetweenSlices(InterpolateInbetweenSlices(r));
            //r = Matrix.EnlargeAndCenter(Matrix.Convertto1D(r), StructureSet.padsize, r[0].GetLength(0), r[0].GetLength(1), r.GetLength(0));
            return Matrix.Normalize(r);
        }

        private static float[,] AverageSlices(float[,] A, float[,] B)
        {
            float[,] C = Matrix.Zeroes(A.GetLength(0), A.GetLength(1));
            for (int j = 0; j < A.GetLength(1); j++)
                for (int i = 0; i < A.GetLength(0); i++)
                    C[i, j] = (A[i, j] + B[i, j]) / 2;
            return C;
        }

        public static float[][,] InterpolateInbetweenSlices(float[][,] A)
        {
            int newlength = (A.GetLength(0) * 2) - 1;
            float[][,] C = new float[newlength][,];
            //fill in original slices:

            for (int i = 0; i < A.GetLength(0); i++)
                C[2*i] = (float[,])A[i].Clone();

            for (int i = 1; i < A.GetLength(0); i++)
                C[(2 * i) - 1] = AverageSlices(A[i - 1], A[i]);

                        
            return Matrix.ThresholdEq(C,0.4f);
        }


        private float[][,] ExtractDoseData()
        {
            openDicom.DataStructure.DataSet.DataSet m = DF.DataSet;
            openDicom.DataStructure.DataSet.Sequence alldata = m.GetJointSubsequences(); //concatenated long list of all data            
            

            foreach (openDicom.DataStructure.DataSet.DataElement data in alldata)
            {
                //try extracting pixel data, if not found in the PixelData tag
                if (data.Tag.Element == "0010" && data.Tag.Group == "7FE0")
                {
                    foreach (object b in data.Value)
                    {
                        pixeldata = (byte[])b;
                    }
                }

                //Finds rows
                if (data.Tag.Element == "0010" && data.Tag.Group == "0028")
                {
                    foreach (object c in data.Value)
                    {
                        rows = (ushort)c;
                    }
                }
                //finds columns
                if (data.Tag.Element == "0011" && data.Tag.Group == "0028")
                {
                    foreach (object c in data.Value)
                    {
                        columns = (ushort)c;
                    }
                }
                //bits allocated
                if (data.Tag.Element == "0100" && data.Tag.Group == "0028")
                {
                    foreach (object c in data.Value)
                    {
                        bitsalloc = (ushort)c;
                    }
                }
                //number of frames
                if (data.Tag.Element == "0008" && data.Tag.Group == "0028")
                {
                    foreach (object c in data.Value)
                    {
                        numframes = (long)c;
                    }
                }

                //get offset vector
                if (data.Tag.Element == "0032" && data.Tag.Group == "0020")
                {
                    int t = 0;
                    foreach (object c in data.Value)
                    {
                        doseoffset[t] = (decimal)c;
                        t++;
                    }
                    ZStart = (float)doseoffset[2];
                }                

            }
            //return Byte2Float(pixeldata);
            float[] temp_pixeldata = Byte2Float(pixeldata);
            maxdose = (UInt16)temp_pixeldata.Max();
            return Matrix.MakeJaggedFloat(temp_pixeldata, (int)columns, (int)rows, (int)numframes);
        }

        private float[] Byte2Float(byte[] data)
        {
            float[] output = new float[data.GetLength(0) / 2];
            for (int i = 0; i < data.GetLength(0); i += 2)
            {

                output[i / 2] = BitConverter.ToUInt16(data, i);
            }
            return output;
        }
    }
}
