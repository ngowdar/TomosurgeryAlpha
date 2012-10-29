using System;
using System.Drawing;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using openDicom;
using System.IO;
using System.Resources;
using System.Reflection;
using System.Windows.Data;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;

namespace TomosurgeryAlpha
{
    public class DICOMImageFile : Image
    {
        public static decimal[] GlobalOffset;
        public string s_path;
		public static string s_dictionarypath;
        public static openDicom.Registry.DataElementDictionary dd;
		public byte[][] ba_PixelData;
		public float[] fa_PixelData;
		openDicom.File.DicomFile file;        
		int i_width;
		int i_height;
		int i_bitsallocated;
		int i_stride;
        int i_samples;
		float f_Zposition;
		public float[] fa_OffsetCoords;
        public static BackgroundWorker bw_imagemaker;        
		
		
		public DICOMImageFile(string path)
		{
            if (s_dictionarypath == null)
            {
                CreateDictionaryFile(TomosurgeryAlpha.Properties.Resources.dicomdictionary);
            }
            dd = new openDicom.Registry.DataElementDictionary(s_dictionarypath, openDicom.Registry.DictionaryFileFormat.BinaryFile);
            s_path = path;            
            file = new openDicom.File.DicomFile(path, false);

            //edf = new DicomFile(path);            
            openDicom.Image.PixelData pd = new openDicom.Image.PixelData(file.DataSet);
            ba_PixelData = pd.ToBytesArray();

            i_bitsallocated = pd.BitsAllocated;            
            i_width = pd.Columns;
            i_height = pd.Rows;
            //PlanarConfig = pd.PlanarConfiguration;
            i_samples = pd.SamplesPerPixel;
		}

        public void CreateDictionaryFile(byte[] b)
        {
            //Create path
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, "tempdict.bin");

            //Writing byte array to a file
            FileStream fs = new FileStream(path, FileMode.Create);

            BinaryWriter bw = new BinaryWriter(fs);
            bw.Write(b);
            bw.Close();

            DICOMImageFile.s_dictionarypath = path;
            DICOMImageSet.s_dictionarypath = path;
            s_dictionarypath = path;
            //DICOMRT.dictionarypath = path;
            //DICOMdose.dictionarypath = path;
        }     		
		
		public float[] makeFloatArray(byte[][] data)
		{
            byte[] ba = data[0];
            byte[] bytesArray = new byte[ba.GetLength(0)];
            for (int i = 0; i < bytesArray.GetLength(0); i++)
            //for (int j=0; j<ba.GetLength(1); j++)
            {
                int current_pos = i;
                bytesArray[i] = ba[current_pos];
            }
            float[] f = Byte2Float(bytesArray);
            return f;
		}

		public float[] makeFloatArray(double[,,] d)
		{
            float[] result;
            int a = d.GetLength(0); int b = d.GetLength(1); int c = d.GetLength(2);
            result = new float[a * b];
            for (int i = 0; i < a; i++)
                for (int j = 0; j < b; j++)
                    for (int k = 0; k < c; k++)
                        result[k * a * b + j * a + i] = (float)d[i, j, k];
            return result;
		}
		public float[] makeFloatArray(int[,,] d)
		{
            float[] result;
            int a = d.GetLength(0); int b = d.GetLength(1); int c = d.GetLength(2);
            result = new float[a * b];
            for (int i = 0; i < a; i++)
                for (int j = 0; j < b; j++)
                    for (int k = 0; k < c; k++)
                        result[k * a * b + j * a + i] = (float)d[i, j, k];
            return result;
		}
		
		public byte[] make1DByteArray(byte[][] data )
		{
            byte[] ba = data[0];
            byte[] bytesArray = new byte[ba.GetLength(0)];
            for (int i = 0; i < bytesArray.GetLength(0); i++)
            //for (int j=0; j<ba.GetLength(1); j++)
            {
                int current_pos = i;
                bytesArray[i] = ba[current_pos];
            }
            return bytesArray;
		}
		
		public float[] Byte2Float(byte[] data)
		{
            float[] output = new float[data.GetLength(0) / 2];
            for (int i = 0; i < data.GetLength(0); i += 2)
            {

                output[i / 2] = BitConverter.ToUInt16(data, i);
            }
            return output; 
		}

        public float[] getAbsoluteImgPosition()
        {
            openDicom.DataStructure.DataSet.DataSet m = file.DataSet;
            openDicom.DataStructure.DataSet.Sequence alldata = m.GetJointSubsequences(); //concatenated long list of all data            
            float[] position = new float[3];
            foreach (openDicom.DataStructure.DataSet.DataElement data in alldata)
            {
                //try extracting pixel data, if not found in the PixelData tag
                if (data.Tag.Element == "0032" && data.Tag.Group == "0020")
                {
                    int t = 0;
                    foreach (object b in data.Value)
                    {
                        position[t] = (float)Convert.ToDouble(b);
                        t++;
                    }
                }
            }            
            return position;
        }

        
    }
}
