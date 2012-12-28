using System;
using System.Drawing;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using TomosurgeryAlpha.Properties;
using openDicom;
using System.IO;
using System.Resources;
using System.Reflection;
using System.Windows.Data;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.Threading;
using System.Threading.Tasks;
using openDicom.DataStructure.DataSet;
using openDicom.File;
using openDicom.Image;
using openDicom.Registry;

namespace TomosurgeryAlpha
{
    public class DICOMImageFile : Image
    {
        public static decimal[] GlobalOffset;
        public static string s_dictionarypath;
        public static DataElementDictionary dd;
        public byte[][] ba_PixelData;
        private DicomFile file;
        private int i_bitsallocated;
        private int i_height;
        private int i_samples;
        private int i_width;
        public string s_path;


        public DICOMImageFile(string path)
        {
            if (s_dictionarypath == null)
            {
                CreateDictionaryFile(Resources.dicomdictionary);
            }
            dd = new DataElementDictionary(s_dictionarypath, DictionaryFileFormat.BinaryFile);
            s_path = path;
            file = new DicomFile(path, false);

            //edf = new DicomFile(path);            
            var pd = new PixelData(file.DataSet);
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
            string path = Path.Combine(PathSet.ActiveDirectory, "tempdict.bin");

            //Writing byte array to a file
            var fs = new FileStream(path, FileMode.Create);

            var bw = new BinaryWriter(fs);
            bw.Write(b);
            bw.Close();

            s_dictionarypath = path;
            DICOMImageSet.s_dictionarypath = path;
            s_dictionarypath = path;
            //DICOMRT.dictionarypath = path;
            //DICOMdose.dictionarypath = path;
        }

        public float[] makeFloatArray(byte[][] data)
        {
            byte[] ba = data[0];
            var bytesArray = new byte[ba.GetLength(0)];
            for (int i = 0; i < bytesArray.GetLength(0); i++)
                //for (int j=0; j<ba.GetLength(1); j++)
            {
                int current_pos = i;
                bytesArray[i] = ba[current_pos];
            }
            float[] f = Byte2Float(bytesArray);
            return f;
        }

        public float[] Byte2Float(byte[] data)
        {
            var output = new float[data.GetLength(0)/2];
            for (int i = 0; i < data.GetLength(0); i += 2)
            {
                output[i/2] = BitConverter.ToUInt16(data, i);
            }
            return output;
        }

        public float[] GetAbsoluteImgPosition()
        {
            DataSet m = file.DataSet;
            Sequence alldata = m.GetJointSubsequences(); //concatenated long list of all data            
            var position = new float[3];
            foreach (DataElement data in alldata)
            {
                //try extracting pixel data, if not found in the PixelData tag
                if (data.Tag.Element == "0032" && data.Tag.Group == "0020")
                {
                    int t = 0;
                    foreach (object b in data.Value)
                    {
                        position[t] = (float) Convert.ToDouble(b);
                        t++;
                    }
                }
            }
            return position;
        }
    }
}