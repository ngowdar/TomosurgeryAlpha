using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Media.Imaging;

namespace TomosurgeryAlpha
{
    public class DICOMImageSet : Image
    {
        public DICOMinfo Dinfo;
        public ArrayList dicomArray;
        public static string s_dictionarypath;
        public static openDicom.Registry.DataElementDictionary dd;
        public BackgroundWorker bw_imagemaker;
        public event ProgressChangedEventHandler ImageWorkerProgressChanged;
        public event RunWorkerCompletedEventHandler ImageWorkerCompleted;
        public float[][] f_imagearray;
        public float[] imagePosition;
        public float[] ZIndexArray;
        public int NumberOfImages;
        public DICOMImageSet(string[] paths)
        {
            dd = new openDicom.Registry.DataElementDictionary(s_dictionarypath, openDicom.Registry.DictionaryFileFormat.BinaryFile);                        

            /*For the first image, grab the absolute image position. This
             * only matters for the z-coordinate, because the x and y position
             * of the top-left pixel won't change with each image.
             */            
            
            NumberOfImages = paths.Length;
        }
        #region ImageMaker Background Worker methods
        public void Initialize_imagemaker_bw()
        {
            bw_imagemaker = new BackgroundWorker();
            bw_imagemaker.DoWork += new DoWorkEventHandler(bw_imagemaker_DoWork);
            bw_imagemaker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(bw_imagemaker_RunWorkerCompleted);
            bw_imagemaker.ProgressChanged += new ProgressChangedEventHandler(bw_imagemaker_ProgressChanged);
            bw_imagemaker.WorkerReportsProgress = true;
        }

        void bw_imagemaker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (ImageWorkerProgressChanged != null)
                ImageWorkerProgressChanged.Invoke(null, e);
        }

        void bw_imagemaker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            SetZIndex();
            Normalize();
            Dinfo = new DICOMinfo();
            Dinfo.NumImages = NumberOfImages;
            Dinfo.CreateInfo();
            if (ImageWorkerCompleted != null)
                ImageWorkerCompleted.Invoke(null, e);
        }

        public void Normalize()
        {
            float max=0; float temp;
            for (int i = 0; i < f_imagearray.GetLength(0); i++)
            { 
                temp = f_imagearray[i].Max();
                if (temp > max)
                    max = temp;
            }
            
            for (int i = 0; i < f_imagearray.GetLength(0); i++)
            {                
                f_imagearray[i] = Matrix.ScalarMultiply(ref f_imagearray[i], max);
            }

        }

        private void SetZIndex()
        {
            if (f_imagearray != null)
            {
                ZIndexArray = new float[f_imagearray.GetLength(0)];
                ZIndexArray[0] = imagePosition[2];
                for (int i = 1; i < f_imagearray.GetLength(0); i++)
                    ZIndexArray[i] = imagePosition[2] - i * 2;
            }
        }

        void bw_imagemaker_DoWork(object sender, DoWorkEventArgs e)
        {            
            CreateDICOMsFromFileList((string[])e.Argument);
            //bw_imagemaker.ReportProgress(100);
        }

        public void CreateDICOMsFromFileList(string[] files)
        {
            
            DICOMImageFile tempImg = new DICOMImageFile(files[0]);
            imagePosition = tempImg.getAbsoluteImgPosition();            
            f_imagearray = new float[files.GetLength(0)][];
            f_imagearray[0] = tempImg.makeFloatArray(tempImg.ba_PixelData);
            int ImagesSoFar = 1;
            Parallel.For(1, files.GetLength(0), s =>
            {
                DICOMImageFile i = new DICOMImageFile((string)files[s]);
                float[] f = i.makeFloatArray(i.ba_PixelData);
                //Bitmap b = CreateDICOMImageFromFloat(f, 100, 612);
                f_imagearray[s] = f;
                //imagebitmaps[s] = b;
                ImagesSoFar++;
                bw_imagemaker.ReportProgress(ImagesSoFar);
            });
            //imagearrays.TrimToSize();
        }               
        
        #endregion
    }

    public struct DICOMinfo
    {
        public string Name { get; set; }
        public string Info { get; set; }
        public int NumImages { get; set; }
        public void CreateInfo()
        {
            Info = "DICOM: " + NumImages + " images";
        }
    }
}
