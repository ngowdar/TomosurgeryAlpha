using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Media.Imaging;
using openDicom.Registry;

namespace TomosurgeryAlpha
{
    public class DICOMImageSet : Image
    {
        public static string s_dictionarypath;
        public static DataElementDictionary dd;
        public DICOMinfo Dinfo;
        public int NumberOfImages;
        public float[] ZIndexArray;
        public BackgroundWorker bw_imagemaker;
        public float[][] f_imagearray;
        public float[] imagePosition;

        public DICOMImageSet(string[] paths)
        {
            dd = new DataElementDictionary(s_dictionarypath, DictionaryFileFormat.BinaryFile);

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
            bw_imagemaker.DoWork += bw_imagemaker_DoWork;
            bw_imagemaker.RunWorkerCompleted += bw_imagemaker_RunWorkerCompleted;
            bw_imagemaker.ProgressChanged += bw_imagemaker_ProgressChanged;
            bw_imagemaker.WorkerReportsProgress = true;
        }

        private void bw_imagemaker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            if (ImageWorkerProgressChanged != null)
                ImageWorkerProgressChanged.Invoke(null, e);
        }

        private void bw_imagemaker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
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
            float max = 0;
            float temp;
            for (int i = 0; i < f_imagearray.GetLength(0); i++)
            {
                temp = f_imagearray[i].Max();
                if (temp > max)
                    max = temp;
            }

            for (int i = 0; i < f_imagearray.GetLength(0); i++)
            {
                f_imagearray[i] = Matrix.ScalarMultiply(f_imagearray[i], max);
            }
        }

        private void SetZIndex()
        {
            if (f_imagearray != null)
            {
                ZIndexArray = new float[f_imagearray.GetLength(0)];
                ZIndexArray[0] = imagePosition[2];
                for (int i = 1; i < f_imagearray.GetLength(0); i++)
                    ZIndexArray[i] = imagePosition[2] - i*2;
            }
        }

        private void bw_imagemaker_DoWork(object sender, DoWorkEventArgs e)
        {
            CreateDICOMsFromFileList((string[]) e.Argument);
            //bw_imagemaker.ReportProgress(100);
        }

        public void CreateDICOMsFromFileList(string[] files)
        {
            var tempImg = new DICOMImageFile(files[0]);
            imagePosition = tempImg.GetAbsoluteImgPosition();
            DICOMImageFile.GlobalOffset = new decimal[3]
                                              {
                                                  (decimal) imagePosition[0], (decimal) imagePosition[1],
                                                  (decimal) imagePosition[2]
                                              };
            f_imagearray = new float[files.GetLength(0)][];
            f_imagearray[0] = tempImg.makeFloatArray(tempImg.ba_PixelData);
            int ImagesSoFar = 1;
            Parallel.For(1, files.GetLength(0), s =>
                                                    {
                                                        var i = new DICOMImageFile(files[s]);
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

        public event ProgressChangedEventHandler ImageWorkerProgressChanged;
        public event RunWorkerCompletedEventHandler ImageWorkerCompleted;
    }

    public struct DICOMinfo
    {
        public string Info { get; set; }
        public int NumImages { get; set; }

        public void CreateInfo()
        {
            Info = "DICOM: " + NumImages + " images";
        }
    }
}