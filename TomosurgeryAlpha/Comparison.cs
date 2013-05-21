using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Microsoft.Win32;
using System.Diagnostics;
using System.Windows.Media.Imaging;


namespace TomosurgeryAlpha
{
    public static class Comparison
    {
        public static DICOMDoseFile A;
        public static DICOMDoseFile B;
        public static WriteableBitmap Awb;
        public static WriteableBitmap Bwb;

        public static void LoadA()
        {
            var loadconfig = new OpenFileDialog();
            loadconfig.Title = "Select a dosespace or DICOM-RT dose file";
            //float[] dosespace;
            if (loadconfig.ShowDialog() != false)
            {
                Debug.WriteLine(Path.GetExtension(loadconfig.FileName));
                if (Path.GetExtension(loadconfig.FileName) == ".dcm")
                {
                    A = LoadDICOMDoseFile(loadconfig.FileName);
                    //StatusTxtBox.Text = "DICOM dose successfully loaded";
                }
                else if (Path.GetExtension(loadconfig.FileName) == ".txt")
                {
                    A = LoadDoseSpace(loadconfig.FileName);
                }
            }

        }

        public static void LoadB()
        {
            var loadconfig = new OpenFileDialog();
            loadconfig.Title = "Select a dosespace or DICOM-RT dose file";
            //float[] dosespace;
            if (loadconfig.ShowDialog() != false)
            {
                Debug.WriteLine(Path.GetExtension(loadconfig.FileName));
                if (Path.GetExtension(loadconfig.FileName) == ".dcm")
                {
                    B = LoadDICOMDoseFile(loadconfig.FileName);
                    //StatusTxtBox.Text = "DICOM dose successfully loaded";
                }
                else if (Path.GetExtension(loadconfig.FileName) == ".txt")
                {
                    B = LoadDoseSpace(loadconfig.FileName);
                }
            }
            
        }

        public static DICOMDoseFile LoadDoseSpace(string p)
        {
            return new DICOMDoseFile(p, false);
        }

        public static DICOMDoseFile LoadDICOMDoseFile(string p)
        {
            return new DICOMDoseFile(p, true);
        }

    }
}
