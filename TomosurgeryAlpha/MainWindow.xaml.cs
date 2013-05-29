using System;
using System.Drawing;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Forms;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.IO;
using System.Threading.Tasks;
using System.Threading;
using System.ComponentModel;
using System.Diagnostics;
using Ionic.Utils;
using Binding = System.Windows.Data.Binding;
using Color = System.Windows.Media.Color;
using MessageBox = System.Windows.MessageBox;
using MouseEventArgs = System.Windows.Input.MouseEventArgs;
using OpenFileDialog = Microsoft.Win32.OpenFileDialog;
using Path = System.IO.Path;
using Point = System.Drawing.Point;

namespace TomosurgeryAlpha
{
    /// <summary>
    ///   Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public static int N;
        public static string SessionName;
        public static WriteableBitmap wb_DICOM;
        public static WriteableBitmap wb_DS;
        public static WriteableBitmap wb_DDS;
        public static WriteableBitmap wb_Plan;
        private bool AlignmentOn;
        public int CursorRadius = 7;
        public double DICOM_aspectMultiplier = 1.0;
        public DoseKernel DK;
        public bool HasPreviewBeenCreated = false;
        public bool IsDICOMLoaded = false;
        public bool IsDoseEditable = false;
        public bool IsDoseLoaded = false;
        public bool IsPathSetCreated = false;
        public bool IsSSLoaded = false;
        public Coordinates LGKcoords;
        
        public bool Normalize = false;
        public PathSet PS;
        public bool PlanOptimized = false;
        public double Plan_aspectMultiplier_x = 1.0;
        public double Plan_aspectMultiplier_y = 1.0;
        public StructureSet SS;
        public int[] blue_colormap;
        private PointF[] circle;
        public DICOMDoseFile ddf;
        public int[] green_colormap;
        private float[,] mask;
        public Stopwatch plantimer;
        public int[] red_colormap;
        public DICOMImageSet set;

        public MainWindow()
        {
            InitializeComponent();
            DICOM_aspectMultiplier = 256/DICOM_imgbox.Width;
            SetParameterSliderLimits();
            ColormapTool(Properties.Resources.BWheatmap3);
            circle = new PointF[1];
            CreateCirclePoints();
            ResetMask();
        }

        #region UI Update Methodsq

        private void UpdateTextBlock(string s)
        {
            textBlock2.Text = s + "\n";
        }

        private void UpdateStatusBar(string s)
        {
            textBlock1.Text = s;
        }

        private void UpdateProgressBar(double v)
        {
            progressBar1.Value = v;
        }

        #endregion

        public void InitializeDICOM()
        {
            wb_DICOM = new WriteableBitmap(256, 256, 96, 96, PixelFormats.Bgr32, null);
            DICOM_imgbox.Source = wb_DICOM;
            DICOM_imgbox.Stretch = Stretch.UniformToFill;
            DICOM_imgbox.MouseMove += DICOM_imgbox_MouseMove;
            DICOM_imgbox.MouseLeftButtonDown += DICOM_imgbox_MouseLeftButtonDown;
            DICOM_imgbox.MouseWheel += DICOM_imgbox_MouseWheel;
        }

        public void InitializeDDS()
        {
            wb_DDS = new WriteableBitmap(256, 256, 96, 96, PixelFormats.Bgr32, null);
            DDS_imgbox.Source = wb_DDS;
            DDS_imgbox.Stretch = Stretch.UniformToFill;
            DDS_imgbox.MouseMove += DDS_imgbox_MouseMove;
            DDS_imgbox.MouseWheel += DDS_imgbox_MouseWheel;
            DDS_imgbox.MouseLeftButtonDown += DDS_imgbox_MouseLeftButtonDown;
            DDS_imgbox.MouseRightButtonDown += DDS_imgbox_MouseRightButtonDown;
        }

        public void LOAD_CONFIG_FILE(string path)
        {
            /* STANDARD CONFIG FORMAT
             * 
             * This is the line order:
             * ===========================
             * tumorpath
             * headerpath
             * shotsize (this is a number: 4, 8, 16)
             * rasterwidth
             * stepsize
             * slicethickness
             * dosecalculationthickness
             * workingdirectory
             * 
             */
            string tumorpath;
            string headerpath;
            int which_size;
            string rasterwidth;
            string stepsize;
            string slicethickness;
            int dosecalculationthickness;
            string workdir;
            string dosesavepath;

            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read))
            using (var br = new StreamReader(fs))
            {
                tumorpath = br.ReadLine();
                headerpath = br.ReadLine();
                which_size = Convert.ToInt16(br.ReadLine());
                rasterwidth = br.ReadLine();
                stepsize = br.ReadLine();
                slicethickness = br.ReadLine();
                dosecalculationthickness = Convert.ToInt16(br.ReadLine());
                workdir = br.ReadLine();
                dosesavepath = br.ReadLine();
            }
            txt_rasterwidth.Text = rasterwidth;
            txt_slicethickness.Text = slicethickness;
            txt_stepsize.Text = stepsize;
            PathSet.DCT = dosecalculationthickness;


            //SET WORKING DIRECTORY
            SetWorkingDirectory(workdir);

            //LOAD DOSE KERNEL
            if (which_size == 4)
                Load4mmDefault();
            else
            {
                MessageBox.Show(
                    "You haven't created the Load8mmDefault and Load16mmDefault functions yet. Using 4mm for now");
                Load4mmDefault();
            }

            //LOAD STRUCTURE
            SS = new StructureSet(headerpath, tumorpath);
            if (SS.f_structurearray != null)
            {
                IsSSLoaded = true;
                slider2.Minimum = 0;
                slider2.Maximum = SS.f_structurearray.GetLength(0);
                tabControl1.SelectedIndex = 1;
                slider2.Value = (SS.f_structurearray.GetLength(0)/2);
                DisplayStructure(SS.f_structurearray.GetLength(0)/2);
                AddStructureLoadedToListBox();
                Plan_btn.IsEnabled = true;
            }

            //LOAD DICOM DOSE
            LoadDoseSpace(dosesavepath);
        }

        public void WRITE_CONFIG_FILE(string name)
        {
            /* STANDARD CONFIG FORMAT
             * 
             * This is the line order:
             * ===========================
             * tumorpath
             * headerpath 
             * shotsize (this is a number: 4, 8, 16)
             * rasterwidth
             * stepsize
             * slicethickness
             * dosecalculationthickness
             * workingdirectory
             * DICOMdosepath "dosesavepath"
             */
            string path = Path.Combine(PathSet.ActiveDirectory, name);
            using (var fs = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            using (var sw = new StreamWriter(fs))
            {
                sw.WriteLine(SS.tumorpath);
                sw.WriteLine(SS.headerpath);
                sw.WriteLine(4);
                sw.WriteLine(txt_rasterwidth.Text);
                sw.WriteLine(txt_stepsize.Text);
                sw.WriteLine(txt_slicethickness.Text);
                sw.WriteLine(PathSet.DCT);
                sw.WriteLine(PathSet.ActiveDirectory);
                sw.WriteLine(Analysis.ddf.dosesavepath);
                //TODO: finish.
            }
        }

        public int GetCurrentSlice()
        {
            //slider2.Maximum = PathSet.NumSlices - 1;
            //slider2.Minimum = 0;
            
            int value;
            value = (int) Math.Round(slider2.Value);
            UpdateZPositionInfo(value, (int) slider2.Maximum);
            return value;
        }

        public int GetCurrentZSlice()
        {
            slider2.Maximum = PS.DoseSpace.GetLength(0) - 1;
            slider2.Minimum = 0;
            return (int) Math.Round(slider2.Value);
        }

        private void UpdateZPositionInfo(int pos, int of)
        {
            string s = "" + pos + " of " + of;
            zpos_index_lbl.Content = s;
        }

        private void Initialize_DICOMimghandlers()
        {
            set.Initialize_imagemaker_bw();
            set.ImageWorkerCompleted += DICOMImageSet_ImageWorkerCompleted;
            set.ImageWorkerProgressChanged += DICOMImageSet_ProgressChanged;
        }

        private void LoadDICOM_Menu_Click(object sender, RoutedEventArgs e)
        {
            CreateDictionaryFile(Properties.Resources.dicomdictionary);

            string[] dicomfiles = null;

            //CreateDictionaryFile();            
            var opendicom = new OpenFileDialog();
            opendicom.Title = "Select Image .DCM files (i.e. IMG00000_.dcm)";
            opendicom.Filter = "DICOM image files (.dcm)|*.dcm";
            opendicom.Multiselect = true;
            if (opendicom.ShowDialog() != false)
            {
                dicomfiles = opendicom.FileNames;
                //DICOMImageFile.imagebitmaps = new Bitmap[dicomfiles.GetLength(0)];                
                set = new DICOMImageSet(dicomfiles);
                Initialize_DICOMimghandlers();
                //Initialize and run the imagemaker background worker.                               
                set.bw_imagemaker.RunWorkerAsync(dicomfiles);
            }
        }

        private void CreateDictionaryFile(byte[] b)
        {
            //Create path
            string path = ReturnFullPath("tempdict.bin");

            //Writing byte array to a file
            var fs = new FileStream(path, FileMode.Create);

            var bw = new BinaryWriter(fs);
            bw.Write(b);
            bw.Close();

            DICOMImageFile.s_dictionarypath = path;
            DICOMImageSet.s_dictionarypath = path;
            //DICOMRT.dictionarypath = path;
            //DICOMdose.dictionarypath = path;
        }

        private string ReturnFullPath(string add)
        {
            string output;
            output = "C:\\TomosurgeryWorking\\" + add;
            return output;
        }

        private void sliderA_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            float[,] Aimg;
            sliderA.Minimum = 0;
            sliderA.Maximum = Comparison.A.OriginalDose.GetLength(0) - 1;
            int slice = (int) Math.Round(sliderA.Value);
            /*DS_z_lbl.Content = "Z: " +
                                       (DICOMDoseFile.doseoffset[2] - (decimal)slider2.Minimum -
                                        (decimal)(Math.Round(slider2.Value, 2) - slider2.Minimum));
            DS_index_lbl.Content = "Actual Z: " + slice;*/
            
            Aimg = Comparison.A.OriginalDose[slice];
            DisplayOnA(Aimg);
        }

        private void sliderB_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            float[,] Bimg;
            sliderB.Minimum = 0;
            sliderB.Maximum = Comparison.B.OriginalDose.GetLength(0) - 1;
            int slice = (int)Math.Round(sliderB.Value);
            /*DS_z_lbl.Content = "Z: " +
                                       (DICOMDoseFile.doseoffset[2] - (decimal)slider2.Minimum -
                                        (decimal)(Math.Round(slider2.Value, 2) - slider2.Minimum));
            DS_index_lbl.Content = "Actual Z: " + slice;*/
            
            Bimg = Comparison.B.OriginalDose[slice];
            DisplayOnB(Bimg);
        }

        private void slider2_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            var slice = (int) Math.Round(slider2.Value);
            if (tabControl1.SelectedIndex == 0) //"DICOM" tab
            {
                slider2.Minimum = 0;
                slider2.Maximum = set.f_imagearray.GetLength(0) - 1;
                //WriteableBitmap wb = (WriteableBitmap)DICOM_imgbox.Source;                
                DisplayDICOM(set.f_imagearray[slice]);

                //slicepos_label.Content = "Actual Z: " + (int)DICOMImageSet.ZIndexArray[slice];
                imgnumb_label.Content = "Image: " + slice;
            }
            else if (tabControl1.SelectedIndex == 1) //"Structure/DDS" tab
            {
                if (SS.fj_Tumor != null)
                {
                    slider2.Minimum = 0;
                    slider2.Maximum = SS.fj_Tumor.GetLength(0) - 1;
                    if (slice >= SS.fj_Tumor.GetLength(0) || slice < 0)
                    {
                        slice = (SS.fj_Tumor.GetLength(0)/2);
                    }
                }

                DisplayStructure(slice);
                DDS_z_lbl.Content = "Z: " +
                                    ((decimal) StructureSet.f_global_zoffset -
                                     (decimal) 0.25*((decimal) Math.Round(slider2.Value, 2)));
                DDS_index_lbl.Content = "Actual Z: " + slice;
            }

            else if (tabControl1.SelectedIndex == 3) //"DS" tab
            {
                float[,] img;
                if (dicomdose_rb_btn.IsChecked == true)
                {
                    //slider2.Minimum = Math.Abs((float)StructureSet.f_global_zoffset - (float)DICOMDoseFile.doseoffset[2]);
                    //slider2.Maximum = DICOMDoseFile.DoseJ.GetLength(0)-slider2.Minimum;
                    DS_z_lbl.Content = "Z: " +
                                       (DICOMDoseFile.doseoffset[2] - (decimal) slider2.Minimum -
                                        (decimal) (Math.Round(slider2.Value, 2) - slider2.Minimum));
                    DS_index_lbl.Content = "Actual Z: " + slice;
                    if (origdose_rb_btn.IsChecked == true)
                    {
                        slider2.Minimum = 0;
                        slider2.Maximum = ddf.OriginalDose.GetLength(0) - 1;
                        //if (Normalize)
                          //  img = Matrix.NormalizeBy(DICOMDoseFile.OriginalDose[slice], 65535.0f);
                        //else
                            img = ddf.OriginalDose[slice];
                        if (ShowStruct_chkbox.IsChecked == true)
                            DisplayDoseWithStruct(img);
                        else
                            Display2DFloat(img);
                    }
                    else if (newdose_rb_btn.IsChecked == true)
                    {
                        slider2.Minimum = 0;
                        slider2.Maximum = ddf.Dose.GetLength(0) - 1;
                        //Display2DFloat(DICOMDoseFile.Dose[slice]);
                        if (Normalize)
                            img = Matrix.NormalizeBy(ddf.Dose[slice], 65535.0f);
                        else
                            img = ddf.Dose[slice];
                        if (ShowStruct_chkbox.IsChecked == true)
                            DisplayDoseWithStruct(img);
                        else
                            Display2DFloat(img);
                    }
                }
                else if (PS.DoseSpace != null)
                {
                    slider2.Minimum = 0;
                    slider2.Maximum = PS.DoseSpace.GetLength(0) - 1;
                    //Display2DFloat(PS.DoseSpace[GetCurrentSlice()]);

                    if (Normalize)
                        img = Matrix.NormalizeBy(PS.DoseSpace[GetCurrentSlice()], PS.Max);
                    else
                        img = PS.DoseSpace[GetCurrentZSlice()];
                    if (ShowStruct_chkbox.IsChecked == true)
                        DisplayDoseWithStruct(img);
                    else
                        Display2DFloat(img);
                    
                    DS_z_lbl.Content = slice;
                }
            }

            else if (tabControl1.SelectedIndex == 2) //i.e. the "Plan" tab
            {
                if (PS != null)
                {
                    if (PS.RasterPaths != null)
                    {
                        slider2.Minimum = 0;
                        slider2.Maximum = PathSet.NumSlices - 1;
                        dataGrid1.SelectedIndex = GetCurrentSlice();
                        UpdateSliceLabels();
                        DisplayPlan();
                    }
                    else
                        SilenceTrackingLabels();
                }
            }
        }


        private void LoadStructure_Menu_Click(object sender, RoutedEventArgs e)
        {
            Stream dicomfile = null;

            var opendicom = new OpenFileDialog();
            opendicom.Multiselect = true;
            if (opendicom.ShowDialog() != false)
            {
                if (opendicom.FileNames.GetLength(0) > 1)
                {
                    string h = "";
                    string t = "";
                    foreach (string s in opendicom.FileNames)
                    {
                        var f = new FileInfo(s);
                        if (f.Extension == ".h")
                            h = s;
                        if (f.Extension == ".txt")
                            t = s;
                    }
                    string path = Path.GetDirectoryName(opendicom.FileName);
                    SetWorkingDirectory(path);
                    SS = new StructureSet(h, t);
                    string structname = Path.GetFileNameWithoutExtension(SS.tumorpath);
                    string structinfo = string.Concat("\n Width//Height//Frames: ", StructureSet.SS_dim[0], ", ",
                                                      StructureSet.SS_dim[1], ", ", StructureSet.SS_dim[2]);
                    StatusTxtBox.Text = structname + structinfo;
                }
                else
                {
                    dicomfile = opendicom.OpenFile();
                    //TM = new DICOMRT(opendicom.FileName, 0);                       
                }

                if (SS.f_structurearray != null)
                {
                    IsSSLoaded = true;
                    slider2.Minimum = 0;
                    slider2.Maximum = SS.f_structurearray.GetLength(0);
                    tabControl1.SelectedIndex = 1;
                    slider2.Value = (SS.f_structurearray.GetLength(0)/2);
                    DisplayStructure(SS.f_structurearray.GetLength(0)/2);
                    AddStructureLoadedToListBox();

                    PS = new PathSet(SS.fj_Tumor, Convert.ToInt16(txt_slicethickness.Text),
                             Convert.ToInt16(txt_slicethickness.Text), DK, SS);

                    Binding myBinding = new Binding("Name");
                    myBinding.Source = PathSet.ST;
                    suggthick_lbl.SetBinding(ContentProperty, myBinding);



                    AttachPSHandlers();
                    tabControl1.SelectedIndex = 1;
                    //slider2.Minimum = 0;
                    //slider2.Maximum = PathSet.NumSlices - 1;
                    if (DK != null)
                        PathSet.DK = DK;
                    if (SS != null)
                        PS.SS = SS;
                    //StatusTxtBox.Text = WriteArrayAsList(PS.SlicePositions, "slice positions");

                    Plan_btn.IsEnabled = true;
                }
                else
                    MessageBox.Show("Structure didn't load properly...");
            }
        }

        private void DDStab_IsVisibleChanged(object sender, DependencyPropertyChangedEventArgs e)
        {
            slider2.Minimum = 0;
            if (SS != null)
            {
                slider2.Maximum = SS.f_structurearray.GetLength(0);
                DisplayStructure((int) slider2.Value);
            }
        }

        private void button1_Click(object sender, RoutedEventArgs e)
        {
            if (aligning_helper_label.IsEnabled == false)
            {
                AlignmentOn = true;
                aligning_helper_label.IsEnabled = true;
                aligning_helper_label.Content = "Please click on the upper left post location.";
                LGKcoords = new Coordinates();
            }
            else if (LGKcoords.leftset == false)
            {
                LGKcoords.left = set.ZIndexArray[(int) Math.Round(slider2.Value)];
                LGKcoords.left_slider = (int) Math.Round(slider2.Value);
                LGKcoords.leftset = true;
                aligning_helper_label.Content = "Line up right fiducial...";
            }
            else if (LGKcoords.leftset && LGKcoords.rightset == false)
            {
                LGKcoords.right = set.ZIndexArray[(int) Math.Round(slider2.Value)];
                LGKcoords.right_slider = (int) Math.Round(slider2.Value);
                LGKcoords.rightset = true;
                aligning_helper_label.Content = "LGK Coordinates set.";

                AlignmentOn = false;
                LGKcoords.SetLGK_Zoffset();
                label1.Content = "( , )";
            }
        }


        public void CreateAlignmentMarkings(double x, double y)
        {
            var P = new PointF[7];
            float mid_y = (float) (y + (y + Coordinates.YPostWidth))/2;
            var mid_xright = (float) (x + Coordinates.XPostWidth);
            var mid_x = (float) ((x + (x + Coordinates.XPostWidth))/2);
            P[0] = new PointF((float) x, (float) y);
            P[1] = new PointF((float) x, (float) (y + Coordinates.YPostWidth));
            P[2] = new PointF((float) (x + Coordinates.XPostWidth), (float) y);
            P[3] = new PointF((float) (x + Coordinates.XPostWidth), (float) (y + Coordinates.YPostWidth));
            P[4] = new PointF((float) x, mid_y);
            P[5] = new PointF(mid_xright, mid_y);
            P[6] = new PointF(mid_x, mid_y);
            LGKcoords.align = P;
            LGKcoords.x100 = mid_x;
            LGKcoords.y100 = mid_y;
        }

        private void CreatePaths()
        {
            //PS.ProcessVolume(Convert.ToInt16(txt_slicethickness.Text), Convert.ToInt16(txt_slicethickness.Text));
            if (IsSSLoaded)
            {
                CreatePathSet();
                CreateDataGridColumns();
                if (IsDoseLoaded)
                {
                    CreatePreviewDose(PS);
                }
                PopulateDataGrid(PS);
                tabControl1.SelectedIndex = 3;
                IsPathSetCreated = true;
                HasPreviewBeenCreated = true;
                DisplayPlan();
            }
            else
            {
                MessageBox.Show(
                    "I can't create a plan until you properly load a structure object. Please load a DICOM-RT file first.");
            }
        }

        private void CreatePreviewDose(PathSet PS)
        {
            foreach (RasterPath rp in PS.RasterPaths)
            {
                rp.Calculate2DDoseSpace(DK.midplane, null);
                rp.CreateSliceInfo();
            }
        }

        private void CreateDataGridColumns()
        {
            var slicenum = new DataGridTextColumn();
            var numlines = new DataGridTextColumn();
            var numshots = new DataGridTextColumn();
            var coverage = new DataGridTextColumn();
            slicenum.Binding = new Binding("SliceNumber");
            slicenum.Header = "Slice #";
            numlines.Binding = new Binding("NumberOfLines");
            numlines.Header = "Lines";
            numshots.Binding = new Binding("NumberOfShots");
            numshots.Header = "Shots";
            coverage.Binding = new Binding("Coverage");
            coverage.Header = "Coverage";
            dataGrid1.Columns.Add(slicenum);
            dataGrid1.Columns.Add(numlines);
            dataGrid1.Columns.Add(numshots);
            dataGrid1.Columns.Add(coverage);
        }

        private void PopulateDataGrid(PathSet ps)
        {
            dataGrid1.Items.Clear();
            for (int i = 0; i < ps.RasterPaths.Count; i++)
            {
                var rp = (RasterPath) ps.RasterPaths[i];
                rp.info.SliceNumber = i + 1;
                dataGrid1.Items.Add(rp.info);
            }
        }

        private void RefreshDataGrid()
        {
            for (int i = 0; i < dataGrid1.Items.Count; i++)
            {
                var rp = (RasterPath) PS.RasterPaths[i];
                rp.CreateSliceInfo();
                rp.info.SliceNumber = i + 1;
                dataGrid1.Items[i] = rp.info;
            }
        }

        private string WriteArrayAsList(int[] array, string name)
        {
            string s = "" + name + ": ";
            for (int i = 0; i < array.GetLength(0); i++)
            {
                name += ", " + array[i];
            }
            return s;
        }

        private void CreatePathSet()
        {
            RasterPath.StepSize = Convert.ToInt16(txt_stepsize.Text);
            RasterPath.RasterWidth = Convert.ToInt16(txt_rasterwidth.Text);
            int padding = Convert.ToInt16(txt_slicethickness.Text)/2;
            //double sum1 = Matrix.SumAll(SS.fj_Tumor);
            //double sum2 = SS.f_structurearray[130].Sum();
            PS.ProcessVolume(Convert.ToInt16(txt_slicethickness.Text), Convert.ToInt16(txt_slicethickness.Text));
        }

        private void AttachPSHandlers()
        {
            PS.PathsetWorkerCompleted += PS_1_PathsetWorkerCompleted;
            PS.PathsetWorkerProgressChanged += PS_1_PathsetWorkerProgressChanged;
            RasterPath.SliceWorkerProgressChanged += RasterPath_SliceWorkerProgressChanged;
            RasterPath.ShotWeightsProgressHandler += RasterPath_ShotWeightsProgressHandler;
            PS.OptimizationWorkerCompleted += PS_2_OptimizationWorkerCompleted;
            PS.OptimizationWorkerProgress += PS_2_OptimizationProgress;
            PS.SliceweightWorkerProgress += PS_3_SliceweightWorkerProgress;
            PS.SliceweightWorkerCompleted += PS_3_SliceweightWorkerCompleted;
        }

        private void RasterPath_ShotWeightsProgressHandler(object sender, ProgressChangedEventArgs e)
        {
            //double array contains first value representing which slice, then the rest of the values are the shotweights.
            var ShotWeights = (double[]) e.UserState;
        }

        private void PS_3_SliceweightWorkerProgress(object sender, ProgressChangedEventArgs e)
        {
            UpdateProgressBar(e.ProgressPercentage);
        }

        private void PS_2_OptimizationProgress(object sender, ProgressChangedEventArgs e)
        {
            UpdateProgressBar(e.ProgressPercentage);
        }

        private void PS_3_SliceweightWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            //PS.DoseSpace = Matrix.Normalize(PS.DoseSpace);
            if (plantimer.IsRunning)
                Debug.WriteLine("Slice weighting time mark: " + plantimer.Elapsed);

            PS.RetrieveFinalizedDoseSpace();
            if (plantimer.IsRunning)
            {
                Debug.WriteLine("TOTAL TIME: " + plantimer.Elapsed);
                plantimer.Stop();
            }
            SetUpAnalysis();
            FindBeamOnTime();
            tabControl1.SelectedIndex = 4;


            MessageBox.Show(
                "Optimization complete. You may run an analysis using the analysis tab, or save/export using the buttons shown.");
            UpdateTextBlock2("Optimization complete. See analysis tab for details. Save using buttons below.");
            UpdateStatusBar("Ready.");
            save_plan_btn.IsEnabled = true;
            export_shots_btn.IsEnabled = true;
        }

        private void PS_2_OptimizationWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (plantimer.IsRunning)
                Debug.WriteLine("Current time after temp doses written (PS2): " + plantimer.Elapsed);
            UpdateTextBlock2("Temporary doses written. Optimizing slice weights...");
        }

        private void LoadTest_menu_Click(object sender, RoutedEventArgs e)
        {
            LoadDefaultSetup();
        }

        private void LoadDefaultSetup()
        {
            CreateTumorObject(80);
            IsSSLoaded = true;
            slider2.Minimum = 0;
            slider2.Maximum = SS.f_structurearray.GetLength(0);
            tabControl1.SelectedIndex = 1;
            slider2.Value = (SS.f_structurearray.GetLength(0)/2);
            DisplayStructure(SS.f_structurearray.GetLength(0)/2);
            AddStructureLoadedToListBox();
        }

        private void SetParameterSliderLimits()
        {
            int stepsize = Convert.ToInt16(txt_stepsize.Text);
            int rastwidth = Convert.ToInt16(txt_rasterwidth.Text);
            int slicethick = Convert.ToInt16(txt_slicethickness.Text);
            slider_stepsize.Minimum = 1;
            slider_stepsize.Maximum = 40;
            slider_rasterwidth.Minimum = 1;
            slider_rasterwidth.Maximum = 40;
            slider_slicethickness.Minimum = 1;
            slider_slicethickness.Maximum = 50;
        }


        private void slider_rasterwidth_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            txt_rasterwidth.Text = slider_rasterwidth.Value.ToString();
            RasterPath.RasterWidth = (int) Math.Round(slider_rasterwidth.Value);
            //PathSet.RasterWidth = RasterPath.RasterWidth;
            if (PS != null)
            {
                if (PS.RasterPaths != null && !PlanOptimized)
                {
                    RasterPath rp = GetCurrentSliceRP();
                    rp.Lines =  rp.LineSpacer(rp.boundaries[0], rp.boundaries[1], rp.LineSidePadding, rp.WhichSlice);
                    UpdateCurrentSlicePaths();
                    RefreshDataGrid();
                }
                else if (PlanOptimized)
                {
                    redwarn_lbl.Content = "Plan already exists! Please click RePlan to see changes";
                    redwarn_lbl.IsEnabled = true;
                }
            }
        }

        private void slider_stepsize_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            txt_stepsize.Text = slider_stepsize.Value.ToString();

            //Adjust padding appropriately with stepsize
            //Find proportion of max first
            //if (slider_edgepad != null)
            //{
            //    double edgeratio = Convert.ToDouble(txt_edgepadding.Text)/slider_edgepad.Value;
            //    double sideratio = Convert.ToDouble(txt_sidepad.Text)/slider_sidepad.Value;
            //    slider_edgepad.Minimum = 1.0;
            //    slider_sidepad.Minimum = 1.0;
            //    //slider_edgepad.Maximum = slider_stepsize.Value/2;
            //    //slider_sidepad.Maximum = slider_stepsize.Value/2;
            //    slider_edgepad.Maximum = 20;
            //    slider_sidepad.Maximum = 20;

            //    if (edgeratio <= 1.0 && edgeratio > 0)
            //        slider_edgepad.Value = (edgeratio*slider_edgepad.Maximum);
            //    else
            //        slider_edgepad.Value = slider_edgepad.Maximum/2;

            //    if (sideratio <= 1.0 && sideratio > 0)
            //        slider_sidepad.Value = (sideratio*slider_sidepad.Maximum);
            //    else
            //        slider_sidepad.Value = slider_sidepad.Maximum/2;

            //    txt_sidepad.Text = slider_sidepad.Value.ToString();
            //    txt_edgepadding.Text = slider_edgepad.Value.ToString();
            //}
            if (PS != null)
            {
                if (PS.RasterPaths != null && !PlanOptimized)
                {
                    UpdateCurrentSlicePaths();
                    RefreshDataGrid();
                }
                else if (PlanOptimized)
                {
                    redwarn_lbl.Content = "Plan already exists! Please click RePlan to see changes";
                    redwarn_lbl.IsEnabled = true;
                }
            }
        }

        private void slider_slicethickness_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            txt_slicethickness.Text = slider_slicethickness.Value.ToString();
            /*if (PS != null)
            {
                if (!PlanOptimized)
                {
                    txt_rasterwidth.IsEnabled = false;
                    txt_stepsize.IsEnabled = false;
                    slider_rasterwidth.IsEnabled = false;
                    slider_stepsize.IsEnabled = false;
                    //recalc_btn.IsEnabled = true;
                    redwarn_lbl.IsEnabled = true;
                    redwarn_lbl.Content = "Please click re-plan to see new slice thickness.";
                    UpdateCurrentSlicePaths();
                    RefreshDataGrid();
                }
                else if (PlanOptimized)
                {
                    redwarn_lbl.Content = "Plan already exists! Please click RePlan to see changes";
                    redwarn_lbl.IsEnabled = true;
                }
            }*/
            //recalc_btn.IsEnabled = true;
        }

        private RasterPath GetCurrentSliceRP()
        {
            if (PS.RasterPaths != null)
                return (RasterPath)PS.RasterPaths[GetCurrentSlice()];
            else
                return null;
        }

       
        private void ShiftRasterLines(int d)
        {
            RasterPath rp = GetCurrentSliceRP();
            int[] bounds = rp.boundaries;
            //First check to see if a new line needs to be added
            int start = rp.boundaries[0];
            int end = rp.boundaries[1];
            

            for (int i = 0; i < rp.Lines.GetLength(0); i++)
            {
                
                /*if (rp.Lines[i] + d >= end)
                    rp.Lines[i] = start + (d - 1);
                else if (rp.Lines[i] + d <= start)
                    rp.Lines[i] = end + (d + 1);
                else*/
                rp.Lines[i] = rp.Lines[i] + d;
                //rp.Lines.ToArray();
            }
            if (rp.Lines[0] - start > RasterPath.RasterWidth)
            {
                rp.Lines[rp.Lines.GetLength(0) - 1] = rp.Lines[0] - RasterPath.RasterWidth;
                int[] longer = new int[rp.Lines.GetLength(0) + 1];
                longer[0] = rp.Lines[0] - RasterPath.RasterWidth;
                for (int i = 0; i < rp.Lines.GetLength(0); i++)
                    longer[i + 1] = rp.Lines[i];
                rp.Lines = (int[]) longer.Clone();
                Array.Sort(rp.Lines);

            }
            if (end - rp.Lines[rp.Lines.GetLength(0) - 1] > RasterPath.RasterWidth)
            {
                int[] longer = new int[rp.Lines.GetLength(0)+1];
                longer[longer.GetLength(0) - 1] = rp.Lines[rp.Lines.GetLength(0) - 1] + RasterPath.RasterWidth-1;
                for (int i = 0; i < rp.Lines.GetLength(0); i++)
                    longer[i] = rp.Lines[i];
                rp.Lines = (int[])longer.Clone();
                Array.Sort(rp.Lines);
            }
            if (rp.Lines[0] < start)
            {
                int[] shrink = new int[rp.Lines.GetLength(0) - 1];
                for (int i = 0; i < shrink.GetLength(0); i++)
                    shrink[i] = rp.Lines[i + 1];
                rp.Lines = (int[])shrink.Clone();
            }
            if (rp.Lines[rp.Lines.GetLength(0) - 1] >= end)
            {   
                int[] shrink = new int[rp.Lines.GetLength(0)-1];
                for (int i = 0; i < shrink.GetLength(0); i++)
                    shrink[i] = rp.Lines[i];
                rp.Lines = (int[])shrink.Clone();
            }

            Array.Sort(rp.Lines);
            rp.FindAllShotPoints();
            rp.Calculate2DDoseSpace(DK.midplane, null);
            rp.CreateSliceInfo();
            DisplayPlan();
        }

        private void UpdateCurrentSlicePaths()
        {
            RasterPath rp = GetCurrentSliceRP();
            int stepsize = Convert.ToInt16(txt_stepsize.Text);
            int rasterwidth = (int)Math.Round(Convert.ToDouble(txt_rasterwidth.Text));
            int edgepad = (int)Math.Round(Convert.ToDouble(txt_edgepadding.Text));
            int sidepad = Convert.ToInt16(txt_sidepad.Text);
            rp.ChangeParamsUpdatePoints(stepsize, rasterwidth, edgepad, sidepad);
            rp.Calculate2DDoseSpace(DK.midplane, null);
            rp.CreateSliceInfo();
            DisplayPlan();
        }

        private void plan_dose_btn_Click(object sender, RoutedEventArgs e)
        {
            //RasterPath rp = GetCurrentSliceRP();
            //rp.Calculate2DDoseSpace(DK.midplane);
            //UpdateCurrentSlicePaths();
        }

        private void txt_rasterwidth_TextChanged(object sender, TextChangedEventArgs e)
        {
        }

        private void txt_stepsize_TextChanged(object sender, TextChangedEventArgs e)
        {
        }

        private void txt_slicethickness_TextChanged(object sender, TextChangedEventArgs e)
        {
        }

        private void recalc_btn_Click(object sender, RoutedEventArgs e)
        {
            if (PS != null)
            {
                PS.RecalculateSlices(Convert.ToInt16(txt_slicethickness.Text), Convert.ToInt16(txt_slicethickness.Text));
                TurnOnSliders();
                //recalc_btn.IsEnabled = false;
                UpdateCurrentSlicePaths();
                DisplayPlan();
            }
        }

        private void TurnOnSliders()
        {
            txt_rasterwidth.IsEnabled = true;
            txt_stepsize.IsEnabled = true;
            slider_rasterwidth.IsEnabled = true;
            slider_stepsize.IsEnabled = true;
            slider_slicethickness.IsEnabled = true;
        }

        private void LoadDose_menu_Click(object sender, RoutedEventArgs e)
        {
            var ofd = new OpenFileDialog();
            if (ofd.ShowDialog() != false)
            {
                DK = new DoseKernel(ofd.FileName);
                plan_dpRB.IsEnabled = true;
                IsDoseLoaded = true;
                RasterPath.doseN = DoseKernel.N;
                AddDoseLoadedToListBox();
            }
        }

        private void AddDoseLoadedToListBox()
        {
            listBox1.Items.Add(DK.DKI.Info);
        }

        private void AddStructureLoadedToListBox()
        {
            listBox1.Items.Add(SS.SI.Info);
        }

        private void AddDICOMLoadedToListBox()
        {
            listBox1.Items.Add(set.Dinfo.Info);
        }

        private void plan_dpRB_Checked(object sender, RoutedEventArgs e)
        {
            if (HasPreviewBeenCreated)
            {
                for (int i = 0; i < PathSet.NumSlices; i++)
                    ((RasterPath) PS.RasterPaths[i]).Calculate2DDoseSpace(DK.midplane, null);
            }
            else
                UpdateCurrentSlicePaths();
            //UpdateCurrentSlicePaths();
            RefreshDataGrid();
            DisplayPlan();
        }

        private void SilenceTrackingLabels()
        {
            planX_lbl.Content = "";
            planY_lbl.Content = "";
            planZ_lbl.Content = "";
        }

        private void UpdateTrackingLabels(MouseEventArgs e)
        {
            double aspectx = (StructureSet.BIG_dim[0]/plan_imgbox.Width);
            double aspecty = (StructureSet.BIG_dim[1]/plan_imgbox.Height);
            planX_lbl.Content = "X: " + Math.Round(e.GetPosition(plan_imgbox).X * aspectx,2);
            planY_lbl.Content = "Y: " + Math.Round(e.GetPosition(plan_imgbox).Y * aspecty,2);
            if (PS != null)
                if (PS.SlicePositions != null)
                planZ_lbl.Content = "Z: " + PS.SlicePositions[GetCurrentSlice()];
            else
                planZ_lbl.Content = "Z: unknown";
        }

        private void UpdateSliceLabels()
        {
            planlines_lbl.Content = "Lines: " + (GetCurrentSliceRP()).NumOfLines;
            planshots_lbl.Content = "Shots: " + (GetCurrentSliceRP()).NumOfLines;
            planslice_lbl.Content = "Slice #: " + (GetCurrentSlice() + 1);
            plansize_lbl.Content = "Size: " + StructureSet.BIG_dim[0] + " x " + StructureSet.BIG_dim[1];
            slider_edgepad.Value = (GetCurrentSliceRP()).LineEdgePadding;
            slider_sidepad.Value = (GetCurrentSliceRP()).LineSidePadding;
            txt_edgepadding.Text = slider_edgepad.Value.ToString();
            txt_sidepad.Text = slider_sidepad.Value.ToString();
        }

        private void PlanViewtab_GotFocus(object sender, RoutedEventArgs e)
        {
            //CreatePaths();
        }

        private void plan_btn_Click(object sender, RoutedEventArgs e)
        {
            plantimer = new Stopwatch();
            plantimer.Start();

            PathSet.StepSize = Convert.ToInt16(txt_stepsize.Text);
            PathSet.RasterWidth = Convert.ToInt16(txt_rasterwidth.Text);
            PathSet.LineSidePadding = Convert.ToInt16(txt_edgepadding.Text);
            PathSet.LineEdgePadding = Convert.ToInt16(txt_edgepadding.Text);
            

            if (PlanOptimized)
            {
                redwarn_lbl.Content = "";
                PlanOptimized = false;
                save_plan_btn.IsEnabled = false;
                export_shots_btn.IsEnabled = false;
                listBox2.Items.Clear();
                dataGrid1.Items.Clear();
            }
            CreatePaths();
            if (IsPathSetCreated)
                Opt_btn.IsEnabled = true;
            //AdjustCursorSize();
            SetWorkingDirectory();
            Debug.WriteLine("Current Plantimer: " + plantimer.Elapsed);
            tabControl1.SelectedIndex = 2;
            //plantimer.Stop();

        }

        private void UpdateTextBlock2(string update)
        {
            textBlock2.Text = update;
        }


        private void Opt_btn_Click(object sender, RoutedEventArgs e)
        {
            UpdateStatusBar("Running optimization...this may take some time.");
            SetWorkingDirectory();
            Analysis.SessionName = NameTxtBox.Text;
            if (PS.ShotsWeighted == false)
                PS.Ps1ShotOptimizeWorker.RunWorkerAsync();
            else
            {
                //Calculate the slicedoses so far...
                PS.CreateDoseMatrix(DK, PS.folderpath); //<- This calls second background worker "PS_Two"
            }
            //UpdateTextBlock2("Optimizing shot weights...please wait.");
        }

        private void SetWorkingDirectory()
        {
            string path;
            if (PathSet.ActiveDirectory != null)
            {
                PS.folderpath = PathSet.ActiveDirectory;
                PathSet.SessionName = NameTxtBox.Text;
            }
            else
            {
                path = LoadConfig();
                PathSet.ActiveDirectory = path;
            }
        }

        private string LoadConfig()
        {
            string path = Directory.GetCurrentDirectory();
            path = Path.Combine(path, "config.ini");
            if (File.Exists(path))
            {
                using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read))
                using (var br = new StreamReader(fs))
                {
                    string s = br.ReadLine();
                    PathSet.ActiveDirectory = s;
                    return s;
                }
            }
            else
            {
                string _folderName = "c:\\dinoch";

                _folderName = (Directory.Exists(_folderName)) ? _folderName : "";
                var dlg1 = new FolderBrowserDialogEx
                               {
                                   Description = "Select a folder for temporary files and dosefile:",
                                   ShowNewFolderButton = true,
                                   ShowEditBox = true,
                                   //NewStyle = false,
                                   SelectedPath = _folderName,
                                   ShowFullPathInEditBox = false,
                               };
                dlg1.RootFolder = Environment.SpecialFolder.MyComputer;

                DialogResult result = dlg1.ShowDialog();

                if (result == System.Windows.Forms.DialogResult.OK)
                {
                    _folderName = dlg1.SelectedPath;
                }
                PathSet.ActiveDirectory = _folderName;
                string configfile = Path.Combine(_folderName, "config.ini");
                string configpath = Path.Combine(Directory.GetCurrentDirectory(), "config.ini");
                using (var fs = new FileStream(configfile, FileMode.Create, FileAccess.Write))
                using (var bw = new StreamWriter(fs))
                {
                    bw.Write(_folderName);
                }
                return _folderName;
            }
        }

        private void SetWorkingDirectory(string fullpath)
        {
            PathSet.ActiveDirectory = fullpath;
        }

        private void dataGrid1_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (dataGrid1.SelectedIndex >= 0)
                PopulateShotList((RasterPath) PS.RasterPaths[dataGrid1.SelectedIndex]);
        }

        private void PopulateShotList(RasterPath rp)
        {
            listBox2.Items.Clear();
            //PointF[] pf = rp.ReturnSinglePoints();
            if (rp.shot_points != null)
            {
                foreach (PointF[] pf in rp.shot_points)
                {
                    string s = "X: " + pf[0].X + " Y: [" + pf[0].Y;
                    if (pf.GetLength(0) > 1)
                    {
                        for (int i = 1; i < pf.GetLength(0); i++)
                        {
                            s += ", " + (pf[i]).Y;
                        }
                    }
                    s += "]";
                    listBox2.Items.Add(s);
                }
                
            }
        }

        private void listBox2_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            /*if (dataGrid1.SelectedIndex >= 0)
            {
                var rp = (RasterPath) PS.RasterPaths[dataGrid1.SelectedIndex];
                if (listBox2.SelectedIndex >= 0)
                {
                    PointF p = (rp.ReturnSinglePoints())[listBox2.SelectedIndex];
                    if (tabControl1.SelectedIndex == 2)
                        DisplaySingleShotDS(p);
                    else
                    {
                        tabControl1.SelectedIndex = 3;
                        DisplaySingleShot(p, rp.ShotWeights[listBox2.SelectedIndex]);
                    }
                }
            }*/
        }
        private void treadmill_right_btn_Click(object sender, RoutedEventArgs e)
        {
            ShiftShotLocations(1);
        }

        public void ShiftShotLocations(int shift)
        {
            if (listBox2.SelectedIndex >= 0)
            {
                var rp = GetCurrentSliceRP();
                PointF[] line = rp.shot_points[listBox2.SelectedIndex];
                int[] linebounds = rp.LineBoundaries(rp.Lines[listBox2.SelectedIndex]);
                PointF[] newline = TreadmillShotLocations(line, linebounds, shift);
                rp.shot_points[listBox2.SelectedIndex] = newline;
                if (newline.GetLength(0) != line.GetLength(0))
                {
                    int diff = newline.GetLength(0) - line.GetLength(0);
                    rp.NumOfShots += diff;
                    //Need to re-init weight array
                }
                //rp.FindAllShotPoints();
                rp.Calculate2DDoseSpace(DK.midplane, rp.InitWeightArray());
                rp.CreateSliceInfo();
                DisplayPlan();
            }
        }

        private void treadmill_left_btn_Click(object sender, RoutedEventArgs e)
        {
            ShiftShotLocations(-1);
        }

        private PointF[] TreadmillShotLocations(PointF[] line, int[] bounds, int shift)
        {
            float loc = (int)Math.Round(line[0].X);
            int edgepad = Convert.ToInt16(txt_edgepadding.Text);
            int rasterwidth = Convert.ToInt16(txt_rasterwidth.Text);
            PointF[] shiftline = new PointF[line.GetLength(0)];
            int[] shifty = new int[line.GetLength(0)];
            for (int i = 0; i < line.GetLength(0); i++)
            {
                int y = (int)line[i].Y + shift;
                shifty[i] = y;
            }
            if (shifty[0] - bounds[0] >= rasterwidth)
            {
                int[] newshifty = new int[shifty.GetLength(0) + 1];
                newshifty[0] = shifty[0] - rasterwidth;
                for (int i = 0; i < shifty.GetLength(0); i++)
                    newshifty[i + 1] = shifty[0];
                shifty = (int[]) newshifty.Clone();
                Array.Sort(shifty);
            }
            if (bounds[1] - shifty[shifty.GetLength(0)-1] >= rasterwidth)
            {
                int[] newshifty = new int[shifty.GetLength(0) + 1];
                newshifty[newshifty.GetLength(0) - 1] = shifty[shifty.GetLength(0) - 1] + rasterwidth;
                for (int i = 0; i < shifty.GetLength(0); i++)
                    newshifty[i] = shifty[i];
                shifty = (int[]) newshifty.Clone();
                Array.Sort(shifty);
            }
            if (shifty[0] < bounds[0])
            {
                int[] shrink = new int[shifty.GetLength(0) - 1];
                for (int i = 0; i < shrink.GetLength(0); i++)
                    shrink[i] = shifty[i+1];
                shifty = (int[]) shrink.Clone();
            }
            if (shifty[shifty.GetLength(0)-1] >= bounds[1])
            {
                int[] shrink = new int[shifty.GetLength(0) - 1];
                for (int i = 0; i < shrink.GetLength(0); i++)
                {
                    shrink[i] = shifty[i];
                }
                shifty = (int[]) shrink.Clone();
            }

            shiftline = new PointF[shifty.GetLength(0)];
            Array.Sort(shifty);
            for (int i = 0; i < shiftline.GetLength(0); i++)
            {
                shiftline[i] = new PointF(loc, shifty[i]);
            }

            return shiftline;
        }



        public ArrayList Line;
        public bool _isDrawing = false;

        private void plan_imgbox_MouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            //DrawCirclePoints(e.GetPosition(plan_imgbox).X, e.GetPosition(plan_imgbox).Y);
            _isDrawing = true;
            Line = new ArrayList();
            System.Windows.Point p = e.GetPosition(plan_imgbox);
            Line.Add(p);
        }

        public float[,] modslice;

        private void plan_imgbox_MouseLeftButtonUp(object sender, MouseButtonEventArgs e)
        {
            _isDrawing = false;
            ModifyPrioritySlice(modslice);

        }

        private void plan_imgbox_MouseMove(object sender, MouseEventArgs e)
        {
            if (_isDrawing)
            {
                System.Windows.Point p = e.GetPosition(plan_imgbox);
                Line.Add(p);
                Line.TrimToSize();
                modslice = DrawCurrentLine(Line);
            }
            UpdateTrackingLabels(e);
        }

        private void plan_imgbox_MouseEnter(object sender, MouseEventArgs e)
        {
            plan_imgbox.MouseMove += plan_imgbox_MouseMove_1;
            cursor_ellipse.IsEnabled = true;
            UpdateTrackingLabels(e);
        }

        private void plan_imgbox_MouseLeave(object sender, MouseEventArgs e)
        {
            plan_imgbox.ReleaseMouseCapture();
            plan_imgbox.MouseMove -= plan_imgbox_MouseMove_1;
            cursor_ellipse.IsEnabled = false;

            SilenceTrackingLabels();
        }

        private void plan_imgbox_MouseMove_1(object sender, MouseEventArgs e)
        {
            System.Windows.Point p = e.GetPosition(plan_imgbox);
            var aspectPoint = new System.Windows.Point(p.X*Plan_aspectMultiplier_x, p.Y*Plan_aspectMultiplier_y);
            var halfpoint = (int) ((15*Plan_aspectMultiplier_x - 1)/2);
            cursor_ellipse.Arrange(new Rect(p.X - halfpoint, p.Y - halfpoint, (int) (15*Plan_aspectMultiplier_x),
                                            (int) (15*Plan_aspectMultiplier_y)));
            if (e.LeftButton == MouseButtonState.Pressed)
                DrawCirclePoints(e.GetPosition(plan_imgbox).X, e.GetPosition(plan_imgbox).Y);
            //canvas1.Arrange(new Rect(p.X, p.Y, 15, 15));
            UpdateTrackingLabels(e);
            //if (IsMouseCaptured)
            //    ForceCursor = true;
        }

        private void plan_imgbox_IsMouseDirectlyOverChanged(object sender, DependencyPropertyChangedEventArgs e)
        {
            if (IsMouseDirectlyOver)
                Mouse.Capture(plan_imgbox, CaptureMode.Element);
            else if (!IsMouseDirectlyOver)
                plan_imgbox.ReleaseMouseCapture();
        }

        private float[,] DrawCurrentLine(ArrayList line)
        {
            RasterPath rp = GetCurrentSliceRP();
            float[,] temp_mask = Matrix.Zeroes(rp.slice.GetLength(0), rp.slice.GetLength(1));
            double aspect_x = rp.slice.GetLength(0) / plan_imgbox.Height;
            double aspect_y = rp.slice.GetLength(1) / plan_imgbox.Width;
            unsafe
            {
                foreach (Point p in line)
                {
                    foreach (PointF pf in circle)
                    {
                        var pBackBuffer = (int) wb_Plan.BackBuffer;
                        var rel_y = (int) pf.Y;
                        var rel_x = (int) pf.X;
                        var global_x = ((int) (p.X*aspect_x) - CursorRadius) + rel_x;
                        var global_y = ((int) (p.Y*aspect_y) - CursorRadius) + rel_y;

                        if (temp_mask[global_x, global_y] == 1.0f)
                            continue;
                        else
                        {
                            temp_mask[global_x, global_y] = 1;
                            pBackBuffer += global_y * wb_Plan.BackBufferStride;
                            pBackBuffer += global_x * 4;
                            
                            int color_data = 0 << 16; // R
                            color_data |= 255 << 8; // G
                            color_data |= 0 << 0; // B

                            // Assign the color data to the pixel.
                            *((int*)pBackBuffer) = color_data;
                        }
                    }
                }
                try
                {
                    // Specify the area of the bitmap that changed.
                    wb_Plan.AddDirtyRect(new Int32Rect(0, 0, (int)wb_Plan.Width, (int)wb_Plan.Height));
                }
                catch (Exception ex)
                {
                    string s = ex.ToString();
                }
                finally
                {
                    wb_Plan.Unlock();
                }
            }
            return temp_mask;

        }

        private void ModifyPrioritySlice(float[,] temp_mask)
        {
            RasterPath rp = GetCurrentSliceRP();
            rp.mask = (float[,]) temp_mask.Clone();
            DisplayPlan();
        }


        private void DrawCirclePoints(double x, double y)
        {
            double Plan_aspectMultiplier = wb_Plan.PixelHeight / plan_imgbox.Height;
            //AdjustCursorSize();
            int i = ((int)(x * Plan_aspectMultiplier)) - CursorRadius;
            int j = ((int)(y * Plan_aspectMultiplier)) - CursorRadius;

            ResetMask();
            //i and j is the single center of the circle, adjusted for the aspect ratio.
            //u and v refer to each of the individual pixels in the circle cursor mask.
            if (circle != null)
                foreach (PointF p in circle)
                {
                    var u = (int)p.Y;
                    var v = (int)p.X;
                    mask[i + u, j + v] = 1;
                }

            wb_Plan.Lock();

            unsafe
            {
                foreach (PointF p in circle)
                {
                    // Get a pointer to the back buffer.
                    var pBackBuffer = (int)wb_Plan.BackBuffer;
                    var u = (int)p.Y;
                    var v = (int)p.X;
                    // Find the address of the pixel to draw.
                    pBackBuffer += (j + u) * wb_Plan.BackBufferStride;
                    pBackBuffer += (i + v) * 4;
                    // Compute the pixel's color.

                    int color_data = 255 << 16; // R
                    color_data |= 0 << 8; // G
                    color_data |= 0 << 0; // B

                    // Assign the color data to the pixel.
                    *((int*)pBackBuffer) = color_data;
                }
                try
                {
                    // Specify the area of the bitmap that changed.
                    wb_Plan.AddDirtyRect(new Int32Rect(0, 0, (int)wb_Plan.Width, (int)wb_Plan.Height));
                }
                catch (Exception ex)
                {
                    string s = ex.ToString();
                }
                finally
                {
                    wb_Plan.Unlock();
                }
            }
        }


        private void canvas1_MouseEnter(object sender, MouseEventArgs e)
        {
        }

        private void canvas1_MouseLeave(object sender, MouseEventArgs e)
        {
        }

        private void canvas1_MouseMove(object sender, MouseEventArgs e)
        {
        }

        private void cursor_ellipse_MouseEnter(object sender, MouseEventArgs e)
        {
        }

        private void cursor_ellipse_MouseMove(object sender, MouseEventArgs e)
        {
        }

        private void AddAnalysisTab()
        {
            var analysis = new TabItem();
            analysis.Header = "Analysis";
            tabControl1.Items.Add(analysis);
        }

        private void canvas1_IsMouseDirectlyOverChanged(object sender, DependencyPropertyChangedEventArgs e)
        {
        }

        private void tabControl1_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
        }

        private void Select4mmDefault_Click(object sender, RoutedEventArgs e)
        {
            Load4mmDefault();
        }

        private void Select8mmDefault_Click(object sender, RoutedEventArgs e)
        {
            Load8mmDefault();
        }

        private void Load4mmDefault()
        {
            DK = new DoseKernel(4);
            plan_dpRB.IsEnabled = true;
            IsDoseLoaded = true;
            RasterPath.doseN = DoseKernel.N;
            RasterPath.N = DoseKernel.N;
            AddDoseLoadedToListBox();

            //Set parameters to '4mm' defaults
            slider_rasterwidth.Value = 20;
            slider_slicethickness.Value = 20;
            slider_stepsize.Value = 20;
            txt_rasterwidth.Text = "20";
            txt_slicethickness.Text = "20";
            txt_stepsize.Text = "20";
            txt_edgepadding.Text = "10";
            slider_edgepad.Maximum = 20;
            slider_sidepad.Minimum = 1.0;
            slider_edgepad.Value = 10;
            txt_sidepad.Text = "10";
            slider_sidepad.Maximum = 20;
            slider_sidepad.Minimum = 1.0;
            slider_sidepad.Value = 10;
            RasterPath.ComparisonKernelSize = 60;
            RasterPath.ShotSize = 4;
        }

        private void Load8mmDefault()
        {
            //Load 8mm default kernel file
            DK = new DoseKernel(8);
            plan_dpRB.IsEnabled = true;
            IsDoseLoaded = true;
            RasterPath.doseN = DoseKernel.N;
            RasterPath.N = DoseKernel.N;


            RasterPath.ComparisonKernelSize = 100;
            AddDoseLoadedToListBox();

            //Set parameters to default values for '8mm'
            slider_rasterwidth.Value = 36;
            slider_slicethickness.Value = 36;
            slider_stepsize.Value = 38;
            txt_rasterwidth.Text = "36";
            txt_slicethickness.Text = "36";
            txt_stepsize.Text = "38";

            txt_edgepadding.Text = "20";
            slider_edgepad.Maximum = 40;
            slider_sidepad.Minimum = 1.0;
            slider_edgepad.Value = 20;
            txt_sidepad.Text = "20";
            slider_sidepad.Maximum = 40;
            slider_sidepad.Minimum = 1.0;
            slider_sidepad.Value = 20;
            RasterPath.ShotSize = 8;
        }

        private void CalcSaveDose_btn_Click(object sender, RoutedEventArgs e)
        {
            //CreateDoseMatrix(PS.ActiveDirectory);
        }

        private void CreateDoseMatrix(string folderpath)
        {
            //Moved this to the end of the Optimize step of the 2nd button click.
            //PS.CreateDoseMatrix(DK, path); //Calls the PS_InitialDose background worker

            Display2DFloat(PS.DoseSpace[PS.DoseSpace.GetLength(0)/2]);
            if (PS.DoseSpace != null)
            {
                Analysis_datagrid.IsEnabled = true;
            }
        }

        private void SetUpAnalysis()
        {
            PS.GetTotalShotNumber();
            Analysis_datagrid.IsEnabled = true;
            if (PS != null && SS != null)
                Analysis.RunAnalysis(PS, SS, PathSet.RxDose);
            else if (Analysis.ddf != null)
                Analysis.RunAnalysis(Analysis.ddf, SS.fj_Tumor, SS.fj_CS, 0.5);
            Analysis_datagrid.ItemsSource = GetAnalysisInfo();

            Analysis.AddAnalysisReport("Analysis.txt");
        }

        private List<AnalysisInfo> GetAnalysisInfo()
        {
            return Analysis.AIList;
        }

        private void LoadDS_menu_Click_1(object sender, RoutedEventArgs e)
        {
            var loadDS = new OpenFileDialog();
            loadDS.Title = "Select the Dosespace text file";
            loadDS.Filter = "Text files (.txt)|*.txt";
            if (loadDS.ShowDialog() != false)
            {
                PS = new PathSet(loadDS.FileName);
            }
        }

        private void run_analysis_btn_Click(object sender, RoutedEventArgs e)
        {
            SetUpAnalysis();
            FindBeamOnTime();
        }

        private void clearhistory_btn_Click(object sender, RoutedEventArgs e)
        {
            Analysis_datagrid.Items.Clear();
        }

        private void save_plan_btn_Click(object sender, RoutedEventArgs e)
        {
            //Calculate the slicedoses so far...
            UpdateStatusBar("Writing the dose to file...this will take a while.");
            UpdateProgressBar(10);
            PS.CreateDoseMatrix(DK, PS.folderpath);
            UpdateProgressBar(33);
            
            PS.WriteDoseSpaceToFile(string.Concat(SessionName, "_FinalDoseSpace.txt"));
            UpdateProgressBar(50);
            WriteParametersToFile();
            UpdateProgressBar(75);
            ExportShotList();
            UpdateProgressBar(100);
        }

        private void WriteParametersToFile()
        {
            string subfolder = PathSet.ActiveDirectory;
            string filename = string.Concat(SessionName, "_PlanParameters.txt");
            string path = Path.Combine(subfolder, filename);
            using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
            using (var bw = new StreamWriter(fs))
            {
                bw.WriteLine("Raster Width: " + RasterPath.RasterWidth);
                bw.WriteLine("Shot Step Size: " + RasterPath.StepSize);
                bw.WriteLine("Slice Thickness: " + PathSet.ST);
                bw.WriteLine("Prescription dose: " + PathSet.RxDose);
                bw.WriteLine("Tolerance dose: " + PathSet.ToleranceDose);
            }
        }



        public void ExportShotList()
        {
            string subfolder = PathSet.ActiveDirectory;
            string path = Path.Combine(subfolder, "ShotList.txt");
            using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
            using (var bw = new StreamWriter(fs))
            {
                for (int i = 0; i < PS.RasterPaths.Count; i++)
                {
                    double sliceweight = PS.SliceWeights[i];
                    var rp = (RasterPath) PS.RasterPaths[i];
                    for (int shot = 0; shot < rp.shots.GetLength(0); shot++)
                    {
                        double weight = rp.ShotWeights[shot]*sliceweight;
                        PointF original_p = rp.shots[shot];
                        PointF p = ConvertXYCoordinates_MillimeterRes(original_p);
                        float z = ConvertZCoordinate_MillimeterRes(PS.SlicePositions[i]);
                        string sep = ", ";
                        bw.WriteLine(p.X + sep + p.Y + sep + z + sep + Math.Round(weight, 2) + ";");
                    }
                }
            }
        }

        public PointF ConvertXYCoordinates_MillimeterRes(PointF p)
        {
            float centerx = StructureSet.centralpoint[0];
            float centery = StructureSet.centralpoint[1];
            
            float x = ((p.X - centerx)/4.0f) + 100.0f;
            float y = ((p.Y - centery) / 4.0f) + 100.0f;

            return new PointF(x, y);
        }

        public float ConvertZCoordinate_MillimeterRes(int z)
        {
            float centerz = StructureSet.centralpoint[2];
            float newz = ((z - centerz)/4.0f) + 100.0f;
            return newz;
        }

        public double BeamOnTime(double doserate, double MaxGy)
        {
            /*
             * To calculate beam-on time, take max gray level (ex: 40 Gy). Find out how many minutes it takes to put down 
             * the max level by dividing by the doserate at the center of the isocenter. That number of minutes corresponds
             * to a dose of 1.0, the max normalized dose. Then the beam-on time of subsequent shots should just be the weight x that value.
             */
            var weights = (double[]) PS.SliceWeights.Clone();
            double MinutesToMax = MaxGy/doserate;
            double BeamOnTime = 0;
            for (int i = 0; i < PS.RasterPaths.Count; i++)
            {
                double sliceweight = weights[i];
                var rp = (RasterPath) PS.RasterPaths[i];
                for (int shot = 0; shot < rp.shots.GetLength(0); shot++)
                {
                    double weight = rp.ShotWeights[shot]*sliceweight;
                    BeamOnTime += (weight*MinutesToMax);
                }
            }
            PS.FindDoseContributionsToReferencePoint(weights, doserate, MaxGy);
            return BeamOnTime;
        }

        private void export_shots_btn_Click(object sender, RoutedEventArgs e)
        {
            ExportShotList();
            UpdateStatusBar("Shots written to file.");
        }

        private void MenuItem_Click_1(object sender, RoutedEventArgs e)
        {
        }

        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
        }

        private void Button_Click_2(object sender, RoutedEventArgs e)
        {
            SetWorkingDirectory();
            Load4mmDefault();
            LoadDefaultSetup();
        }

        private void Set_Direc_Click(object sender, RoutedEventArgs e)
        {
            //SetWorkingDirectory();
            ChooseDirectory();
        }

        private void ChooseDirectory()
        {
            string _folderName = "c:\\dinoch";

            _folderName = (Directory.Exists(_folderName)) ? _folderName : "";
            var dlg1 = new FolderBrowserDialogEx
                           {
                               Description = "Select a folder for temporary files and dosefile:",
                               ShowNewFolderButton = true,
                               ShowEditBox = true,
                               //NewStyle = false,
                               SelectedPath = _folderName,
                               ShowFullPathInEditBox = false,
                           };
            dlg1.RootFolder = Environment.SpecialFolder.MyComputer;

            DialogResult result = dlg1.ShowDialog();

            if (result == System.Windows.Forms.DialogResult.OK)
            {
                _folderName = dlg1.SelectedPath;
            }
            PathSet.ActiveDirectory = _folderName;
            string configfile = Path.Combine(_folderName, "config.ini");
            string configpath = Path.Combine(Directory.GetCurrentDirectory(), "config.ini");
            using (var fs = new FileStream(configpath, FileMode.Create, FileAccess.Write))
            using (var bw = new StreamWriter(fs))
            {
                bw.Write(_folderName);
            }
        }


        private void ViewIso_chkbox_Checked(object sender, RoutedEventArgs e)
        {
            if (dicomdose_rb_btn.IsChecked == true)
            {
                Display2DFloat(ddf.Dose[GetCurrentSlice()]);
            }
            else if (plandose_rb_btn.IsChecked == true)
                Display2DFloat(PS.DoseSpace[GetCurrentSlice()]);
            else
                MessageBox.Show("Pick a radio button.");
        }

        private void GPU_chkbox_Checked(object sender, RoutedEventArgs e)
        {
            GPU.GPUenabled = true;
        }

        private void GPU_chkbox_Unchecked(object sender, RoutedEventArgs e)
        {
            GPU.GPUenabled = false;
        }

        private void LoadConfig_Menu_Click(object sender, RoutedEventArgs e)
        {
            var loadconfig = new OpenFileDialog();
            loadconfig.Title = "Select the config file";
            loadconfig.Filter = "Text files (.txt)|*.txt";
            if (loadconfig.ShowDialog() != false)
            {
                LOAD_CONFIG_FILE(loadconfig.FileName);
            }
        }

        private void Button_Click_3(object sender, RoutedEventArgs e)
        {
        }

        private void Load_dosespace_btn_Click(object sender, RoutedEventArgs e)
        {
            var loadconfig = new OpenFileDialog();
            loadconfig.Title = "Select a dosespace or DICOM-RT dose file";
            //float[] dosespace;
            if (loadconfig.ShowDialog() != false)
            {
                Debug.WriteLine(Path.GetExtension(loadconfig.FileName));
                if (Path.GetExtension(loadconfig.FileName) == ".dcm")
                {
                    LoadDICOMdose(loadconfig.FileName);
                    StatusTxtBox.Text = "DICOM dose successfully loaded";
                }
                else if (Path.GetExtension(loadconfig.FileName) == ".txt")
                {
                    LoadDoseSpace(loadconfig.FileName);
                }
            }
        }

        private void LoadDoseSpace(string p)
        {
            Analysis.ddf = new DICOMDoseFile(p, false);
            dicomdose_rb_btn.IsEnabled = true;
        }

        private void LoadDICOMdose(string p)
        {
            Analysis.ddf = new DICOMDoseFile(p, true);
            dicomdose_rb_btn.IsEnabled = true;
        }

        private void RadioButton_Checked_1(object sender, RoutedEventArgs e)
        {
            origdose_rb_btn.IsEnabled = true;
            newdose_rb_btn.IsEnabled = true;
            origdose_rb_btn.IsChecked = true;
            newdose_rb_btn.IsChecked = false;
            if (ddf.OriginalDose != null)
            {
                slider2.Minimum = 0;
                slider2.Maximum = ddf.OriginalDose.GetLength(0) - 1;
                slider2.Value = (slider2.Maximum/2);
                if (origdose_rb_btn.IsChecked == true)
                {
                    if (ShowStruct_chkbox.IsChecked == true)
                        DisplayDoseWithStruct(ddf.OriginalDose[GetCurrentSlice()]);
                    else
                        Display2DFloat(ddf.OriginalDose[GetCurrentSlice()]);
                }
                else if (newdose_rb_btn.IsChecked == true)
                {
                    if (ShowStruct_chkbox.IsChecked == true)
                        DisplayDoseWithStruct(ddf.Dose[GetCurrentSlice()]);
                    else
                        Display2DFloat(ddf.Dose[GetCurrentSlice()]);
                }
            }
        }

        public void DisplayDose(float[,] d)
        {
            DisplayImageSizeInfo(d.GetLength(0), d.GetLength(1));
            //d = Matrix.Normalize(d);
            //Display2DFloat(d);
            
        }

        

        private void plandose_rb_btn_Copy_Checked(object sender, RoutedEventArgs e)
        {
            dicomdose_rb_btn.IsChecked = false;
            slider2.Minimum = 0;
            if (PS != null)
            {
                slider2.Maximum = PS.DoseSpace.GetLength(0) - 1;
                Display2DFloat(PS.DoseSpace[GetCurrentSlice()]);
            }
        }

        private void SaveConfig_Menu_Click(object sender, RoutedEventArgs e)
        {
            WRITE_CONFIG_FILE("CONFIG_1.txt");
        }

        private void run_DDF_analysis_btn_Click(object sender, RoutedEventArgs e)
        {
            //Analysis.RunAnalysis(Analysis.ddf,StructureSet.originalTumor, 0.5);
            StatusTxtBox.Text = "Running Analysis...";
            Analysis.AnalyzeDICOMdose(Analysis.ddf, SS);
            StatusTxtBox.Text = "Analysis completed.";
            Analysis_datagrid.ItemsSource = GetAnalysisInfo();
        }

        private void origdose_rb_btn_Checked(object sender, RoutedEventArgs e)
        {
            newdose_rb_btn.IsChecked = false;
            slider2.Minimum = 0;
            slider2.Maximum = ddf.OriginalDose.GetLength(0) - 1;
            Display2DFloat(ddf.OriginalDose[GetCurrentSlice()]);
        }

        private void newdose_rb_btn_Checked(object sender, RoutedEventArgs e)
        {
            origdose_rb_btn.IsChecked = false;
            slider2.Minimum = 0;
            slider2.Maximum = ddf.Dose.GetLength(0) - 1;
            Display2DFloat(ddf.Dose[GetCurrentSlice()]);
        }

        private void Normalized_chkbox_Checked(object sender, RoutedEventArgs e)
        {
            if (Normalize == false)
                Normalize = true;
        }

        private void ShowStruct_chkbox_Unchecked(object sender, RoutedEventArgs e)
        {
            if (Normalize == true)
                Normalize = false;
        }

        private void Shift_down_btn_Click(object sender, RoutedEventArgs e)
        {
            ShiftRasterLines(-1);
        }

        private void Shift_up_btn_Click(object sender, RoutedEventArgs e)
        {
            ShiftRasterLines(1);
        }

        private void BeamOnTime_btn_Click(object sender, RoutedEventArgs e)
        {
            FindBeamOnTime();
        }

        private void FindBeamOnTime()
        {
            double doserate = Convert.ToDouble(DoseRate_txtbox.Text);
            Analysis.doserate = doserate;
            double maxdose = Convert.ToDouble(MaxDose_txtbox.Text);
            Analysis.maxdose = maxdose;
            double time = BeamOnTime(Convert.ToDouble(DoseRate_txtbox.Text), Convert.ToDouble(MaxDose_txtbox.Text));
            Analysis.beamontime = time;
            //Analysis.AddLineToReport("For max dose " + maxdose + " at " + doserate + " Gy/min, TOTAL BEAM-ON TIME ==> " +time);
            Time_lbl.Content = "" + Math.Round(time, 2) + " mins";
        }

        private void slider_sidepad_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            txt_sidepad.Text = slider_sidepad.Value.ToString();
            if (PS != null)
            {
                if (PS.RasterPaths != null && !PlanOptimized)
                {
                    //UpdateCurrentSlicePaths();
                    RefreshDataGrid();
                    RasterPath rp = GetCurrentSliceRP();
                    int[] bounds = rp.boundaries;
                    int start = rp.boundaries[0];
                    int end = rp.boundaries[1];
                    rp.Lines = rp.LineSpacer(start, end, RasterPath.StepSize/2, rp.WhichSlice);
                    //rp.Lines.ToArray();


                    rp.FindAllShotPoints();
                    rp.Calculate2DDoseSpace(DK.midplane, null);
                    rp.CreateSliceInfo();
                    DisplayPlan();
                }
                else if (PlanOptimized)
                {
                    redwarn_lbl.Content = "Plan already exists! Please click RePlan to see changes";
                    redwarn_lbl.IsEnabled = true;
                }
            }
        }

        private void slider_edgepad_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            txt_edgepadding.Text = slider_edgepad.Value.ToString();
            if (PS != null)
            {
                if (PS.RasterPaths != null && !PlanOptimized)
                {
                    UpdateCurrentSlicePaths();
                    RefreshDataGrid();
                }
                else if (PlanOptimized)
                {
                    redwarn_lbl.Content = "Plan already exists! Please click RePlan to see changes";
                    redwarn_lbl.IsEnabled = true;
                }
            }
        }

        private void SetRasterConfig_btn_Click(object sender, RoutedEventArgs e)
        {
            RasterPath rp = GetCurrentSliceRP();
            rp.LineSidePadding = Convert.ToInt16(txt_sidepad.Text);
            rp.LineEdgePadding = Convert.ToInt16(txt_edgepadding.Text);
            rp.FindAllShotPoints();
            rp.Calculate2DDoseSpace(DK.midplane, null);
            rp.CreateSliceInfo();
            DisplayPlan();
        }

        #region Pathset BackgroundWorkers

        private void RasterPath_SliceWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            var ShotWeights = (double[]) e.UserState;
        }


        private void PS_1_PathsetWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            UpdateProgressBar(e.ProgressPercentage);
            UpdateStatusBar("Calculating optimal shot weights and slice weights...please wait");
            //UpdateTextBlock("Progress: " + e.ProgressPercentage);
        }

        private void PS_1_PathsetWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (plantimer.IsRunning)
            {
                Debug.WriteLine("Current time after planning: " + plantimer.Elapsed);
            }
            UpdateStatusBar("Writing dose file...");
            UpdateTextBlock2("Finished shot weighting. Now writing dose, please wait...");
            RefreshDataGrid();

            foreach (RasterPath rp in PS.RasterPaths)
            {
                rp.Calculate2DDoseSpace(RasterPath.dosemidplane, null);
                rp.CreateSliceInfo();
            }


            //Auto view the first slice in dose form
            plan_dpRB.IsChecked = true;
            dataGrid1.SelectedIndex = 0;

            //Adjust controls to reflect the "Optimized" status.
            PlanOptimized = true;
            Plan_btn.Content = "Re-plan";
            Opt_btn.Content = "...Continue!";
            //Opt_btn.IsEnabled = false;
            //save_plan_btn.IsEnabled = false;
            //export_shots_btn.IsEnabled = true;
            //CalcSaveDose_btn.IsEnabled = true;
            DisplayPlan();
        }

        #endregion

        #region CreatingTestFiles

        public void CreateTumorObject(int radius)
        {
            N = 161;
            int length = radius + N;
            StructureSet.size = length;
            StructureSet.BIG_dim = new int[3] {length, length, length};

            var tumor = new float[length][];
            float[] slice;
            float dist;
            for (int z = 0; z < length; z++)
            {
                int x2;
                int y2;
                int z2;
                slice = new float[length*length];
                z2 = (int) Math.Pow(z - (length - 1)/2, 2);
                for (int i = 0; i < length; i++)
                {
                    x2 = (int) Math.Pow(i - (length - 1)/2, 2);
                    for (int j = 0; j < length; j++)
                    {
                        y2 = (int) Math.Pow(j - (length - 1)/2, 2);
                        dist = (float) Math.Sqrt(x2 + y2 + z2);
                        if (dist < radius/2)
                            slice[i + j*length] = 1;
                        else
                            slice[i + j*length] = 0;
                    }
                }
                tumor[z] = slice;
            }
            SS = new StructureSet(tumor);
            SS.SI.Size = new int[3] {length, length, length};
            SS.SI.CreateTestInfo(radius);
            SS.fj_Tumor = Convert1DJaggedto2DJagged(tumor, length, length);
            Plan_btn.IsEnabled = true;
        }

        public float[][,] Convert1DJaggedto2DJagged(float[][] f, int x, int y)
        {
            var F = new float[f.GetLength(0)][,];

            for (int k = 0; k < f.GetLength(0); k++)
            {
                var temp = new float[x,y];
                for (int j = 0; j < y; j++)
                    for (int i = 0; i < x; i++)
                    {
                        temp[i, j] = f[k][i + j*x];
                    }
                F[k] = temp;
            }
            return F;
        }

        #endregion

        #region Background Workers

        private void DICOMImageSet_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            double imgnum = e.ProgressPercentage;
            double percent = (imgnum/set.NumberOfImages)*100;
            UpdateStatusBar("Currently processing image " + imgnum + " of " + set.NumberOfImages + "...");
            UpdateProgressBar(percent);
        }

        private void DICOMImageSet_ImageWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            //dicom_radiobtn.IsEnabled = true;
            //dicom_radiobtn.IsChecked = true;
            UpdateStatusBar("Total: " + set.NumberOfImages + " frames loaded.");
            slider2.Minimum = 0;
            slider2.Maximum = set.f_imagearray.GetLength(0) - 1;
            slider2.Value = 0;
            ColormapTool(Properties.Resources.BWheatmap3);
            DICOM_imgbox.Stretch = Stretch.Uniform;
            //Bitmap temp = (Bitmap)DICOMImage.imagebitmaps[0];
            //DICOM_imgbox.Source = DoseHeatPoints.getBitmapSource((Bitmap)DICOMImage.imagebitmaps[0], temp.Width);
            var wb = new WriteableBitmap(256, 256, 96, 96, PixelFormats.Bgr32, null);
            DICOM_imgbox.Source = wb;
            IsDICOMLoaded = true;
            DisplayDICOM(set.f_imagearray[0]);
            AddDICOMLoadedToListBox();
        }

        #endregion

        #region UI Element Handlers

        private void DDS_imgbox_MouseRightButtonDown(object sender, MouseButtonEventArgs e)
        {
            ErasePixel(ref wb_DDS, DDS_imgbox, e);
        }

        private void DDS_imgbox_MouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            DrawPixel(ref wb_DDS, DDS_imgbox, e);
        }

        private void DS_imgbox_MouseWheel(object sender, MouseWheelEventArgs e)
        {
            throw new NotImplementedException();
        }

        private void DS_imgbox_MouseMove(object sender, MouseEventArgs e)
        {
            double aspectMultiplier;
            double[] img;
            if (plandose_rb_btn.IsChecked == true)
            {
                //Calculate an aspect multiplier based on the size of the matrix.
                double aspectmult = (DS_imgbox.Width/PS.DoseSpace[0].GetLength(0));

                img = new double[3];
                img[0] = Math.Round(e.GetPosition(this.DS_imgbox).X/aspectmult, 2);
                img[1] = Math.Round(e.GetPosition(this.DS_imgbox).Y/aspectmult, 2);
                img[2] = GetCurrentSlice();
                UpdateDSLabels(img);
                UpdateStatusWindow(img);
            }
            else if (dicomdose_rb_btn.IsChecked == true)
            {
                if (ddf.Dose != null)
                {
                    double aspectmultx = (DS_imgbox.Width/ddf.Dose[0].GetLength(0));
                    double aspectmulty = (DS_imgbox.Height/ddf.Dose[0].GetLength(1));
                    img = new double[3];
                    img[0] = Math.Round(e.GetPosition(this.DS_imgbox).X/aspectmultx, 2);
                    img[1] = Math.Round(e.GetPosition(this.DS_imgbox).Y/aspectmulty, 2);
                    img[2] = GetCurrentSlice();
                    UpdateDICOMDoseLabels(img);
                    UpdateStatusWindow(img);
                }
            }
        }

        private void UpdateDICOMDoseLabels(double[] point)
        {
            decimal xstart = Math.Round(DICOMDoseFile.doseoffset[0] + (decimal) point[0], 2);
            decimal ystart = Math.Round(DICOMDoseFile.doseoffset[1] + (decimal) point[1], 2);
            decimal zstart = Math.Round(DICOMDoseFile.doseoffset[2] + (decimal) point[2], 2);
            DS_x_lbl.Content = "X: " + xstart;
            DS_y_lbl.Content = "Y: " + ystart;
            DS_z_lbl.Content = "Z: " + zstart;
        }

        private void UpdateStatusWindow(double[] point)
        {
            string dicomoffset = "Not set yet.";
            string doseoffset = "Not set yet.";
            string structoffset = "Not set yet.";

            if (DICOMImageFile.GlobalOffset != null)
                dicomoffset = "DICOM ref. pos: " + Math.Round(DICOMImageFile.GlobalOffset[0], 1) + ", " +
                              Math.Round(DICOMImageFile.GlobalOffset[1], 1) + ", " +
                              Math.Round(DICOMImageFile.GlobalOffset[2], 1);
            if (DICOMDoseFile.doseoffset != null)
                doseoffset = "Dose ref. pos: " + Math.Round(DICOMDoseFile.doseoffset[0], 1) + ", " +
                             Math.Round(DICOMDoseFile.doseoffset[1], 1) + ", " +
                             Math.Round(DICOMDoseFile.doseoffset[2], 1);
            if (StructureSet.f_SSoffset != null)
                structoffset = "Struct ref pos: " + Math.Round(StructureSet.f_SSoffset[0], 1) + ", " +
                               Math.Round(StructureSet.f_SSoffset[1], 1) + ", " +
                               Math.Round(StructureSet.f_SSoffset[2], 1);

            string imgtxt = "Cursor: " + Math.Round(point[0], 1) + ", " + Math.Round(point[1], 1) + ", " +
                            Math.Round(point[2], 1);
            string dicomtxt = "DICOM: " + Math.Round(point[0], 1) + ", " + Math.Round(point[1], 1) + ", " +
                              Math.Round(point[2], 1);
            string chosenoffset = "";
            if (tabControl1.SelectedIndex == 2) //"DS" tab
                chosenoffset = doseoffset;
            else if (tabControl1.SelectedIndex == 1) //"Structure"
                chosenoffset = structoffset;
            else if (tabControl1.SelectedIndex == 0) //"DICOM"
                chosenoffset = dicomoffset;
            StatusTxtBox.Text = imgtxt + "\n" + dicomtxt + "\n" + "\n" + chosenoffset;
        }

        private void UpdateDSLabels(double[] point)
        {
            double xstart;
            double ystart;
            double zstart;
            if (DICOMDoseFile.doseoffset != null)
            {
                xstart = Math.Round((double) DICOMDoseFile.doseoffset[0] + point[0], 1);
                ystart = Math.Round((double) DICOMDoseFile.doseoffset[1] + point[1], 1);
                zstart = Math.Round((double) DICOMDoseFile.doseoffset[2] + point[2], 1);
            }
            else
            {
                xstart = point[0];
                ystart = point[1];
                zstart = point[2];
            }
            DS_x_lbl.Content = "X: " + point[0];
            DS_y_lbl.Content = "Y: " + point[1];
            DS_z_lbl.Content = "Z: " + point[2];
            UpdateStatusWindow(new double[3] {xstart, ystart, zstart});
        }

        private void UpdateDDSLabels(double[] point)
        {
            double xstart;
            double ystart;
            double zstart;
            if (StructureSet.f_global_xoffset != null)
            {
                xstart = Math.Round(StructureSet.f_global_xoffset + point[0], 1);
                ystart = Math.Round(StructureSet.f_global_yoffset + point[1], 1);
                zstart = Math.Round(StructureSet.f_global_zoffset + point[2], 1);
            }
            else
            {
                xstart = point[0];
                ystart = point[1];
                zstart = point[2];
            }
            DDS_x_lbl.Content = "X: " + xstart;
            DDS_y_lbl.Content = "Y: " + ystart;
            DDS_z_lbl.Content = "Z: " + zstart;
            UpdateStatusWindow(new double[3] {xstart, ystart, zstart});
        }

        private void DDS_imgbox_MouseWheel(object sender, MouseWheelEventArgs e)
        {
            throw new NotImplementedException();
        }

        private void DDS_imgbox_MouseMove(object sender, MouseEventArgs e)
        {
            double[] img;
            if (SS != null)
            {
                img = new double[3];
                double aspectmult = (SS.f_structurearray.GetLength(0)/DDS_imgbox.Width);
                img[0] = Math.Round(e.GetPosition(this.DDS_imgbox).X*aspectmult, 2);
                img[1] = Math.Round(e.GetPosition(this.DDS_imgbox).Y*aspectmult, 2);
                img[2] = Math.Round(slider2.Value, 2);
                UpdateDDSLabels(img);
            }
            //if (e.LeftButton == MouseButtonState.Pressed)
            //    DrawPixel(ref wb_DDS,DDS_imgbox, e);
            //else if (e.RightButton == MouseButtonState.Pressed)
            //    ErasePixel(ref wb_DDS, DDS_imgbox, e);
        }

        private void DICOM_imgbox_MouseWheel(object sender, MouseWheelEventArgs e)
        {
            double val = slider2.Value + e.Delta;
            if (val > slider2.Maximum)
                val = slider2.Maximum;
            else if (val < slider2.Minimum)
                val = slider2.Minimum;
            slider2.Value = val;
        }

        private void DICOM_imgbox_MouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            throw new NotImplementedException();
        }

        private void DICOM_imgbox_MouseMove(object sender, MouseEventArgs e)
        {
            throw new NotImplementedException();
        }

        private void DICOM_imgbox_MouseEnter(object sender, MouseEventArgs e)
        {
            tracking_label.IsEnabled = true;
            tracking_label.Content = "Image Pos: (" + Math.Round(e.GetPosition(this.DICOM_imgbox).X, 2) + ", " +
                                     Math.Round(e.GetPosition(this.DICOM_imgbox).Y, 2) + ")";
            if (LGKcoords.finished)
            {
                label1.IsEnabled = true;
                var x =
                    (decimal) Math.Round((e.GetPosition(this.DICOM_imgbox).X*DICOM_aspectMultiplier) - LGKcoords.X, 2);
                decimal y = 256 -
                            (decimal)
                            Math.Round((e.GetPosition(this.DICOM_imgbox).Y*DICOM_aspectMultiplier) + LGKcoords.Y, 2);
                //decimal z = (decimal)Math.Round((int)DICOMImage.img_zindex[(int)sliderbar.Value] - LGKcoords.Z, 2);
                var z =
                    (decimal)
                    Math.Round((int) set.ZIndexArray[(int) Math.Round(95*slider2.Value/slider2.Maximum)] - LGKcoords.Z,
                               2);
                label1.Content = "LGK Frame: <" + x + ", " + y + ", " + z + ">";
            }
        }

        private void DICOM_imgbox_MouseMove_1(object sender, MouseEventArgs e)
        {
            //Get absolute mouse position, multiply it by aspect multiplier (usually < 1, since the image is probably 256x256)
            var img = new decimal[3];
            var dicom = new decimal[3];

            img[0] = (decimal) Math.Round(e.GetPosition(this.DICOM_imgbox).X*DICOM_aspectMultiplier, 2);
            img[1] = (decimal) Math.Round(e.GetPosition(this.DICOM_imgbox).Y*DICOM_aspectMultiplier, 2);
            img[2] = (decimal) Math.Round(slider2.Value, 2);

            //The DICOM position is the mouse position plus the top left pixel of the DICOM reference coordinates.
            //The Z position is retrieved by the sliderbar value
            if (set != null)
            {
                dicom[0] = img[0] + Convert.ToDecimal(set.imagePosition[0]);
                dicom[1] = img[1] + Convert.ToDecimal(set.imagePosition[1]);
                dicom[2] = Convert.ToDecimal(set.imagePosition[2]) - img[2]*2;
                tracking_label.Content = "Image Pos: (" + img[0] + ", " + img[1] + ", " + img[2] + ")";
                UpdateStatusWindow((double[]) dicom.Clone());
            }
            if (LGKcoords.finished)
            {
                decimal[] lgk = LGKcoords.Image2LGKCoordinates(img);
                lgk[2] = (decimal) Math.Round(LGKcoords.getLGK_Z_forDICOM(slider2.Value), 2);
                label1.Content = "LGK Frame: <" + lgk[0] + ", " + lgk[1] + ", " + lgk[2] + ">";
                DICOMCoordLabel.Content = "DICOM Location: (" + (Math.Round(dicom[0], 2)) + ", " +
                                          (Math.Round(dicom[1], 2)) + ", " + (Math.Round(dicom[2], 2)) + ")";
            }
        }

        private void DICOM_imgbox_MouseLeave(object sender, MouseEventArgs e)
        {
            tracking_label.IsEnabled = false;
            label1.IsEnabled = false;
        }

        private void DICOM_imgbox_MouseWheel_1(object sender, MouseWheelEventArgs e)
        {
            double value = (int) slider2.Value + e.Delta;
            if (e.Delta > 0)
                value = slider2.Value + slider2.LargeChange;
            else if (e.Delta < 0)
                value = slider2.Value - slider2.LargeChange;
            if (value < slider2.Minimum)
                value = slider2.Minimum;
            else if (value > slider2.Maximum)
                value = slider2.Maximum;
            slider2.Value = value;
        }

        private void DICOM_imgbox_MouseLeftButtonDown_1(object sender, MouseButtonEventArgs e)
        {
            if (aligning_helper_label.IsEnabled)
            {
                LGKcoords.SetLGK_XYoffset(e.GetPosition(DICOM_imgbox).X*DICOM_aspectMultiplier,
                                          e.GetPosition(DICOM_imgbox).Y*DICOM_aspectMultiplier);
                CreateAlignmentMarkings(e.GetPosition(DICOM_imgbox).X*DICOM_aspectMultiplier,
                                        e.GetPosition(DICOM_imgbox).Y*DICOM_aspectMultiplier);
                aligning_helper_label.Content = "Line up left fiducial...";
                button1.Content = "Lined up";
            }
        }

        #endregion

        #region Display Methods

        private double FindDistanceBetweenPoints(Point p1, Point p2)
        {
            double distance = Math.Round(Math.Sqrt(Math.Pow((p2.X - p1.X), 2) + Math.Pow((p2.Y - p1.Y), 2)), 1);
            return distance;
        }

        private double FindDistanceBetweenPoints(PointF p1, PointF p2)
        {
            double distance = Math.Round(Math.Sqrt(Math.Pow((p2.X - p1.X), 2) + Math.Pow((p2.Y - p1.Y), 2)), 1);
            return distance;
        }

        public void DrawPixel(ref WriteableBitmap writeableBitmap, System.Windows.Controls.Image i, MouseEventArgs e)
        {
            var column = (int) e.GetPosition(i).X;
            var row = (int) e.GetPosition(i).Y;
            // Reserve the back buffer for updates.
            writeableBitmap.Lock();

            unsafe
            {
                // Get a pointer to the back buffer.
                var pBackBuffer = (int) writeableBitmap.BackBuffer;

                // Find the address of the pixel to draw.
                pBackBuffer += row*writeableBitmap.BackBufferStride;
                pBackBuffer += column*4;

                // Compute the pixel's color.
                int color_data = 255 << 16; // R
                color_data |= 128 << 8; // G
                color_data |= 255 << 0; // B

                // Assign the color data to the pixel.
                *((int*) pBackBuffer) = color_data;
            }

            // Specify the area of the bitmap that changed.
            writeableBitmap.AddDirtyRect(new Int32Rect(column, row, 1, 1));

            // Release the back buffer and make it available for display.
            writeableBitmap.Unlock();
        }

        private void CreateCirclePoints()
        {
            var c = new int[CursorRadius*2 + 1,CursorRadius*2 + 1];
            var points = new ArrayList();
            for (int y = 0; y < c.GetLength(0); y++)
                for (int x = 0; x < c.GetLength(1); x++)
                {
                    if (FindDistanceBetweenPoints(new Point(7, 7), new Point(x, y)) <= CursorRadius)
                    {
                        c[x, y] = 1;
                        points.Add(new PointF(x, y));
                    }
                    else
                        c[x, y] = 0;
                }
            circle = (PointF[]) points.ToArray(typeof (PointF));
        }

        private PointF[] CreateShotMarker(PointF center, int diameter) //diameter should be odd!!!
        {
            var al = new ArrayList();
            int midpoint = (diameter - 1)/2;
            var c = new int[diameter,diameter];
            for (int y = 0; y < c.GetLength(0); y++)
                for (int x = 0; x < c.GetLength(1); x++)
                {
                    if (FindDistanceBetweenPoints(new PointF(midpoint, midpoint), new PointF(x, y)) < midpoint)
                    {
                        al.Add(new PointF(center.X - midpoint + x, center.Y - midpoint + y));
                    }
                }
            return (PointF[]) al.ToArray(typeof (PointF));
        }

        private void AdjustCursorSize()
        {
            cursor_ellipse.Height = (int) (cursor_ellipse.Height/Plan_aspectMultiplier_y);
            cursor_ellipse.Width = (int) (cursor_ellipse.Width/Plan_aspectMultiplier_x);
        }

        private void ResetMask()
        {
            PathSet.mask = Matrix.Zeroes((int) DICOM_imgbox.Width, (int) DICOM_imgbox.Height);
        }

        public void ErasePixel(ref WriteableBitmap writeableBitmap, System.Windows.Controls.Image i, MouseEventArgs e)
        {
            var column = (int) e.GetPosition(i).X;
            var row = (int) e.GetPosition(i).Y;
        }

        public void ColormapTool(Bitmap b)
        {
            red_colormap = new int[256];
            blue_colormap = new int[256];
            green_colormap = new int[256];
            for (int i = 0; i < 256; i++)
            {
                red_colormap[i] = b.GetPixel(i, 0).R;
                blue_colormap[i] = b.GetPixel(i, 0).B;
                green_colormap[i] = b.GetPixel(i, 0).G;
            }
        }

        public void DisplayImageSizeInfo(int x, int y)
        {
            string s = "" + x + " x " + y;
            structure_size_lbl.Content = s;
        }

        public void SetCoverageLabelColor(double cov)
        {
            byte r = 0;
            byte g = 0;

            if (cov > 0.8)
            {
                if (cov > 0.9)
                {
                    r = 255;
                    g = (byte) Math.Round(((cov - 0.9)/0.1)*255);
                }
                else
                    r = (byte) Math.Round(((cov - 0.8)/0.1)*255);
            }
            var scb = new SolidColorBrush(Color.FromRgb((byte) (255 - r), (byte) (0 + g), 0));
            RealtimeCoverage_lbl.Foreground = scb;
            RealtimeCoverage_lbl.Content = "Max Potential Coverage: " + Math.Round(cov, 2);
        }

        public void DisplayPlan()
        {
            //UpdateSliceLabels();
            slider2.Minimum = 0;
            slider2.Maximum = PS.SlicePositions.GetLength(0) - 1;
            var rp = (RasterPath) PS.RasterPaths[GetCurrentSlice()];
            int[] lastlines = null;
            int[] nextlines = null;
            if (rp.WhichSlice > 0)
                lastlines = ((RasterPath) PS.RasterPaths[rp.WhichSlice - 1]).Lines;
            if (rp.WhichSlice < PathSet.NumSlices-1)
                nextlines = ((RasterPath) PS.RasterPaths[rp.WhichSlice + 1]).Lines;
            
            double cov = rp.info.Coverage;
            SetCoverageLabelColor(cov);
            sliceindicator_lbl.Content = "Slice: " + rp.WhichSlice;


            int zpos = PS.SlicePositions[GetCurrentSlice()];
            PointF[] points = rp.ReturnSinglePoints();
            //float[,] img = SS.fj_Tumor[zpos];
            float[,] img = rp.slice;
            float[,] p_img = rp.priorityslice;

            var doseimg = new float[img.GetLength(0),img.GetLength(1)];
            //float[] img = SS.f_structurearray[zpos]; 
            bool dp = false;
            if (plan_dpRB.IsChecked == true && rp.dosespace != null)
            {
                if (rp.dosespace != null)
                {
                    doseimg = Convert1Dto2D(rp.dosespace, img.GetLength(0), img.GetLength(1));
                    dp = true;
                }
            }

            //float max = img.Max();
            int sizex = img.GetLength(0);
            int sizey = img.GetLength(1);
            DisplayImageSizeInfo(sizex, sizey);
            Plan_aspectMultiplier_x = img.GetLength(0)/plan_imgbox.Width;
            Plan_aspectMultiplier_y = img.GetLength(1)/plan_imgbox.Height;
            //AdjustCursorSize();
            //img = Matrix.Normalize(img);
            Matrix.Normalize(ref img);
            //float sum = img.Sum();
            wb_Plan = new WriteableBitmap(sizey, sizex, 96, 96, PixelFormats.Bgr32, null);
            plan_imgbox.Source = wb_Plan;
            wb_Plan.Lock();

            unsafe
            {
                //First draw the structure slice as a background
                for (int y = 0; y < sizex; y++)
                    for (int x = 0; x < sizey; x++)
                    {
                        var pBackBuffer = (int) wb_Plan.BackBuffer;
                        pBackBuffer += y*wb_Plan.BackBufferStride;
                        pBackBuffer += x*4;
                        int value;
                        int alpha = 0;
                        int red = 0;
                        int extra = 0;

                        
                        int pval = (int) Math.Round(p_img[y, x]*255);
                        value = (int) Math.Round(img[y, x]*255);
                        int mod;
                        if (rp.mask != null)
                            mod = (int)Math.Round(rp.mask[y, x]);
                        else
                            mod = 0;


                        if (mod > 0)
                            alpha = 0;
                        else if (value > 0)
                            if (pval > 0)
                                alpha = 255;
                            else
                            {
                                alpha = 200;
                            }
                        else
                            alpha = 0;
                        int color_data = 0;
                        double multiplier;
                        if (dp)
                        {
                            var dose = (int) Math.Round(doseimg[y, x]);

                            if (dose >= 0.5)
                            {
                                multiplier = 0.1 + ((dose - 0.5)*1.5);
                                red = (int) Math.Round(multiplier*255);
                            }

                            color_data = red << 16;
                            color_data |= alpha + extra << 8;
                            color_data |= alpha + extra << 0;
                        }
                        else if (!dp)
                        {
                            //if (value > 0)
                            //    alpha = 255;
                            //else
                            //    alpha = 0;
                            color_data = 0 << 16;
                            color_data |= alpha << 8;
                            color_data |= alpha << 0;
                        }


                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }


                //Overlay the planned points
                foreach (PointF pf in points)
                {
                    PointF[] circlepoints = CreateShotMarker(pf, 3);
                    foreach (PointF p in circlepoints)
                    {
                        var y = (int) Math.Round(p.X);
                        var x = (int) Math.Round(p.Y);

                        var pBackBuffer = (int) wb_Plan.BackBuffer;
                        pBackBuffer += y*wb_Plan.BackBufferStride;
                        pBackBuffer += x*4;

                        int color_data = 255 << 16; // R
                        color_data |= 255 << 8; // G
                        color_data |= 0 << 0; // B

                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }
                }
                if (Adj_Slices_RB_btn.IsChecked == true)
                {
                    //Overlay previous lines
                    if (lastlines != null)
                    {
                        for (int i = 0; i < lastlines.GetLength(0); i++)
                        {
                            for (int j = 0; j < 40; j++)
                            {
                                var y = lastlines[i];
                                var x = j;
                                var pBackBuffer = (int) wb_Plan.BackBuffer;
                                pBackBuffer += y*wb_Plan.BackBufferStride;
                                pBackBuffer += x*4;

                                int color_data = 0 << 16;
                                color_data |= 0 << 8;
                                color_data |= 255 << 0;

                                *((int*) pBackBuffer) = color_data;

                            }
                        }
                    }
                    if (nextlines != null)
                    {
                        for (int i = 0; i < nextlines.GetLength(0); i++)
                        {
                            for (int j = 0; j < 40; j++)
                            {
                                var y = nextlines[i];
                                var x = j;
                                var pBackBuffer = (int) wb_Plan.BackBuffer;
                                pBackBuffer += y*wb_Plan.BackBufferStride;
                                pBackBuffer += x*4;

                                int color_data = 0 << 16;
                                color_data |= 128 << 8;
                                color_data |= 255 << 0;

                                *((int*) pBackBuffer) = color_data;

                            }
                        }
                    }
                }
            }
            try
            {
                wb_Plan.AddDirtyRect(new Int32Rect(0, 0, sizey, sizex));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_Plan.Unlock();
            }
        }

        private float[,] Convert1Dto2D(float[] f, int x, int y)
        {
            var output = new float[x,y];
            for (int i = 0; i < y; i++)
                for (int j = 0; j < x; j++)
                    output[j, i] = f[i*x + j];
            return output;
        }

        private void DisplayStructure(int slice)
        {
            float[,] structureimg;
            if (combined_rb.IsChecked == true)
                structureimg = SS.fj_Combined[slice];
                    //structureimg = Convert1Dto2D(SS.f_structurearray[slice], SS.fj_Tumor[0].GetLength(0), SS.fj_Tumor[0].GetLength(1));
            else if (tumor_rb.IsChecked == true)
                structureimg = SS.fj_Tumor[slice];
            else if (cs_rb.IsChecked == true)
                structureimg = SS.fj_CS[slice];
            else
                structureimg = new float[StructureSet.BIG_dim[0],StructureSet.BIG_dim[1]];
            int sizex = structureimg.GetLength(0);
            int sizey = structureimg.GetLength(1);
            StructureSet.BIG_dim = StructureSet.SS_dim;
            DisplayImageSizeInfo(structureimg.GetLength(0), structureimg.GetLength(1));

            wb_DDS = new WriteableBitmap(sizex, sizey, 96, 96, PixelFormats.Bgr32, null);
            //int max = (int)structureimg.Max();
            DDS_imgbox.Source = wb_DDS;
            wb_DDS.Lock();
            unsafe
            {
                for (int y = 0; y < sizey; y++)
                    for (int x = 0; x < sizex; x++)
                    {
                        var pBackBuffer = (int) wb_DDS.BackBuffer;
                        pBackBuffer += y*wb_DDS.BackBufferStride;
                        pBackBuffer += x*4;
                        var value = (int) structureimg[x, y];
                        int alpha;
                        int color_data = 0;

                        if (value > 0)
                        {
                            alpha = 255;
                            if (value == 1)
                            {
                                color_data = red_colormap[alpha] << 16; // R
                                color_data |= green_colormap[alpha] << 8; // G
                                color_data |= blue_colormap[alpha] << 0; // B
                            }
                            else if (value > 1)
                            {
                                color_data = red_colormap[alpha] << 16;
                                color_data |= 0 << 8;
                                color_data |= 0 << 0;
                            }
                        }
                        else
                        {
                            alpha = 0;
                            color_data = 0 << 16;
                            color_data |= 0 << 8;
                            color_data |= 0 << 0;
                        }


                        //int color_data = alpha << 16;
                        //color_data |= alpha << 8;
                        //color_data |= alpha << 0;

                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }
            }
            try
            {
                wb_DDS.AddDirtyRect(new Int32Rect(0, 0, sizex, sizey));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_DDS.Unlock();
            }
        }

        public void DisplayDICOM(float[] image)
        {
            DisplayImageSizeInfo(256, 256);
            double max = image.Max();
            wb_DICOM = new WriteableBitmap(256, 256, 96, 96, PixelFormats.Bgr32, null);
            DICOM_imgbox.Source = wb_DICOM;
            int startpt = 100;
            int endpt = 512;
            wb_DICOM.Lock();

            unsafe
            {
                for (int y = 0; y < 256; y++)
                    for (int x = 0; x < 256; x++)
                    {
                        var pBackBuffer = (int) wb_DICOM.BackBuffer;
                        pBackBuffer += y*wb_DICOM.BackBufferStride;
                        pBackBuffer += x*4;
                        double value = (int) image[x + y*256];

                        var alph = (int) Math.Round(255*(value/max));
                        if (alph < 0)
                            alph = 0;
                        if (alph >= 255)
                            alph = 255;

                        int color_data = red_colormap[alph] << 16; // R
                        color_data |= green_colormap[alph] << 8; // G
                        color_data |= blue_colormap[alph] << 0; // B
                        //int color_data = value << 16;
                        //color_data |= value << 8;
                        //color_data |= value << 0;

                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }
            }
            try
            {
                wb_DICOM.AddDirtyRect(new Int32Rect(0, 0, 256, 256));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_DICOM.Unlock();
            }
        }

        private void DisplayWindowCenteredAboutPoint(float[,] img, PointF p)
        {
            int imgsize = img.GetLength(0);
            int startx = (int) p.X - (imgsize/2);
            int starty = (int) p.Y - (imgsize/2);
            int endx = (int) p.X + (imgsize/2);
            int endy = (int) p.Y + (imgsize/2);

            if (startx < 0)
                startx = 0;
            if (starty < 0)
                starty = 0;
            if (endx >= PS.X)
                endx = PS.X - 1;
            if (endy >= PS.Y)
                endy = PS.Y - 1;

            var temp = new float[endx - startx + 1,endy - starty + 1];
            for (int j = 0; j < temp.GetLength(1); j++)
                for (int i = 0; i < temp.GetLength(0); i++)
                    temp[i, j] = img[i + startx, j + starty];

            Display2DFloat(temp);
        }

        private void DisplayOnA(float[,] f)
        {
            DisplayInfoUnderA(f.GetLength(0), f.GetLength(1));
            double default_rx_dose = 0.5;
            bool viewiso = Aiso_chkbox.IsChecked == true;


            Comparison.Awb = new WriteableBitmap(f.GetLength(1), f.GetLength(0), 96, 96, PixelFormats.Bgr32, null);
            imgA.Source = Comparison.Awb;
            Comparison.Awb.Lock();

            unsafe
            {
                for (int j = 0; j < f.GetLength(0); j++)
                    for (int i = 0; i < f.GetLength(1); i++)
                    {
                        var pBackBuffer = (int)Comparison.Awb.BackBuffer;
                        pBackBuffer += j * Comparison.Awb.BackBufferStride;
                        pBackBuffer += i * 4;
                        int value;
                        int alpha;
                        alpha = (int)Math.Round((f[j, i]) * 255);
                        int color_data = 0;
                        if (viewiso && alpha > (default_rx_dose * 255))
                        {
                            color_data = alpha << 16;
                            color_data |= 0 << 8;
                            color_data |= 0 << 0;
                        }
                        else
                        {
                            color_data = alpha << 16;
                            color_data |= alpha << 8;
                            color_data |= alpha << 0;
                        }

                        // Assign the color data to the pixel.
                        *((int*)pBackBuffer) = color_data;
                    }
            }
            try
            {
                Comparison.Awb.AddDirtyRect(new Int32Rect(0, 0, f.GetLength(1), f.GetLength(0)));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                Comparison.Awb.Unlock();
            }
        }

        private void DisplayOnB(float[,] f)
        {
            DisplayInfoUnderB(f.GetLength(0), f.GetLength(1));
            double default_rx_dose = 0.5;
            bool viewiso = Biso_chkbox.IsChecked == true;


            Comparison.Bwb = new WriteableBitmap(f.GetLength(1), f.GetLength(0), 96, 96, PixelFormats.Bgr32, null);
            imgB.Source = Comparison.Bwb;
            Comparison.Bwb.Lock();

            unsafe
            {
                for (int j = 0; j < f.GetLength(0); j++)
                    for (int i = 0; i < f.GetLength(1); i++)
                    {
                        var pBackBuffer = (int)Comparison.Bwb.BackBuffer;
                        pBackBuffer += j * Comparison.Bwb.BackBufferStride;
                        pBackBuffer += i * 4;
                        int value;
                        int alpha;
                        alpha = (int)Math.Round((f[j, i]) * 255);
                        int color_data = 0;
                        if (viewiso && alpha > (default_rx_dose * 255))
                        {
                            color_data = alpha << 16;
                            color_data |= 0 << 8;
                            color_data |= 0 << 0;
                        }
                        else
                        {
                            color_data = alpha << 16;
                            color_data |= alpha << 8;
                            color_data |= alpha << 0;
                        }

                        // Assign the color data to the pixel.
                        *((int*)pBackBuffer) = color_data;
                    }
            }
            try
            {
                Comparison.Bwb.AddDirtyRect(new Int32Rect(0, 0, f.GetLength(1), f.GetLength(0)));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                Comparison.Bwb.Unlock();
            }

        }

        private void DisplayInfoUnderA(int x, int y)
        {
            string s = "Image A: " + x + " x " + y;
        }

        private void DisplayInfoUnderB(int x, int y)
        {
            string s = "Image B: " + x + " x " + y;
        }

        private void Display2DFloat(float[,] f)
        {
            DisplayImageSizeInfo(f.GetLength(0), f.GetLength(1));

            //f = Matrix.Normalize(f);
            bool viewiso = false;
            bool viewcov = false;
            if (ViewIso_chkbox.IsChecked == true)
            {
                viewiso = true;
                viewcov = false;
            }
            if (ViewCoverage_chkbox.IsChecked == true)
            {
                viewcov = true;
                viewiso = false;
            }
            //int maxvalue = 255; //change if increase range
            double default_rx_dose = 0.5;

            //int iso = (int)Math.Round(maxvalue*default_rx_dose);
            wb_DS = new WriteableBitmap(f.GetLength(1), f.GetLength(0), 96, 96, PixelFormats.Bgr32, null);
            DS_imgbox.Source = wb_DS;
            wb_DS.Lock();

            unsafe
            {
                for (int j = 0; j < f.GetLength(0); j++)
                    for (int i = 0; i < f.GetLength(1); i++)
                    {
                        var pBackBuffer = (int) wb_DS.BackBuffer;
                        pBackBuffer += j*wb_DS.BackBufferStride;
                        pBackBuffer += i*4;
                        int value;
                        int alpha;
                        alpha = (int) Math.Round((f[j, i])*255);
                        int color_data = 0;
                        if (viewiso && alpha > (default_rx_dose*255))
                        {
                            color_data = alpha << 16;
                            color_data |= 0 << 8;
                            color_data |= 0 << 0;
                        }
                        else
                        {
                            color_data = alpha << 16;
                            color_data |= alpha << 8;
                            color_data |= alpha << 0;
                        }

                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }
            }
            try
            {
                wb_DS.AddDirtyRect(new Int32Rect(0, 0, f.GetLength(1), f.GetLength(0)));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_DS.Unlock();
            }
        }

        private void DisplayDoseWithStruct(float[,] f)
        {
            int slice = (int) slider2.Value;
            float[,] structureimg = SS.fj_Combined[slice];

            DisplayImageSizeInfo(f.GetLength(0), f.GetLength(1));

            //f = Matrix.Normalize(f);
            bool viewiso = false;
            bool viewcov = false;
            if (ViewIso_chkbox.IsChecked == true)
            {
                viewiso = true;
                viewcov = false;
            }
            if (ViewCoverage_chkbox.IsChecked == true)
            {
                viewcov = true;
                viewiso = false;
            }
            //int maxvalue = 255; //change if increase range
            double default_rx_dose = PathSet.RxDose;

            //int iso = (int)Math.Round(maxvalue*default_rx_dose);
            wb_DS = new WriteableBitmap(f.GetLength(1), f.GetLength(0), 96, 96, PixelFormats.Bgr32, null);
            DS_imgbox.Source = wb_DS;
            wb_DS.Lock();

            unsafe
            {
                for (int j = 0; j < f.GetLength(0); j++)
                    for (int i = 0; i < f.GetLength(1); i++)
                    {
                        var pBackBuffer = (int)wb_DS.BackBuffer;
                        pBackBuffer += j * wb_DS.BackBufferStride;
                        pBackBuffer += i * 4;
                        int value;
                        int alpha;
                        int beta = 0;
                        alpha = (int)Math.Round((f[j, i]) * 255);
                        if (structureimg[j, i] > 0 && alpha >= (default_rx_dose * 255))
                            beta = 255;
                        else if (structureimg[j, i] > 0 && alpha < (default_rx_dose * 255))
                            beta = 128;
                        int color_data = 0;
                        if (viewiso && alpha > (default_rx_dose * 255))
                        {
                            color_data = alpha << 16;
                            color_data |= 0 << 8;
                            color_data |= beta << 0;
                        }
                        else
                        {
                            color_data = alpha << 16;
                            color_data |= alpha << 8;
                            color_data |= beta << 0;
                        }

                        // Assign the color data to the pixel.
                        *((int*)pBackBuffer) = color_data;
                    }
            }
            try
            {
                wb_DS.AddDirtyRect(new Int32Rect(0, 0, f.GetLength(1), f.GetLength(0)));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_DS.Unlock();
            }

        }

        private void DisplaySingleShot(PointF p, double weight)
        {
            UpdateSliceLabels();
            var rp = (RasterPath) PS.RasterPaths[GetCurrentSlice()];
            int zpos = PS.SlicePositions[GetCurrentSlice()];
            PointF[] points = rp.ReturnSinglePoints();
            float[] img = SS.f_structurearray[zpos];
            bool dp = false;
            if (plan_dpRB.IsChecked == true && rp.dosespace != null)
            {
                img = rp.dosespace;
                dp = true;
            }
            float max = img.Max();
            int sizex = SS.fj_Tumor[0].GetLength(0);
            int sizey = SS.fj_Tumor[0].GetLength(1);

            float sum = img.Sum();
            wb_Plan = new WriteableBitmap(sizey, sizex, 96, 96, PixelFormats.Bgr32, null);
            plan_imgbox.Source = wb_Plan;
            wb_Plan.Lock();

            unsafe
            {
                //First draw the structure slice as a background
                for (int y = 0; y < sizex; y++)
                    for (int x = 0; x < sizey; x++)
                    {
                        var pBackBuffer = (int) wb_Plan.BackBuffer;
                        pBackBuffer += y*wb_Plan.BackBufferStride;
                        pBackBuffer += x*4;
                        int value;
                        int alpha;
                        alpha = (int) Math.Round((img[y + x*sizex])*255);
                        if (!dp)
                        {
                            if (alpha > 0)
                                alpha = 255;
                            else
                                alpha = 0;
                        }

                        //int color_data = red_colormap[alph] << 16; // R
                        //color_data |= green_colormap[alph] << 8;   // G
                        //color_data |= blue_colormap[alph] << 0;   // B
                        int color_data = alpha << 16;
                        color_data |= alpha << 8;
                        color_data |= alpha << 0;

                        // Assign the color data to the pixel.
                        *((int*) pBackBuffer) = color_data;
                    }

                //Overlay the planned points
                foreach (PointF pf in points)
                {
                    var y = (int) Math.Round(pf.X);
                    var x = (int) Math.Round(pf.Y);

                    var pBackBuffer = (int) wb_Plan.BackBuffer;
                    pBackBuffer += y*wb_Plan.BackBufferStride;
                    pBackBuffer += x*4;

                    int color_data = 255 << 16; // R
                    color_data |= 0 << 8; // G
                    color_data |= 0 << 0; // B

                    // Assign the color data to the pixel.
                    *((int*) pBackBuffer) = color_data;
                }
            }
            try
            {
                wb_Plan.AddDirtyRect(new Int32Rect(0, 0, sizey, sizex));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb_Plan.Unlock();
            }
        }

        private void DisplaySingleShotDS(PointF p)
        {
            var rp = (RasterPath) PS.RasterPaths[dataGrid1.SelectedIndex];
            //float[,] f = RasterPath.GetMultiplied_DS_Subset(rp.dosespace, p.X, p.Y, RasterPath.dosemidplane);
            float[,] f = Matrix.Subset(rp.dosespace, RasterPath.Y, RasterPath.X, (int) p.Y, (int) p.X, N);
            //f = Matrix.Normalize(f);
            //DisplayWindowCenteredAboutPoint(f, p);
            Display2DFloat(f);
        }

        #endregion

        private void loadA_btn_Click(object sender, RoutedEventArgs e)
        {
            Comparison.LoadA();
            DisplayOnA(Comparison.A.OriginalDose[0]);

        }

        private void loadB_btn_Click(object sender, RoutedEventArgs e)
        {
            Comparison.LoadB();
            DisplayOnB(Comparison.B.OriginalDose[0]);
        }

        private void SaveFinalReport_btn_Click(object sender, RoutedEventArgs e)
        {
            string savepath = Analysis.CreateFinalReportFile();
            ExportShotList();
            Savepath_loc_txt.Content = "Report saved to: " + savepath;

        }

        

        
        
    }

    
}