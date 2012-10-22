using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace TomosurgeryAlpha
{
    /// <summary>
    /// Interaction logic for RestrictSliceWeight.xaml
    /// </summary>
    public partial class RestrictSliceWeight : Window
    {
        public static double[] SliceWeight; //0 is number, 1 is weight.

        public RestrictSliceWeight()
        {
            InitializeComponent();
        }

        private void Reweightslices_btn_Click(object sender, RoutedEventArgs e)
        {
            SliceWeight = new double[2]{Convert.ToDouble(Slicenum_box.Text), Convert.ToDouble(Sliceweight_box.Text)};
            return;

        }
    }
}
