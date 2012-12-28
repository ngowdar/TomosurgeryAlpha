using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

namespace TomosurgeryAlpha
{
    public struct Coordinates
    {
        public const int XPostWidth = 190; //This is what I found the distance between the L and R fiducials to be.
        public const int YPostWidth = 120; //This is what I found the distance between the ant and post fiducials to be.
        public PointF[] align;
        public bool finished;
        public bool leftset;
        public bool rightset;
        public double X { get; set; } //This is the X-offset value (or, location of the LGK origin in the image space)
        public double Y { get; set; }
        public double Z { get; set; }
        public double left { get; set; } //The DICOM coordinate when the left z-fiducial is aligned.
        public double right { get; set; } //The DICOM coordinate when the right z-fiducial is aligned.
        public double left_slider { get; set; }
        public double right_slider { get; set; }

        //The DICOM coordinates of the LGK (100,100,100) point.
        public double x100 { get; set; }
        public double y100 { get; set; }
        public double z100 { get; set; }

        public double slider_100 { get; set; }

        //The center point should be equal to (100,100,100) in LGK coordinates. So, this method finds the offset from this point.
        public void SetLGK_XYoffset(double x, double y)
        {
            /* since the midpoint in DICOM space = (100,100,100) in LGK space, then (100-midpoint) should be coordinate for LGK origin in DICOM space.
             * Cursor location minus the origin location (in image coordinates) should give LGK coordinates.
             * However, must remember that y-axis is flipped (image y-positive is down, LGK y-positive is up).
             * */

            X = ((x + x + XPostWidth)/2) - 100;
            Y = ((y + y + YPostWidth)/2) + 100;
            finished = false;
        }

        public void SetLGK_Zoffset()
        {
            z100 = ((left + right)/2);
            Z = z100 - 100;
            slider_100 = (left_slider + right_slider)/2;
            finished = true;
        }

        public double getLGK_Z_forDICOM(double slidervalue)
        {
            //return 100 + (slidervalue - slider_100) * 2.0;
            return slidervalue*2 + Z + 100;
        }

/*
        public double getLGK_Z_forDose(double slidervalue)
        {
            //return 100 + (slidervalue - slider_100) * 2.0;
            double d = getLGK_Z_forDICOM(20) + slidervalue;
            //return Convert.ToDouble(DICOMdose.doseoffset[2]) - (slidervalue + Z + 100);
            return d;
        }
*/

        public double getLGK_X(double x)
        {
            return x - X;
        }

        public double getLGK_Y(double y)
        {
            return Y - y;
        }

        public decimal[] Image2LGKCoordinates(decimal[] dd)
        {
            var d = new decimal[3];
            d[0] = Math.Round((decimal) getLGK_X(Convert.ToDouble(dd[0])), 2);
            d[1] = Math.Round((decimal) getLGK_Y(Convert.ToDouble(dd[1])), 2);
            d[2] = Math.Round((decimal) getLGK_Z_forDICOM(Convert.ToDouble(dd[2])), 2);
            return d;
        }

/*
        public decimal[] LGK2ImageCoordinates(decimal[] dd)
        {
            var d = new decimal[3];
            d[0] = Math.Round((dd[0] + (decimal) X), 2);
            d[1] = Math.Round((decimal) Y - dd[1], 2);
            d[2] = Math.Round(((dd[2] - 100)) + (decimal) slider_100, 2);
            return d;
        }
*/
    }
}