using System;
using System.IO;
using System.Resources;
using System.Reflection;
using System.Windows.Controls;
using System.Windows;
using System.Drawing;
using System.Windows.Data;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.Windows.Input;
using System.ComponentModel;
using System.Threading;
using System.Threading.Tasks;

namespace TomosurgeryAlpha
{
    public abstract class Image
    {
        public virtual void DrawPixel(ref WriteableBitmap writeableBitmap, int x, int y)
        {
            int column = x;
            int row = y;

            // Reserve the back buffer for updates.
            writeableBitmap.Lock();
            
            unsafe
            {
                // Get a pointer to the back buffer.
                int pBackBuffer = (int)writeableBitmap.BackBuffer;

                // Find the address of the pixel to draw.
                pBackBuffer += row * writeableBitmap.BackBufferStride;
                pBackBuffer += column * 4;

                // Compute the pixel's color.
                int color_data = 255 << 16; // R
                color_data |= 255 << 8;   // G
                color_data |= 255 << 0;   // B

                // Assign the color data to the pixel.
                *((int*)pBackBuffer) = color_data;
            }

            // Specify the area of the bitmap that changed.
            writeableBitmap.AddDirtyRect(new Int32Rect(column, row, 1, 1));

            // Release the back buffer and make it available for display.
            writeableBitmap.Unlock();
        }

        public virtual WriteableBitmap DisplayFloatArray(ref WriteableBitmap wb, float[] d)
        {
            //NEEDS TO BE NORMALIZED BEFOREHAND
            wb = new WriteableBitmap((int)wb.Width, (int)wb.Height, 96, 96, PixelFormats.Bgr32, null);
            wb.Lock();

            unsafe
            {
                for (int y = 0; y < wb.Height; y++)
                    for (int x = 0; x < wb.Width; x++)
                    {
                        int pBackBuffer = (int)wb.BackBuffer;
                        pBackBuffer += y * wb.BackBufferStride;
                        pBackBuffer += x * 4;
                        int value = (int)(d[(y * x) + x]);
                        int color_data = value << 16;
                        color_data |= value << 8;
                        color_data |= value << 0;
                        // Assign the color data to the pixel.
                        *((int*)pBackBuffer) = color_data;
                    }
            }
            try
            {
                wb.AddDirtyRect(new Int32Rect(0, 0, d.GetLength(0), d.GetLength(1)));
            }
            catch (Exception ex)
            {
                string s = ex.ToString();
            }
            finally
            {
                wb.Unlock();
            }
            return wb;
        }
    }
}
