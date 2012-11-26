using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.ComponentModel;
using System.Drawing;
using System.Diagnostics;

namespace TomosurgeryAlpha
{
    public static class Matrix
    {
        public static float[] ScalarMultiply(float[] f, float m)
        {
            float[] F = new float[f.GetLength(0)];
            Parallel.For(0, f.GetLength(0), (i) =>
                {
                    F[i] = f[i] * m;
                });
            return F;
        }

        public static float[][] MakeJaggedFloat(float[, ,] f)
        {
            int[] size = GetSize(f);
            float[][] ff = new float[size[2]][];
            for (int k = 0; k < size[2]; k++)
                for (int i = 0; i < size[0]; i++)
                    for (int j = 0; j < size[1]; j++)
                        ff[k][j*size[0]+i] = f[i,j,k];
            return ff;
        }

        public static float[][,] MakeJaggedFloat(float[] f, int x, int y, int z)
        {
            int[] size = new int[3]{x, y, z};
            float[][,] ff = new float[size[2]][,];
            for (int k = 0; k < z; k++)
            {
                float[,] temp = new float[x, y];
                for (int j = 0; j < y; j++)
                    for (int i = 0; i < x; i++)
                        temp[i, j] = f[k * x * y + j * x + i];
                ff[k] = temp;
            }
            return ff;
        }

        public static float[] Zero1DFloat(int xyz)
        {
            return (float[])Array.CreateInstance(typeof(float), xyz);
        }

        public static float[][] ZeroJaggedFloat(int z, int xy)
        {
            float[][] f = new float[z][];
            float[] temp = Zero1DFloat(xy);
            for (int i = 0; i < z; i++)
                f[i] = (float[])temp.Clone();
            return f;
        }

        #region GetSize (overloaded)

        //3-dimensional arrays
        public static int[] GetSize(float[, ,] f)
        {
            return new int[3]{f.GetLength(0),f.GetLength(1),f.GetLength(2)};
        }
        public static int[] GetSize(int[, ,] f)
        {
            return new int[3] { f.GetLength(0), f.GetLength(1), f.GetLength(2) };
        }
        public static int[] GetSize(double[, ,] f)
        {
            return new int[3] { f.GetLength(0), f.GetLength(1), f.GetLength(2) };
        }

        //2-dimensional arrays
        public static int[] GetSize(float[,] f)
        {
            return new int[2] { f.GetLength(0), f.GetLength(1)};
        }
        public static int[] GetSize(int[,] f)
        {
            return new int[2] { f.GetLength(0), f.GetLength(1) };
        }
        public static int[] GetSize(double[,] f)
        {
            return new int[2] { f.GetLength(0), f.GetLength(1) };
        }

        //3-D jagged arrays (1D array of 2D arrays)
        public static int[] GetSize(float[][,] f)
        {
            return new int[3] { f[0].GetLength(0), f[0].GetLength(1), f.GetLength(0) };
        }
        public static int[] GetSize(int[][,] f)
        {
            return new int[3] { f[0].GetLength(0), f[0].GetLength(1), f.GetLength(0) };
        }
        public static int[] GetSize(double[][,] f)
        {
            return new int[3] { f[0].GetLength(0), f[0].GetLength(1), f.GetLength(0) };
        }
        #endregion
        #region SumAll
        public static int SumAll(float[, ,] f)
        {
            int[] size = GetSize(f);
            float sum = 0;
            for (int k = 0; k < size[2]; k++)
                for (int i = 0; i < size[0]; i++)
                    for (int j = 0; j < size[1]; j++)
                    {                        
                        sum += f[i, j, k];
                    }
            return (int)sum;
        }
        public static int SumAll(double[, ,] f)
        {
            int[] size = GetSize(f);            
            double sum = 0;
            for (int k = 0; k < size[2]; k++)
                for (int i = 0; i < size[0]; i++)
                    for (int j = 0; j < size[1]; j++)
                    {
                        sum += f[i, j, k];
                    }
            return (int)sum;
        }

        public static int SumAll(int[, ,] f)
        {
            int[] size = GetSize(f);
            float sum = 0;
            for (int k = 0; k < size[2]; k++)
                for (int i = 0; i < size[0]; i++)
                    for (int j = 0; j < size[1]; j++)
                    {
                        sum += f[i, j, k];
                    }
            return (int)sum;
        }
        public static int SumAll(float[,] f)
        {
            int[] size = GetSize(f);
            float sum = 0;            
                for (int i = 0; i < size[0]; i++)
                    for (int j = 0; j < size[1]; j++)
                    {
                        sum += f[i, j];
                    }
            return (int)sum;
        }
        public static int SumAll(double[,] f)
        {
            int[] size = GetSize(f);
            double sum = 0;
            for (int i = 0; i < size[0]; i++)
                for (int j = 0; j < size[1]; j++)
                {
                    sum += f[i, j];
                }
            return (int)sum;
        }
        public static int SumAll(int[,] f)
        {
            int[] size = GetSize(f);
            float sum = 0;
            for (int i = 0; i < size[0]; i++)
                for (int j = 0; j < size[1]; j++)
                {
                    sum += f[i, j];
                }
            return (int)sum;
        }
        #endregion


        public static float[] makeFloatArray(byte[][] data)
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

        public static byte[] makeByteArray(byte[][] data)
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

        
        
        public static float[] Byte2Float(byte[] data)
        {
            float[] output = new float[data.GetLength(0)/2];            
            for (int i = 0; i < data.GetLength(0); i += 2)
            {  
                
                output[i / 2] = BitConverter.ToUInt16(data, i);                
            }
            return output;            
        }

        internal static void AddSubsetJagged(float[][,] dd, int startx, int starty, int startz, int addxsize, int addysize)
        {
            
            
        }

        //public static float[][,] EnlargeAndCenter(float[][,] dd, int enlargement, int startx, int starty, int startz)
        //{
        //    LinMatrix m = new LinMatrix(dd);
        //    m.SetSizeofM(m.X + enlargement, m.Y + enlargement, m.Z + enlargement);
        //    //m.Add(m.Convertto1D(dd), enlargement / 2, enlargement / 2, enlargement / 2, dd[0].GetLength(0), dd[1].GetLength(1));
        //    m.Add_NoGPU(m.Convertto1D(dd), enlargement / 2, enlargement / 2, enlargement / 2, dd[0].GetLength(0), dd[1].GetLength(1));
            
        //    return m.ConvertToJagged();
        //}

        public static float[] Grab1DSlice(float[] d, int x, int y, int z)
        {
            float[] slice = new float[x*y];
            for (int j = 0; j < y; j++)
                for (int i = 0; i < x; i++)
                    slice[z*(x*y)+(j*x) + i] = d[z*(x*y) + (j*x) + i];
            return slice;
        }

        public static float[] Convertto1D(float[][,] d)
        {
            float[] result;
            int z = d.GetLength(0); int x = d[0].GetLength(0); int y = d[0].GetLength(1);
            result = new float[z * x * y];
            for (int i = 0; i < z; i++)
                for (int j = 0; j < y; j++)
                    for (int k = 0; k < x; k++)
                        result[(i * x * y) + (j * x) + k] = (float)d[i][k, j];
            return result;
        }

        //Added 10/17/2012
        public static float[][,] EnlargeAndCenter(float[] dd, int padsize, int sizex, int sizey, int sizez)
        {
            float[][,] output = new float[sizez+(2*padsize)][,];
            //Set everything to zeroes.
            for (int i = 0; i < output.GetLength(0); i++)
                output[i] = Matrix.Zeroes(sizex + (2 * padsize), sizey + (2 * padsize));

            //Starting at z=padsize, start filling shit in.
            for (int k = 0; k < sizez; k++)
            {
                float[,] s = new float[sizex+padsize+padsize,sizey+padsize+padsize];
                for (int j = 0; j < sizey; j++)
                    for (int i = 0; i < sizex; i++)
                    {
                        s[padsize + i, padsize + j] = dd[(k * sizex * sizey) + (j * sizex) + i];
                    }
                output[k + padsize] = s;
            }
            return output;
        }

        public static void WriteFloatArray2BMP(float[,] temp, string p)
        {
            string path = System.IO.Path.Combine(PathSet.ActiveDirectory, p);
            float[,] temp2 = (float[,])Matrix.Normalize(temp).Clone();
            int color = 0;
            Bitmap b = new Bitmap(temp.GetLength(0), temp.GetLength(1));
            for (int j = 0; j < temp.GetLength(1); j++)
                for (int i = 0; i < temp.GetLength(0); i++)
                {
                    color = (int)(temp2[i, j] * 255);
                    b.SetPixel(i, j, Color.FromArgb(color, color, color));
                }
            b.Save(path);
        }

        public static float[,] TransposeMatrix(float[,] d)
        {
            float[,] transpose = new float[d.GetLength(1), d.GetLength(0)];
            Parallel.For(0, d.GetLength(0), (i) =>
                {
                    for (int j = 0; j < d.GetLength(1); j++)
                    {
                        transpose[j, i] = d[i, j];
                    }
                });
            return transpose;
        }

        public static float[][,] TransposeMatrix(float[][,] d)
        {
            float[][,] output = new float[d.GetLength(0)][,];
            for (int i = 0; i < d.GetLength(0); i++)
                output[i] = TransposeMatrix(d[i]);
            return output;
        }


        public static float[,] Add(float[,] A, float[,] B)
        {
            float[,] sum = new float[A.GetLength(0), A.GetLength(1)];
            Parallel.For(0, A.GetLength(0), (x) =>
            {
                for (int y = 0; y < A.GetLength(1); y++)
                    sum[x, y] = A[x, y] + B[x, y];
            });
            return sum;
        }

        public static float[,] ThresholdEq(float[,] d, float th)
        {
            float[,] t = new float[d.GetLength(0), d.GetLength(1)];
            for (int i = 0; i < d.GetLength(0); i++)
                for (int j = 0; j < d.GetLength(1); j++)
                {
                    if (d[i, j] >= th)
                        t[i, j] = 1;
                    else
                        t[i, j] = 0;
                }
            return t;
        }

        public static float[][,] ThresholdEq(float[][,] d, float th)
        {
            float[][,] t = new float[d.GetLength(0)][,];
            for (int k = 0; k < d.GetLength(0); k++)
                t[k] = ThresholdEq(d[k], 0.4f);
            return t;
        }




        internal static float[,] ThresholdEq(float[,] d, int th)
        {            
                float[,] t = new float[d.GetLength(0), d.GetLength(1)];
                Parallel.For(0, d.GetLength(0), (i) =>
                    {
                        for (int j = 0; j < d.GetLength(1); j++)
                        {
                            if (d[i, j] >= th)
                                t[i, j] = 1;
                            else t[i, j] = 0;
                        }
                    });

                return t;            
        }

        public static float FindMax(float[,] d)
        {
            float max = 0;
            for (int i = 0; i < d.GetLength(0); i++)
                for (int j = 0; j < d.GetLength(1); j++)
                    if (d[i, j] > max)
                        max = d[i, j];
            return max;
        }

        public static float[,] Normalize(float[,] d)
        {
            float max = FindMax(d);
            float divisor;
            if (max <= 0)
                divisor = 1.0f;
            else
                divisor = 1 / FindMax(d);            
            
            return ScalarMultiply(d, divisor);
        }

        public static float[][,] Normalize(float[][,] d)
        {
            float max = 0;
            for (int k = 0; k < d.GetLength(0); k++)
            {                
                float tempmax = FindMax(d[k]);
                if (tempmax > max)
                    max = tempmax;
            }

            float divisor;
            if (max <= 0)
                divisor = 1.0f;
            else
                divisor = 1 / max;     

            for (int k = 0; k < d.GetLength(0); k++)
            {               
                d[k] = ScalarMultiply(d[k], divisor);
            }
            return d;
        }

        public static float[][,] NormalizeBy(float[][,] d, float max)
        {
            float divisor;
            if (max <= 0)
                divisor = 1.0f;
            else
                divisor = 1 / max;

            for (int k = 0; k < d.GetLength(0); k++)
            {
                d[k] = ScalarMultiply(d[k], divisor);
            }
            return d;
        }

        public static float[,] ScalarMultiply(float[,] a, float scalar)
        {
            float[,] product = new float[a.GetLength(0), a.GetLength(1)];
            Parallel.For(0, a.GetLength(0), (x) =>
            {
                for (int y = 0; y < a.GetLength(1); y++)
                    product[x, y] = a[x, y] * scalar;
            });
            return product;
        }

        internal static float[][,] ScalarMultiply(float[][,] a, float scalar)
        {
            float[][,] product = new float[a.GetLength(0)][,];
            Parallel.For(0, a.GetLength(0), (k) =>
            {
                product[k] = new float[a[0].GetLength(0), a[0].GetLength(1)];
                for (int y = 0; y < a[0].GetLength(1); y++)
                    for (int x = 0; x < a[0].GetLength(0); x++)
                        product[k][x, y] = a[k][x, y] * scalar;
            });
            return product;
        }

        //public static float[,] Add(float[,] A, float[,] B)
        //{
        //    float[,] sum = new float[A.GetLength(0), A.GetLength(1)];
        //    Parallel.For(0, A.GetLength(0), (x) =>
        //    {
        //        for (int y = 0; y < A.GetLength(1); y++)
        //            sum[x, y] = A[x, y] + B[x, y];
        //    });
        //    return sum;
        //}

        

        internal static float[,] MultiplyElements(float[,] A, float[,] B)
        {
            float[,] product = new float[A.GetLength(0), A.GetLength(1)];
            if (A.GetLength(0) == B.GetLength(0) && A.GetLength(1) == B.GetLength(1))
            {
                Parallel.For(0, A.GetLength(0), (x) =>
                {
                    for (int y = 0; y < A.GetLength(1); y++)
                        product[x, y] = A[x, y] * B[x, y];
                });
                return product;
            }
            else
            {
                System.Windows.MessageBox.Show("MultiplyElements: Matrices aren't same length!!");
                return product;
            }
        }

        internal static float[][,] MultiplyElements(float[][,] A, float[][,] B)
        {
            float[][,] product = new float[A.GetLength(0)][,];
            if (A.GetLength(0) != B.GetLength(0))
            {
                System.Windows.MessageBox.Show("MultiplyElements: Matrices aren't same length!!");
            }
            else
            {
                //Loop through each z-slice and call the 2D version of the function.
                Parallel.For(0, A.GetLength(0), (k) =>
                    {
                        product[k] = MultiplyElements(A[k], B[k]);
                    });
            }
            return product;
        }

        internal static float[,] Subset(float[,] A, int centerx, int centery, int subsetsize)
        {
            float[,] Window;
            int halfwindow = (subsetsize / 2);
            int startx = centerx - halfwindow; int starty = centery - halfwindow;
            int endx = centerx + halfwindow; int endy = centery + halfwindow;
            bool xfits = false; bool yfits = false;

            if (endx >= A.GetLength(0))
            {
                endx = A.GetLength(0)-1;
                xfits = false;
            }
            else if (endx < A.GetLength(0))
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= A.GetLength(1))
            {
                endy = (A.GetLength(1) - 1);
                yfits = false;
            }
            else if (endy < A.GetLength(1))
            { yfits = true; }

            if (starty < 0)
            {
                starty = 0;
                yfits = false;
            }
            else if (starty >= 0)
                yfits = true;

            Window = new float[endx - startx, endy - starty];
            //Parallel.For(0, Window.GetLength(1), (j) =>
            //{
            //    for (int i = 0; i < subsetsize; i++)
            //        Window[i, j] = A[i+startx,j+starty];
            //}); 

            Parallel.For(0, endy - starty, (j) =>
                {
                    for (int i = 0; i < endx - startx; i++)
                        Window[i, j] = A[i + startx, j + starty];
                }); 
            return Window;
        }

        internal static float[,] Subset(float[] A, int Ax, int Ay, int centerx, int centery, int subsetsize)
        {
            float[,] Window;
            //int Asize = (int)Math.Sqrt(A.GetLength(0));
            int halfwindow = (subsetsize / 2);
            int startx = centerx - halfwindow; int starty = centery - halfwindow;
            int endx = centerx + halfwindow; int endy = centery + halfwindow;
            bool xfits = false; bool yfits = false;

            if (endx >= Ax)
            {
                endx = (Ax - 1);
                xfits = false;
            }
            else if (endx < Ax)
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= Ay)
            {
                endy = (Ay - 1);
                yfits = false;
            }
            else if (endy < Ay)
            { yfits = true; }

            if (starty < 0)
            {
                starty = 0;
                yfits = false;
            }
            else if (starty >= 0)
                yfits = true;
            //Parallel.For(0, Window.GetLength(1), (j) =>
            //{
            //    for (int i = 0; i < subsetsize; i++)
            //        Window[i, j] = A[i + startx + (j + starty)*Asize];
            //});
            Window = new float[endx - startx, endy - starty];
            Parallel.For(0, endy - starty, (j) =>
            {
                for (int i = 0; i < endx - startx; i++)
                    Window[i, j] = A[(i + startx) + (j + starty)*Ax];
            });
            return Window;
        }

        internal static float[,] MultiplySubset(float[,] A, float[,] b, int centerx, int centery)
        {
            float[,] AA = new float[A.GetLength(0), A.GetLength(1)];
            int halfB = (b.GetLength(0) - 1) / 2;
            int startx = centerx - halfB; int starty = centery - halfB;
            int endx = centerx + halfB; int endy = centery + halfB;
            bool xfits = false; bool yfits = false;

            if (endx >= A.GetLength(0))
            {
                endx = (AA.GetLength(0) - 1);
                xfits = false;
            }
            else if (endx < A.GetLength(0))
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= A.GetLength(1))
            {
                endy = (AA.GetLength(1) - 1);
                yfits = false;
            }
            else if (endy < AA.GetLength(1))
            { yfits = true; }

            if (starty < 0)
            {
                starty = 0;
                yfits = false;
            }
            else if (starty >= 0)
                yfits = true;
            
            //First, fill in all the values in AA, equal to A.
            Parallel.For(0, AA.GetLength(1), (j) =>
            {
                for (int i = 0; i < AA.GetLength(0); i++)
                    AA[i, j] = A[i, j];
            });
            //Then, go back through and replace the subset values with the newly multiplied values.
            Parallel.For(0, endy - starty, (j) =>
                {
                    for (int i = 0; i < endx - startx; i++)
                        AA[i+startx, j+starty] = A[i + startx, j + starty] * b[i, j];
                }); 
            return AA;
        }

        internal static float[,] MultiplySubset(float[] A, float[,] b, int centerx, int centery, int sizex, int sizey)
        {
            //Check if b will fit inside A
            int Asize = (int)Math.Sqrt(A.GetLength(0));
            float[,] window = new float[sizex, sizey];
            int halfB = (b.GetLength(0) - 1) / 2;
            int halfsize = sizex / 2;
            int startx = centerx - halfsize; int starty = centery - halfsize;
            int endx = centerx + halfsize; int endy = centery + halfsize;
            int startdose = halfB - halfsize;
            bool xfits = false; bool yfits = false;

            if (endx >= Asize)
            {
                endx = (window.GetLength(0) - 1);
                xfits = false;
            }
            else if (endx < Asize)
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= Asize)
            {
                endy = (window.GetLength(1) - 1);
                yfits = false;
            }
            else if (endy < window.GetLength(1))
            { yfits = true; }

            if (starty < 0)
            {
                starty = 0;
                yfits = false;
            }
            else if (starty >= 0)
                yfits = true;
            //for (int i = 0; i < endx - startx; i++)
            //    for (int j = 0; j < endy - starty; j++)
            //        AA[i + startx, j + starty] = A[(i + startx)+(j + starty)*sizex] * b[i, j];
            Parallel.For(0, sizey, (j) =>
            {
                for (int i = 0; i < sizex; i++)
                    window[i, j] = A[(i+startx) + ((j+starty)* Asize)];
            });

            Parallel.For(0, endy - starty, (j) =>
                {
                    for (int i = 0; i < endx - startx; i++)
                        window[i, j] = A[(i + startx) + (j + starty) * Asize] * b[startdose+i, startdose+j];
                });

            return window;
        }

        internal static float[] Normalize(float[] img)
        {            
            float[] n = new float[img.GetLength(0)];
            float max = img.Max();
            Parallel.For(0, img.GetLength(0), (i) =>
                {
                    n[i] = img[i] / max;
                });
            return n;
        }

        internal static float[,] Zeroes(int p1, int p2)
        {
            float[,] output = new float[p1, p2];
            for (int i = 0; i < p1; i++)
                for (int j = 0; j < p2; j++)
                    output[i, j] = 0;
            return output;
        }

        internal static double SumAll(float[][,] p)
        {
            double sum = 0;
            for (int k = 0; k < p.GetLength(0); k++)
                for (int j = 0; j < p[0].GetLength(1); j++)
                    for (int i = 0; i < p[0].GetLength(0); i++)
                    {
                        sum += p[k][i, j];
                    }
            return sum;
        }

        public static double[, ,] LinearInterpSlicesX(double[, ,] d)
        {
            int xx = 2 * d.GetLength(0);
            double[, ,] temp = new double[xx - 1, d.GetLength(1), d.GetLength(2)];
            //Interpolating along x            
            for (int x = 0; x < d.GetLength(0); x++)
            {
                for (int i = 0; i < d.GetLength(1); i++)
                    for (int j = 0; j < d.GetLength(2); j++)
                    {
                        if (x == d.GetLength(0) - 1) //i.e. if it is the last slice
                        {
                            temp[2 * x, i, j] = d[x, i, j];
                        }
                        else
                        {
                            temp[2 * x, i, j] = d[x, i, j];
                            temp[(2 * x + 1), i, j] = (d[x, i, j] + d[(x + 1), i, j]) / 2;
                        }
                    }
                //int p = 100 * (x / d.GetLength(0));
                //MatrixInterp.ReportProgress(p);
            }
            return temp;
        }

        public static double[, ,] LinearInterpSlicesY(double[, ,] d)
        {
            int yy = 2 * d.GetLength(1);
            double[, ,] temp = new double[d.GetLength(0), yy - 1, d.GetLength(2)];
            //Interpolating along z            
            for (int y = 0; y < d.GetLength(1); y++)
            {
                for (int i = 0; i < d.GetLength(0); i++)
                    for (int j = 0; j < d.GetLength(2); j++)
                    {
                        //if (y == 0)
                        //    temp[i, 0, j] = d[i, 0, j];

                        if (y == d.GetLength(1) - 1) //i.e. if it is the last slice
                        {
                            temp[i, 2 * y, j] = d[i, y, j];
                        }
                        else
                        {
                            temp[i, (2 * y), j] = d[i, y, j];
                            temp[i, (2 * y + 1), j] = (d[i, y, j] + d[i, (y + 1), j]) / 2;
                        }
                    }
                //int p = 100 * (y / d.GetLength(1));
                //MatrixInterp.ReportProgress(p);
            }
            return temp;
        }

        public static float[,] LinearlyInterpolateSlices(float[,] d)
        {
            //Covnert d matrix into a float[][] array;            

            float[] row;
            float[][] rows = new float[d.GetLength(0)][];
            for (int j = 0; j < d.GetLength(0); j++)
            {
                row = new float[d.GetLength(1)];
                for (int i = 0; i < d.GetLength(1); i++)
                    row[i] = d[j, i];
                rows[j] = InterpolateSingleVector(row);
            }
            return InterpolateVectors(rows);
        }

        public static float[] InterpolateSingleVector(float[] f)
        {
            float[] output = new float[(2*f.GetLength(0)) - 1];
            //Fill in original elements.
            for (int i = 0; i < f.GetLength(0); i++)
                output[2*i] = f[i];
            //Average to find the others.
            for (int i = 1; i < f.GetLength(0); i++)
            {
                output[(2 * i) - 1] = (output[(2 * i) - 2] + output[2 * i]) / 2;
            }            
            return output;

        }

        public static float[,] InterpolateVectors(float[][] f)
        {
            float[,] A = new float[(2 * f.GetLength(0)) - 1, f[0].GetLength(0)];
            //Fill in values into positions within A
            for (int i = 0; i < f.GetLength(0); i++)
            {
                float[] vector = f[i];
                for (int j = 0; j < vector.GetLength(0); j++)
                    A[2 * i, j] = vector[j];
            }
            for (int j = 0; j < A.GetLength(1); j++)
                {
                    for (int i = 1; i < f.GetLength(0); i++)
                        A[(2 * i) - 1, j] = (A[(2 * i) - 2, j] + A[(2 * i), j]) / 2;
                }
            return A;
        }

        public static void WriteArrayAsList(string prefix, float[] f)
        {
            string output = prefix + ": [" + Math.Round(f[0], 2);
            for (int i = 1; i < f.GetLength(0); i++)
                output += ", " + Math.Round(f[i], 2);
            output += "]";
            Debug.WriteLine(output);
        }



        internal static float[,] Subtract(float[,] slice, float[,] p)
        {
            float[,] sub = new float[slice.GetLength(0), slice.GetLength(1)];
            for (int j = 0; j < slice.GetLength(1); j++)
                for (int i = 0; i < slice.GetLength(0); i++)
                    sub[i, j] = slice[i, j] - p[i, j];
            return sub;
        }

        public static float[][,] PrepareDDSMask(float[][,] tumor, int iterations)
        {
            float[][,] mask = new float[tumor.GetLength(0)][,];
            for (int k = 0; k < tumor.GetLength(0); k++)
            {
                float[,] slice = new float[tumor[0].GetLength(0), tumor[0].GetLength(1)];
                for (int i = 0; i < iterations; i++)
                    slice = DilateSlice(tumor[k]);
                mask[k] = Matrix.Subtract(slice, tumor[k]);
            }
            return mask;
        }

        public static float[,] DilateSlice(float[,] p)
        {
            float[,] output = (float[,])p.Clone();
            float top;
            float right;
            float left; float bottom;
            for (int j = 1; j < p.GetLength(1) - 1; j++)
                for (int i = 1; i < p.GetLength(0) - 1; i++)
                {
                    if (p[i, j] > 0)
                        continue;
                    else
                    {
                        top = p[i, j - 1];
                        left = p[i - 1, j];
                        right = p[i + 1, j];
                        bottom = p[i, j + 1];


                        float sum = top + left + right + bottom;
                        if (sum > 0)
                            output[i, j] = 1;
                    }
                }
            return output;
        }

        internal static float FindMax(float[][,] d)
        {
            float max = 0;
            for (int k = 0; k < d.GetLength(0); k++)
                   for (int j = 0; j < d[0].GetLength(1); j++)
                       for (int i = 0; i < d[0].GetLength(0); i++)
                       {
                           float temp = d[k][i, j];
                           if (temp > max)
                               max = temp;
                       }
            return max;
        }

        

        internal static float[][,] LinearlyCombine(float[][,] fj_Tumor, float[][,] fj_CS, int p)
        {
            float[][,] output = new float[fj_Tumor.GetLength(0)][,];
            for (int k = 0; k < fj_Tumor.GetLength(0); k++)
            {
                float[,] temp = Matrix.Zeroes(fj_Tumor[0].GetLength(0), fj_Tumor[0].GetLength(1));
                for (int j = 0; j < temp.GetLength(1); j++)
                    for (int i = 0; i < temp.GetLength(0); i++)
                    {
                        temp[i, j] = fj_Tumor[k][i, j] + p * fj_CS[k][i, j];
                    }
                output[k] = temp;
            }
            return output;
        }
    }
}
