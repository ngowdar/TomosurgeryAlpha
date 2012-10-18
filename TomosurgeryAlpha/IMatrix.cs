using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.ComponentModel;

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
            int[] size = new int[3]{f.GetLength(0),f.GetLength(1),f.GetLength(2)};
            float[][,] ff = new float[size[2]][,];
            for (int k = 0; k < z; k++)
                for (int j = 0; j < y; j++)
                    for (int i = 0; i < x; i++)                    
                        ff[k][i,j] = f[k*x*y + j*x + i];
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
            float divisor = 1 / FindMax(d);
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
            for (int k = 0; k < d.GetLength(0); k++)
            {
                float divisor = 1 / max;
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

        internal static float[,] Subset(float[] A, int centerx, int centery, int subsetsize)
        {
            float[,] Window;
            int Asize = (int)Math.Sqrt(A.GetLength(0));
            int halfwindow = (subsetsize / 2);
            int startx = centerx - halfwindow; int starty = centery - halfwindow;
            int endx = centerx + halfwindow; int endy = centery + halfwindow;
            bool xfits = false; bool yfits = false;

            if (endx >= Asize)
            {
                endx = (Asize - 1);
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
                endy = (Asize - 1);
                yfits = false;
            }
            else if (endy < Asize)
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
                    Window[i, j] = A[i + startx + (j + starty)*Asize];
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

        
    }
}
