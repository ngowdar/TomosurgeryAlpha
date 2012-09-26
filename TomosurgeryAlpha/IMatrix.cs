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
        public static float[] ScalarMultiply(ref float[] f, float m)
        {
            for (int i = 0; i < f.GetLength(0); i++)
                f[i] = f[i] * m;
            return f;
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

        public static float[][,] EnlargeAndCenter(float[][,] dd, int enlargement, int startx, int starty, int startz)
        {
            LinMatrix m = new LinMatrix(dd);
            m.SetSizeofM(m.X + enlargement, m.Y + enlargement, m.Z + enlargement);
            //m.Add(m.Convertto1D(dd), enlargement / 2, enlargement / 2, enlargement / 2, dd[0].GetLength(0), dd[1].GetLength(1));
            m.Add_NoGPU(m.Convertto1D(dd), enlargement / 2, enlargement / 2, enlargement / 2, dd[0].GetLength(0), dd[1].GetLength(1));
            
            return m.ConvertToJagged();
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

        internal static object Subset(float[,] slice, int startx, int starty, int p, int p_2)
        {
            throw new NotImplementedException();
        }

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
                for (int k = 0; k < A.GetLength(0); k++)
                    product[k] = MultiplyElements(A[k], B[k]);
            }
            return product;
        }

        internal static object Subset(float[] ds, int startx, int starty, int p, int p_2)
        {
            throw new NotImplementedException();
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
            else if (endx < AA.GetLength(0))
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= AA.GetLength(1))
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
            for (int i = 0; i < endx - startx; i++)
                for (int j = 0; j < endy - starty; j++)
                    AA[i + startx, j + starty] = A[i + startx,j + starty] * b[i, j];
            return AA;
        }

        internal static float[,] MultiplySubset(float[] A, float[,] b, int centerx, int centery, int sizex, int sizey)
        {
            //Check if b will fit inside A
            float[,] AA = new float[sizex, sizey];
            int halfB = (b.GetLength(0) - 1) / 2;
            int startx = centerx - halfB; int starty = centery - halfB;
            int endx = centerx + halfB; int endy = centery + halfB;
            bool xfits = false; bool yfits = false;

            if (endx >= A.GetLength(0))
            {
                endx = (AA.GetLength(0) - 1);
                xfits = false;
            }
            else if (endx < AA.GetLength(0))
            { xfits = true; }

            if (startx < 0)
            {
                startx = 0;
                xfits = false;
            }
            else if (startx >= 0)
            { xfits = true; }

            if (endy >= AA.GetLength(1))
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
            for (int i = 0; i < endx - startx; i++)
                for (int j = 0; j < endy - starty; j++)
                    AA[i + startx, j + starty] = A[(i + startx)+(j + starty)*sizex] * b[i, j];
            return AA;
        }

        internal static float[] Normalize(float[] img)
        {
            float[] n = new float[img.GetLength(0)];
            float max = img.Max();
            for (int i = 0; i < img.GetLength(0); i++)
                n[i] = img[i] / max;
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
