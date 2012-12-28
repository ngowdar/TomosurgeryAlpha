using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenCLTemplate;

namespace TomosurgeryAlpha
{
    public class LinMatrix
    {
        public float[] M;

        private CLCalc.Program.Kernel Matrix;
        public int X;
        public int Y;
        public int Z;
        public int stride;

        public LinMatrix(int x, int y, int z)
        {
            SetSizeofM(x, y, z);
            GetStride();
        }

        public LinMatrix(float[][,] d)
        {
            Z = d.GetLength(0);
            X = d[0].GetLength(0);
            Y = d[0].GetLength(1);
            M = Convertto1D(d);
            GetStride();
        }

        public void SetSizeofM(int x, int y, int z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public void GetStride()
        {
            stride = X*Y;
        }

        public float[] Convertto1D(float[][,] d)
        {
            float[] result;
            int z = d.GetLength(0);
            int x = d[0].GetLength(0);
            int y = d[0].GetLength(1);
            result = new float[z*x*y];
            for (int i = 0; i < z; i++)
                for (int j = 0; j < y; j++)
                    for (int k = 0; k < x; k++)
                        result[(i*x*y) + (j*x) + k] = d[i][k, j];
            return result;
        }

        public float[] GrabLinearSlice(int z)
        {
            var slice = new float[stride];
            for (int i = 0; i < stride; i++)
                slice[i] = M[i + (z*stride)];
            return slice;
        }

        public float[,] Grab2DSlice(int z)
        {
            var slice = new float[X,Y];
            for (int j = 0; j < Y; j++)
                for (int i = 0; i < X; i++)
                    slice[i, j] = M[z*X*Y + j*X + i];
            return slice;
        }

        public float[][,] ConvertToJagged()
        {
            var P = new float[Z][,];
            for (int i = 0; i < Z; i++)
                P[i] = Grab2DSlice(i);
            return P;
        }

        public float[] Add_NoGPU(float[] B, int startx, int starty, int startz, int bx, int by)
        {
            var start = new float[3] {startx, starty, startz}; //where to start the B matrix within the bigger one.
            var size = new float[6] {X, Y, Z, B.GetLength(0), bx, by};
            var csResult = new float[X*Y*Z];

            int widthA = X;
            int heightA = Y;
            int depthA = Z;
            int widthB = bx;
            int heightB = by;
            int depthB = B.GetLength(0)/(bx*by);
            int startA = startz*widthA*heightA + starty*widthA + startx;

            int endx = startx + widthB;
            int endy = starty + heightB;
            int endz = startz + depthB;

            if (endx >= widthA)
                endx = widthA;
            if (endy >= heightA)
                endy = heightA;
            if (endz >= depthA)
                endz = depthA;
            for (int id = 0; id < B.GetLength(0); id++)
            {
                int Ax = id%widthA;
                int Ay = ((id%(widthA*heightA)) - Ax)/widthA;
                int Az = (id - Ax - (Ay*widthA))/(widthA*heightA);

                if (Ax >= startx && Ay >= starty && Az >= startz)
                {
                    if (Ax < endx && Ay < endy && Az < endz)
                    {
                        int Bx = Ax - startx;
                        int By = Ay - starty;
                        int Bz = Az - startz;
                        int B_id = Bz*(heightB*widthB) + By*(widthB) + Bx;
                        csResult[id] = M[id] + B[B_id];
                    }
                }
                else
                    csResult[id] = M[id];
            }
            return csResult;
        }

        public float[] Add(float[] B, int startx, int starty, int startz, int Bx, int By)
        {
            var start = new float[3] {startx, starty, startz};
            var size = new float[6] {X, Y, Z, B.GetLength(0), Bx, By};
            var csResult = new float[X*Y*Z];
            /*In order to avoid OOM, will start at the startz value, and pass
            in one 2D matrix at a time.
             * */
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                var source = new CLSource();
                CLCalc.Program.Compile(new[] {source.AddSubset_3string});
                Matrix = new CLCalc.Program.Kernel("AddSubset3");
            }
            var dev_A = new CLCalc.Program.Variable(M);
            var dev_B = new CLCalc.Program.Variable(B);
            var dev_start = new CLCalc.Program.Variable(start);
            var dev_size = new CLCalc.Program.Variable(size);
            var dev_Result = new CLCalc.Program.Variable(csResult);
            var args = new CLCalc.Program.Variable[5] {dev_A, dev_B, dev_start, dev_size, dev_Result};
            Matrix.Execute(args, new[] {X*Y*Z});
            dev_Result.ReadFromDeviceTo(csResult);
            M = csResult;
            return csResult;
        }

        public float[][,] Convertto3D(float[] oned)
        {
            var r = new float[Z][,];
            var temp = new float[X,Y];
            for (int k = 0; k < Z; k++)
            {
                temp = new float[X,Y];
                for (int j = 0; j < Y; j++)
                    for (int i = 0; i < X; i++)
                        temp[i, j] = oned[k*X*Y + j*X + i];
                r[k] = temp;
            }
            return r;
        }

        public float[,] Convertto2D(float[] oned)
        {
            var r = new float[X,Y];
            for (int j = 0; j < Y; j++)
                for (int i = 0; i < X; i++)
                    r[i, j] = oned[j*X + i];
            return r;
        }
    }

    public class CLSource
    {
        public string AddSubset_2string = @"

  __kernel void
AddSubset2(
__global	float * A,
__global	float * B,
__global	float * start_coords,
__global	float * size,
__global	float * Result
)
 {
 	//Vector element index
 	int id = get_global_id(0); 		

	int widthA = size[0];
 	int heightA = size[1];
    
 	int widthB = size[3];
 	int heightB = size[4];
    
 	//Starting position within the A array, corresponding to the first pixel of the B matrix.
	int start = start_coords[1]*widthA + start_coords[0];


	int endx = start_coords[0] + widthB;
 	int endy = start_coords[1] + heightB;
 	


//Choose end-coordinates based on size of B.
 	if (endx >= widthA)
 		endx = widthA;
 	if (endy >= heightA)
 		endy = heightA;
 	
 	int Ax = id % widthA;
 	int Ay = (id-Ax)/widthA;
 	

if (Ax >= start_coords[0] && Ay >= start_coords[1])
{
    if (Ax < endx && Ay < endy)
    {
        int Bx = Ax - start_coords[0];
        int By = Ay - start_coords[1];
        
        int B_id = By*(widthB)+Bx;
        Result[id] = A[id]+B[B_id];
    }
}
else
    Result[id] = A[id];

 }";
        public string AddSubset_3string = @"

  __kernel void
AddSubset3(
__global	float * A,
__global	float * B,
__global	float * start_coords,
__global	float * size,
__global	float * Result
)
 {
 	//Vector element index
 	int id = get_global_id(0); 		

	int widthA = size[0];
 	int heightA = size[1];
    int depthA = size[2];
 	int widthB = size[3];
 	int heightB = size[4];
    int depthB = size[5];
 	//Starting position within the A array, corresponding to the first pixel of the B matrix.
	int start = start_coords[2]*widthA*heightA + start_coords[1]*widthA + start_coords[0];


	int endx = start_coords[0] + widthB;
 	int endy = start_coords[1] + heightB;
 	int endz = start_coords[2] + depthB;


//Choose end-coordinates based on size of B.
 	if (endx >= widthA)
 		endx = widthA;
 	if (endy >= heightA)
 		endy = heightA;
 	if (endz >= depthA)
 		endz = depthA;
 	
 	int Ax = id % widthA;
 	int Ay = ((id % (widthA*heightA))-Ax)/widthA;
 	int Az = (id-Ax-(Ay*widthA))/(widthA*heightA);	

if (Ax >= start_coords[0] && Ay >= start_coords[1] && Az >= start_coords[2])
{
    if (Ax < endx && Ay < endy && Az < endz)
    {
        int Bx = Ax - start_coords[0];
        int By = Ay - start_coords[1];
        int Bz = Az - start_coords[2];
        int B_id = Bz*(heightB*widthB)+By*(widthB)+Bx;
        Result[id] = A[id]+B[B_id];
    }
}
else
    Result[id] = A[id];

 }";
        public string MM2D = @"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

 __kernel void
MultiplyElements2_CL(
__global	float * A,
__global	float * B,
__global	float * start_coords,
__global	float * size,
__global	float * Result)
 {
 	//Vector element index
 	int id = get_global_id(0); 	

 	int widthA = size[0];
 	int heightA = size[1];    
 	int widthB = size[2];
 	int heightB = size[3];
    

    //Starting position within the A array, corresponding to the first pixel of the B matrix.
	int start = start_coords[1]*widthA + start_coords[0];


	int endx = start_coords[0] + widthB;
 	int endy = start_coords[1] + heightB; 	


    //Choose end-coordinates based on size of B.
 	if (endx >= widthA)
 		endx = widthA;
 	if (endy >= heightA)
 		endy = heightA; 	

 	int Ax = id % widthA;
 	int Ay = (id-Ax)/widthA; 	

if (Ax >= start_coords[0] && Ay >= start_coords[1])
{
    if (Ax < endx && Ay < endy)
    {
        int Bx = Ax - start_coords[0];
        int By = Ay - start_coords[1];
        
        int B_id = By*(widthB)+Bx;
        Result[id] = A[id]*B[B_id];
    }
}
else
    Result[id] = A[id];
}
";
        public string MM3D = @"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

 __kernel void
MultiplyElements3_CL(
__global	float * A,
__global	float * B,
__global	float * start_coords,
__global	float * size,
__global	float * Result)
 {
 	//Vector element index
 	int id = get_global_id(0);
 	//int idy = get_global_id(1);
 	//In order to use two global ids, need to 
 	//set the group size = to one row.

 	int widthA = size[0];
 	int heightA = size[1];
    int depthA = size[2];
 	int widthB = size[3];
 	int heightB = size[4];
    int depthB = size[5];

    //Starting position within the A array, corresponding to the first pixel of the B matrix.
	int start = start_coords[2]*widthA*heightA + start_coords[1]*widthA + start_coords[0];


	int endx = start_coords[0] + widthB;
 	int endy = start_coords[1] + heightB;
 	int endz = start_coords[2] + depthB;


    //Choose end-coordinates based on size of B.
 	if (endx >= widthA)
 		endx = widthA;
 	if (endy >= heightA)
 		endy = heightA;
 	if (endz >= depthA)
 		endz = depthA;	

 	int Ax = id % widthA;
 	int Ay = ((id % (widthA*heightA))-Ax)/widthA;
 	int Az = (id-Ax-(Ay*widthA))/(widthA*heightA);

if (Ax >= start_coords[0] && Ay >= start_coords[1] && Az >= start_coords[2])
{
    if (Ax < endx && Ay < endy && Az < endz)
    {
        int Bx = Ax - start_coords[0];
        int By = Ay - start_coords[1];
        int Bz = Az - start_coords[2];
        int B_id = Bz*(heightB*widthB)+By*(widthB)+Bx;
        Result[id] = A[id]*B[B_id];
    }
}
else
    Result[id] = A[id];
}
";
    }
}