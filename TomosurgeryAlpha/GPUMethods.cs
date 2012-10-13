using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenCLTemplate;

namespace TomosurgeryAlpha
{
    public static class GPU
    {
        public static bool GPUenabled = false;

        public static float[,] ScalarMultiply(float[,] A, float scalar)
        {
            
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                CL_Sourcecode src = new CL_Sourcecode();
                CLCalc.Program.Compile(new string[] { src.ScalarMultiply2DFloat }); //puts the code into the kernel
                CLCalc.Program.Kernel ScalarMultiply_kernel = new CLCalc.Program.Kernel("ScalarMultiply");


                //Convert parameter variables into 1D arrays
                float[] A_linear = ConvertTo1D(A);
                float[] scalar_linear = new float[1] { scalar };
                float[] result_linear = new float[A_linear.GetLength(0)];

                //Add variables to the OpenCL program
                CLCalc.Program.Variable dev_A_matrix = new CLCalc.Program.Variable(A_linear);
                CLCalc.Program.Variable dev_scalar = new CLCalc.Program.Variable(scalar_linear);
                CLCalc.Program.Variable dev_result = new CLCalc.Program.Variable(result_linear);
                CLCalc.Program.Variable[] args = new CLCalc.Program.Variable[3] { dev_result, dev_A_matrix, dev_scalar };

                //Run the program
                if (ScalarMultiply_kernel != null)
                    ScalarMultiply_kernel.Execute(args, A_linear.GetLength(0));
                dev_result.ReadFromDeviceTo(result_linear);

                return BackTo2D(result_linear, A.GetLength(0));
            }
            else
            {
                Debug.WriteLine("OpenCL code not working!");
                return null;
            }
            

        }

        private static float[,] BackTo2D(float[] result_linear, int size)
        {
            float[,] result = new float[size, size];
            for (int j = 0; j < size; j++)
                for (int i = 0; i < size; i++)
                    result[i, j] = result_linear[j * size + i];
            return result;
        }
        private static float[][,] BackTo3D(float[] result_linear, int size, int size2, int size3)
        {
            float[][,] result = new float[size3][,];
            for (int k = 0; k < size3; k++)
            {
                float[,] slice = new float[size, size2];
                for (int j = 0; j < size2; j++)
                    for (int i = 0; i < size; i++)
                        slice[i, j] = result_linear[k*size*size2+ j * size + i];
                result[k] = slice;
            }
            
            return result;
        }

        private static float[] ConvertTo1D(float[,] A)
        {            
            int size = A.GetLength(0); int size2 = A.GetLength(1);
            float[] result = new float[size*size2];
            for (int j = 0; j < size2; j++)
                for (int i = 0; i < size; i++)
                    result[j * size + i] = (float)A[i, j];
            return result;
        }
        private static float[] ConvertTo1D(float[][,] A)
        {
            int size = A[0].GetLength(0); int size2 = A[0].GetLength(1); int size3 = A.GetLength(0);
            float[] result = new float[size * size2 * size3];
            for (int k = 0; k < size3; k++)
            {
                for (int j = 0; j < size2; j++)
                    for (int i = 0; i < size; i++)
                        result[k*size*size2 + j * size + i] = (float)A[k][i, j];
            }
            return result;
        }

        public static float[][,] ScalarMultiply(float[][,] A, float scalar)
        {
            float[][,] wtf = new float[A.GetLength(0)][,];
            return wtf;
        }


        internal static float[] ScalarMultiply(float[] A, float scalar)
        {
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                CL_Sourcecode src = new CL_Sourcecode();
                CLCalc.Program.Compile(new string[] { src.ScalarMultiply2DFloat }); //puts the code into the kernel
                CLCalc.Program.Kernel ScalarMultiply_kernel = new CLCalc.Program.Kernel("ScalarMultiply");


                //Convert parameter variables into 1D arrays
                float[] A_linear = (float[])A.Clone();
                float[] scalar_linear = new float[1] { scalar };
                float[] result_linear = new float[A_linear.GetLength(0)];

                //Add variables to the OpenCL program
                CLCalc.Program.Variable dev_A_matrix = new CLCalc.Program.Variable(A_linear);
                CLCalc.Program.Variable dev_scalar = new CLCalc.Program.Variable(scalar_linear);
                CLCalc.Program.Variable dev_result = new CLCalc.Program.Variable(result_linear);
                CLCalc.Program.Variable[] args = new CLCalc.Program.Variable[3] { dev_A_matrix, dev_scalar, dev_result };

                //Run the program
                if (ScalarMultiply_kernel != null)
                    ScalarMultiply_kernel.Execute(args, A_linear.GetLength(0));
                dev_result.ReadFromDeviceTo(result_linear);

                return result_linear;
            }
            else
            {
                Debug.WriteLine("OpenCL code not working!");
                return null;
            }

        }
        
        public static float[][,] AddToDoseSpace(float[] slicedose, float[][,] output_ds, int which_z, float weight)
        {
            //int TranslateZBy = SlicePositions[which_slice] - (DoseCalculationThickness / 2); // <- NEED TO CHANGE APPROPRIATELY
            //int startz = SlicePositions[which_slice] - DoseCalculationThickness/2;

            //NEED TO ADD BOUNDARY CONDITION in case slice position trims off some of the dose slab.
            //Then change slicethickness in for loop limit to another variable based on size.
            //Parallel.For(0, DoseCalculationThickness, (z) =>
            //{
            //    int current_z = which_z + z;
            //    if (output_ds[current_z] == null)
            //        output_ds[current_z] = Matrix.Zeroes(volume[0].GetLength(0), volume[0].GetLength(1));
            //    output_ds[current_z] = Matrix.Add(output_ds[current_z], GrabSlice(slicedose, z, volume[0].GetLength(0), volume[0].GetLength(1)));
            //});

            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                CL_Sourcecode src = new CL_Sourcecode();
                CLCalc.Program.Compile(new string[] { src.Add3DFloat }); //puts the code into the kernel
                CLCalc.Program.Kernel Add_kernel = new CLCalc.Program.Kernel("Add");


                //Convert parameter variables into 1D arrays
                float[] A_linear = ConvertTo1D(output_ds);
                float[] size = new float[3] { output_ds[0].GetLength(0), output_ds[0].GetLength(1), output_ds.GetLength(0) };
                float[] which_z_linear = new float[1] { which_z };
                float[] weight_linear = new float[1] { weight };
                //float[] result_linear = (float[])A_linear.Clone();

                //Add variables to the OpenCL program
                CLCalc.Program.Variable dev_A_matrix = new CLCalc.Program.Variable(A_linear);
                CLCalc.Program.Variable dev_slicedose = new CLCalc.Program.Variable(slicedose);
                CLCalc.Program.Variable dev_weight = new CLCalc.Program.Variable(weight_linear);
                CLCalc.Program.Variable dev_size = new CLCalc.Program.Variable(size);
                CLCalc.Program.Variable dev_whichz = new CLCalc.Program.Variable(which_z_linear);
                //CLCalc.Program.Variable dev_result = new CLCalc.Program.Variable(result_linear);
                CLCalc.Program.Variable[] args = new CLCalc.Program.Variable[5] { dev_A_matrix, dev_slicedose, dev_weight, dev_size, dev_whichz };

                //Run the program
                if (Add_kernel != null)
                    Add_kernel.Execute(args, slicedose.GetLength(0));
                dev_A_matrix.ReadFromDeviceTo(A_linear);
                return BackTo3D(A_linear, (int)size[0],(int)size[1],(int)size[2]);
            }
            else
            {
                Debug.WriteLine("OpenCL code not working!");
                return null;
            }
            
        }
    }
}
