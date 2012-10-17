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
        public static float[][,] BackTo3D(float[] result_linear, int size, int size2, int size3)
        {
            int z = result_linear.GetLength(0) / (size*size2);
            float[][,] result = new float[z][,];
            for (int k = 0; k < z; k++)
            {
                float[,] slice = new float[size, size2];
                for (int j = 0; j < size2; j++)
                    for (int i = 0; i < size; i++)
                        slice[i, j] = result_linear[(k * size * size2) + (j * size) + i];
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
        public static float[] ConvertTo1D(float[][,] A)
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

        public static float[] WeightOriginalDS(int[] SlicePositions, double[] weights, int[] size, int DCT, string folderpath)
        {
            
            string filename = "OriginalDS.txt";
            string path = System.IO.Path.Combine(folderpath, filename);
            float[] originalds = Matrix.Normalize(PathSet.ReadDoseSpaceFromFile(path));
            Debug.WriteLine("oDS sum: " + originalds.Sum());
            Debug.WriteLine("oDS sum normalized: " + Matrix.Normalize(originalds).Sum());
            float[] wDS = new float[originalds.GetLength(0)];

            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                CL_Sourcecode src = new CL_Sourcecode();
                CLCalc.Program.Compile(new string[] { src.WeightOriginalDS }); //puts the code into the kernel               
                CLCalc.Program.Kernel AddSliceDose_kernel = new CLCalc.Program.Kernel("AddWeight2OriginalDS");

                //Convert parameter variables into 1D arrays 
                float[] positions = new float[SlicePositions.GetLength(0)];
                float[] w = new float[SlicePositions.GetLength(0)];
                for (int n = 0; n < SlicePositions.GetLength(0); n++)
                {
                    positions[n] = (float)SlicePositions[n];
                    w[n] = (float)weights[n];
                }                
                float[] param = new float[3] { SlicePositions.GetLength(0), size[0], PathSet.DCT };

                //Add variables to the OpenCL program
                CLCalc.Program.Variable dev_wDS = new CLCalc.Program.Variable(wDS);
                CLCalc.Program.Variable dev_oDS = new CLCalc.Program.Variable(originalds);
                CLCalc.Program.Variable dev_positions = new CLCalc.Program.Variable(positions);
                CLCalc.Program.Variable dev_weights = new CLCalc.Program.Variable(w);
                CLCalc.Program.Variable dev_param = new CLCalc.Program.Variable(param);

                //CLCalc.Program.Variable dev_result = new CLCalc.Program.Variable(result_linear);
                CLCalc.Program.Variable[] args = new CLCalc.Program.Variable[5] { dev_wDS, dev_oDS, dev_positions, dev_weights, dev_param };

                //Run the program
                if (AddSliceDose_kernel != null)
                    AddSliceDose_kernel.Execute(args, new int[2]{PathSet.DCT, size[0]});
                dev_wDS.ReadFromDeviceTo(wDS);
            }     
            return wDS;
            

        }


        public static float[] PrepareDoseSpace(float[] output_ds, int[] SlicePositions, double[] weights, int[] size, int DCT, string folderpath)
        {            
            float[] slicedose;
            string filename = "";
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.Unknown)
                CLCalc.InitCL();
            if (CLCalc.CLAcceleration == CLCalc.CLAccelerationType.UsingCL)
            {
                CL_Sourcecode src = new CL_Sourcecode();
                CLCalc.Program.Compile(new string[] { src.AddSliceDose2DoseSpace }); //puts the code into the kernel

                for (int s = 0; s < SlicePositions.GetLength(0); s++)
                {
                    
                    filename = String.Concat("slice_", s);
                    string path = System.IO.Path.Combine(folderpath, filename);
                    slicedose = Matrix.Normalize(PathSet.ReadSliceDoseFromFile(path));

                    
                    
                    CLCalc.Program.Kernel AddSliceDose_kernel = new CLCalc.Program.Kernel("AddSliceDose2DoseSpace");
                    
                    //Convert parameter variables into 1D arrays                 
                    float position = (float)SlicePositions[s] - (DCT / 2);                    
                    float[] param = new float[4] { position, (float)weights[s], size[0], s };                    

                    //Add variables to the OpenCL program
                    CLCalc.Program.Variable dev_wDS = new CLCalc.Program.Variable(output_ds);
                    CLCalc.Program.Variable dev_slicedose = new CLCalc.Program.Variable(slicedose);                    
                    CLCalc.Program.Variable dev_params = new CLCalc.Program.Variable(param);
                    
                    //CLCalc.Program.Variable dev_result = new CLCalc.Program.Variable(result_linear);
                    CLCalc.Program.Variable[] args = new CLCalc.Program.Variable[3] { dev_wDS, dev_slicedose, dev_params};

                    //Run the program
                    if (AddSliceDose_kernel != null)
                        AddSliceDose_kernel.Execute(args, slicedose.GetLength(0));
                    dev_wDS.ReadFromDeviceTo(output_ds);                    
                }
                return output_ds;
            }
            else
            {
                Debug.WriteLine("OpenCL code not working!");
                return null;
            }
            
        }
    }
}
