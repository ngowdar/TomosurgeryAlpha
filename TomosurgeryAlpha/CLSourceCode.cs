﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TomosurgeryAlpha
{
    public class CL_Sourcecode
    {
        public string ScalarMultiply2DFloat = @"
 __kernel void
ScalarMultiply(
__global	float * A,
__global	float * c,
__global	float * Result
)
 {
 	//Vector element index
 	int id = get_global_id(0);

 	
 	Result[id] = A[id]*c[0];
 }
";
        public string Add3DFloat = @"
 __kernel void
Add(
__global	float * A,
__global	float * B,
__global	float * weight,
__global	float * size,
__global    float * start_z
)
 {
 	//Vector element index
 	int B_loc = get_global_id(0); 		
    int A_loc = start_z[0]*(size[0]*size[1])+B_loc;    
    A[A_loc] = A[A_loc]+(B[B_loc]*weight[0]);
 }

";

        public string AddSliceDose2DoseSpace = @"
__kernel void
AddSliceDose2DoseSpace(
__global	float * weightedDS,
__global	float * slicedose,
__global    float * params)
 {
 	//Vector element index, number of workers = slicedose size.
 	int id = get_global_id(0); 	
    int startz = (int)params[0];
    float weight = params[1];
    int size = (int)params[2];
    int slicenum = (int)params[3];
 	
 	
    weightedDS[(startz*size*size)+id] += slicedose[id]*weight;
  }
";
    }
}
