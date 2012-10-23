using System;
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
        //NOT FINISHED YET.
        public string WeightOriginalDS = @"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

__kernel void
AddWeight2OriginalDS(
__global    float * wDS,
__global    float * originalDS,
__global    float * positions,
__global    float * weights,
__global    float * params) //size, dosecalculationthickness
 {
 	//Vector element index, number of workers = slicedose size.
 	int id_z = get_global_id(0); 
    int id_x = get_global_id(1); //represents which slice of size:dosecalculationthickness
    int real_x;
    
    int startz;
    float weight;
    int index;
    int numslices = (int)params[0];
    int size = (int)params[1];
    int dosecalcthick = (int)params[2];
 	
 	//Loop through slices, and then loop through the length of one slicedose array
    for (int s = 0; s < numslices; s++)
    {
        startz = positions[s] - (dosecalcthick / 2);
        weight = weights[s];
        int extra = 0;
        
        if (id_x < (size-1))
        {
            real_x = id_x;
            extra = 1;
        }
        else if (id_x >= (size-1))
        {
            real_x = id_x-(size-1);
            extra = 0;
        }


        index = ((startz + id_z) * size * size) + (real_x * size);

        //each worker handles half a single row (241 elements), 
        

        for (int x = 0; x < (size - 1) / 2; x++)
        {
            wDS[index+(2*x)+extra] += weight*originalDS[index+(2*x)+extra];
        }       
if (id_x >= (size-1))
            wDS[index+size] += weight*originalDS[index + size];

    }
}

";
    }
}
