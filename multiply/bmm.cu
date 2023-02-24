//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// TILEX and TILEY are used to set the number of threads in a CUDA block 
#define TILEX 32
#define TILEY 16

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

dim3 getDimGrid(const int m, const int n) {
	dim3 dimGrid(n/TILEX,n/TILEY);
	return dimGrid;
}
dim3 getDimBlock(const int m, const int n) {
	dim3 dimBlock(TILEX,TILEY);
	return dimBlock;
}
__global__ void kernelFunc(float* ad, float* bd, float* cd, const int m, const int n) {

	// write your GPU kernel function here



	if (TILEX > TILEY){
    __shared__ float As[TILEY][TILEY];
  	__shared__ float Bs[TILEY ][TILEX];
 
  	int x = threadIdx.x + blockIdx.x * blockDim.x; 
  	int y = threadIdx.y + blockIdx.y * blockDim.y; 

  	float Pvalue = 0;
  	for (int m1 = 0; m1 < n/TILEY; m1++){
    		if(threadIdx.x<TILEY){	 
    		As[threadIdx.y][threadIdx.x] = ad[y * n + (m1 * TILEY + threadIdx.x)]; 
    		}

    		Bs[threadIdx.y][threadIdx.x] = bd[(m1 * TILEY + threadIdx.y) * n + x];
   	 	__syncthreads();


    		for (int k = 0; k < TILEY; k++)
      		Pvalue += As[threadIdx.y][k] * Bs[k][threadIdx.x];
   		__syncthreads();
  	}
  
  // write back to the global memory
  	cd[y * n + x] = Pvalue;   
  } 
    


// TILEX < TILEY :
    
   if (TILEX < TILEY){
    __shared__ float As[TILEY][TILEX];
  	__shared__ float Bs[TILEX ][TILEX];
 
  	int x = threadIdx.x + blockIdx.x * blockDim.x; 
  	int y = threadIdx.y + blockIdx.y * blockDim.y; 

  	float Pvalue = 0;
  	for (int m1 = 0; m1 < n/TILEX; m1++){
    			 
    		As[threadIdx.y][threadIdx.x] = ad[y * n + (m1 * TILEX + threadIdx.x)]; 
    		
      if(threadIdx.y<TILEX){
    	  Bs[threadIdx.y][threadIdx.x] = bd[(m1 * TILEX + threadIdx.y) * n + x];
      }
   	 	__syncthreads();


    		for (int k = 0; k < TILEX; k++)
      		Pvalue += As[threadIdx.y][k] * Bs[k][threadIdx.x];
   		__syncthreads();
  	}
  
  // write back to the global memory
  	cd[y * n + x] = Pvalue;   
  }






  // TILEX == TILEY :
    
   if (TILEX == TILEY){
    __shared__ float As[TILEY][TILEX];
  	__shared__ float Bs[TILEY][TILEX];
 
  	int x = threadIdx.x + blockIdx.x * blockDim.x; 
  	int y = threadIdx.y + blockIdx.y * blockDim.y; 

  	float Pvalue = 0;
  	for (int m1 = 0; m1 < n/TILEY; m1++){
    			 
    		As[threadIdx.y][threadIdx.x] = ad[y * n + (m1 * TILEY + threadIdx.x)]; 
    	  Bs[threadIdx.y][threadIdx.x] = bd[(m1 * TILEX + threadIdx.y) * n + x];
      
   	 	__syncthreads();

    		for (int k = 0; k < TILEX; k++)
      		Pvalue += As[threadIdx.y][k] * Bs[k][threadIdx.x];
   		__syncthreads();
  	}
  
  // write back to the global memory
  	cd[y * n + x] = Pvalue;   
  }





}
