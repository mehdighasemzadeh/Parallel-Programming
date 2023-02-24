//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "fft.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

__global__ void bit_reversal(float* input , float* output,  unsigned int N)
{   
    int thread_id = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*16384));
 
    float log2N = log2f(N);
    int n = 0;
    int  x = thread_id;
    for (int i = 0; i < log2N; i++)
    {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }

      output[thread_id] = input[n];
    __syncthreads();
}












__global__ void fft(float *x_r , float *x_i ,float *rev_x_r, float *rev_x_i, int N, int stage, int butterfly_width, int step){

	const double twoPIdivN = 2 * PI / N;
	int thread_id = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*16384));
	//======= temp variable========= 
	float wn_r , wn_i;
	float temp1_r , temp1_i ;
	float temp2_r , temp2_i ;
  
	//=======bit r ===============
	/*
	if(stage==1){
	  int r = bit_reversal(thread_id, N);
	  rev_x_r[thread_id] = x_r[r];
	  __syncthreads();
	}
  
  
	if(stage==1){
	  int r = bit_reversal(thread_id, N);
	  rev_x_i[thread_id] = x_i[r];
	  __syncthreads();
	}
	*/
  
  
  
	//============== fft parametr ===============
	int pos = thread_id / butterfly_width * step;
	int j = thread_id % butterfly_width;
	int res = pos + j;
	if (res < N){
		
	  //Wn = e^(-j*2*PI/N) converted with euler's formula(real and imaginary parts)
	  wn_r =  cos(twoPIdivN * j * N / step);
	  wn_i = -sin(twoPIdivN * j * N / step);
  
	  temp1_r = rev_x_r[res];
	  temp1_i = rev_x_i[res];
	  temp2_r = rev_x_r[res + butterfly_width] * wn_r - rev_x_i[res + butterfly_width] * wn_i;
	  temp2_i = rev_x_i[res + butterfly_width] * wn_r + rev_x_r[res + butterfly_width] * wn_i;
  
	  rev_x_r[res]                   = temp1_r + temp2_r;
	  rev_x_i[res]                   = temp1_i + temp2_i;
	  rev_x_r[res + butterfly_width] = temp1_r - temp2_r;
	  rev_x_i[res + butterfly_width] = temp1_i - temp2_i;
		  
	  __syncthreads();
	}
	//===========end if ====================
  
}







void fft_caller(float *x_r , float *x_i , float *rev_x_r , float *rev_x_i , int N)
{

	dim3 blocks(16384,4);
	dim3 threadsPerBlock(1024,1);
 /*
    if (N>=1024 && N< 67108864 ){
      dim3 blocks(N/1024,1);
	  dim3 threadsPerBlock(1024,1);
    }

    if (N==67108864){
     	dim3 blocks(32768,2);
	    dim3 threadsPerBlock(1024,1);
    }
	*/

  	float stages = log2f(N);
  	int butterfly_width, step;

  	if (N > 1)
  	{
    	for (int stage = 1; stage <= stages; stage++)
    {   
        //printf("%d ", stage);
        step = 1 << stage;
        butterfly_width = step >> 1;
        fft<<<blocks, threadsPerBlock>>>(x_r , x_i , rev_x_r, rev_x_i, N, stage, butterfly_width, step);
    }
 
  }
}

  










//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, /*float* X_r_d, float* X_i_d,*/ const unsigned int N, const unsigned int M)
{
	// In this function, both inputs and outputs are on GPU.
	// No need for cudaMalloc, cudaMemcpy or cudaFree.
	
	// set thread count
	dim3 blocks(16384,4);
	dim3 threadsPerBlock(1024,1);
	/*
    if (N>=1024 && N< 67108864 ){
      dim3 blocks(N/1024,1);
	  dim3 threadsPerBlock(1024,1);
    }

    if (N==67108864){
     	dim3 blocks(32768,2);
	    dim3 threadsPerBlock(1024,1);
    }

	*/

	float* rev_x_r;
	float* rev_x_i;


	cudaMalloc((void**)&rev_x_r, N * sizeof(float));
	bit_reversal<<<blocks, threadsPerBlock>>>(x_r_d , rev_x_r, N);
	cudaMemcpy(x_r_d, rev_x_r , N * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaFree(rev_x_r);


	cudaMalloc((void**)&rev_x_i, N * sizeof(float));
	bit_reversal<<<blocks, threadsPerBlock>>>(x_i_d , rev_x_i, N);
	cudaMemcpy(x_i_d, rev_x_i , N * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaFree(rev_x_i);


	fft_caller(rev_x_r , rev_x_i , x_r_d ,x_i_d , N);
	
	



	
}
