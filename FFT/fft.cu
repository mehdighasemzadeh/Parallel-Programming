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

__device__ int bit_reversal(unsigned int N , unsigned int M , int i)
{   
    int thread_id = i;
    int log2N = M;
    int n = 0;
    int  x = thread_id;
    for (int i = 0; i < log2N; i++)
    {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n ; 
}


__global__ void bit_reversal_helper(float* inputr,float* inputi , unsigned int N , unsigned int M ){
  int thread_id = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*32768));
  int j ;
  float temp1 ;
  float temp2 ; 
  int i =  thread_id ; 
  j = bit_reversal(N, M,i);
    if (i<j){
      temp1 = inputr[i] ; 
      inputr[i] = inputr[j];
      inputr[j] = temp1 ;

      temp2 = inputi[i] ; 
      inputi[i] = inputi[j];
      inputi[j] = temp2 ;

    }
}





__global__ void fft(float *rev_x_r, float *rev_x_i, int N, int stage, int butterfly_width, int step){

	const double twoPIdivN = 2 * PI / N;
	int thread_id = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(512*32768));
	//======= temp variable========= 
	float wn_r , wn_i;
	float temp1_r , temp1_i ;
	float temp2_r , temp2_i ;
  
	//============== fft parametr ===============
	int pos = (thread_id / butterfly_width) * step;
	int j = thread_id % butterfly_width;
	int res = pos + j;
	if (res < N){
		
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
		  
	  //__syncthreads();
	}
	//===========end if ====================
  
}













__global__ void fft_SM( float *rev_x_rg, float *rev_x_ig, int N){

	const double twoPIdivN = 2 * PI / N;
	int thread_id = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(512*32768));
  __shared__ float rev_x_r[1024];
  __shared__ float rev_x_i[1024];
  rev_x_r[2*tx] = rev_x_rg[2*thread_id];
  rev_x_r[2*tx+1] = rev_x_rg[2*thread_id+1];
  rev_x_i[2*tx] = rev_x_ig[2*thread_id];
  rev_x_i[2*tx+1] = rev_x_ig[2*thread_id+1];
  __syncthreads();

	//======= temp variable=========
  float wn_r , wn_i;
	float temp1_r , temp1_i ;
	float temp2_r , temp2_i ;

  int step;
  int stages = 10;
  int butterfly_width;
  for (int stage = 1; stage <= stages; stage++)
  {   
    step = 1 << stage;
    butterfly_width = step >> 1; 

  
  
	//============== fft parametr ===============
	int pos = (tx / butterfly_width) * step;
	int j = tx % butterfly_width;
	int res = pos + j;
	if (res < N){
		
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
	}

  rev_x_rg[2*thread_id] = rev_x_r[2*tx];
  rev_x_rg[2*thread_id+1] = rev_x_r[2*tx+1];
  rev_x_ig[2*thread_id] = rev_x_i[2*tx];
  rev_x_ig[2*thread_id+1] = rev_x_i[2*tx+1];
	//===========end if ====================


  
}







void fft_helper(float *rev_x_r , float *rev_x_i , int N , unsigned int M)
{
	  dim3 blocks;
	  dim3 threadsPerBlock;

    if(N<=1024){
    blocks.x = 1 ;
	  blocks.y = 1 ;
    threadsPerBlock.x = N/2;
    threadsPerBlock.y = 1 ;
    }
    if (N>1024 && N< 67108864 ){
    blocks.x = N/1024 ;
	  blocks.y = 1 ;
    threadsPerBlock.x = 512;
    threadsPerBlock.y = 1 ;
    }

    if (N==67108864){
    blocks.x = 32768 ;
	  blocks.y = 2 ;
    threadsPerBlock.x = 512;
    threadsPerBlock.y = 1 ;
    }
	
    //======== 10 frist stage =======
    fft_SM<<<blocks, threadsPerBlock>>>(rev_x_r, rev_x_i, N);
   //============stage 11 to end================================
    int stages = M;
  	int butterfly_width, step;

  	if (N > 1)
  	{
    	for (int stage = 11; stage <= stages; stage++)
    {   
        step = 1 << stage;
        butterfly_width = step >> 1;
        fft<<<blocks, threadsPerBlock>>>(rev_x_r, rev_x_i, N, stage, butterfly_width, step);
    }
 
  }

}

  










//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, /*float* X_r_d, float* X_i_d,*/ const unsigned int N, const unsigned int M)
{
	// In this function, both inputs and outputs are on GPU.
	// No need for cudaMalloc, cudaMemcpy or cudaFree.
	
	// set thread count
    dim3 blocks;
	  dim3 threadsPerBlock;
    if(N<1024){
    blocks.x = 1 ;
	  blocks.y = 1 ;
    threadsPerBlock.x = N;
    threadsPerBlock.y = 1 ;
  }
    if (N>=1024 && N< 67108864 ){
    blocks.x = N/1024 ;
	  blocks.y = 1 ;
    threadsPerBlock.x = 1024;
    threadsPerBlock.y = 1 ;
    }

    if (N==67108864){
    blocks.x = 32768 ;
	  blocks.y = 2 ;
    threadsPerBlock.x = 1024;
    threadsPerBlock.y = 1 ;
    }

	

bit_reversal_helper<<<blocks, threadsPerBlock>>>(x_r_d , x_i_d , N , M);
fft_helper(x_r_d ,x_i_d , N ,M);
	
	



	
}
