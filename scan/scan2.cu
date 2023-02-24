// ONLY MODIFY THIS FILE

#include "scan2.h"
#include "gpuerrors.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

__global__ void scanl(float* ad, float* cd ,  float* cqd){
		int n =1024;
		__shared__ float As[1024];
  	__shared__ float Cs[1024];
		int offset = 1;
		int x = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*32768));

		As[threadIdx.x] = ad[x];
		// reduce part---------
		for (int d=512 ; d > 0 ; d = d/2){
			__syncthreads() ;
			if ( threadIdx.x < d){
				int ai = offset * (2*  threadIdx.x +1) -1;
				int bi = offset * (2*  threadIdx.x +2) -1;
				As[bi] += As [ai];
			}
			offset *= 2;
			}
		// end reduce part
		
		
		
		if ( threadIdx.x == 0){
				Cs[threadIdx.x + n - 1] = As [threadIdx.x + n - 1] ;
				As [threadIdx.x + n - 1] = 0 ; 
		}


		for (int d = 1; d < n; d *= 2){
				offset /= 2 ;
				__syncthreads() ;
				if ( threadIdx.x < d){
					int ai = offset * (2* threadIdx.x +1) -1;
					int bi = offset * (2* threadIdx.x +2) -1;
					float t = As[ai ];
					As[ai] = As[bi ];
					As[bi] += t;
				}
		}

		__syncthreads () ;

		if (threadIdx.x >=1){
				Cs[threadIdx.x - 1] = As[threadIdx.x] ;
		}
		__syncthreads () ;
		cd[x] = Cs[threadIdx.x];
		if(threadIdx.x == 0 ){
		cqd[blockIdx.x + 32768*blockIdx.y ] = Cs[threadIdx.x + 1023];
		}

			
}




__global__ void add(float* cd ,  float* cqd){
		int x = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*32768));
		if (blockIdx.x + 32768*blockIdx.y>0){
		cd[x] = cd[x] + cqd[blockIdx.x + 32768*blockIdx.y -1];
		}
 		
}



__global__ void scan_s(float* ad, float* cd ,const int n ){
		int x = threadIdx.x;
		int i ;
		cd[x] = ad[x];
		for(i=1 ; i<n ;i++){
				cd[i+x] = cd[i+x -1] + ad[i];
		}

		
}




__global__ void add_l(float* cd ,  const float back){
		int x = (threadIdx.x + blockIdx.x * blockDim.x) + (blockIdx.y*(1024*32768));
		cd[x] += back ; 
		}





__global__ void scan128(float* ad, float* cd){
		int n =128;
		__shared__ float As[128];
  	__shared__ float Cs[128];
		int offset = 1;
		int x = threadIdx.x + blockIdx.x * blockDim.x;

		As[threadIdx.x] = ad[x];
		// reduce part---------
		for (int d=64 ; d > 0 ; d = d/2){
			__syncthreads() ;
			if ( threadIdx.x < d){
				int ai = offset * (2*  threadIdx.x +1) -1;
				int bi = offset * (2*  threadIdx.x +2) -1;
				As[bi] += As [ai];
			}
			offset *= 2;
			}
		// end reduce part
		
		
		
		if ( threadIdx.x == 0){
				Cs[threadIdx.x + n - 1] = As [threadIdx.x + n - 1] ;
				As [threadIdx.x + n - 1] = 0 ; 
		}


		for (int d = 1; d < n; d *= 2){
				offset /= 2 ;
				__syncthreads() ;
				if ( threadIdx.x < d){
					int ai = offset * (2* threadIdx.x +1) -1;
					int bi = offset * (2* threadIdx.x +2) -1;
					float t = As[ai ];
					As[ai] = As[bi ];
					As[bi] += t;
				}
		}

		__syncthreads () ;

		if (threadIdx.x >=1){
				Cs[threadIdx.x - 1] = As[threadIdx.x] ;
		}
		__syncthreads () ;
		cd[x] = Cs[threadIdx.x];			
}









__global__ void scan64(float* ad, float* cd){
		int n =64;
		__shared__ float As[64];
  	__shared__ float Cs[64];
		int offset = 1;
		int x = threadIdx.x + blockIdx.x * blockDim.x;

		As[threadIdx.x] = ad[x];
		// reduce part---------
		for (int d=32 ; d > 0 ; d = d/2){
			__syncthreads() ;
			if ( threadIdx.x < d){
				int ai = offset * (2*  threadIdx.x +1) -1;
				int bi = offset * (2*  threadIdx.x +2) -1;
				As[bi] += As [ai];
			}
			offset *= 2;
			}
		// end reduce part
		
		
		
		if ( threadIdx.x == 0){
				Cs[threadIdx.x + n - 1] = As [threadIdx.x + n - 1] ;
				As [threadIdx.x + n - 1] = 0 ; 
		}


		for (int d = 1; d < n; d *= 2){
				offset /= 2 ;
				__syncthreads() ;
				if ( threadIdx.x < d){
					int ai = offset * (2* threadIdx.x +1) -1;
					int bi = offset * (2* threadIdx.x +2) -1;
					float t = As[ai ];
					As[ai] = As[bi ];
					As[bi] += t;
				}
		}

		__syncthreads () ;

		if (threadIdx.x >=1){
				Cs[threadIdx.x - 1] = As[threadIdx.x] ;
		}
		__syncthreads () ;
		cd[x] = Cs[threadIdx.x];			
}









__global__ void scan32(float* ad, float* cd){
		int n =32;
		__shared__ float As[32];
  		__shared__ float Cs[32];
		int offset = 1;
		int x = threadIdx.x + blockIdx.x * blockDim.x;

		As[threadIdx.x] = ad[x];
		// reduce part---------
		for (int d=16 ; d > 0 ; d = d/2){
			__syncthreads() ;
			if ( threadIdx.x < d){
				int ai = offset * (2*  threadIdx.x +1) -1;
				int bi = offset * (2*  threadIdx.x +2) -1;
				As[bi] += As [ai];
			}
			offset *= 2;
			}
		// end reduce part
		
		
		
		if ( threadIdx.x == 0){
				Cs[threadIdx.x + n - 1] = As [threadIdx.x + n - 1] ;
				As [threadIdx.x + n - 1] = 0 ; 
		}


		for (int d = 1; d < n; d *= 2){
				offset /= 2 ;
				__syncthreads() ;
				if ( threadIdx.x < d){
					int ai = offset * (2* threadIdx.x +1) -1;
					int bi = offset * (2* threadIdx.x +2) -1;
					float t = As[ai ];
					As[ai] = As[bi ];
					As[bi] += t;
				}
		}

		__syncthreads () ;

		if (threadIdx.x >=1){
				Cs[threadIdx.x - 1] = As[threadIdx.x] ;
		}
		__syncthreads () ;
		cd[x] = Cs[threadIdx.x];			
}












void scan20(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 1024 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(1024);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<1,1024>>>(cqd,cqd,ad);
	add<<<1024,1024>>>(cd,cqd);

	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);


}




void scan21(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 2048* sizeof(float));
	cudaMalloc((void**)&cqd1, 2 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(2048);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<2,1024>>>(cqd,cqd,cqd1);
	scan_s<<<1,1>>>(cqd1,cqd1,2);
	add<<<2,1024>>>(cqd,cqd1);
	add<<<blocks,1024>>>(cd,cqd);


	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}


void scan22(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 4096* sizeof(float));
	cudaMalloc((void**)&cqd1, 4 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(4096);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<4,1024>>>(cqd,cqd,cqd1);
	scan_s<<<1,1>>>(cqd1,cqd1,4);
	add<<<4,1024>>>(cqd,cqd1);
	add<<<blocks,1024>>>(cd,cqd);


	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}





void scan23(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 8192 * sizeof(float));
	cudaMalloc((void**)&cqd1, 8 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(8192);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<8,1024>>>(cqd,cqd,cqd1);
	scan_s<<<1,1>>>(cqd1,cqd1,8);
	add<<<8,1024>>>(cqd,cqd1);
	add<<<blocks,1024>>>(cd,cqd);


	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}






void scan24(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 16384 * sizeof(float));
	cudaMalloc((void**)&cqd1, 16 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(16384);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<16,1024>>>(cqd,cqd,cqd1);
	scan_s<<<1,1>>>(cqd1,cqd1,16);
	add<<<16,1024>>>(cqd,cqd1);
	add<<<blocks,1024>>>(cd,cqd);


	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}





void scan25(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 32768 * sizeof(float));
	cudaMalloc((void**)&cqd1, 32 * sizeof(float));
  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);


	dim3 blocks(32768);
	dim3 threadsPerBlock(1024);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<32,1024>>>(cqd,cqd,cqd1);
	scan32<<<1,32>>>(cqd1,cqd1);
	add<<<32,1024>>>(cqd,cqd1);
	add<<<blocks,1024>>>(cd,cqd);


	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}






void scan26(float* a, float* c,int n) {

	float* ad;
	float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 65536 * sizeof(float));
	cudaMalloc((void**)&cqd1, 64 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);
  


	dim3 blocks(32768,2);
	dim3 threadsPerBlock(1024,1);
	scanl<<<blocks,threadsPerBlock>>>(ad,cd,cqd);
	scanl<<<64,1024>>>(cqd,cqd,cqd1);
	scan64<<<1,64>>>(cqd1,cqd1);
	add<<<64,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(cd,cqd);

	cudaMemcpy(c, cd, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}





void scan27(float* a, float* c,int n) {

	float* ad;
	//float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, n * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice);
  


	dim3 blocks(32768,4);
	dim3 threadsPerBlock(1024,1);
	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);

	cudaMemcpy(c, ad, n * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(ad);
  	//cudaFree(cd);
	cudaFree(cqd);
	cudaFree(cqd1);



}




void scan28(float* a, float* c,int n) {



	float* ad;
	//float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, (n/2) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a, (n/2) * sizeof(float), cudaMemcpyHostToDevice);
  


	dim3 blocks(32768,4);
	dim3 threadsPerBlock(1024,1);
	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);

	cudaMemcpy(c, ad, (n/2) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);





//-----------------------------------part 2 ------------------------------------


	cudaMalloc((void**)&ad, (n/2) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a+n/2, (n/2) * sizeof(float), cudaMemcpyHostToDevice);
  



	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);
	add_l<<<blocks,threadsPerBlock>>>(ad,c[n/2-1]);

	cudaMemcpy(c+n/2, ad, (n/2) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);




}





void scan29(float* a, float* c,int n) {

	float* ad;
	//float* cd;
	float* cqd;
	float* cqd1;

  	cudaMalloc((void**)&ad, (n/4) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a, (n/4) * sizeof(float), cudaMemcpyHostToDevice);
  


	dim3 blocks(32768,4);
	dim3 threadsPerBlock(1024,1);
	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);

	cudaMemcpy(c, ad, (n/4) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);





//-----------------------------------part 2 ------------------------------------


	cudaMalloc((void**)&ad, (n/4) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a+n/4, (n/4) * sizeof(float), cudaMemcpyHostToDevice);
  



	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);
	add_l<<<blocks,threadsPerBlock>>>(ad,c[n/4-1]);

	cudaMemcpy(c+n/4, ad, (n/4) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);








//-----------------------------------part 3 ------------------------------------


	cudaMalloc((void**)&ad, (n/4) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a+n/2, (n/4) * sizeof(float), cudaMemcpyHostToDevice);
  



	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);
	add_l<<<blocks,threadsPerBlock>>>(ad,c[n/2-1]);

	cudaMemcpy(c+n/2, ad, (n/4) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);







//-----------------------------------part 4 ------------------------------------


	cudaMalloc((void**)&ad, (n/4) * sizeof(float));
  	//cudaMalloc((void**)&cd, n * sizeof(float));
	cudaMalloc((void**)&cqd, 131072 * sizeof(float));
	cudaMalloc((void**)&cqd1, 128 * sizeof(float));

  	cudaMemcpy(ad, a+(n/4)*3, (n/4) * sizeof(float), cudaMemcpyHostToDevice);
  



	scanl<<<blocks,threadsPerBlock>>>(ad,ad,cqd);
	scanl<<<128,1024>>>(cqd,cqd,cqd1);
	scan128<<<1,128>>>(cqd1,cqd1);
	add<<<128,1024>>>(cqd,cqd1);
	add<<<blocks,threadsPerBlock>>>(ad,cqd);
	add_l<<<blocks,threadsPerBlock>>>(ad,c[(n/4)*3-1]);

	cudaMemcpy(c+(n/4)*3, ad, (n/4) * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(ad);
	cudaFree(cqd);
	cudaFree(cqd1);



}










void gpuKernel(float* a, float* c,int n) {

 switch (n) {
  case 1048576:
    scan20(a,c,n);
    break;

  case 2097152:
    scan21(a,c,n);
    break;

  case 4194304:
    scan22(a,c,n);
    break;

  case 8388608:
    scan23(a,c,n);
    break;

  case 16777216:
    scan24(a,c,n);
    break;

  case 33554432:
    scan25(a,c,n);
   break;

  case 67108864:
    scan26(a,c,n);
    break;

  case 134217728:
    scan27(a,c,n);
    break;


  case 268435456:
    scan28(a,c,n);
    break;


  case 536870912:
    scan29(a,c,n);
    break;
  }
 	

}