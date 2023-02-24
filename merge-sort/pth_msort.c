#include "pth_msort.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>



//========threads data========
typedef struct {
    unsigned int size ; 
    int* input;
    int* output;
    int low ; 
    int mid ; 
    int high ; 
    }inputs1;

//========threads data for last merge ============
typedef struct { 
    int* input;
    int x1 ; 
    int x2 ; 
    int y1 ;
    int y2 ; 
    int* output;
  }inputsp3;




 //==============binarySearch start=====================
 int binarySearch(int arr[], int l, int r, int x)
{
	while (l <= r) {
		int m = l + (r - l) / 2;

		// Check if x is present at mid
		if (arr[m] <= x && arr[m+1] >= x  ) 
			return m;

		// If x greater, ignore left half
		if (arr[m] < x)
			l = m + 1;

		// If x is smaller, ignore right half
		else
			r = m - 1;
	}

	// if we reach here, then element was
	// not present
	return -1;
}
 //============== end binarySearch =====================
















//===============merge code===================================

void merge(int *out, int *in, int low, int mid, int high)
{
	int ai = low; 
  int bi = mid + 1;
  int i = low;

	while (ai <= mid && bi <= high){
		if (out[ai] <= out[bi]){
			in[i] = out[ai];
			ai++;
		}else{
			in[i] = out[bi];
			bi++;
		} i++;
	}
	
	if (ai > mid){
		memcpy(&in[i], &out[bi], (high-bi+1)*sizeof(int));
	}else{
		memcpy(&in[i], &out[ai], (mid-ai+1)*sizeof(int));
	}
	memcpy(&out[low], &in[low], (high-low+1)*sizeof(int));
}




void mergesortHelper(int *out, int *in, int low, int high){
	if (low == high) return;
	int mid = (low + high)/2;
	
	mergesortHelper(out, in, low, mid);
	mergesortHelper(out, in, mid+1, high);
	merge(out, in, low, mid, high);
}


//=====================function for last merge=======================
void merge3(int *arr, int l1, int l2 , int l3 , int l4 , int* c_out){
    int* c = c_out;
    int ptr1 = l1 ;
    int ptr2 = l3 ;
    int k = 0 ;
    while (ptr1 < l2 && ptr2 < l4) {
        if (arr[ptr1] < arr[ptr2]) {
            c[k] = arr[ptr1];
            k++;
            ptr1++;
        }
        else  {
            c[k] = arr[ptr2];
            k++;
            ptr2++;
        }
    }
    while (ptr1 < l2) {
        c[k] = arr[ptr1];
        k++;
        ptr1++;
    }
    while (ptr2 < l4) {
        c[k] = arr[ptr2];
        k++;
        ptr2++;
    }
}


















//==================threads function===================
void* mergesort_t1(void* input) {
   inputs1* s ; 
   s = (inputs1*) input;
   mergesortHelper( s->input, s->output , s->low , s->high);
   return NULL;
}


void* mergesort_t21(void* input) {
   inputs1* s ; 
   s = (inputs1*) input;
   merge( s->input , s->output , s->low , s->mid , s->high );
   return NULL;
}


void* mergesort_t31(void* input) {
   inputsp3* s ; 
   s = (inputsp3*) input;
   merge3( s->input, s->x1, s->x2, s->y1 ,s->y2, s->output);
   return NULL;
}





//==========main f==============================
void mergeSortParallel (const int* values, unsigned int N, int* sorted) {
  int* arr = (int*) values;
  int i ;
  int j ;
  int temp ;
  //creart pthread_t handles
  pthread_t handle1, handle2, handle3, handle4;

  //==========start step1===========================

  //------thread1-----------
  inputs1 input1 = {N , arr , sorted ,0,0,N/4-1};
  pthread_create( &handle1, NULL, mergesort_t1, (void*) (&input1) );

  //------thread2-----------
  inputs1 input2 = {N , arr , sorted ,N/4,0,N/2-1};
  pthread_create( &handle2, NULL, mergesort_t1, (void*) (&input2) );

  //------thread3-----------
  inputs1 input3 = {N , arr , sorted ,N/2,0,3* (N/4) -1};
  pthread_create( &handle3, NULL, mergesort_t1, (void*) (&input3) );

  //------thread4-----------
  inputs1 input4 = {N , arr , sorted ,3* (N/4) , 0 ,N-1};
  pthread_create( &handle4, NULL, mergesort_t1, (void*) (&input4) );

  pthread_join(handle1, NULL);
  pthread_join(handle2, NULL);
  pthread_join(handle3, NULL);
  pthread_join(handle4, NULL);
 
  //===================end step 1 ====================================

  //=========serial part for step1=====================================
  /*
  mergesortHelper(input1.input, input1.output, input1.low, input1.high);
  mergesortHelper(input2.input, input2.output, input2.low, input2.high);
  mergesortHelper(input3.input, input3.output, input3.low, input3.high);
  mergesortHelper(input4.input, input4.output, input4.low, input4.high);
  */
  //=========end serial part for step1=====================================



  //======================start step2=======================================
  //==============thread 1==================================================

  inputs1 input21 = {N , arr , sorted ,0 , N/4-1 ,N/2-1};
  pthread_create( &handle1, NULL, mergesort_t21, (void*) (&input21) );

  //==============thread 2==================================================
  inputs1 input22 = {N , arr , sorted ,N/2 ,3*(N/4)-1 ,N-1};
  pthread_create( &handle2, NULL, mergesort_t21, (void*) (&input22) );

  pthread_join(handle1, NULL);
  pthread_join(handle2, NULL);

  //------------------serial code for step2----------------------------------
  //merge(input21.input, input21.output,input21.low,input21.mid ,input21.high);
  //merge(input22.input, input22.output,input22.low,input22.mid ,input22.high);
  //==============end step2====================================================


 //=================start step 3 =================

    int* a1 = (int*) malloc ( sizeof(int) * 5 );
    int* a2 = (int*) malloc ( sizeof(int) * 5 );
    int result;


    //=============start search==========================
    result = binarySearch(arr + N/2 , 0, N/2, arr[N/8]);
    a2[0] = result + N/2;
    if (result==-1 && arr[N/8] < arr[N/2]){
        a2[0] = N/2 ;
    }
    if (result==-1 && arr[N/8] > arr[N]){
        a2[0] = N-1 ;
    }



    result = binarySearch(arr + N/2 , 0, N/2, arr[3*N/8]);
    a2[1] =  result + N/2;
    if (result==-1 && arr[3*N/8] < arr[N/2]){
        a2[1] = N/2 ;
    }
    if (result==-1 && arr[3*N/8] > arr[N]){
        a2[1] = N-1 ;
    }



    result = binarySearch(arr  , 0, N/2, arr[3*N/4]);
    a1[1] = result ; 
    if (result==-1 && arr[3*N/4] < arr[0]){
        a1[1] = 0 ;
    }
    if (result==-1 && arr[3*N/4] > arr[N/2-1]){
        a1[1] = N/2-1 ;
    }
    //=============end search==========================

    
    //==================sort subarray=================
    a1[0] = 0;
    a1[2] = N/8-1;
    a1[3] = 3*N/8-1;
    a1[4] = N/2 - 1 ;
    
    a2[2] = N/2;
    a2[3] = 3*N/4-1;
    a2[4] = N -1 ;

    for(i=0;i<5;i++){
        for(j=0;j<4;j++){
            if (a1[j] > a1 [j+1]){
                temp = a1[j];
                a1[j] = a1[j+1];
                a1[j+1] = temp; 
            }

        }
    }



    for(i=0;i<5;i++){
        for(j=0;j<4;j++){
            if (a2[j] > a2[j+1]){
                temp = a2[j];
                a2[j] = a2[j+1];
                a2[j+1] = temp; 
            }

        }
    }


    //=======================finish sort=================================
 
    //========thread #1========================
    inputsp3 input31 = {arr , a1[0] , a1[1]+1 ,a2[0],a2[1]+1,sorted}; 
    //merge3(arr , a1[0],   a1[1]+1  , a2[0]    , a2[1]+1  ,  sorted );   
    pthread_create( &handle1, NULL, mergesort_t31, (void*) (&input31) );


    //========thread #2========================
    inputsp3 input32 = {arr ,a1[1]+1 ,  a1[2]+1 ,a2[1]+1,a2[2]+1,sorted + a1[1] - a1[0] + a2[1] - a2[0] + 2}; 
    //merge3(arr , a1[1]+1, a1[2]+1  , a2[1]+1  , a2[2]+1   , sorted + a1[1] - a1[0] + a2[1] - a2[0] + 2   );  
    pthread_create( &handle2, NULL, mergesort_t31, (void*) (&input32) );



    //========thread #3========================
    inputsp3 input33 = {arr ,a1[2]+1 ,   a1[3]+1 ,a2[2]+1,a2[3]+1,sorted + a1[2] - a1[1] + a2[2] - a2[1] + a1[1] - a1[0] + a2[1] - a2[0] +2}; 
    //merge3(arr , a1[2]+1, a1[3]+1  , a2[2]+1  , a2[3]+1   , sorted + a1[2] - a1[1] + a2[2] - a2[1] + a1[1] - a1[0] + a2[1] - a2[0] +2  ); 
    pthread_create( &handle3, NULL, mergesort_t31, (void*) (&input33) );



    //========thread #4========================
    inputsp3 input34 = {arr ,a1[3]+1 , a1[4]+1 ,a2[3]+1, a2[4]+1,sorted + a1[3] - a1[2] + a2[3] - a2[2] + a1[2] - a1[1] + a2[2] - a2[1] + a1[1] - a1[0] + a2[1] - a2[0] +2 }; 
    // merge3(arr , a1[3]+1, a1[4]+1  , a2[3]+1  , a2[4]+1   , sorted + a1[3] - a1[2] + a2[3] - a2[2] + a1[2] - a1[1] + a2[2] - a2[1] + a1[1] - a1[0] + a2[1] - a2[0] +2  );
    pthread_create( &handle4, NULL, mergesort_t31, (void*) (&input34) );

    pthread_join(handle1, NULL);
    pthread_join(handle2, NULL);
    pthread_join(handle3, NULL);
    pthread_join(handle4, NULL);







}