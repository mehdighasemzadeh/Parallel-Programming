# Parallel-Programming

This repository contains:
  * Implementing parallel matrix multiplication, parallel scan algorithm, and parallel Fast Fourier Transform (FFT) using cuda.
  * Parallelizing merge sort algorithm using pthread library.
 
## Parallel Matrix Multiplication
Run 

```
cd multiply
nvcc -O2 bmm_main.cu bmm.cu -o bmm
./bmm M
```
2^M determines matrix-size


## Parallel Scan Algorithm
Run 

```
cd scan
nvcc scan2.cu scan2_main.cu -o scan2
./scan2 M
```
2^M determines array-size 


## FFT
Run 

```
cd FFT
nvcc -O2 fft_main.cu fft.cu -o fft
./fft M
```
2^M determines array-size

## Merge-Sort
four thread are uesed for this Parallelizing

Run 

```
cd merge-sort
gcc -O2 pth_msort_test.c pth_msort.c -lpthread -lm
./a.out M
```
2^M determines array-size

## All codes are executable on google colab, and an [ipynb file](Parallel.ipynb) is provided for that.

