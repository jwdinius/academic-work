/*
  Copyright (C) 2009-2012 Fraunhofer SCAI, Schloss Birlinghoven, 53754 Sankt Augustin, Germany;
  all rights reserved unless otherwise stated.
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  MA 02111-1307 USA
*/

#include <cstdio>
#include <sys/time.h>

#include "cublas.h"

#include "QR.hpp"
#include "GPUTimer.hpp"

/** Index function to address the two-dimensional arrays
    Q and R

    Matrices are stored in column-major order (like Fortran).

    i is the row, j is the column (index starts at 1)
    ld is the number of elements for each column
*/

#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))

/* ---------------------------------------------------------------------- */

/*  Tuning can be done by different block sizes. */

#define BLOCK1 64

// 8800 GT:    128 x 1
// C1060:      128 x 1
#define BLOCK1X    64
#define BLOCK1Y    8

// 8800 GT:      64 x 4
// C1060:        64 x 8
#define BLOCK2X   512
#define BLOCK2Y   1
/* ---------------------------------------------------------------------- */

/** Kernel for matrix-vector multiplication

    R(k,k:n) = matmulv( Q(1:m,k:n), Q(1:m) )

    Same as this BLAS-2 call:

    call sgemv('T', m, n-k+1, 1.0, Q(1,k), M, Q(1,k), 1, 0.0, R(k,k), N)

    The threads in x-dimension are used for parallelization of
    the dot products, the threads in y-dimension compute different
    elements of the result vector.

    Each thread (t1,t2)  will be responsible for BLOCK1X columns and BLOCK1Y
    rows of the matrix Q.
*/

__global__ void mult(float* Q, float* R, int m, int n, int k)
{
  __shared__ float RS[BLOCK1Y][BLOCK1X];
  __shared__ float QK[BLOCK1Y];

  int tid1 = threadIdx.x;
  int tid2 = threadIdx.y;

  int i = blockIdx.x * BLOCK1Y + tid2 + k;

  float S = 0.0f;

  if (i < k or i > n) return;

  for (int j = tid1+1; j <= m; j+=BLOCK1X) {
    if (tid1 == 0) QK[tid2] = Q[IDX2F(j,k,m)];
    __syncthreads();
    S += QK[tid2] * Q[IDX2F(j,i,m)];
  }

  // thread writes result in shared array RS

  RS[tid2][tid1] = S;

  int NT = BLOCK1X;

  while (NT > 1) {
    // first half of threads sums up
    __syncthreads();
    NT = NT >> 1 ;
    if (tid1 < NT) {
      RS[tid2][tid1] += RS[tid2][tid1+NT];
    }
  }

  // now thread 0 writes the result

  if (tid1 == 0) {
    R[IDX2F(k,i,n)] = RS[tid2][0];
  }
}

/* ---------------------------------------------------------------------- */

/** This kernel scales the row k of the matrix R

    R(k,k:n) = R(k,k:n) * S
*/

__global__ void scaleR(float* Q, float* R, int m, int n, int k, float S)
{
  int i = blockIdx.x * BLOCK1 + threadIdx.x + k;

  if (i >= k and i <= n) {
     R[IDX2F(k,i,n)] *= S;
  }
}

/* ---------------------------------------------------------------------- */

/** This kernel scales the column k of the matrix Q.

    Q(1:m,k) = Q(1:m,k) * S
*/

__global__ void scaleQ(float* Q, float* R, int m, int n, int k, float S)
{
  int i = blockIdx.x * BLOCK1 + threadIdx.x + 1;

  if (i <= m) {
     Q[IDX2F(i,k,m)] *= S;
  }
}

/* ---------------------------------------------------------------------- */

/** This kernel updates the matrix Q by a product of two vectors.

    Q(1:m,k+1:n) -= R(k,k+1:n) * Q(1:m,k)

    same as this BLAS-2 call:

    call sger(M, N-K, -1.0, Q(1,K), 1, R(K,K+1), N, Q(1,K+1), M)

    Each thread (t1,t2)  will be responsible for BLOCK2X columns and BLOCK2Y
    rows of the matrix Q.
*/

__global__ void update(float* Q, float* R, int m, int n, int k)
{
  __shared__ float RK[BLOCK2Y];
  __shared__ float QK[BLOCK2X];

  int tid1 = threadIdx.x;
  int tid2 = threadIdx.y;

  int j = blockIdx.y * BLOCK2Y + tid2 + k + 1;

  if (j < k+1 or j > n) return;

  if (tid1 == 0) {
    RK[tid2] = R[IDX2F(k,j,n)];
  }

  for (int i = tid1 + 1; i <= m; i += BLOCK2X ) {

    if (tid2 == 0) {
       QK[tid1] = Q[IDX2F(i,k,m)];
    }

    __syncthreads();

    Q[IDX2F(i,j,m)] -= QK[tid1] * RK[tid2];
  }
}

/* ---------------------------------------------------------------------- */

/**  QR factorization of a matrix

     @param[in]      m is number of rows for Q and R
     @param[in]      n is number of columns for Q and R
     @param[in,out]  Q is a matrix of size m x n, column major order
     @param[out]     R is a matrix of size m x n, column major order

     @returns 0 if successful

     Q(in) = Q(out) * R, where Q(out) is orthonormal and R upper-triangular
*/

int QR(float* Q, float* R, unsigned int m, unsigned int n) 
{
  float* QGPU;   // Q on GPU
  float* RGPU;   // R on GPU
  
  GPUTimer timeMult;   // timer for matrix-vector multiplication
  GPUTimer timeScale;  // timer for the both scaling operations
  GPUTimer timeUpdate; // timer for matrix update with vector * transposed vector

#if defined(TIMING) and defined(LIST)
  FILE* f;
  char filename[256];
#endif

#ifdef DEBUG
  printf("allocate QGPU (%d, %d)\n", m, n);
#endif

  CUDA_SAFE_CALL(cudaMalloc((void**) &QGPU, m * n * sizeof(float)));

#ifdef DEBUG
  printf("allocate RGPU (%d, %d)\n", n, n);
#endif

  CUDA_SAFE_CALL(cudaMalloc((void**) &RGPU, n * n * sizeof(float)));

#ifdef DEBUG
  printf("transfer arrays to GPU\n");
#endif

  CUDA_SAFE_CALL(cudaMemcpy(QGPU, Q, m * n * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(RGPU, R, n * n * sizeof(float), cudaMemcpyHostToDevice));

  dim3 dimGrid;
  dim3 dimBlock;
 
#ifdef DEBUG
  printf("Run main loop\n");
#endif

#if defined(TIMING) and defined(LIST)
  sprintf(filename, "list.gpu.%d_%d", m, n);
  f = fopen(filename, "w");
#endif

  for (unsigned int k = 1; k <= n; k++) {

#ifdef TIMING
    timeMult.start();
#endif 

    dimGrid  = dim3((n - k + BLOCK1Y) / BLOCK1Y, 1, 1);
    dimBlock = dim3(BLOCK1X, BLOCK1Y, 1);

    // CUDA_SAFE_CALL(cudaGetLastError());

    mult<<<dimGrid,dimBlock>>>(QGPU, RGPU, m, n, k);

    // CUDA_SAFE_CALL(cudaThreadSynchronize());

#ifdef TIMING
    timeMult.stop();
#endif

    float S;

    CUDA_SAFE_CALL(cudaMemcpy(&S, &RGPU[IDX2F(k,k,n)], sizeof(float), cudaMemcpyDeviceToHost));

    S = sqrt(S);

    S = 1.0 / S;

#ifdef TIMING
    timeScale.start();
#endif 

    dimGrid  = dim3((m + BLOCK1 -1)/ BLOCK1, 1, 1);
    dimBlock = dim3(BLOCK1, 1, 1);
    scaleQ<<<dimGrid,dimBlock>>>(QGPU, RGPU, m, n, k, S);

    dimGrid  = dim3((n - k + BLOCK1)/ BLOCK1, 1, 1);
    dimBlock = dim3(BLOCK1, 1, 1);
    scaleR<<<dimGrid,dimBlock>>>(QGPU, RGPU, m, n, k, S);

#ifdef TIMING
    timeScale.stop();
    timeUpdate.start();
#endif

    dimGrid = dim3(1, (n - k + BLOCK2Y)/ BLOCK2Y, 1);
    dimBlock = dim3(BLOCK2X, BLOCK2Y, 1);
    update<<<dimGrid,dimBlock>>>(QGPU, RGPU, m, n, k);

#ifdef TIMING
    timeUpdate.stop();
#endif
  
#if defined(TIMING) and defined(LIST)
   fprintf(f, "%d %f %f %f\n", n + 1 - k, timeMult.last(), timeScale.last(), timeUpdate.last());
#endif

  }
  
#if defined(DEBUG) or defined(TIMING)
   printf("Detailed timings: mult = %g, scale = %g, update = %g\n",
           timeMult.get(), timeScale.get(), timeUpdate.get());
#endif

#if defined(TIMING) and defined(LIST)
   fclose(f);
#endif

  // transfer data back from GPU to CPU

  CUDA_SAFE_CALL(cudaMemcpy(Q, QGPU, m * n * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(R, RGPU, n * n * sizeof(float), cudaMemcpyDeviceToHost));

  // free device memory

  CUDA_SAFE_CALL(cudaFree(QGPU));
  CUDA_SAFE_CALL(cudaFree(RGPU));

  return EXIT_SUCCESS; 
}

/* ---------------------------------------------------------------------- */

__global__ void dummy_kernel()
{
}

/* ---------------------------------------------------------------------- */

/** Function to initialize the GPU. Takes the commandline arguments to allow to
    choose the GPU.
*/

void initGPU(int argc, char **argv)
{
  // call a dummy kernel that initializes GPU

  int deviceCount;

  CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));

#ifdef DEBUG
  printf ("initGPU: %d devices available\n", deviceCount);
#endif 

  int device = 0;   // default device

#ifdef DEBUG
  printf ( "Number of arguments = %d\n", argc);
#endif

  if (argc > 1) {
     sscanf(argv[1], "%d", &device);
  }

#ifdef DEBUG
  printf ("initGPU: try to use device %d\n", device);
#endif

  CUDA_SAFE_CALL(cudaSetDevice(device));

  cudaDeviceProp deviceProp;

  CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, device));

  if (deviceProp.major < 1) {
     fprintf(stderr, "cutil error: device %d does not support CUDA.\n", device);
     exit(-1);
  }

  fprintf(stderr, "Using device %d: %s\n", device, deviceProp.name);

  dim3 dimGrid(16, 16, 1);
  dim3 dimBlock(16, 16, 1);

  dummy_kernel<<<dimGrid, dimBlock>>>();

#ifdef DEBUG
  printf ("initGPU done, did already run a small dummy kernel\n");
#endif
}

/* ---------------------------------------------------------------------- */

void freeGPU()
{
  // nothing to be done 
}
