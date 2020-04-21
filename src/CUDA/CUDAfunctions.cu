/*
 * CUDAfunctions.cu
 *
 *  Created on: Nov 27, 2014
 *      Author: marcel
 */

#include "CUDAfunctions.hpp"
#include <stdio.h>
#include <cuda.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <cmath>



__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
				__double_as_longlong(val + __longlong_as_double(assumed)));
		// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
		} while (assumed != old);

	return __longlong_as_double(old);
}

__global__ void modulus(cuDoubleComplex* initialpptr, long int dimx,long int dimy,long int dimz, double* value){
	long int xmax=dimx-1;
	long int ymax=dimy-1;
	long int zmax=dimz-1;
	double partial=0.0;

	long int l = threadIdx.x;
	long int m = blockIdx.x;

        char* pointer = (char*)initialpptr; 
	size_t pitch = dimx*sizeof(cuDoubleComplex);
	size_t slicePitch = pitch * dimy;
	char* slice = pointer + m * slicePitch;
	cuDoubleComplex* initial=(cuDoubleComplex*)(slice + l * pitch);


	if(l==0 || l ==ymax){
		partial+=cuCreal(cuCmul(initial[0],cuConj(initial[0])))/4.0;
		for(long int k = 1; k < xmax; ++k){
			partial+=cuCreal(cuCmul(initial[k],cuConj(initial[k])))/2.0;
		}
		partial+=cuCreal(cuCmul(initial[xmax],cuConj(initial[xmax])))/4.0;
		if(m==0 || m ==zmax){
			partial+=cuCreal(cuCmul(initial[0],cuConj(initial[0])))/8.0;
			for(long int k = 1; k < xmax; ++k){
				partial+=cuCreal(cuCmul(initial[k],cuConj(initial[k])))/4.0;
			}
			partial+=cuCreal(cuCmul(initial[xmax],cuConj(initial[xmax])))/8.0;
		}
	}
	else if(m==0 || m ==zmax){
		partial+=cuCreal(cuCmul(initial[0],cuConj(initial[0])))/4.0;
		for(long int k = 1; k < xmax; ++k){
			partial+=cuCreal(cuCmul(initial[k],cuConj(initial[k])))/2.0;
		}
		partial+=cuCreal(cuCmul(initial[xmax],cuConj(initial[xmax])))/4.0;
	}

	partial+=cuCreal(cuCmul(initial[0],cuConj(initial[0])))/2.0;
	for(long int k = 1; k < xmax; ++k){
		partial+=cuCreal(cuCmul(initial[k],cuConj(initial[k])));
	}
	partial+=cuCreal(cuCmul(initial[xmax],cuConj(initial[xmax])))/2.0;


	atomicAdd(value,partial);
}


__global__ void dividebyscalar3D(cuDoubleComplex* initialpptr, long int dimx,long int dimy,long int dimz, double value){
	long int xmax=dimx-1;

	long int l = threadIdx.x;
	long int m = blockIdx.x;

	char* pointer = (char*)initialpptr;
	size_t pitch = dimx*sizeof(cuDoubleComplex);
	size_t slicePitch = pitch * dimy;
	char* slice = pointer + m * slicePitch;
	cuDoubleComplex* initial=(cuDoubleComplex*)(slice + l * pitch);
	cuDoubleComplex complexvalue = make_cuDoubleComplex(value,0.0);

	for(long int k = 0; k <= xmax; ++k){
		initial[k]=cuCdiv(initial[k],complexvalue);
	}

}

__global__ void matrixmulti3D(cuDoubleComplex* initialpptr, cuDoubleComplex* potentialpptr, long int dimx,long int dimy,long int dimz){
	long int xmax=dimx-1;


	long int l = threadIdx.x;
	long int m = blockIdx.x;

	char* pointer = (char*)initialpptr;
	size_t pitch = dimx*sizeof(cuDoubleComplex);
	size_t slicePitch = pitch * dimy;
	char* slice = pointer + m * slicePitch;
	cuDoubleComplex* initial=(cuDoubleComplex*)(slice + l * pitch);

	char* pointertwo = (char*)potentialpptr;
	size_t pitchtwo = dimx*sizeof(cuDoubleComplex);
	size_t slicePitchtwo = pitchtwo * dimy;
	char* slicetwo = pointertwo + m * slicePitchtwo;
	cuDoubleComplex* potential=(cuDoubleComplex*)(slicetwo + l * pitchtwo);

	for(long int k = 0; k <= xmax; ++k){
		initial[k]=cuCmul(initial[k],potential[k]);
	}

}

__host__ void getmodulus3D(cuDoubleComplex* initialpptr, long int dimx,long int dimy,long int dimz, double xstep, double ystep, double zstep, double* modulusout){
	double* modulusinside;
        cudaMalloc(&modulusinside,sizeof(double));
	double zero = 0.0;
	cudaMemcpy(modulusinside, &zero, sizeof(double),cudaMemcpyHostToDevice);
	modulus<<<dimz,dimy>>>(initialpptr,dimx,dimy,dimz,modulusinside);
	double value;
	cudaMemcpy(&value, modulusinside, sizeof(double),cudaMemcpyDeviceToHost);
        value=sqrt(value*(dimx-1)*(dimy-1)*(dimy-1)*xstep*ystep*zstep);
        *modulusout=value;
        cudaFree(modulusinside);
}

 void groundlevel_periodicFFT_CUDA(long int timesteps,std::complex<double>* initial,std::complex<double>* potential,std::complex<double>* fftmulti,
		double effectivemass, long int dimx,long int dimy,long int dimz, double xstep, double ystep, double zstep){
	double modulusvalue;

	size_t sizeComplex = dimx*dimy*dimz*sizeof(cuDoubleComplex);

	cuDoubleComplex* initialPtr;
	cuDoubleComplex* potentialPtr;
	cuDoubleComplex* fftmultiPtr;

	cudaMalloc(&initialPtr, sizeComplex);
	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to allocate\n");
	return;
	}

	cudaMalloc(&potentialPtr, sizeComplex);
	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to allocate\n");
	return;
	}

	cudaMalloc(&fftmultiPtr, sizeComplex);
	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to allocate\n");
	return;
	}

        cudaMemcpy(initialPtr,initial,sizeComplex,cudaMemcpyHostToDevice);
	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to copy to device\n");
	return;
	}	
        cudaMemcpy(potentialPtr,potential,sizeComplex,cudaMemcpyHostToDevice);
        	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to copy to device\n");
	return;
	}
        cudaMemcpy(fftmultiPtr,fftmulti,sizeComplex,cudaMemcpyHostToDevice);
        	if (cudaGetLastError() != cudaSuccess){
	 fprintf(stderr, "Cuda error: Failed to copy to device\n");
	return;
	}


        
	getmodulus3D(initialPtr,dimx,dimy,dimz,xstep,ystep,zstep,&modulusvalue);  
	dividebyscalar3D<<<dimz,dimy>>>(initialPtr,dimx,dimy,dimz,modulusvalue);
                
	cufftHandle plan;
	int n[3] = {dimx,dimy,dimz};

        //int inembed[] = {dimx,dimy,initialPitchedPtr.pitch/sizeof(cufftComplex)}; // Input size with pitch
        //int onembed[] = {dimx,dimy,initialPitchedPtr.pitch/sizeof(cufftComplex)}; // Output size with pitch

	/* Create a 3D FFT plan. */
	if (cufftPlanMany(&plan, 3, n,
	 NULL, 1, 1, // *inembed, istride, idist
	 NULL, 1, 1, // *onembed, ostride, odist
	 CUFFT_Z2Z, 1) != CUFFT_SUCCESS){
	 fprintf(stderr, "CUFFT error: Plan creation failed");
	return;
	}

     
        
        
	for(long int i = 0 ; i< timesteps; i++){
		printf("step: %i...\n",i);
		matrixmulti3D<<<dimz,dimy>>>(initialPtr,potentialPtr,dimx,dimy,dimz);
		cufftExecZ2Z(plan,(cufftDoubleComplex*)initialPtr, (cufftDoubleComplex*)initialPtr, CUFFT_FORWARD);
		matrixmulti3D<<<dimz,dimy>>>(initialPtr,fftmultiPtr,dimx,dimy,dimz);
		cufftExecZ2Z(plan,(cufftDoubleComplex*)initialPtr , (cufftDoubleComplex*)initialPtr, CUFFT_INVERSE);
		matrixmulti3D<<<dimz,dimy>>>(initialPtr,potentialPtr,dimx,dimy,dimz);
		double tonorm;
                getmodulus3D(initialPtr,dimx,dimy,dimz,xstep,ystep,zstep,&tonorm);
		dividebyscalar3D<<<dimz,dimy>>>(initialPtr,dimx,dimy,dimz,tonorm);
	}
        
                
        cudaMemcpy(initial,initialPtr,sizeComplex,cudaMemcpyDeviceToHost);
	cudaFree(initialPtr);
	cudaFree(potentialPtr);
	cudaFree(fftmultiPtr);

}


