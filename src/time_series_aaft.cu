/* For Brian */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Timer.h"

#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>

#include <helper_cuda.h>
#include <cuComplex.h>
#include <cufft.h>
#include <curand.h>

#define SIGSIZE 440
#define SIGDIM 10000
#define NBLK 256
#define TIMESLOT 439
#define pi 3.14159

// 'aaft' is output.
void amplitudeAdjustedFourierTransform(double *aaft, const double *data, const int row, const int maxColumns) {
	
	return;
}

__device__ float angle_trans(const cuComplex& z){
	return atan2(cuCimagf(z), cuCrealf(z));
}

__global__ void fft_polar_angle(cufftComplex *data, float *angle, float *mag, int data_size){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if(idx >= data_size){
		return;
	}
	//abs of fft
	mag[idx] = cuCabsf(data[idx]);
	//angle of fft
	angle[idx] = angle_trans(data[idx]);
	return;
}

// do p(2:N)=[p1 -flipud(p1)];
__global__ void odd_surr_trans(float *angle, float *ran, int data_size, int sig_size, int half_sig_size){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if(idx >= data_size){
		return;
	}
	

	int data_col = idx/sig_size;
	int data_idx = idx%sig_size;
	// p(1) is not necessary for changing
	if(data_idx ==0){
		return;
	}
	
	int half_idx;
	//p(2: 2+half-1)
	if(data_idx <= half_sig_size){
		half_idx = (data_idx-1) + data_col*half_sig_size;
		angle[idx] = 2*pi*ran[half_idx];
			
	// -flipup(p1)	
	}else{
		int diff = data_idx - half_sig_size;
		int reverse_data_idx = half_sig_size- diff;
		half_idx = reverse_data_idx + data_col*half_sig_size;
		angle[idx] = -2*pi*ran[half_idx];
	}

	return;

}

__global__ void even_surr_trans(float *angle, float *mag, float *ran, int data_size, int sig_size, int half_sig_size){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if(idx >= data_size){
		return;
	}
	
	int data_col = idx/sig_size;
	int data_idx = idx%sig_size;
	// angle part
	// p(2:N)=[p1' p(h+1) -flipud(p1)'];
	int half_idx;
	// 0 nothing
	// 1->half_sig_size-1
	if(data_idx == 0 || data_idx == half_sig_size){
		angle[idx] = angle[idx];
	}
	else if(1<= data_idx < half_sig_size){
		half_idx = (data_idx-1) + data_col*half_sig_size;
		angle[idx] = 2*pi*ran[half_idx];
	}
	// half_sig_size
	// half_sig_size+1->data_size-1
	if(data_idx > half_sig_size){
		int diff = data_idx - half_sig_size+1;
		int reverse_data_idx = half_sig_size- diff;
		
		half_idx = reverse_data_idx + data_col*half_sig_size;
		angle[idx] = -2*pi*ran[half_idx];

		// magnitude part
		// m=[flipud(m(2:h))];
		diff = data_idx - (half_sig_size);
		reverse_data_idx = (half_sig_size) - diff;
		int mag_idx = reverse_data_idx + data_col*sig_size;
		mag[idx] = mag[mag_idx];
	}

	return;

}
// s(:,j)=m.*exp(i*p);
__global__ void i_mul_trans(cufftComplex *result, const float *mag, const float *angle, int data_size){
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= data_size){
		return;
	}
	float mag_val = mag[idx];
	float angle_val = angle[idx];
	result[idx].x = mag_val*cosf(angle_val);
	result[idx].y = mag_val*sinf(angle_val);

	return;
}

__global__ void get_real_trans(float *result, const cufftComplex *data_list, const int data_size, const int sig_size){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;

	if (idx >= data_size){
		return;
	}
	result[idx] = data_list[idx].x/sig_size;
	return;
}

__global__ void real2cufft_trans(cufftComplex *result, const float *input, const int data_size){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;

	if (idx >= data_size){
		return;
	}
	result[idx].x = input[idx];
	result[idx].y = 0;
	return;
}

void phaseran(float *result, const int data_num, const int time_size){
	int data_size = data_num*time_size;
	int mem_size = sizeof(cufftComplex)*data_size;
	
	cufftComplex *d_signal;
	checkCudaErrors(cudaMalloc((void **) &d_signal, mem_size));
	float *d_input;
	checkCudaErrors(cudaMalloc(&d_input, sizeof(float)*data_size));
	checkCudaErrors(cudaMemcpy(d_input, result, sizeof(float)*data_size, cudaMemcpyHostToDevice));
	
	real2cufft_trans<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_signal, d_input, data_size);
	
	//cufft
	cufftHandle plan_r, plan;
	
	if (cufftPlan1d(&plan_r, time_size, CUFFT_R2C, data_num) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}
	if (cufftPlan1d(&plan, time_size, CUFFT_C2C, data_num) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}

	//forward transform
	// printf("---Transform fft--- \n");
	cufftExecC2C(plan, d_signal, d_signal, CUFFT_FORWARD);
	checkCudaErrors(cudaFree(d_input));
	
	//do angle implement in matlab
	float *d_angle, *d_mag;
	checkCudaErrors(cudaMalloc(&d_angle, sizeof(float)*data_size));
	checkCudaErrors(cudaMalloc(&d_mag, sizeof(float)*data_size));

	fft_polar_angle<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_signal, d_angle, d_mag, data_size);

	checkCudaErrors(cudaFree(d_signal));

	// start parallel surrogate
	int half_col_size = time_size/2;
	int half_size = half_col_size*data_num;
	float *d_ran_series;
	
	curandGenerator_t gen;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, rand()%10000);
	
	if(time_size%2==0){
		//assign half minus 1
		int half_minus_one_size = (half_col_size-1)*data_num;
		checkCudaErrors(cudaMalloc(&d_ran_series, sizeof(float)*half_minus_one_size));
		curandGenerateUniform(gen, d_ran_series, half_minus_one_size);

		even_surr_trans<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_angle, d_mag, d_ran_series, data_size, time_size, half_col_size);
	
	}else{		
		
		//assign half 
		checkCudaErrors(cudaMalloc(&d_ran_series, sizeof(float)*half_size));
		//random generator
		
		curandGenerateUniform(gen, d_ran_series, half_size);
		
		// do column vector trans p(2:N)=[p1 -flipud(p1)];
		odd_surr_trans<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_angle, d_ran_series, data_size, time_size, half_col_size);
	}
	checkCudaErrors(cudaFree(d_ran_series));
	
	// multiply with m.*exp(i*p) = m*cos(p) + m*i*sin(p)
	cufftComplex *d_i_mul;
	checkCudaErrors(cudaMalloc((void **) &d_i_mul, mem_size));
	
	i_mul_trans<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_i_mul, d_mag, d_angle, data_size);

	// backward transform
	// printf("---Inverse fft transform --- \n");
	cufftExecC2C(plan, d_i_mul, d_i_mul, 
							   CUFFT_INVERSE);
	float *d_result;
	checkCudaErrors(cudaMalloc(&d_result, sizeof(float)*data_size));

	get_real_trans<<<(data_size+NBLK-1)/NBLK, NBLK>>>(d_result, d_i_mul, data_size, time_size);
	
	checkCudaErrors(cudaMemcpy(result, d_result, sizeof(float)*data_size, cudaMemcpyDeviceToHost));
	
	cufftDestroy(plan);
	
	checkCudaErrors(cudaFree(d_angle));
	checkCudaErrors(cudaFree(d_mag));
	checkCudaErrors(cudaFree(d_result));
	checkCudaErrors(cudaFree(d_i_mul));
	cudaDeviceReset();
	return;
}

int main(int argc, char **argv)
{	
	float *result = (float *)malloc(sizeof(float)*SIGSIZE*SIGDIM);	
	Timer phaseran_timer;
	phaseran_timer.Start();
	for(int i = 0; i <1 ; i++){

		for(int i = 0; i<SIGSIZE*SIGDIM;i++){
			result[i] = i;
		}
		
		phaseran(result, SIGDIM, SIGSIZE);
		
	}
	phaseran_timer.Pause();
	printf_timer(phaseran_timer);
	free(result);
	return 0;
}	