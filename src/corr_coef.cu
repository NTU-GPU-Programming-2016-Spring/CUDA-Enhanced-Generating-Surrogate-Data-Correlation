#include <cstdio>
#include <cstdlib>
#include <thrust/count.h>
#include <thrust/execution_policy.h>

#define REPEAT_TIMES 100000

struct is_less_than
{
  float num;

  is_less_than(float x) {
    num = x;
  }

  __host__ __device__
  bool operator() (float x)
  {
    if (x < num)
      return true;
    else
      return false;
  }
};

__device__ float calSum(const float *x, int length)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i];
  }
  return sum;
}

__device__ float calAvg(const float *x, int length)
{
  return calSum(x, length) / length;
}

__device__ float calSquareSum(const float *x, int length)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i] * x[i];
  }
  return sum;
}

__device__ float calMultiplySum(const float *x, const float *y, int length)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i] * y[i];
  }
  return sum;
}

__device__ float calStd(const float *x, int length)
{
  const float x_square_sum = calSquareSum(x, length);
  const float x_avg = calAvg(x, length);

  return sqrtf((x_square_sum - length * x_avg * x_avg) / (length - 1));
}

__device__ float calCorrCoef(const float *x, const float *y, const int length)
{
  const float xy_sum = calMultiplySum(x, y, length);
  const float x_avg = calAvg(x, length);
  const float y_avg = calAvg(y, length);
  const float x_std = calStd(x, length);
  const float y_std = calStd(y, length);

  return (xy_sum - length * x_avg * y_avg) / ((length - 1) * x_std * y_std);
}

__device__ float calFisherTransform(const float x, const int time_size)
{
  // z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
  return 0.5 * logf((1+x) / (1-x)) / rsqrt((float)time_size/2.34 - 3);
}

__global__ void calculateInterSubjectCorrelation(float *isc_array, float *all_corr_coef_matrix, const float *all_data_matrix, const int subject_size, const int time_size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= REPEAT_TIMES)
    return;

  const float *data_matrix = all_data_matrix + idx * subject_size * time_size;
  float *corr_coef_matrix = all_corr_coef_matrix + idx * subject_size * subject_size;
  float isc_sum = 0;

  for (int i = 0; i < subject_size; i++) {
    const float *data_x = data_matrix + i * time_size;
    corr_coef_matrix[i * subject_size + i] = 0.0;

    for (int j = i + 1; j < subject_size; j++) {
      const float *data_y = data_matrix + j * time_size;
      const float coef = calCorrCoef(data_x, data_y, time_size);

      corr_coef_matrix[i * subject_size + j] = coef;
      corr_coef_matrix[j * subject_size + i] = 0.0;

      isc_sum += calFisherTransform(coef, time_size);
    }
  }

  const int isc_pairs = subject_size * (subject_size - 1) / 2;
  isc_array[idx] = isc_sum / isc_pairs;
}

int main(int argc, char **argv)
{
  int subject_size = 8, time_size = 440;

  // TODO padding
  float *h_data_matrix, *h_coef_matrix, *h_isc_array;
  float *d_data_matrix, *d_coef_matrix, *d_isc_array;

  h_data_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * subject_size * time_size);
  h_coef_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * subject_size * subject_size);
  h_isc_array = (float *)malloc(sizeof(float) * REPEAT_TIMES);

  for(int i = 0; i < REPEAT_TIMES; i++) {
    int start = i * subject_size * time_size;
    for(int j = 0; j < subject_size; j++) {
      for(int k = 0; k < time_size; k++) {
        h_data_matrix[start + j * time_size + k] = rand();
      }
    }
  }

  cudaMalloc(&d_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size);
  cudaMalloc(&d_coef_matrix, sizeof(float) * REPEAT_TIMES * subject_size * subject_size);
  cudaMalloc(&d_isc_array, sizeof(float) * REPEAT_TIMES);

  cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size, cudaMemcpyHostToDevice);

  int blocksize = 32;
  int nblock = REPEAT_TIMES/blocksize + REPEAT_TIMES%blocksize==0?0:1;

  calculateInterSubjectCorrelation<<<nblock, blocksize>>>(d_isc_array, d_coef_matrix, d_data_matrix, subject_size, time_size);
  cudaDeviceSynchronize();

  cudaMemcpy(h_coef_matrix, d_coef_matrix, sizeof(float) * REPEAT_TIMES * subject_size * subject_size, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_isc_array, d_isc_array, sizeof(float) * REPEAT_TIMES, cudaMemcpyDeviceToHost);

  int result = thrust::count_if(thrust::host, h_isc_array, h_isc_array + REPEAT_TIMES, is_less_than(h_isc_array[0]));

  // {
  //   const int idx = rand() % REPEAT_TIMES;

  //   for(int i = 0; i < subject_size; i++) {
  //     for(int j = 0; j < time_size; j++) {
  //       printf("%f ", h_data_matrix[idx * subject_size * time_size + i * time_size + j]);
  //     }
  //     printf("\n");
  //   }
  //   printf("--------------------\n");


  //   for (int i = 0; i < subject_size; i ++ ) {
  //     for (int j = 0; j < subject_size; j ++) {
  //       printf("%f ", h_coef_matrix[idx * subject_size * subject_size + i * subject_size + j]);
  //     }
  //     printf("\n");
  //   }
  //   printf("--------------------\n");

  //   printf("%f\n", h_isc_array[idx]);

  //   printf("--------------------\n");
  //   printf("%d\n", result);
  //   for(int i = 0; i < REPEAT_TIMES; i++) {
  //     printf("%f ", h_isc_array[i]);
  //   }
  //   printf("\n");

  // }


  cudaFree(d_data_matrix);
  cudaFree(d_coef_matrix);
  cudaFree(d_isc_array);
  free(h_data_matrix);
  free(h_coef_matrix);
  free(h_isc_array);

  return 0;
}
