#include <cstdio>
#include <cstdlib>

#define REPEAT_TIMES 100000

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

__global__ void calculateCorrelationCoefficientMatrix(float *all_corr_coef_matrix, const float *all_data_matrix, const int subject_size, const int time_size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= REPEAT_TIMES)
    return;

  const float *data_matrix = all_data_matrix + idx * subject_size * time_size;
  float *corr_coef_matrix = all_corr_coef_matrix + idx * subject_size * subject_size;

  for (int i = 0; i < subject_size; i++) {
    const float *data_x = data_matrix + i * time_size;
    corr_coef_matrix[i * subject_size + i] = 1.0;

    for (int j = i + 1; j < subject_size; j++) {
      const float *data_y = data_matrix + j * time_size;
      const float coef = calCorrCoef(data_x, data_y, time_size);

      corr_coef_matrix[i * subject_size + j] = coef;
      corr_coef_matrix[j * subject_size + i] = coef;
    }
  }
}

int main(int argc, char **argv)
{
  int subject_size = 3, time_size = 5;

  // TODO padding
  float *h_data_matrix, *h_coef_matrix;

  h_data_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * subject_size * time_size);
  h_coef_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * subject_size * subject_size);

  for(int i = 0; i < REPEAT_TIMES; i++) {
    int start = i * subject_size * time_size;
    for(int j = 0; j < subject_size; j++) {
      for(int k = 0; k < time_size; k++) {
        h_data_matrix[start + j * time_size + k] = rand();
      }
    }
  }

  float *d_data_matrix, *d_coef_matrix;

  cudaMalloc(&d_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size);
  cudaMalloc(&d_coef_matrix, sizeof(float) * REPEAT_TIMES * subject_size * subject_size);

  cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size, cudaMemcpyHostToDevice);

  int blocksize = 32;
  int nblock = REPEAT_TIMES/blocksize + REPEAT_TIMES%blocksize==0?0:1;

  calculateCorrelationCoefficientMatrix<<<nblock, blocksize>>>(d_coef_matrix, d_data_matrix, subject_size, time_size);

  cudaDeviceSynchronize();

  cudaMemcpy(h_coef_matrix, d_coef_matrix, sizeof(float) * REPEAT_TIMES * subject_size * subject_size, cudaMemcpyDeviceToHost);

  for(int i = 0; i < subject_size; i++) {
    for(int j = 0; j < time_size; j++) {
      printf("%f ", h_data_matrix[i * time_size + j]);
    }
    printf("\n");
  }
  printf("--------------------\n");
  for (int i = 0; i < subject_size; i ++ ) {
    for (int j = 0; j < subject_size; j ++) {
      printf("%f ", h_coef_matrix[i * subject_size + j]);
    }
    printf("\n");
  }

  cudaFree(d_data_matrix);
  cudaFree(d_coef_matrix);
  free(h_data_matrix);
  free(h_coef_matrix);

  return 0;
}
