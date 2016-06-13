#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <thrust/count.h>
#include <thrust/execution_policy.h>

#define REPEAT_TIMES 10000

#define cudaCheckErrors(msg) \
  do { \
    cudaError_t __err = cudaGetLastError(); \
    if (__err != cudaSuccess) { \
      fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
          msg, cudaGetErrorString(__err), \
          __FILE__, __LINE__); \
      fprintf(stderr, "*** FAILED - ABORTING\n"); \
      exit(1); \
    } \
  } while (0)

struct is_bigger_than
{
  float num;

  is_less_than(float x) {
    num = x;
  }

  __host__ __device__
  bool operator() (float x)
  {
    if (x >= num)
      return true;
    else
      return false;
  }
};

__device__ float calSum(const float *x, const int length, const int step)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step];
  }
  return sum;
}

__device__ float calAvg(const float *x, const int length, const int step)
{
  return calSum(x, length, step) / length;
}

__device__ float calSquareSum(const float *x, const int length, const int step)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step] * x[i*step];
  }
  return sum;
}

__device__ float calMultiplySum(const float *x, const float *y, const int length, const int step)
{
  float sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step] * y[i*step];
  }
  return sum;
}

__device__ float calStd(const float *x, const int length, const int step)
{
  const float x_square_sum = calSquareSum(x, length, step);
  const float x_avg = calAvg(x, length, step);

  return sqrtf((x_square_sum - length * x_avg * x_avg) / (length - 1));
}

__device__ float calCorrCoef(const float *x, const float *y, const int length, const int step)
{
  const float xy_sum = calMultiplySum(x, y, length, step);
  const float x_avg = calAvg(x, length, step);
  const float y_avg = calAvg(y, length, step);
  const float x_std = calStd(x, length, step);
  const float y_std = calStd(y, length, step);

  return (xy_sum - length * x_avg * y_avg) / ((length - 1) * x_std * y_std);
}

__device__ float calFisherTransform(const float x, const int time_size)
{
  // z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
  return 0.5 * logf((1+x) / (1-x)) / rsqrt((float)time_size/2.34 - 3);
}

__global__ void calculateCorrelationCoefficientMatrixAndFisherTransform(float *all_corr_coef_matrix, const float *all_data_matrix, const int subject_size, const int time_size)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int matrix_works = subject_size * (subject_size - 1) / 2;
  if (idx >= REPEAT_TIMES * matrix_works)
    return;

  const int n_matrix = idx / matrix_works;
  int remain_works = idx % matrix_works;
  int x = 1, y = 1;

  for(int i = 0; i < subject_size; i ++) {
    const int row_works = subject_size - i - 1;
    if (remain_works < row_works) {
      x = i;
      y = i + 1 + remain_works;
      break;
    }

    remain_works -= row_works;
  }

  const float *data_matrix = all_data_matrix + n_matrix * time_size * subject_size;
  const float coef = calCorrCoef(data_matrix + x, data_matrix + y, time_size, subject_size);
  const float zvalue = calFisherTransform(coef, time_size);

  all_corr_coef_matrix[idx] = zvalue;
}


__global__ void calculateInterSubjectCorrelation(float *isc_array, const float *all_corr_coef_matrix, const int subject_size)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= REPEAT_TIMES)
    return;

  const int matrix_works = subject_size * (subject_size - 1) / 2;
  const float *corr_coef_matrix = all_corr_coef_matrix + idx * matrix_works;

  float sum = 0;

  for (int i = 0; i < matrix_works; i++)
    sum += corr_coef_matrix[i];

  isc_array[idx] = sum / matrix_works;
}

__global__ void rearrangeMatrixPosition(float *data_matrix, const float *source_matrix, const int subject_size, const int time_size)
{
  // 1st subject_size, 2nd REPEAT_TIMES, 3rd time_size
  // to
  // 1st REPEAT_TIMES, 2nd time_size, 3rd subject_size

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= REPEAT_TIMES * time_size * subject_size)
    return;

  const int subject_idx = idx / (REPEAT_TIMES * time_size);
  const int repeat_idx = (idx % (REPEAT_TIMES * time_size)) / time_size;
  const int time_idx = (idx % (REPEAT_TIMES * time_size)) % time_size;

  const int data_idx = repeat_idx * time_size * subject_size + time_idx * subject_size + subject_idx;
  
  data_matrix[data_idx] = source_matrix[idx];
}

void printMatrix(const float *data, const int first, const int second, const int third)
{
  float *tmp = (float *)malloc(sizeof(float) * first * second * third);
  cudaMemcpy(tmp, data, sizeof(float) * first * second * third, cudaMemcpyDeviceToHost);

  printf("%% 1st:%d 2nd:%d 3rd:%d\n", first, second, third);
  for (int i = 0; i < first; i++ ) {
    for (int j = 0; j < second; j++ ) {
      printf("%% ");
      for (int k = 0; k < third; k++ ) {
        printf("%f ", tmp[i*second*third + j*third + k]);
      }
      printf("\n");
    }
    printf("%%\n");
  }

  free(tmp);
}

void transposeMatrix(float *data_matrix, const int subject_size, const int time_size)
{
  float *tmp;
  cudaMalloc(&tmp, sizeof(float) * subject_size * REPEAT_TIMES * time_size);
  cudaMemcpy(tmp, data_matrix, sizeof(float) * subject_size * REPEAT_TIMES * time_size, cudaMemcpyDeviceToDevice);

  int total_works = REPEAT_TIMES * subject_size * time_size;
  int blocksize = 128;
  int nblock = total_works/blocksize + (total_works%blocksize==0?0:1);

  rearrangeMatrixPosition<<<nblock, blocksize>>>(data_matrix, tmp, subject_size, time_size);
  cudaDeviceSynchronize();

  cudaFree(tmp);
}

int main(int argc, char **argv)
{
  // srand(time(NULL));
  std::clock_t start;

  const int subject_size = 8, time_size = 440;
  const int matrix_works = subject_size * (subject_size - 1) / 2;
  const int total_works = REPEAT_TIMES * matrix_works;

  float *h_data_matrix, *h_coef_matrix, *h_isc_array;
  float *d_data_matrix, *d_coef_matrix, *d_isc_array;

  h_data_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * time_size * subject_size);
  h_coef_matrix = (float *)malloc(sizeof(float) * REPEAT_TIMES * matrix_works);
  h_isc_array = (float *)malloc(sizeof(float) * REPEAT_TIMES);

  start = std::clock();
  for(int i = 0; i < REPEAT_TIMES; i++)
    for(int j = 0; j < time_size; j++)
      for(int k = 0; k < subject_size; k++)
        h_data_matrix[i * time_size * subject_size + j * subject_size + k] = rand();
  printf("%% Generating data: %fs\n", (std::clock() - start) / (float) CLOCKS_PER_SEC);

  cudaMalloc(&d_data_matrix, sizeof(float) * REPEAT_TIMES * time_size * subject_size);
  cudaMalloc(&d_coef_matrix, sizeof(float) * REPEAT_TIMES * matrix_works);
  cudaMalloc(&d_isc_array, sizeof(float) * REPEAT_TIMES);
  cudaCheckErrors("cudaMalloc");

  cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size, cudaMemcpyHostToDevice);
  cudaMemset(d_coef_matrix, 0, sizeof(float) * REPEAT_TIMES * matrix_works);
  cudaCheckErrors("cudaMemcpy and cudaMemset");

  start = std::clock();
  transposeMatrix(d_data_matrix, subject_size, time_size);
  cudaCheckErrors("transposeMatrix");
  printf("%% transposeMatrix: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  int blocksize = 128;
  int nblock = total_works/blocksize + (total_works%blocksize==0?0:1);

  start = std::clock();
  calculateCorrelationCoefficientMatrixAndFisherTransform<<<nblock, blocksize>>>(d_coef_matrix, d_data_matrix, subject_size, time_size);
  cudaDeviceSynchronize();
  cudaCheckErrors("calculateCorrelationCoefficientMatrixAndFisherTransform");
  printf("%% calculateCorrelationCoefficientMatrixAndFisherTransform: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  nblock = REPEAT_TIMES/blocksize + (REPEAT_TIMES%blocksize==0?0:1);
  start = std::clock();
  calculateInterSubjectCorrelation<<<nblock, blocksize>>>(d_isc_array, d_coef_matrix, subject_size);
  cudaDeviceSynchronize();
  cudaCheckErrors("calculateInterSubjectCorrelation");
  printf("%% calculateInterSubjectCorrelation: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  cudaMemcpy(h_coef_matrix, d_coef_matrix, sizeof(float) * REPEAT_TIMES * matrix_works, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_isc_array, d_isc_array, sizeof(float) * REPEAT_TIMES, cudaMemcpyDeviceToHost);

  int result = thrust::count_if(thrust::host, h_isc_array, h_isc_array + REPEAT_TIMES, is_bigger_than(h_isc_array[0]));

  {
    cudaMemcpy(h_data_matrix, d_data_matrix, sizeof(float) * REPEAT_TIMES * subject_size * time_size, cudaMemcpyDeviceToHost);
    const int idx = rand() % REPEAT_TIMES;
    printf("%% idx: %d\n", idx);

    printf("data = [ ");
    for(int i = 0; i < time_size; i++) {
      for(int j = 0; j < subject_size; j++) {
        printf("%f ", h_data_matrix[idx * time_size * subject_size + i * subject_size + j]);
      }
      printf(";");
    }
    printf("];\n");
    printf("tmp=tril(corrcoef(data),-1);\n");
    printf("rr=tmp(find(tmp));\n");
    printf("z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));\n");
    printf("zm=mean(z)\n");
    printf("z'\n");
    printf("exit;\n");

    printf("%% coef_matrix: ");
    for (int i = 0; i < matrix_works; i++ )
        printf("%f ", h_coef_matrix[idx * matrix_works + i]);
    printf("\n");

    printf("%% --------------------\n");

    printf("%% mean: %f\n", h_isc_array[idx]);

    printf("%% --------------------\n");
    printf("%% %d\n", result);
  }

  cudaFree(d_data_matrix);
  cudaFree(d_coef_matrix);
  cudaFree(d_isc_array);
  free(h_data_matrix);
  free(h_coef_matrix);
  free(h_isc_array);

  cudaCheckErrors("cudaFree");

  return 0;
}

