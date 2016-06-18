#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <thrust/count.h>
#include <thrust/execution_policy.h>

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

__device__ double calSum(const double *x, const int length, const int step)
{
  double sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step];
  }
  return sum;
}

__device__ double calAvg(const double *x, const int length, const int step)
{
  return calSum(x, length, step) / length;
}

__device__ double calSquareSum(const double *x, const int length, const int step)
{
  double sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step] * x[i*step];
  }
  return sum;
}

__device__ double calMultiplySum(const double *x, const double *y, const int length, const int step)
{
  double sum = 0;
  for(int i = 0; i < length; i++) {
    sum += x[i*step] * y[i*step];
  }
  return sum;
}

__device__ double calStd(const double *x, const int length, const int step)
{
  const double x_square_sum = calSquareSum(x, length, step);
  const double x_avg = calAvg(x, length, step);

  return sqrt((x_square_sum - length * x_avg * x_avg) / (length - 1));
}

__device__ double calCorrCoef(const double *x, const double *y, const int length, const int step)
{
  const double xy_sum = calMultiplySum(x, y, length, step);
  const double x_avg = calAvg(x, length, step);
  const double y_avg = calAvg(y, length, step);
  const double x_std = calStd(x, length, step);
  const double y_std = calStd(y, length, step);

  return (xy_sum - length * x_avg * y_avg) / ((length - 1) * x_std * y_std);
}

__device__ double calFisherTransform(const double x, const int time_size)
{
  // z=0.5.*log((1+rr)./(1-rr));
  return 0.5 * log((1+x) / (1-x));
}

__device__ double calInverseFisherTransform(const double x)
{
  // zm= (exp(2.*zm)-1)./(exp(2.*zm)+1);
  return (exp(2*x) - 1) / (exp(2*x) + 1);
}

__global__ void calculateCorrelationCoefficientMatrix(double *all_corr_coef_matrix, const double *all_data_matrix, const int subject_size, const int time_size, const int repeat_times)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int matrix_works = subject_size * (subject_size - 1) / 2;
  if (idx >= repeat_times * matrix_works)
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

  const double *data_matrix = all_data_matrix + n_matrix * time_size * subject_size;
  const double coef = calCorrCoef(data_matrix + x, data_matrix + y, time_size, subject_size);
  const double zvalue = calFisherTransform(coef, time_size);

  all_corr_coef_matrix[idx] = zvalue;
}


__global__ void calculateInterSubjectCorrelation(double *isc_array, const double *all_corr_coef_matrix, const int subject_size, const int repeat_times)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= repeat_times)
    return;

  const int matrix_works = subject_size * (subject_size - 1) / 2;
  const double *corr_coef_matrix = all_corr_coef_matrix + idx * matrix_works;

  double sum = 0;

  for (int i = 0; i < matrix_works; i++)
    sum += corr_coef_matrix[i];

  const double mean = sum / matrix_works;
  isc_array[idx] = calInverseFisherTransform(mean);
}

__global__ void rearrangeMatrixPosition(double *data_matrix, const double *source_matrix, const int subject_size, const int time_size, const int repeat_times)
{
  // 1st subject_size, 2nd repeat_times, 3rd time_size
  // to
  // 1st repeat_times, 2nd time_size, 3rd subject_size

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= repeat_times * time_size * subject_size)
    return;

  const int subject_idx = idx / (repeat_times * time_size);
  const int repeat_idx = (idx % (repeat_times * time_size)) / time_size;
  const int time_idx = (idx % (repeat_times * time_size)) % time_size;

  const int data_idx = repeat_idx * time_size * subject_size + time_idx * subject_size + subject_idx;

  data_matrix[data_idx] = source_matrix[idx];
}

void printMatrix(const double *data, const int first, const int second, const int third)
{
  double *tmp = (double *)malloc(sizeof(double) * first * second * third);
  cudaMemcpy(tmp, data, sizeof(double) * first * second * third, cudaMemcpyDeviceToHost);

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

void correlationCoefficient(double *d_isc_array, const double *d_aaft_matrix, const int subject_size, const int time_size, const int repeat_times)
{
  // std::clock_t start;

  const int blocksize = 128;
  int total_works, nblock;

  double *d_data_matrix, *d_coef_matrix;

  cudaMalloc(&d_data_matrix, sizeof(double) * repeat_times * time_size * subject_size);
  cudaMalloc(&d_coef_matrix, sizeof(double) * repeat_times * subject_size * (subject_size - 1) / 2);
  cudaCheckErrors("cudaMalloc");

  // start = std::clock();
  total_works = repeat_times * subject_size * time_size;
  nblock = total_works/blocksize + (total_works%blocksize==0?0:1);
  rearrangeMatrixPosition<<<nblock, blocksize>>>(d_data_matrix, d_aaft_matrix, subject_size, time_size, repeat_times);
  cudaDeviceSynchronize();
  cudaCheckErrors("rearrangeMatrixPosition");
  // printf("%% transposeMatrix: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  // start = std::clock();
  total_works = repeat_times * subject_size * (subject_size - 1) / 2;
  nblock = total_works/blocksize + (total_works%blocksize==0?0:1);
  calculateCorrelationCoefficientMatrix<<<nblock, blocksize>>>(d_coef_matrix, d_data_matrix, subject_size, time_size, repeat_times);
  cudaDeviceSynchronize();
  cudaCheckErrors("calculateCorrelationCoefficientMatrix");
  // printf("%% calculateCorrelationCoefficientMatrix: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  // start = std::clock();
  nblock = repeat_times/blocksize + (repeat_times%blocksize==0?0:1);
  calculateInterSubjectCorrelation<<<nblock, blocksize>>>(d_isc_array, d_coef_matrix, subject_size, repeat_times);
  cudaDeviceSynchronize();
  cudaCheckErrors("calculateInterSubjectCorrelation");
  // printf("%% calculateInterSubjectCorrelation: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

  // {
  //   printMatrix(d_aaft_matrix, subject_size, repeat_times, time_size);
  //   printMatrix(d_data_matrix, repeat_times, time_size, subject_size);

  //   double *h_data_matrix, *h_coef_matrix, *h_isc_array;
  //   h_data_matrix = (double *)malloc(sizeof(double) * repeat_times * subject_size * time_size);
  //   h_coef_matrix = (double *)malloc(sizeof(double) * total_works);
  //   h_isc_array = (double *)malloc(sizeof(double) * repeat_times);

  //   cudaMemcpy(h_data_matrix, d_data_matrix, sizeof(double) * repeat_times * subject_size * time_size, cudaMemcpyDeviceToHost);
  //   cudaMemcpy(h_coef_matrix, d_coef_matrix, sizeof(double) * total_works, cudaMemcpyDeviceToHost);
  //   cudaMemcpy(h_isc_array, d_isc_array, sizeof(double) * repeat_times, cudaMemcpyDeviceToHost);

  //   const int idx = rand() % repeat_times;
  //   printf("%% idx: %d\n", idx);

  //   printf("data = [ ");
  //   for(int i = 0; i < time_size; i++) {
  //     for(int j = 0; j < subject_size; j++) {
  //       printf("%f ", h_data_matrix[idx * time_size * subject_size + i * subject_size + j]);
  //     }
  //     printf(";");
  //   }
  //   printf("];\n");
  //   printf("tmp=tril(corrcoef(data),-1);\n");
  //   printf("rr=tmp(find(tmp));\n");
  //   printf("z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));\n");
  //   printf("zm=mean(z)\n");
  //   printf("z'\n");
  //   printf("exit;\n");

  //   printf("%% coef_matrix: ");
  //   const int matrix_works = subject_size * (subject_size - 1) / 2;
  //   for (int i = 0; i < matrix_works; i++ )
  //       printf("%f ", h_coef_matrix[idx * matrix_works + i]);
  //   printf("\n");

  //   printf("%% --------------------\n");

  //   printf("%% mean: %f\n", h_isc_array[idx]);
  // }

  cudaFree(d_data_matrix);
  cudaFree(d_coef_matrix);
}

// int main(int argc, char **argv)
// {
//   srand(time(NULL));
//   std::clock_t start;
// 
//   const int subject_size = 8, time_size = 440, repeat_times = 10000;
// 
//   double *h_data_matrix;
//   double *d_data_matrix, *d_isc_array;
// 
//   h_data_matrix = (double *)malloc(sizeof(double) * repeat_times * time_size * subject_size);
// 
//   start = std::clock();
//   for(int i = 0; i < repeat_times; i++)
//     for(int j = 0; j < time_size; j++)
//       for(int k = 0; k < subject_size; k++)
//         h_data_matrix[i * time_size * subject_size + j * subject_size + k] = rand();
//   printf("%% Generating data: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);
// 
//   cudaMalloc(&d_data_matrix, sizeof(double) * repeat_times * time_size * subject_size);
//   cudaMalloc(&d_isc_array, sizeof(double) * repeat_times);
//   cudaCheckErrors("cudaMalloc");
// 
//   cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(double) * repeat_times * subject_size * time_size, cudaMemcpyHostToDevice);
//   cudaCheckErrors("cudaMemcpy");
// 
//   start = std::clock();
//   correlationCoefficient(d_isc_array, d_data_matrix, subject_size, time_size, repeat_times);
//   printf("%% correlationCoefficient: %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);
// 
//   free(h_data_matrix);
//   cudaFree(d_data_matrix);
//   cudaFree(d_isc_array);
//   cudaCheckErrors("cudaFree");
// 
//   return 0;
// }

