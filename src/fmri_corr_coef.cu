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

__device__ double calSum(const double *x, const int length)
{
  double sum = 0;
  for(int i = 0; i < length; i++)
    sum += x[i];
  return sum;
}

__device__ double calAvg(const double *x, const int length)
{
  return calSum(x, length) / length;
}

__device__ double calSquareSum(const double *x, const int length)
{
  double sum = 0;
  for(int i = 0; i < length; i++)
    sum += x[i] * x[i];
  return sum;
}

__device__ double calMultiplySum(const double *x, const double *y, const int length)
{
  double sum = 0;
  for(int i = 0; i < length; i++)
    sum += x[i] * y[i];
  return sum;
}

__device__ double calStd(const double *x, const int length)
{
  const double x_square_sum = calSquareSum(x, length);
  const double x_avg = calAvg(x, length);

  return sqrt((x_square_sum - length * x_avg * x_avg) / (length - 1));
}

__device__ double calCorrCoef(const double *x, const double *y, const int length)
{
  const double xy_sum = calMultiplySum(x, y, length);
  const double x_avg = calAvg(x, length);
  const double y_avg = calAvg(y, length);
  const double x_std = calStd(x, length);
  const double y_std = calStd(y, length);

  return (xy_sum - length * x_avg * y_avg) / ((length - 1) * x_std * y_std);
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

  const double *data_matrix = all_data_matrix + n_matrix * subject_size * time_size;
  const double *data_x = data_matrix + x * time_size;
  const double *data_y = data_matrix + y * time_size;
  const double coef = calCorrCoef(data_x, data_y, time_size);

  all_corr_coef_matrix[idx] = coef;
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

  isc_array[idx] = sum / matrix_works;
}

// void printMatrix(const double *data, const int first, const int second, const int third)
// {
//   double *tmp = (double *)malloc(sizeof(double) * first * second * third);
//   cudaMemcpy(tmp, data, sizeof(double) * first * second * third, cudaMemcpyDeviceToHost);
// 
//   printf("%% 1st:%d 2nd:%d 3rd:%d\n", first, second, third);
//   for (int i = 0; i < first; i++ ) {
//     for (int j = 0; j < second; j++ ) {
//       printf("%% ");
//       for (int k = 0; k < third; k++ ) {
//         printf("%f ", tmp[i*second*third + j*third + k]);
//       }
//       printf("\n");
//     }
//     printf("%%\n");
//   }
// 
//   free(tmp);
// }

void correlationCoefficient(double *d_isc_array, const double *d_aaft_matrix, const int subject_size, const int time_size, const int repeat_times)
{
  // std::clock_t start;

  const int blocksize = 128;
  int total_works, nblock;

  double *d_coef_matrix;

  cudaMalloc(&d_coef_matrix, sizeof(double) * repeat_times * subject_size * (subject_size - 1) / 2);
  cudaCheckErrors("cudaMalloc");

  // start = std::clock();
  total_works = repeat_times * subject_size * (subject_size - 1) / 2;
  nblock = total_works/blocksize + (total_works%blocksize==0?0:1);
  calculateCorrelationCoefficientMatrix<<<nblock, blocksize>>>(d_coef_matrix, d_aaft_matrix, subject_size, time_size, repeat_times);
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
  //   printMatrix(d_aaft_matrix, repeat_times, subject_size, time_size);

  //   double *h_data_matrix, *h_isc_array;
  //   h_data_matrix = (double *)malloc(sizeof(double) * repeat_times * subject_size * time_size);
  //   h_isc_array = (double *)malloc(sizeof(double) * repeat_times);

  //   cudaMemcpy(h_data_matrix, d_aaft_matrix, sizeof(double) * repeat_times * subject_size * time_size, cudaMemcpyDeviceToHost);
  //   cudaMemcpy(h_isc_array, d_isc_array, sizeof(double) * repeat_times, cudaMemcpyDeviceToHost);

  //   const int idx = rand() % repeat_times;
  //   printf("%% idx: %d\n", idx);

  //   printf("data = [ ");
  //   for(int j = 0; j < time_size; j++) {
  //     for(int i = 0; i < subject_size; i++) {
  //       printf("%f ", h_data_matrix[idx * subject_size * time_size + i * time_size + j]);
  //     }
  //     printf(";");
  //   }
  //   printf("];\n");
  //   printf("tmp=tril(corrcoef(data),-1);\n");
  //   printf("rr=tmp(find(tmp));\n");
  //   printf("zm=mean(rr)\n");
  //   printf("exit;\n");

  //   printf("%% --------------------\n");
  //   printf("%% mean: %f\n", h_isc_array[idx]);
  // }

  cudaFree(d_coef_matrix);
}

int main(int argc, char **argv)
{
  srand(time(NULL));
  std::clock_t start;

  const int subject_size = 8, time_size = 440, repeat_times = 25000;

  double *h_data_matrix;
  double *d_data_matrix, *d_isc_array;

  for (int repeat = 0; repeat < 10; repeat ++ ) {
    h_data_matrix = (double *)malloc(sizeof(double) * repeat_times * time_size * subject_size);

    start = std::clock();
    for(int i = 0; i < repeat_times * subject_size * time_size; i++)
      h_data_matrix[i] = rand();
    printf("%% %2d Generating data: %fs\n", repeat, (std::clock() - start) / (double) CLOCKS_PER_SEC);

    cudaMalloc(&d_data_matrix, sizeof(double) * repeat_times * time_size * subject_size);
    cudaMalloc(&d_isc_array, sizeof(double) * repeat_times);
    cudaCheckErrors("cudaMalloc");

    cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(double) * repeat_times * subject_size * time_size, cudaMemcpyHostToDevice);
    cudaCheckErrors("cudaMemcpy");

    start = std::clock();
    correlationCoefficient(d_isc_array, d_data_matrix, subject_size, time_size, repeat_times);
    printf("%% %2d correlationCoefficient: %fs\n", repeat, (std::clock() - start) / (double) CLOCKS_PER_SEC);

    free(h_data_matrix);
    cudaFree(d_data_matrix);
    cudaFree(d_isc_array);
    cudaCheckErrors("cudaFree");
  }

  return 0;
}

