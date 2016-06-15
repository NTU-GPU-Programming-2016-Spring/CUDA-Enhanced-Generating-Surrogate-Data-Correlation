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

__device__ double calFisherTransform(const double x, const int time_size)
{
  // z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));
  return 0.5 * log((1+x) / (1-x)) / rsqrt((double)time_size/2.34 - 3);
}

__device__ double calInverseFisherTransform(const double x)
{
  // zm= (exp(2.*zm)-1)./(exp(2.*zm)+1);
  return (exp(2*x) - 1) / (exp(2*x) + 1);
}

__global__ void calculateInterSubjectCorrelation(double *isc_array, const double *all_data_matrix, const int subject_size, const int time_size, const int repeat_times)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= repeat_times)
    return;

  const double *data_matrix = all_data_matrix + idx * subject_size * time_size;
  double sum = 0;

  for (int i = 0; i < subject_size; i++) {
    const double *data_x = data_matrix + i * time_size;

    for (int j = i+1; j < subject_size; j++) {
      const double *data_y = data_matrix + j * time_size;
      const double coef = calCorrCoef(data_x, data_y, time_size);
      const double zvalue = calFisherTransform(coef, time_size);
      sum += zvalue;
    }
  }

  const int matrix_works = subject_size * (subject_size - 1) / 2;
  const double zmean = sum / matrix_works;
  isc_array[idx] = calInverseFisherTransform(zmean);
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
  int nblock;

  // start = std::clock();
  nblock = repeat_times/blocksize + (repeat_times%blocksize==0?0:1);
  calculateInterSubjectCorrelation<<<nblock, blocksize>>>(d_isc_array, d_aaft_matrix, subject_size, time_size, repeat_times);
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
  //       printf("%f ", h_data_matrix[idx * time_size * subject_size + i * time_size + j]);
  //     }
  //     printf(";");
  //   }
  //   printf("];\n");
  //   printf("tmp=tril(corrcoef(data),-1);\n");
  //   printf("rr=tmp(find(tmp));\n");
  //   printf("z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(data,1)/2.34-3));\n");
  //   printf("zm=mean(z);\n");
  //   printf("corr_mean=(exp(2.*zm)-1)./(exp(2.*zm)+1)\n");
  //   printf("exit;\n");

  //   printf("%% --------------------\n");
  //   printf("%% mean: %f\n", h_isc_array[idx]);
  // }
}

int main(int argc, char **argv)
{
  srand(time(NULL));
  std::clock_t start;

  const int subject_size = 8, time_size = 440, repeat_times = 25000;

  double *h_data_matrix;
  double *d_data_matrix, *d_isc_array;

  for (int repeat = 0; repeat < 10; repeat++ ) {
    h_data_matrix = (double *)malloc(sizeof(double) * repeat_times * time_size * subject_size);

    start = std::clock();
    for(int i = 0; i < repeat_times * time_size * subject_size; i++)
      h_data_matrix[i] = rand();
    printf("%% %d Generating data: %fs\n", repeat, (std::clock() - start) / (double) CLOCKS_PER_SEC);

    cudaMalloc(&d_data_matrix, sizeof(double) * repeat_times * time_size * subject_size);
    cudaMalloc(&d_isc_array, sizeof(double) * repeat_times);
    cudaCheckErrors("cudaMalloc");

    cudaMemcpy(d_data_matrix, h_data_matrix, sizeof(double) * repeat_times * subject_size * time_size, cudaMemcpyHostToDevice);
    cudaCheckErrors("cudaMemcpy");

    start = std::clock();
    correlationCoefficient(d_isc_array, d_data_matrix, subject_size, time_size, repeat_times);
    printf("%% %d correlationCoefficient: %fs\n", repeat, (std::clock() - start) / (double) CLOCKS_PER_SEC);

    free(h_data_matrix);
    cudaFree(d_data_matrix);
    cudaFree(d_isc_array);
    cudaCheckErrors("cudaFree");
  }

  return 0;
}

