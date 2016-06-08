/* For Kai-Hsiang */

#include <iostream>

// 'coef' is output.
__global__ void correlationCoefficient(double *coef, const double *aaft, const int row, const int maxColumns) {
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int column = 1 + idx % maxColumns;

	/* Code below. */




	/* Code above. */
	return;
}