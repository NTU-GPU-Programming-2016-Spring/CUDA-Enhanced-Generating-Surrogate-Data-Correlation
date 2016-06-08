/* For Brian */

#include <iostream>

// 'aaft' is output.
__global__ void amplitudeAdjustedFourierTransform(double *aaft, const double *data, const int row, const int maxColumns) {
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int column = 1 + idx % maxColumns;

	/* Code below. */


	// if (blockIdx.x == 0) {
	// 	printf("[AAFT] index: %d, value: %f. (row: %d, column: %d)\n", idx, data[idx], row, column);
	// }


	/* Code above. */
	return;
}
