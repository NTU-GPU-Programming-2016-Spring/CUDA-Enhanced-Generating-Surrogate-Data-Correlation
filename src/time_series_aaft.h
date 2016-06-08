#pragma once

__global__ void amplitudeAdjustedFourierTransform(double *aaft, const double *data, const int row, const int maxColumns);