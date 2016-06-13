#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "time_series_aaft.h"
#include "fmri_corr_coef.h"

#define RANDOM_TIMES 100000

// Functions.
std::vector<double> loadBrainData(std::string path, int &rows, int &columns);

// Sub functions.
std::vector<std::string> splitString(std::string str, char delimiter);

int main(int argc, char **argv) {
	// Parameters in command line.
	if (argc <= 2) {
		std::cout << "Have to provide over 2 paths.\n\n";
		return 0;
	}
	std::vector<std::string> csvPath;
	for (int i = 1; i < argc; i++)
		csvPath.push_back(argv[i]);

	// The matrix is expressed by 1D array.
	std::vector< std::vector<double> > data (argc - 1);
	int viewers = data.size();
	double *dataArr;
	// Amount of row and column.
	int rows = 0, columns = 0;
	// GPU variables.
	int threads = 32, blocks;
	double *data_g, *aaft_g, *coef_g;

	// Release memory.
	cudaFree(0);

	// Load CSV each viewer's data.
	for (int i = 0; i < viewers; i++)
		data.at(i) = loadBrainData(csvPath.at(i), rows, columns);

	cudaMalloc(&data_g, sizeof(double) * columns * viewers);
	cudaMalloc(&aaft_g, sizeof(double) * columns * viewers * RANDOM_TIMES);
	cudaMalloc(&coef_g, sizeof(double) * RANDOM_TIMES);

	// Calculate the new time series each row by Amplitude Adjusted Fourier Transform (AAFT), then calculate correlation coefficient.
	// Most of the case, rows are 10,242 times; columns are 450 times.
	for (int i = 0; i < rows; i++) {
		// Concatenate the specific position (row) each viewer's data.
		std::vector<double> subdata;
		for (int j = 0; j < viewers; j++)
			subdata.insert(subdata.end(), data.at(j).begin() + i * columns, data.at(j).begin() + (i + 1) * columns);

		// Convert the vector to array.
		dataArr = &subdata[0];
		cudaMemcpy(data_g, dataArr, sizeof(double) * columns * viewers, cudaMemcpyHostToDevice);
		cudaMemset(aaft_g, 0, sizeof(double) * columns * viewers * RANDOM_TIMES);
		cudaMemset(coef_g, 0, sizeof(double) * RANDOM_TIMES);

		// Kernel - Amplitude Adjusted Fourier Transform (AAFT).
		amplitudeAdjustedFourierTransform(aaft_g, data_g, viewers, RANDOM_TIMES, columns);

		// Kernel - Correlation coefficient.
		blocks = ceil(RANDOM_TIMES / threads);
		correlationCoefficient <<<blocks, threads>>> (coef_g, aaft_g, (i + 1), columns);
	}

	// Release memory.
	cudaFree(data_g);
	cudaFree(aaft_g);
	cudaFree(coef_g);

	// Message.
	std::cout << "Data in CPU: " << " [rows: " << rows << ", columns: " << columns << "], total size: " << data.at(0).size() << ".\n\n";

	return 0;
}

// Load the fMRI data that format is CSV.
std::vector<double> loadBrainData(std::string path, int &rows, int &columns) {
	// CSV file.
	std::ifstream file(path);
	// Line string cache.
	std::string line;
	// 2D matrix expressed by 1D vector.
	std::vector<double> data;

	// Initial.
	rows = 0, columns = 0;

	while (std::getline(file, line)) {
		// Split the string to a vector.
		std::vector<std::string> stringVector = splitString(line, ',');
		// Convert the string vector to double vector.
		std::vector<double> row(stringVector.size());
		std::transform(stringVector.begin(), stringVector.end(), row.begin(), [](const std::string &val) { return std::stod(val); });
		// Update the columns & rows.
		if (! rows) columns = row.size();
		rows ++;
		// Insert the row into data.
		data.insert(data.end(), row.begin(), row.end());
	}
	return data;
}

// Split the string by delimiter to a vector.
std::vector<std::string> splitString(std::string str, char delimiter) {
	std::vector<std::string> result;
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delimiter))
		result.push_back(token);
	return result;
}