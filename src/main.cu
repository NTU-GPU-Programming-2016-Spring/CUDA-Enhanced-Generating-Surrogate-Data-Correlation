#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

#include <helper_cuda.h>

#include "time_series_aaft.h"
#include "fmri_corr_coef.h"

#define RANDOM_TIMES 2000
#define RANDOM_TIMES_UNIT 4096

// Functions.
std::vector<double> loadBrainData(std::string path, int &rows, int &columns);
std::vector<double> repeatVector(std::vector<double> data, int times);

// Sub functions.
std::vector<std::string> splitString(std::string str, char delimiter);

int main(int argc, char **argv) {
	// Parameters in command line.
	if (argc <= 3) {
		std::cout << "Have to provide over 2 paths.\n\n";
		return 0;
	}
	std::vector<std::string> csvPath;
	for (int i = 2; i < argc; i++)
		csvPath.push_back(argv[i]);

	std::clock_t start;
	// The matrix is expressed by 1D array.
	std::vector< std::vector<double> > data (argc - 2);
	int viewers = data.size();
	double *dataArr;
	// Amount of row and column.
	int rows = 0, columns = 0;
	// GPU variables.
	double *data_g, *aaft_g, *coef_g;
	// Output file.
	std::stringstream ss, ssFilename;
	std::ofstream outputFile;
	ssFilename << "debug/" << argv[1] << ".csv";
	outputFile.open(ssFilename.str());

	// Release memory.
	cudaFree(0);

	start = std::clock();
	// Load CSV each viewer's data.
	fprintf(stderr, "Loading CSV data ...\r");
	for (int i = 0; i < viewers; i++)
		data.at(i) = loadBrainData(csvPath.at(i), rows, columns);
	fprintf(stderr, "Loading CSV data ... done with %fs\n", (std::clock() - start) / (double) CLOCKS_PER_SEC);

	start = std::clock();
	cudaMalloc(&data_g, sizeof(double) * columns * viewers * RANDOM_TIMES_UNIT);
	cudaMalloc(&aaft_g, sizeof(double) * columns * viewers * RANDOM_TIMES_UNIT);
	cudaMalloc(&coef_g, sizeof(double) * RANDOM_TIMES);

	// Calculate the new time series each row by Amplitude Adjusted Fourier Transform (AAFT), then calculate correlation coefficient.
	// Most of the case, rows are 10,242 times; columns are 450 times.
	for (int i = 0; i < rows; i++) {
		int executed_times = 0;

		fprintf(stderr, "Processing row ... %5d\r", i);
		while (executed_times < RANDOM_TIMES) {
			const int remain_times = RANDOM_TIMES - executed_times;
			const int random_times = remain_times >= RANDOM_TIMES_UNIT ? RANDOM_TIMES_UNIT : remain_times;

			// Concatenate the specific position (row) each viewer's data.
			std::vector<double> subdata;
			for (int j = 0; j < viewers; j++) {
				std::vector<double> viewer_data(data.at(j).begin() + i * columns, data.at(j).begin() + (i + 1) * columns);
				viewer_data = repeatVector(viewer_data, random_times);
				subdata.insert(subdata.end(), viewer_data.begin(), viewer_data.end());
			}

			// Convert the vector to array.
			dataArr = &subdata[0];
			cudaMemcpy(data_g, dataArr, sizeof(double) * columns * viewers * random_times, cudaMemcpyHostToDevice);
			cudaMemset(aaft_g, 0, sizeof(double) * columns * viewers * random_times);
			cudaMemset(coef_g + executed_times, 0, sizeof(double) * random_times);

			// Kernel - Amplitude Adjusted Fourier Transform (AAFT).
			amplitudeAdjustedFourierTransform(aaft_g, data_g, viewers, random_times, columns);

			// Kernel - Correlation coefficient.
			correlationCoefficient(coef_g + executed_times, aaft_g, viewers, columns, random_times);

			executed_times += random_times;
		}

		double *coef_cpu = (double *)malloc(sizeof(double) * RANDOM_TIMES);
		checkCudaErrors(cudaMemcpy(coef_cpu, coef_g, sizeof(double) * RANDOM_TIMES, cudaMemcpyDeviceToHost));
	
		// Write into file.
		for (int i = 0; i < RANDOM_TIMES; i++)
		    ss << coef_cpu[i] << ((i + 1) != RANDOM_TIMES ? "," : "\n");
	}
	
	outputFile << ss.str();
	outputFile.close();
	
	fprintf(stderr, "Processing row ... %5s with %fs\n", "done", (std::clock() - start) / (double) CLOCKS_PER_SEC);

	// Release memory.
	cudaFree(aaft_g);
	cudaFree(coef_g);
	cudaFree(data_g);

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

// copy a vector multiple times
std::vector<double> repeatVector(std::vector<double> data, int times) {
	std::vector<double> repeat(data.size() * times);
	memcpy(&repeat[0], &data[0], sizeof(double) * data.size());

	size_t num_copied = data.size(), num_total = repeat.size();

	while (num_copied * 2 < num_total) {
		memcpy(&repeat[num_copied], &repeat[0], sizeof(double) * num_copied);
		num_copied *= 2;
	}

	memcpy(&repeat[num_copied], &repeat[0], sizeof(double) * (num_total - num_copied));

	return repeat;
}

