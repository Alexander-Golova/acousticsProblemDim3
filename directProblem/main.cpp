#include "stdafx.h"
#include "taskData.h"
#include "Sources.h"
#include "matrix_utils.h"
#include "exact_solution.h"

using namespace std;

int main()
{
	const Source source;

	vector<vector<vector<complex<double>>>> u(NUMBER_PARTITION_POINTS + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINTS + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINTS + 1, complex<double>())));

	vector<vector<vector<double>>> xi(NUMBER_PARTITION_POINTS + 1,
		vector<vector<double>>(NUMBER_PARTITION_POINTS + 1,
			vector<double>(NUMBER_PARTITION_POINTS + 1, 0.0)));

	GetExactSolution(xi);

	clock_t time = clock();
	clock_t timeBegin = clock();



	return 0;
}

