#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<vector<double>>> & xi)
{
	double sigma = 16;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				xi[i][j][k] = 0.4 * exp(-((i * h - 0.6) * (i * h - 0.6) +
					(j * h - 0.6) * (j * h - 0.6) + (k * h - 0.6) * (k * h - 0.6)) * sigma);
			}
		}
	}
}

void WriteSolutionFile(const vector<vector<vector<double>>> & xi)
{
	ofstream file_xi("exact_xi.txt");
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				file_xi << xi[i][j][k] << " ";
			}
		}
	}
	file_xi.close();
}
