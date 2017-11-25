#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<vector<double>>> & xi)
{
	double sigma = 8;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINT; ++k)
			{
				xi[i][j][k] = 0.2;
			}
		}
	}
}


void WriteSolutionFile(vector<vector<vector<double>>> & xi)
{
	ofstream file_xi("exact_xi.txt");
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINT; ++k)
			{
				file_xi << xi[i][j][k] << " ";
			}
		}
	}
	file_xi.close();
}
