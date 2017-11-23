#include "stdafx.h"
#include "taskData.h"
#include "directProblem_utils.h"
#include "Sources.h"

using namespace std;

void GetSubstantiveMatrix(const vector<vector<vector<vector<vector<vector<complex<double>>>>>>> & a,
	const vector<vector<vector<double>>> & xi, vector<vector<complex<double>>> & substantiveMatrix)
{
	size_t ii, jj;
	complex<double> sumOfTheCoefficients;

	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				ii = i * (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1) + j * (NUMBER_PARTITION_POINTS + 1) + k;
				sumOfTheCoefficients = (0.0, 0.0);
				for (size_t p = 0; p < NUMBER_PARTITION_POINTS; ++p)
				{
					for (size_t q = 0; q < NUMBER_PARTITION_POINTS; ++q)
					{
						for (size_t r = 0; r < NUMBER_PARTITION_POINTS; ++r)
						{
							jj = p * (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1) + q * (NUMBER_PARTITION_POINTS + 1) + r;
							substantiveMatrix[ii][jj] += a[i][j][k][p][q][r] * xi[p][q][r];
						}
					}
				}
				substantiveMatrix[ii][ii] = (1.0, 0.0);
			}
		}
	}
}

void GetRightPartEquation(const Source & source, size_t count, vector<complex<double>> & rightPartEquation)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				ii = i * (NUMBER_PARTITION_POINTS + 1)
					* (NUMBER_PARTITION_POINTS + 1) + j * (NUMBER_PARTITION_POINTS + 1) + k;
				rightPartEquation[ii] = source.Function(source.node[count], i * h, j * h, k * h);
			}
		}
	}
}

void InverseRenumbering(const vector<complex<double>> & numbered_u, vector<vector<vector<complex<double>>>> & u)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				ii = i * (NUMBER_PARTITION_POINTS + 1)
					* (NUMBER_PARTITION_POINTS + 1) + j * (NUMBER_PARTITION_POINTS + 1) + k;
				u[i][j][k] = numbered_u[ii];
			}
		}
	}
}

void GetOverlineU(const Source & source, size_t count,
	const vector<vector<vector<vector<vector<complex<double>>>>>> & overline_a,
	const vector<vector<vector<double>>> & xi, const vector<vector<vector<complex<double>>>> & u,
	vector<vector<complex<double>>> & overline_u)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			overline_u[i][j] = source.Function(source.node[count], i * h, j * h, receiver);
		}
	}
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINTS; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINTS; ++q)
				{
					for (size_t r = 0; r < NUMBER_PARTITION_POINTS; ++r)
					{
						overline_u[i][j] -= overline_a[i][j][p][q][r] * xi[p][q][r] * u[p][q][r];
					}
				}
			}
		}
	}
}
