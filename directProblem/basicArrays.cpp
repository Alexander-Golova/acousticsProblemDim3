#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"

using namespace std;

void GetBasicArrays(vector<vector<vector<vector<vector<vector<complex<double>>>>>>> & a,
	vector<vector<vector<vector<vector<complex<double>>>>>> & overline_a)
{
	// для индексов метода квадратур
	vector<double> index(NUMBER_PARTITION_POINTS + 1);
	for (size_t i = 1; i < NUMBER_PARTITION_POINTS; ++i)
	{
		if (i % 2 != 0)
		{
			index[i] = 4.0 / 3;
		}
		else
		{
			index[i] = 2.0 / 3;
		}
	}
	index[0] = 1.0 / 3;
	index[NUMBER_PARTITION_POINTS] = 1.0 / 3;

	// нахождение массива a
	double dist;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				for (size_t p = 0; p < NUMBER_PARTITION_POINTS; ++p)
				{
					for (size_t q = 0; q < NUMBER_PARTITION_POINTS; ++q)
					{
						for (size_t r = 0; r < NUMBER_PARTITION_POINTS; ++r)
						{
							if ((i != p) || (q != j) || (r != k))
							{
								dist = sqrt((i - p) * (i - p) + (j - q) * (j - q) + (k - r) * (k - r));
								dist *= h;
								a[i][j][k][p][q][r] = index[p] * index[q] * index[r];
								a[i][j][k][p][q][r] *= exp(-dist * (0.0, 1.0) * OMEGA / C_0) / dist;
								a[i][j][k][p][q][r] *= OMEGA * OMEGA * h * h * h;
								a[i][j][k][p][q][r] /= FOUR_PI;
							}
						}
					}
				}
			}
		}
	}

	// нахождение массива overline_a
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
						dist = sqrt((i - p) * (i - p) + (j - q) * (j - q) + (receiver - r) * (receiver - r));
						dist *= h;
						overline_a[i][j][p][q][r] = index[p] * index[q] * index[r];
						overline_a[i][j][p][q][r] *= exp(-dist * (0.0, 1.0) * OMEGA / C_0) / dist;
						overline_a[i][j][p][q][r] *= OMEGA * OMEGA * h * h * h;
						overline_a[i][j][p][q][r] /= FOUR_PI;
					}
				}
			}
		}
	}
}

void WriteBasicArraysFile(const vector<vector<vector<vector<vector<vector<complex<double>>>>>>> & a,
	const vector<vector<vector<vector<vector<complex<double>>>>>> & overline_a)
{
	ofstream f_overline_a("matrix_overline_a.txt");
	f_overline_a << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t p = 0; p <= NUMBER_PARTITION_POINTS; ++p)
			{
				for (size_t q = 0; q <= NUMBER_PARTITION_POINTS; ++q)
				{
					for (size_t r = 0; r <= NUMBER_PARTITION_POINTS; ++r)
					{
						f_overline_a << overline_a[i][j][p][q][r] << " ";
					}
				}
			}
		}
	}
	f_overline_a.close();

	ofstream f_a("matrix_a.txt");
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
			{
				for (size_t p = 0; p <= NUMBER_PARTITION_POINTS; ++p)
				{
					for (size_t q = 0; q <= NUMBER_PARTITION_POINTS; ++q)
					{
						for (size_t r = 0; r <= NUMBER_PARTITION_POINTS; ++r)
						{
							f_a << a[i][j][k][p][q][r] << " ";
						}
					}
				}
			}
		}
	}
	f_a.close();
}
