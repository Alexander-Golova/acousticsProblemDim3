#include "stdafx.h"
#include "directProblemUtils.h"

using namespace std;

void GetSubstantiveMatrix(const vector<complex<float>> & a, const vector<float> & xi, vector<vector<complex<float>>> & substantiveMatrix)
{
	size_t ii, jj;
	size_t coord;
	complex<float> sumOfTheCoefficients;

	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				ii = N_SQUARED * i + N * j + k;
				sumOfTheCoefficients = { 0.0f, 0.0f }; // TODO
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						for (size_t r = 0; r < N; ++r)
						{
							jj = N_SQUARED * p + N * q + r;
							coord = N_FIFTH_DEGREE * i + N_FOURTH_DEGREE * j + N_QUBE * k + N_SQUARED * p + N * q + r;
							substantiveMatrix[ii][jj] += a[coord] * xi[jj]; // TODO
						}
					}
				}
				substantiveMatrix[ii][ii] = { 1.0f, 0.0f };
			}
		}
	}

}

void GetRightPartEquation(const Sources & source, const size_t count, std::vector<std::complex<float>> & rightPartEquation)
{
	size_t ii;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				ii = N_SQUARED * i + N * j  + k;
				rightPartEquation[ii] = source.Function(source.node[count], i * step, j * step, k * step);
			}
		}
	}
}
