#include "stdafx.h"
#include "Sources.h"

using namespace std;

complex<float> Sources::Function(const Point source, const float x, const float y, const float z) const
{
	float dist = sqrtf(pow(x - source.x, 2) + pow(y - source.y, 2) + pow(z - source.z, 2));
	return -exp(I * omega * dist / c_0) * INV_FOUR_PI / dist;
}

void WriteSourceValues(const Sources & source, string name) noexcept
{
	ofstream fileSource(name);
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				for (size_t k = 0; k < N; ++k)
				{
					fileSource << source.Function(source.node[count], i * step, j * step, k * step) << " ";
				}
			}
		}

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				fileSource << source.Function(source.node[count], i * step, j * step, detectorLevel) << " ";
			}
		}
	}
	fileSource.close();
}

void LoadingSourceValues(const Sources & source, vector<complex<float>> & Source_R, vector<complex<float>> & Source_X, string name) noexcept
{
	ifstream fileSource(name);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				for (size_t k = 0; k < N; ++k)
				{
					fileSource >> Source_R[count * N_QUBE + i * N_SQUARED + j * N + k];
				}
			}
		}

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				fileSource >> Source_X[count * N_SQUARED + i * N + j];
			}
		}
	}
	fileSource.close();
}
