#include "stdafx.h"
#include "Detectors.h"

using namespace std;

void SetArrayOverlineA(std::vector<std::complex<float>> & overline_a)
{
	float dist;
	size_t coord;
	vector<float> index(N, 1.0f);
	for (size_t i = 1; i < N - 1; ++i)
	{
		index[i] = 1.0f;
	}
	index[0] = 0.5f;
	index[NUMBER_PARTITION_POINTS] = 0.5f;

	for (size_t j = 0; j < N; ++j)
	{
		for (size_t k = 0; k < N; ++k)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					for (size_t r = 0; r < N; ++r)
					{
						dist = sqrtf(static_cast<float>((j - p) * (j - p) + (k - q) * (k - q) + (detectorLevel - r) * (detectorLevel - r)));
						coord = N_FOURTH_DEGREE * j + N_QUBE * k + N_SQUARED * p + N * q + r;
						overline_a[coord] = -exp(dist * I * omega / c_0) * index[p] * index[q] * index[r] * pow(omega, 2) * pow(step, 2) * INV_FOUR_PI / dist; // TODO
					}
				}
			}
		}
	}
}

void WriteArrayOverlineA(std::vector<std::complex<float>>& overline_a, std::string name)
{
	ofstream f_overline_a(name);
	f_overline_a << fixed << setprecision(6);
	for (size_t i = 0; i < N_FIFTH_DEGREE; ++i)
	{
		f_overline_a << overline_a[i] << " ";
	}
	f_overline_a.close();
}

void LoadingArrayOverlineA(std::vector<std::complex<float>>& overline_a, std::string name)
{
	ifstream f_overline_a(name);
	for (size_t i = 0; i < N_FIFTH_DEGREE; ++i)
	{
		f_overline_a >> overline_a[i];
	}
	f_overline_a.close();
}
