#include "stdafx.h"
#include "Inhomogeneity.h"

using namespace std;

void SetRefractionIndex(vector<float> & xi) noexcept
{
	const float sigma = 16.0f;
	size_t coord;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				coord = N_SQUARED * i + N * j + k;
				xi[coord] = 0.4f * exp(-(pow(i * step - 0.6f, 2) + pow(j * step - 0.6f, 2) + pow(k * step - 0.6f, 2)) * sigma);
			}
		}
	}
}

void WriteRefractionIndex(vector<float>& xi, string name) noexcept
{
	ofstream file_xi(name);
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i < N_QUBE; ++i)
	{
		file_xi << xi[i] << " ";
	}
	file_xi.close();
}

void SetArrayA(vector<complex<float>> & a) noexcept
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

	for (size_t i = 0; i < N; ++i)
	{
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
							dist = sqrtf(static_cast<float>((i - p) * (i - p) + (j - q) * (j - q) + (k - r) * (k - r)));
							coord = N_FIFTH_DEGREE * i + N_FOURTH_DEGREE * j + N_QUBE * k + N_SQUARED * p + N * q + r;
							if (dist > 0.000001f) // DBL_EPSILON
							{								
								a[coord] = -exp(dist * I * omega / c_0) * index[p] * index[q] * index[r] * pow(omega, 2) * pow(step, 2) * INV_FOUR_PI / dist; // TODO
							}
							else
							{
								a[coord] = 0.0f;
							}
						}
					}
				}
			}
		}
	}
}

void WriteArrayA(vector<complex<float>> & a, string name) noexcept
{
	ofstream f_a(name);
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i < N_SIXTH_DEGREE; ++i)
	{
		f_a << a[i] << " ";
	}
	f_a.close();
}

void LoadingArrayA(std::vector<std::complex<float>>& a, std::string name) noexcept
{
	ifstream f_a(name);
	for (size_t i = 0; i < N_SIXTH_DEGREE; ++i)
	{
		f_a >> a[i];
	}
	f_a.close();
}
