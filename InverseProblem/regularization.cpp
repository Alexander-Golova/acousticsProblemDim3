#include "stdafx.h"
#include "regularization.h"

using namespace std;

void GetJacobian(const size_t numberSource, const vector<complex<float>> & a, const vector<complex<float>> & overline_a,
	const vector<complex<float>> & xi, const vector<vector<complex<float>>> & u,
	vector<vector<vector<complex<float>>>> & F_odd, vector<vector<vector<complex<float>>>> & F_even,
	vector<vector<complex<float>>> & F_0, vector<vector<complex<float>>> & F_00)
{
	size_t ii, jj, coord;

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				for (size_t k = 0; k < N; ++k)
				{
					ii = i * N_SQUARED + j * N + k;
					for (size_t p = 0; p < N; ++p)
					{
						for (size_t q = 0; q < N; ++q)
						{
							for (size_t r = 0; r < N; ++r)
							{
								jj = p * N_SQUARED + q * N + r;
								coord = ii * N_QUBE + jj;
								F_odd[count][ii][jj] = a[coord] * u[count][jj]; // TODO
								F_0[ii][jj] = a[coord] * xi[jj]; // TODO
							}
						}
					}
					F_0[ii][ii] += 1.0f;
				}
			}
		}
	}

		for (size_t count = 0; count < numberSource; ++count)
		{
			for (size_t i = 0; i <= N; ++i)
			{
				for (size_t j = 0; j <= N; ++j)
				{
					ii = i * N + j;
					for (size_t p = 0; p < N; ++p)
					{
						for (size_t q = 0; q < N; ++q)
						{
							for (size_t r = 0; r < N; ++r)
							{
								jj = p * N_SQUARED + q * N + r;
								coord = ii * N_QUBE + jj;
								F_even[count][ii][jj] = overline_a[coord] * u[count][jj]; // TODO
								F_00[ii][jj] = overline_a[coord] * xi[coord]; // TODO
							}
						}
					}
				}
			}
		}
}
