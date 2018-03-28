#include "stdafx.h"
#include "initialValue.h"
#include "../directProblem/Sources.h"

using namespace std;

void InitialValueU(const size_t numberSource, vector<vector<complex<float>>> & u, 	vector<complex<float>> & Source_R) noexcept
{
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				for (size_t k = 0; k < N; ++k)
				{
					u[count][i * N_SQUARED + j * N + k] = Source_R[count * N_QUBE + i * N_SQUARED + j * N + k];
				}
			}
		}
	}
}

void InitialValueXi(vector<complex<float>> & xi) noexcept
{
	for (size_t i = 0; i < N_QUBE; ++i)
	{
		xi[i] = static_cast<complex<float>>(0.0f);
	}
}
