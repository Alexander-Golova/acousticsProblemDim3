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

void GetMatrixA(const size_t numberSource, const vector<vector<vector<complex<float>>>>& F_odd,
	const vector<vector<vector<complex<float>>>>& F_even, const vector<vector<complex<float>>>& F_0,
	const vector<vector<complex<float>>>& F_00, vector<vector<vector<complex<float>>>>& A, const float alpha)
{
	// A_0 задачи имеет номер A[numberSource]
	vector<vector<complex<float>>> auxiliaryMatrix(N_QUBE, vector<complex<float>>(N_QUBE));

	MultTransposedMatrix(F_odd[0], F_odd[0], A[numberSource]);
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_odd[count], auxiliaryMatrix);
		AddSquareMatrices(A[numberSource], auxiliaryMatrix);
	}
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_even[count], F_even[count], auxiliaryMatrix);
		AddSquareMatrices(A[numberSource], auxiliaryMatrix);
	}
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_0, A[count]);
		MultTransposedMatrix(F_even[count], F_00, auxiliaryMatrix);
		AddSquareMatrices(A[count], auxiliaryMatrix);
	}

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		A[numberSource][i][i] += alpha;
	}

}

void GetMatrixB(const vector<vector<complex<float>>>& F_0, const vector<vector<complex<float>>>& F_00,
	vector<vector<complex<float>>>& B, const float alpha)
{
	vector<vector<complex<float>>> auxiliaryMatrix(N_QUBE, vector<complex<float>>(N_QUBE));

	MultTransposedMatrix(F_0, F_0, B);
	MultTransposedMatrix(F_00, F_00, auxiliaryMatrix);
	AddSquareMatrices(B, auxiliaryMatrix);

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		B[i][i] += alpha;
	}
}

void GetOperatorF(const size_t numberSource, const vector<complex<float>>& a,
	const vector<complex<float>>& overline_a, const vector<complex<float>>& xi,
	const vector<vector<complex<float>>>& u, const vector<complex<float>>& overline_u,
	const vector<complex<float>>& Source_R, const vector<complex<float>>& Source_X,
	vector<vector<complex<float>>>& F_part_odd, vector<vector<complex<float>>>& F_part_even)
{
	size_t ii, jj, coord;

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				for (size_t k = 0; k <= N; ++k)
				{
					ii = i * N_SQUARED + j * N + k;
					for (size_t p = 0; p < N; ++p)
					{
						for (size_t q = 0; q < N; ++q)
						{
							for (size_t r = 0; r < N; ++r)
							{
								if ((i != p) || (q != j) || (r != k))
								{
									jj = p * N_SQUARED + q * N + r;
									coord = ii * N_QUBE + jj;
									F_part_odd[count][ii] = a[coord] * xi[jj] * u[count][jj]; // TODO
								}
							}
						}
					}
					F_part_odd[count][ii] -= Source_R[count * N_QUBE + ii]; // TODO
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
				ii = i * (N + 1) + j;
				F_part_even[count][ii] = overline_u[count * N_SQUARED + ii] - Source_X[count * N_SQUARED + ii]; // TODO
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						for (size_t r = 0; r < N; ++r)
						{
							jj = p * N_SQUARED + q * N + r;
							coord = ii * N_SQUARED + jj;
							F_part_even[count][ii] += overline_a[coord] * xi[jj] * u[count][jj]; // TODO
						}
					}
				}
			}
		}
	}
}
