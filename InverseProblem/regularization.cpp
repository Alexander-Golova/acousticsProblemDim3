#include "stdafx.h"
#include "regularization.h"

using namespace std;

void GetJacobian(const size_t numberSource, const vector<complex<float>> & a, const vector<complex<float>> & overline_a,
	const vector<complex<float>> & xi, const vector<vector<complex<float>>> & u,
	vector<vector<vector<complex<float>>>> & F_odd, vector<vector<vector<complex<float>>>> & F_even,
	vector<vector<complex<float>>> & F_0, vector<vector<complex<float>>> & F_00) noexcept
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
			for (size_t i = 0; i < N; ++i)
			{
				for (size_t j = 0; j < N; ++j)
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
	const vector<vector<complex<float>>>& F_00, vector<vector<vector<complex<float>>>>& A, const float alpha) noexcept
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
	vector<vector<complex<float>>>& B, const float alpha) noexcept
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
	vector<vector<complex<float>>>& F_part_odd, vector<vector<complex<float>>>& F_part_even) noexcept
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
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				ii = i * N + j;
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

void GetValueDerivedFunction(const size_t numberSource, const vector<complex<float>>& xi, const vector<vector<complex<float>>>& u,
	const vector<vector<vector<complex<float>>>>& F_odd, const vector<vector<vector<complex<float>>>>& F_even,
	const vector<vector<complex<float>>>& F_0, const vector<vector<complex<float>>>& F_00,
	vector<vector<complex<float>>>& F_part_odd, vector<vector<complex<float>>>& F_part_even) noexcept
{
	vector<complex<float>> supportingVector(N_SQUARED);
	vector<complex<float>> supportingVectorSQ(N_QUBE);

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrixVector(F_odd[count], xi, supportingVectorSQ);
		for (size_t i = 0; i < N_QUBE; ++i)
		{
			F_part_odd[count][i] = supportingVectorSQ[i] - F_part_odd[count][i];
		}

		MultMatrixVector(F_0, u[count], supportingVectorSQ);
		for (size_t i = 0; i < N_QUBE; ++i)
		{
			F_part_odd[count][i] += supportingVectorSQ[i];
		}

		MultMatrixVector(F_even[count], xi, supportingVector);
		for (size_t i = 0; i < N_SQUARED; ++i)
		{
			F_part_even[count][i] = supportingVector[i] - F_part_even[count][i];
		}
		MultMatrixVector(F_00, u[count], supportingVector);
		for (size_t i = 0; i < N_SQUARED; ++i)
		{
			F_part_even[count][i] += supportingVector[i];
		}
	}
}

// b_0 задачи находится в b[numberSource]
void Getb(const size_t numberSource,
	const vector<vector<vector<complex<float>>>>& F_odd, const vector<vector<vector<complex<float>>>>& F_even,
	const vector<vector<complex<float>>>& F_0, const vector<vector<complex<float>>>& F_00,
	const vector<vector<complex<float>>>& F_part_odd, const vector<vector<complex<float>>>& F_part_even,
	vector<vector<complex<float>>> & b_right) noexcept
{
	vector<complex<float>> supportingVectorSQ(N_QUBE);

	MultTransposedMatrixVector(F_odd[0], F_part_odd[0], b_right[numberSource]);
	MultTransposedMatrixVector(F_even[0], F_part_even[0], supportingVectorSQ);
	AddVectors(b_right[numberSource], supportingVectorSQ);
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(F_odd[count], F_part_odd[count], supportingVectorSQ);
		AddVectors(b_right[numberSource], supportingVectorSQ);
		MultTransposedMatrixVector(F_even[count], F_part_even[count], supportingVectorSQ);
		AddVectors(b_right[numberSource], supportingVectorSQ);
	}

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(F_0, F_part_odd[count], b_right[count]);
		MultTransposedMatrixVector(F_00, F_part_even[count], supportingVectorSQ);
		AddVectors(b_right[count], supportingVectorSQ);
	}
}

void GetXi(const size_t numberSource,
	vector<vector<vector<complex<float>>>> & A,
	const vector<vector<complex<float>>> & inverseMatrixB,
	vector<vector<complex<float>>> & b_right,
	vector<complex<float>>& xi) noexcept
{
	vector<vector<complex<float>>> auxiliaryMatrix(N_QUBE, vector<complex<float>>(N_QUBE));
	vector<vector<complex<float>>> secondAuxiliaryMatrix(N_QUBE, vector<complex<float>>(N_QUBE));
	vector<complex<float>> supportingVectorSQ(N_QUBE);
	vector<complex<float>> secondSupportingVectorSQ(N_QUBE);

	//для левой части уравнения с xi все складываем в A_00 -> A[numberSource]
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrix(A[count], inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A[count], secondAuxiliaryMatrix);
		SubSquareMatrices(A[numberSource], secondAuxiliaryMatrix);
	}

	//для правой части уравнения с xi все складываем в b0 - b[numberSource]
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrixVector(inverseMatrixB, b_right[count], supportingVectorSQ);
		MultMatrixVector(A[count], supportingVectorSQ, secondSupportingVectorSQ);
		SubVectors(b_right[numberSource], secondSupportingVectorSQ);
	}

	// находим xi
	SolveSlauGaussa(A[numberSource], b_right[numberSource], xi);
}

void GetU(const size_t numberSource,
	const vector<vector<vector<complex<float>>>> & A,
	const vector<vector<complex<float>>> & inverseMatrixB,
	vector<vector<complex<float>>> & b_right,
	const vector<complex<float>>& xi,
	vector<vector<complex<float>>>& u) noexcept
{
	vector<complex<float>> supportingVectorSQ(N_QUBE, complex<float>());

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(A[count], xi, supportingVectorSQ);
		SubVectors(b_right[count], supportingVectorSQ);
		MultMatrixVector(inverseMatrixB, b_right[count], u[count]);
	}
}

void ProjectionXi(vector<complex<float>> & xi) noexcept
{
	for (size_t i = 0; i < N_QUBE; ++i)
	{
		xi[i] = real(xi[i]);
		if (real(xi[i]) <= 0.0f)
		{
			xi[i] = { 0.0f, 0.0f };
		}
	}
}

void PrintXi(const vector<complex<float>> & xi, size_t iteration) noexcept
{
	ofstream f_xi("approximate_xi_" + to_string(iteration + 1) + ".txt");
	f_xi << fixed << setprecision(6);
	for (size_t i = 0; i < N_QUBE; ++i)
	{
		f_xi << real(xi[i]) << " ";
	}
	f_xi.close();
}
