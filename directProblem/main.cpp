#include "stdafx.h"
#include "taskData.h"
#include "Sources.h"
#include "matrix_utils.h"
#include "exact_solution.h"
#include "basicArrays.h"
#include "basicFunctions.h"
#include "directProblem_utils.h"

using namespace std;

int main()
{
	const Source source;

	vector<vector<vector<complex<double>>>> u(NUMBER_PARTITION_POINTS + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINTS + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINTS + 1, complex<double>())));

	vector<vector<vector<double>>> xi(NUMBER_PARTITION_POINTS + 1,
		vector<vector<double>>(NUMBER_PARTITION_POINTS + 1,
			vector<double>(NUMBER_PARTITION_POINTS + 1, 0.0)));

	GetExactSolution(xi);

	clock_t time = clock();
	clock_t timeBegin = clock();

	vector<vector<vector<vector<vector<vector<complex<double>>>>>>> a(NUMBER_PARTITION_POINTS + 1,
		vector<vector<vector<vector<vector<complex<double>>>>>>(NUMBER_PARTITION_POINTS + 1,
			vector<vector<vector<vector<complex<double>>>>>(NUMBER_PARTITION_POINTS + 1,
				vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINTS + 1,
					vector<vector<complex<double>>>(NUMBER_PARTITION_POINTS + 1,
						vector<complex<double>>(NUMBER_PARTITION_POINTS + 1, complex<double>()))))));

	vector<vector<vector<vector<vector<complex<double>>>>>> overline_a(NUMBER_PARTITION_POINTS + 1,
		vector<vector<vector<vector<complex<double>>>>>(NUMBER_PARTITION_POINTS + 1,
			vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINTS + 1,
				vector<vector<complex<double>>>(NUMBER_PARTITION_POINTS + 1,
					vector<complex<double>>(NUMBER_PARTITION_POINTS + 1, complex<double>())))));

	GetBasicArrays(a, overline_a);
	Lasting("Time calculation of basic matrices", time);

	WriteBasicArraysFile(a, overline_a);
	Lasting("Download time major arrays", time);

	WriteSourceValues(source);
	Lasting("The computation time of the source function", time);

	// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEquation[ii]
	vector<complex<double>> rightPartEquation(N_QUBE, complex<double>());
	vector<complex<double>> numbered_u(N_QUBE);
	vector<vector<complex<double>>> substantiveMatrix(N_QUBE,
		vector<complex<double>>(N_QUBE, complex<double>()));
	vector<vector<complex<double>>> overline_u(NUMBER_PARTITION_POINTS + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINTS + 1, complex<double>()));

	GetSubstantiveMatrix(a, xi, substantiveMatrix);
	Lasting("The computation time of the matrix inside the squared", time);

	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);

	for (size_t count = 0; count < source.numberSource; ++count)
	{
		GetRightPartEquation(source, count, rightPartEquation);
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
		InverseRenumbering(numbered_u, u);
		Lasting("Finding the acoustic pressure in R", time);

		GetOverlineU(source, count, overline_a, xi, u, overline_u);
		for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
			{
				file_overline_u << overline_u[i][j] << " ";
			}
		}
		Lasting("Finding the acoustic pressure in X", time);
	}
	file_overline_u.close();

	Lasting("The total time of the program", timeBegin);
	return 0;
}