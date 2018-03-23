#include "stdafx.h"
#include "BasicDataProblem.h"
//#include "CDetectors.h"
#include "Sources.h"
//#include "CInhomogeneity.h"
#include "Inhomogeneity.h"

using namespace std;


int main()
{
	// задаём источники
	Sources source;
	
	// память для индекса рефракции xi
	vector<float> xi(N_QUBE);

	// определем индекс рефракции на неоднородности и печатаем её значение в файл
	SetRefractionIndex(xi);
	WriteRefractionIndex(xi, string("refraction_index.txt"));

	// объявляем массив a
	vector<complex<float>> a(N_SIXTH_DEGREE);
	//определяем массив a на неоднородности
	SetArrayA(a);
	// пишем в файл массив А
	WriteArrayA(a, string("matrix_a.txt"));

	// объявляем массив overline_a
	vector<complex<float>> overline_a(N_FIFTH_DEGREE);
	//определяем массив overline_a на детекторе

	// пишем в файл массив overline_a

	/*
	// задаём детекторы
	CDetectors detectors;

	// задаём уровень детектора
	detectors.SetLevel(1.1f);

	// задаём шаг для интегрирования
	detectors.SetStep(inhomogeneity.GetStep());

	// определяем массив overline_А на детекторе
	detectors.SetArrayA();

	// печатаем массив overline_А
	detectors.WriteArrayA(string("matrix_overline_a"));

	// пробегаемся по всем источникам
	for (size_t count = 0; count < source.GetNumber(); ++count)
	{
		// 
	}

	// получаем поле в детекторе
//	detectors.SetFields();



	/*

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
	*/
	return 0;
}