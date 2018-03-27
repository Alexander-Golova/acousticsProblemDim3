#include "stdafx.h"
#include "BasicDataProblem.h"
#include "Detectors.h"
#include "Sources.h"
#include "Inhomogeneity.h"
#include "directProblemUtils.h"
#include "matrix_utils.h"
#include <boost/timer.hpp>

using namespace std;


int main()
{
	boost::timer t;
	t.restart();

	// задаём источники
	Sources source;
	
	// память для индекса рефракции xi
	vector<float> xi(N_QUBE);

	// определем индекс рефракции на неоднородности и печатаем её значение в файл
	SetRefractionIndex(xi);
	WriteRefractionIndex(xi, string("refraction_index.txt"));

	// на неоднородности
	// объявляем массив a
	vector<complex<float>> a(N_SIXTH_DEGREE);
	// определяем массив a 
	SetArrayA(a);
	// пишем в файл массив А
	WriteArrayA(a, string("matrix_a.txt"));

	// на детекторе
	// объявляем массив overline_a
	vector<complex<float>> overline_a(N_FIFTH_DEGREE);
	// определяем массив overline_a на детекторе
	SetArrayOverlineA(overline_a);
	// пишем в файл массив overline_a
	WriteArrayOverlineA(overline_a, string("matrix_overline_a.txt"));

	// печатаем значение акустического поля на неоднородности и на детекотре в файл
	WriteSourceValues(source, string("Source.txt"));

	// Для получения полей на детекторе u получаем СЛАУ
	// выделяем память для основной матрицы, правой части и вектора неизвестных
	vector<vector<complex<float>>> substantiveMatrix(N_QUBE, vector<complex<float>>(N_QUBE));
	vector<complex<float>> rightPartEquation(N_QUBE);
	// память для поля в детекторе. Значение этого поля передаём в прямую задачу
	//vector<complex<float>> DetectorsFields(N_QUBE);
	
	vector<complex<float>> u(N_QUBE);
	GetSubstantiveMatrix(a, xi, substantiveMatrix);

	// память для значения поля в детекторе N_SQUARED
	vector<complex<float>> overline_u(N_SQUARED);

	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);

	// пробегаемся по всем источникам
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
		// substantiveMatrix[ii][jj] * u[jj] = rightPartEquation[ii]
		GetRightPartEquation(source, count, rightPartEquation);
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, u);
		GetOverlineU(source, count, overline_a, xi, u, overline_u);
		for (size_t i = 0; i < N_SQUARED; ++i)
		{
			file_overline_u << overline_u[i] << " ";
		}
	}
	file_overline_u.close();

	double duration = t.elapsed();
	cout << "The total time of the program " << duration << endl;

	return 0;
}