﻿#include "stdafx.h"
#include "../directProblem/BasicDataProblem.h"
#include "../directProblem/Detectors.h"
#include "../directProblem/Sources.h"
#include "../directProblem/Inhomogeneity.h"
#include "../directProblem/matrix_utils.h"
#include "initialValue.h"
#include "regularization.h"
#include <boost/timer.hpp>

using namespace std;


int main()
{
	size_t numberOfIterations;
	cout << "Enter the number of iterations ";
	cin >> numberOfIterations;

	float alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	float multiplier;
	cout << "Enter q ";
	cin >> multiplier;

	boost::timer t;
	t.restart();
	boost::timer total;
	total.restart();

	const Sources source;

	// выделение памяти и загрузка данных прямой задачи
	vector<complex<float>> a(N_SIXTH_DEGREE);
	LoadingArrayA(a, string("matrix_a.txt"));

	vector<complex<float>> overline_a(N_FIFTH_DEGREE);
	LoadingArrayOverlineA(overline_a, string("matrix_overline_a.txt"));

	// выделение памяти для значений источников и их загрузка
	vector<complex<float>> Source_R(source.numberSource * N_QUBE);
	vector<complex<float>> Source_X(source.numberSource * N_SQUARED);
	LoadingSourceValues(source, Source_R, Source_X, string("Source.txt"));

	// выделение памяти для поля в приемниках и загрузка
	vector<complex<float>> overline_u(source.numberSource * N_SQUARED);
	LoadingOverlineU(source, overline_u, string("matrix_overline_u.txt"));

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<vector<complex<float>>>> F_odd(source.numberSource,
		vector<vector<complex<float>>>(N_QUBE, vector<complex<float>>(N_QUBE)));

	vector<vector<complex<float>>> F_0(N_QUBE, vector<complex<float>>(N_QUBE));

	vector<vector<complex<float>>> F_00(N_SQUARED, vector<complex<float>>(N_QUBE));

	vector<vector<vector<complex<float>>>> F_even(source.numberSource + 1,
		vector<vector<complex<float>>>(N_SQUARED, vector<complex<float>>(N_QUBE)));

	//выделение памяти под массивы A и B
	vector<vector<vector<complex<float>>>> A(source.numberSource + 1,
		vector<vector<complex<float>>>(N_QUBE, vector<complex<float>>(N_QUBE)));

	vector<vector<complex<float>>> B(N_QUBE, vector<complex<float>>(N_QUBE));

	vector<vector<complex<float>>> inverseMatrixB(N_QUBE, vector<complex<float>>(N_QUBE));

	// память для хранения значений основного оператора
	vector<vector<complex<float>>> F_part_odd(source.numberSource, vector<complex<float>>(N_QUBE));
	vector<vector<complex<float>>> F_part_even(source.numberSource, vector<complex<float>>(N_SQUARED));

	// память для b_0, b_1,...
	vector<vector<complex<float>>> b_right(source.numberSource + 1, vector<complex<float>>(N_QUBE));

	// память для u^(1), u^(2), u^(3)
	vector<vector<complex<float>>> u(source.numberSource, vector<complex<float>>(N_QUBE));

	// память для xi
	vector<complex<float>> xi(N_QUBE);

	// Начало вычислительной части
	InitialValueU(source.numberSource, u, Source_R);
	InitialValueXi(xi);

	double duration = t.elapsed();
	cout << "Time allocation and download " << duration << endl;
	t.restart();

	// начало основных итераций
	for (size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;

		//строим левую часть СЛАУ основного метода Ньютона
		// находим якобианы
		GetJacobian(source.numberSource, a, overline_a, xi, u, F_odd, F_even, F_0, F_00);
		duration = t.elapsed();
		cout << "The counting time of the Jacobian matrices " << duration << endl;
		t.restart();

		// находим матрицы А и В
		GetMatrixA(source.numberSource, F_odd, F_even, F_0, F_00, A, alpha);
		duration = t.elapsed();
		cout << "Counting time of matrices A " << duration << endl;
		t.restart();

		GetMatrixB(F_0, F_00, B, alpha);
		duration = t.elapsed();
		cout << "Counting time of matrices B " << duration << endl;
		t.restart();

		cout << "Calculate the left side" << endl;

		//строим правую часть СЛАУ основаного метода Ньютона
		// находим значения основного оператора F
		GetOperatorF(source.numberSource, a, overline_a, xi, u, overline_u, Source_R, Source_X, F_part_odd, F_part_even);
		duration = t.elapsed();
		cout << "Counting time of the main matrix " << duration << endl;
		t.restart();

		GetValueDerivedFunction(source.numberSource, xi, u, F_odd, F_even, F_0, F_00, F_part_odd, F_part_even);

		Getb(source.numberSource, F_odd, F_even, F_0, F_00, F_part_odd, F_part_even, b_right);
		duration = t.elapsed();
		cout << "Calculation time right side " << duration << endl;
		t.restart();

		InvertMatrix(B, inverseMatrixB);
		GetXi(source.numberSource, A, inverseMatrixB, b_right, xi);
		duration = t.elapsed();
		cout << "Time of xi " << duration << endl;
		t.restart();

		GetU(source.numberSource, A, inverseMatrixB, b_right, xi, u);
		duration = t.elapsed();
		cout << "Time of u " << duration << endl;
		t.restart();

		alpha = alpha * multiplier;

		ProjectionXi(xi);

		PrintXi(xi, iteration);

	}

	double durationTotal = total.elapsed();
	cout << "The total time of the program " << durationTotal << endl;

	return 0;
}
