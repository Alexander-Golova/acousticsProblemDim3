#include "stdafx.h"
#include "../directProblem/basicFunctions.h"
#include "../directProblem/BasicDataProblem.h"
#include "../directProblem/Detectors.h"
#include "../directProblem/Sources.h"
#include "../directProblem/Inhomogeneity.h"
#include "../directProblem/matrix_utils.h"
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



    return 0;
}

