#pragma once
#include "BasicDataProblem.h"
#include <vector>
#include <complex>

// задаём структуру источников
struct Sources
{
	// количество источников
	const size_t numberSource = 4;
	// координаты источников
	const std::vector<Point> node = { { 0.0f, 0.0f, -0.1f }, { 0.0f, 1.0f, -0.1f }, { 1.0f, 0.0f, -0.1f }, { 1.0f, 1.0f, -0.1f } };
	std::complex<float> Function(const Point source, const float x, const float y, const float z) const;
};


// печать значений источника 
void WriteSourceValues(const Sources & source, std::string name);

// загрузка
void LoadingSourceValues(const Sources & source, std::vector<std::complex<float>> & Source_R, std::vector<std::complex<float>> & Source_X, std::string name);
