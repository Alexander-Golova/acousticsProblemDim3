#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// источники

struct Source
{
	// количество источников
	const size_t numberSource = 3;
	// координаты источников
	const std::vector<Point> node = {
		{ 0.0, 0.0, 1.1 }, { 1.0, 0.0, 1.1 }, { 0.0, 1.0, 1.1 } };
	std::complex<double> Function(const Point source, const double x, const double y, const double z) const;
};

// печать значений источника в файл "Source.txt"
void WriteSourceValues(const Source & source);
