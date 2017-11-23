#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

void GetBasicArrays(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>>> & a,
	std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>> & overline_a);

void WriteBasicArraysFile(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>>> & a,
	std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>> & overline_a);
