#pragma once
#include "../directProblem/BasicDataProblem.h"
#include "../directProblem/matrix_utils.h"
#include <vector>
#include <complex>


void GetJacobian(const size_t numberSource,
	const std::vector<std::complex<float>> & a,
	const std::vector<std::complex<float>> & overline_a,
	const std::vector<std::complex<float>> & xi,
	const std::vector<std::vector<std::complex<float>>> & u,
	std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	std::vector<std::vector<std::vector<std::complex<float>>>> & F_even,
	std::vector<std::vector<std::complex<float>>> & F_0,
	std::vector<std::vector<std::complex<float>>> & F_00);
