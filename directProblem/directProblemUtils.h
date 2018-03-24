#pragma once
#include "stdafx.h"
#include "BasicDataProblem.h"
#include "Sources.h"

void GetSubstantiveMatrix(const std::vector<std::complex<float>> & a, const std::vector<float> & xi,
	std::vector<std::vector<std::complex<float>>> & substantiveMatrix);

void GetRightPartEquation(const Sources & source, const size_t count, std::vector<std::complex<float>> & rightPartEquation);

void GetOverlineU(const Sources & source, size_t count, const std::vector<std::complex<float>> & overline_a,
	const std::vector<float> & xi, const std::vector<std::complex<float>> & u, std::vector<std::complex<float>> & overline_u);
