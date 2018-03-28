#pragma once
#include "../directProblem/BasicDataProblem.h"
#include "../directProblem/Sources.h"
#include <vector>
#include <complex>

void InitialValueU(const size_t numberSource, std::vector<std::vector<std::complex<float>>> & u, std::vector<std::complex<float>> & Source_R) noexcept;

void InitialValueXi(std::vector<std::complex<float>> & xi) noexcept;

