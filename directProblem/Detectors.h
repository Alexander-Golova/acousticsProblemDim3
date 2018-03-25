#pragma once
#include "stdafx.h"
#include "BasicDataProblem.h"
#include "Sources.h"

void SetArrayOverlineA(std::vector<std::complex<float>> & overline_a);

void WriteArrayOverlineA(std::vector<std::complex<float>> & overline_a, std::string name);

void LoadingArrayOverlineA(std::vector<std::complex<float>> & overline_a, std::string name);

void LoadingOverlineU(const Sources & source, std::vector<std::complex<float>> & overline_u, std::string name);
