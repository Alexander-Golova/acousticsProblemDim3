#pragma once
#include "stdafx.h"
#include "BasicDataProblem.h"


void SetRefractionIndex(std::vector<float> & xi) noexcept;

void WriteRefractionIndex(std::vector<float> & xi, std::string name) noexcept;

void SetArrayA(std::vector<std::complex<float>> & a) noexcept;

void WriteArrayA(std::vector<std::complex<float>> & a, std::string name) noexcept;

void LoadingArrayA(std::vector<std::complex<float>> & a, std::string name) noexcept;
