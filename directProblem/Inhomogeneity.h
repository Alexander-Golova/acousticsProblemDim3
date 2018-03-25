#pragma once
#include "stdafx.h"
#include "BasicDataProblem.h"


void SetRefractionIndex(std::vector<float> & xi);

void WriteRefractionIndex(std::vector<float> & xi, std::string name);

void SetArrayA(std::vector<std::complex<float>> & a);

void WriteArrayA(std::vector<std::complex<float>> & a, std::string name);

void LoadingArrayA(std::vector<std::complex<float>> & a, std::string name);
