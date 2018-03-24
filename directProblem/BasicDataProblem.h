#pragma once
#include "stdafx.h"

const std::complex<float> I = { 0.0f, 1.0f };

const size_t NUMBER_PARTITION_POINTS = 20;

const float inverseN = 1.0f / NUMBER_PARTITION_POINTS;

//шаг
const float step = inverseN;

const size_t N = NUMBER_PARTITION_POINTS + 1;

const size_t N_SQUARED = N * N;

const size_t N_QUBE = N * N * N;

const size_t N_FOURTH_DEGREE = N * N * N * N ;

const size_t N_FIFTH_DEGREE = N * N * N * N * N;

const size_t N_SIXTH_DEGREE = N * N * N * N * N * N;

const float c_0 = 1.0f; // TODO

const float omega = 100.0f; // TODO

const float FOUR_PI = 12.566370614359f; // TODO

const float INV_FOUR_PI = 0.079577471545f;

// уровень приёмника
const float detectorLevel = 1.1f;

struct Point
{
	float x;
	float y;
	float z;
};

std::complex<float> G(const Point source, const Point point);
