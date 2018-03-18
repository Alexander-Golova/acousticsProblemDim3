#pragma once
#include "stdafx.h"
#include "BasicDataProblem.h"
#include "CSourses.h"

class CInhomogeneity
{
public:
	CInhomogeneity();
	~CInhomogeneity();

	void SetCoordinate(Point node);
	Point GetCoordinate() const;
	Point GetStep() const;
	void SetRefractionIndex();
	void WriteRefractionIndex(std::string name) const;
	void SetCoordinateSourses(const std::vector<Point>& coordinates);
	void WriteSourses(std::string name);
	void SetArrayA();
	void WriteArrayA(std::string name) const;

private:
	Point m_coordinate;
	Point m_step;
	std::array<std::array<std::array<float, N>, N>, N > m_refractionIndex;

	std::array<std::complex<float>, N_QUBE> m_fields;

	std::vector<Point> m_coordinateSourses;
	
	std::array<std::array<std::array<std::array<std::array<std::array<std::complex<float>, N>, N>, N>, N>, N>, N> m_a;

	std::array<std::array<std::complex<float>, N_QUBE>, N_QUBE> m_substantiveMatrix;
	std::array<float, N> m_index;
};
