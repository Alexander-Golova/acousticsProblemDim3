#pragma once
#include "BasicDataProblem.h"

class CDetectors
{
public:
	CDetectors();
	~CDetectors();
	void SetLevel(float level);
	void SetStep(Point node);

	void SetFields();
	void WriteFields(std::string name);

	void SetArrayA();
	void WriteArrayA(std::string name) const;

private:
	std::vector<std::complex<float>> m_fields;

	//std::array<std::complex<float>, N_SQUARED> m_fields;
	//std::array<std::array<std::array<std::array<std::array<std::complex<float>, N>, N>, N>, N>, N> m_a;

	std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>>> m_a;

	float m_z;
	Point m_step;
	std::vector<float> m_index;
};
