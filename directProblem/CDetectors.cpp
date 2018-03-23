#include "stdafx.h"
#include "CDetectors.h"

using namespace std;

CDetectors::CDetectors()
{
	/*
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		m_fields[i] = { 0.0f, 0.0f };
		for (size_t j = 0; j < N_QUBE; ++j)
		{
			m_a[i][j] = { 0.0f, 0.0f };
		}
	}
	*/
	// находим индексы метода квадратур
	for (size_t i = 1; i < NUMBER_PARTITION_POINTS; ++i)
	{
		if (i % 2 != 0)
		{
			m_index[i] = 4.0f / 3;
		}
		else
		{
			m_index[i] = 2.0f / 3;
		}
	}
	m_index[0] = 1.0f / 3;
	m_index[NUMBER_PARTITION_POINTS] = 1.0f / 3;
}


CDetectors::~CDetectors()
{
}

void CDetectors::SetLevel(float level)
{
	m_z = level;
}

void CDetectors::SetStep(Point node)
{
	m_step = node;
}


void CDetectors::SetFields()
{
}

void CDetectors::WriteFields(string name)
{
	ofstream file_fields_detectors(name);
	file_fields_detectors << fixed << setprecision(6);
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		file_fields_detectors << m_fields[i] << " ";
	}
	file_fields_detectors.close();
}

void CDetectors::SetArrayA()
{
	float dist;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINTS; ++p) // TODO
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINTS; ++q) // TODO
				{
					for (size_t r = 0; r < NUMBER_PARTITION_POINTS; ++r) // TODO
					{
						dist = sqrtf((i - p) * (i - p) * m_step.x * m_step.x + (j - q) * (j - q) * m_step.y * m_step.y + (m_z - r) * (m_z - r)* m_step.z * m_step.z); // TODO про шаг
						m_a[i][j][p][q][r] = m_index[p] * m_index[q] * m_index[r];
						m_a[i][j][p][q][r] *= exp(-dist * I * omega / c_0) / dist;
						m_a[i][j][p][q][r] *= omega * omega * m_step.x * m_step.y * m_step.z;
						m_a[i][j][p][q][r] *= INV_FOUR_PI;
					}
				}
			}
		}
	}
}

void CDetectors::WriteArrayA(std::string name) const
{
	ofstream f_a(name);
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					for (size_t r = 0; r < N; ++r)
					{
						f_a << m_a[i][j][p][q][r] << " ";
					}
				}
			}
		}
	}
	f_a.close();
}
