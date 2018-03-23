#include "stdafx.h"
#include "CInhomogeneity.h"
#include "CSourses.h"

using namespace std;

CInhomogeneity::CInhomogeneity()
{
	vector<vector<vector<float>>> xi(N, vector<vector<float>>(N, vector<float>(N, 0.0f)));
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


CInhomogeneity::~CInhomogeneity()
{
}

void CInhomogeneity::SetCoordinate(Point node)
{
	m_coordinate = node;
	m_step = { node.x * inverseN, node.y * inverseN, node.z * inverseN };
}

Point CInhomogeneity::GetCoordinate() const
{
	return m_coordinate;
}

Point CInhomogeneity::GetStep() const
{
	return m_step;
}

void CInhomogeneity::SetRefractionIndex()
{
	float sigma = 16.0f;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				m_refractionIndex[i][j][k] = 0.4f * exp(-((i * m_step.x - 0.6f) * (i * m_step.x - 0.6f) +
					(j * m_step.y - 0.6f) * (j * m_step.y - 0.6f) + (k * m_step.z - 0.6f) * (k * m_step.z - 0.6f)) * sigma);
			}
		}
	}
}

void CInhomogeneity::WriteRefractionIndex(string name) const
{
	ofstream file_xi(name);
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				file_xi << m_refractionIndex[i][j][k] << " ";
			}
		}
	}
	file_xi.close();
}

void CInhomogeneity::SetCoordinateSourses(const vector<Point>& coordinates)
{
	copy(begin(coordinates), end(coordinates), begin(m_coordinateSourses));
}

void CInhomogeneity::WriteSourses(string name)
{
	ofstream fileSource(name);
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < static_cast<size_t>(m_coordinateSourses.size()); ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				for (size_t k = 0; k < N; ++k)
				{
					fileSource << G(m_coordinateSourses[count], { i * m_step.x, j * m_step.y, k * m_step.z }) << " ";
				}
			}
		}
	}
	fileSource.close();
}

void CInhomogeneity::SetArrayA()
{
	float dist;
	complex<float> temp;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						for (size_t r = 0; r < N; ++r)
						{
							if ((i != p) || (q != j) || (r != k))
							{
								dist = sqrtf((i - p) * (i - p) * m_step.x * m_step.x + (j - q) * (j - q) * m_step.y * m_step.y + (k - r) * (k - r) * m_step.z * m_step.z);
								temp = exp(-dist * I * omega / c_0) / dist;
								temp *= m_index[p] * m_index[q] * m_index[r];
								temp *= omega * omega * m_step.x * m_step.y * m_step.z;
								temp *= INV_FOUR_PI;
								//m_a_re[i][j][k][p][q][r] = temp.real();
								//m_a_im[i][j][k][p][q][r] = temp.imag();
							}
						}
					}
				}
			}
		}
	}
}

void CInhomogeneity::WriteArrayA(std::string name) const
{
	ofstream f_a(name);
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t k = 0; k < N; ++k)
			{
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						for (size_t r = 0; r < N; ++r)
						{
							//f_a << m_a_re[i][j][k][p][q][r] << " ";
							//f_a << m_a_im[i][j][k][p][q][r] << " ";
						}
					}
				}
			}
		}
	}
	f_a.close();
}
