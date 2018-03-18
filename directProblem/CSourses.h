#pragma once
#include "BasicDataProblem.h"

class CSourses
{
public:
	CSourses();
	~CSourses();
	size_t GetNumber() const;

	void SetCoordinate(Point node);
	std::vector<Point> GetCoordinates() const;
	std::complex<float> Function(const Point source, const Point point) const;

private:
	int m_number = 0;
	std::vector<Point> m_sourceCoordinates;
};
