#include "stdafx.h"
#include "CSourses.h"


using namespace std;

CSourses::CSourses()
{
}

CSourses::~CSourses()
{
}

size_t CSourses::GetNumber() const
{
	return m_number;
}

void CSourses::SetCoordinate(Point node)
{
	m_sourceCoordinates.push_back(node);
	++m_number;
}

std::vector<Point> CSourses::GetCoordinates() const
{
	return m_sourceCoordinates;
}

std::complex<float> CSourses::Function(const Point source, const Point point) const
{
	return G(source, point);
}
