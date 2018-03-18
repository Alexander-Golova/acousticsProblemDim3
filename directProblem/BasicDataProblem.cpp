#include "stdafx.h"
#include "BasicDataProblem.h"

using namespace std;

complex<float> G(const Point source, const Point point)
{
	float dist = sqrt(pow(point.x - source.x, 2) + pow(point.y - source.y, 2) + pow(point.z - source.z, 2));
	return -exp(I * omega * dist / c_0) * INV_FOUR_PI / dist;
}
