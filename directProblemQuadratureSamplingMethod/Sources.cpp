#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<double> Source::Function(const Point source, const double x, const double y, const double z) const
{
	double dist = sqrtf(pow(x - source.x, 2) + pow(y - source.y, 2) + pow(z - source.z, 2));
	return -exp((0.0, 1.0) * OMEGA * dist / C_0) / (4 * PI * dist);
}

void WriteSourceValues(const Source & source)
{
	ofstream fileSource("Source.txt");
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
			{
				for (size_t k = 0; k <= NUMBER_PARTITION_POINT; ++k)
				{
					fileSource << source.Function(source.node[count], i * h, j * h, k * h) << " ";
				}
			}
		}
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t k = 0; k <= NUMBER_PARTITION_POINT; ++k)
			{
				fileSource << source.Function(source.node[count], receiver, j * h, k * h) << " ";
			}
		}
	}
	fileSource.close();
}
