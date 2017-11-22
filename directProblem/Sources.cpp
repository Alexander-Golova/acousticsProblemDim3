#include "stdafx.h"
#include "Sources.h"

using namespace std;

complex<double> Source::Function(const Point source, const double x, const double y, const double z) const
{
	float dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2) + pow(z - source.z, 2));
	return -exp(( 0.0, 1.0 ) * OMEGA * dist / C_0) / (4 * PI * dist);
}

void WriteSourceValues(const Source & source)
{
	ofstream fileSource("Source.txt");
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
			{
				for (size_t k = 0; k <= NUMBER_PARTITION_POINTS; ++k)
				{
					fileSource << source.Function(source.node[count], i * h, j * h, k * h) << " ";
				}
			}
		}

		for (size_t i = 0; i <= NUMBER_PARTITION_POINTS; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINTS; ++j)
			{
				fileSource << source.Function(source.node[count], i * h, j * h, receiver) << " ";
			}
		}
	}
	fileSource.close();
}
