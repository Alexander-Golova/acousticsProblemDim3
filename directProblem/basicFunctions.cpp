#include "stdafx.h"
#include "taskData.h"
#include "basicFunctions.h"

using namespace std;

std::complex<double> G(double x_1, double x_2, double x_3, double y_1, double y_2, double y_3)
{
	double dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2) + pow(x_3 - y_3, 2));
	return -OMEGA * OMEGA * exp((0.0, 1.0) * OMEGA * dist / C_0) / (FOUR_PI * dist);
}

void Lasting(const string & st, clock_t & time)
{
	clock_t timeFinish = clock();
	double d = static_cast<double>(timeFinish - time) / CLOCKS_PER_SEC;
	cout << st << " " << d << endl;
	time = clock();
}