#pragma once

// функция Грина
std::complex<double> G(double x_1, double x_2, double x_3, double y_1, double y_2, double y_3);

//печать времени
void Lasting(const std::string & st, clock_t & time);