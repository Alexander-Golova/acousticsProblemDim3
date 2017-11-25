#pragma once

// функция Грина
std::complex<double> G(const double x_1, const double x_2, const double x_3, const double y_1, const double y_2, const double y_3);

//печать времени
void Lasting(const std::string & st, clock_t & time);