#pragma once

// функция Грина
std::complex<float> G(float x_1, float x_2, float x_3, float y_1, float y_2, float y_3);

//печать времени
void Lasting(const std::string & st, clock_t & time);