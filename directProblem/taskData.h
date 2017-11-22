#pragma once

struct Point
{
	double x;
	double y;
	double z;
};

const double PI = 3.14159265359;

// задание характеристик поля
const double OMEGA = 1.0;
const double C_0 = 1.0;

// координаты приемников: z = receiver
const double receiver = -0.1;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINTS = 10;
const size_t N_SQUARED = (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1);
const size_t N_QUBE = (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1);

// размер куба в котором находится неоднородность
const double DOMAIN_IN_HOMOGENEITY = 1.0;

// шаг по сетке
const double h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINTS;
