#pragma once

struct Point
{
	double x;
	double y;
	double z;
};

const double PI = 3.14159265359;

// ������� ������������� ����
const double OMEGA = 1.0;
const double C_0 = 1.0;

// ���������� ����������: z = receiver
const double receiver = -0.1;

// ���������� ����������� �� ������� ���������
const size_t NUMBER_PARTITION_POINTS = 10;
const size_t N_SQUARED = (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1);
const size_t N_QUBE = (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1) * (NUMBER_PARTITION_POINTS + 1);

// ������ ���� � ������� ��������� ��������������
const double DOMAIN_IN_HOMOGENEITY = 1.0;

// ��� �� �����
const double h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINTS;
