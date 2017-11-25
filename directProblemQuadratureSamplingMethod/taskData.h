#pragma once

struct Point
{
	double x;
	double y;
	double z;
};

// ������� ������������� ����
const double OMEGA = 1.0;
const double C_0 = 1.0;

// ���������� ����������
const double receiver = 1.2;

// ���������� ����������� �� ������� ���������
const size_t NUMBER_PARTITION_POINT = 50;
const size_t N_SQUARED = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);

// ������ �������� � ������� ��������� ��������������
const double DOMAIN_IN_HOMOGENEITY = 1.0;

// ��� �� �����
const double h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINT;

