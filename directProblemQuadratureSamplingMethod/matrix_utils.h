#pragma once

// �������� �������
void GetNullMatrix(std::vector<std::vector<std::complex<double>>> & matrix);

// ������� ������� ������� ������
void SolveSlauGaussa(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & rhs, std::vector<std::complex<double>> & exactSolution);

//���������� �������� ������� ������� �������-������
void InvertMatrix(std::vector<std::vector<std::complex<double>>> matrix,
	std::vector<std::vector<std::complex<double>>> & invertedMatrix);

// �������� ���� ���������� ������ - ��������� ������������ � ������ �������
void AddSquareMatrices(std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs);

// ��������� ���� ���������� ������ - ��������� ������������ � ������ �������
void SubSquareMatrices(std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs);

// �������� ���� �������� - ��������� ������������ � ������ ������
void AddVectors(std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// ��������� ���� �������� - ��������� ������������ � ������ ������
void SubVectors(std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// ��������� ������� �� ������
std::complex<double> MultVectorVector(const std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// ��������� ������� �� ������
void MultMatrixVector(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & vect, std::vector<std::complex<double>> & result);

// ��������� ����������������� ������� �� ������
void MultTransposedMatrixVector(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & vect, std::vector<std::complex<double>> & result);

// ��������� ������
void MultMatrix(const std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs,
	std::vector<std::vector<std::complex<double>>> & result);

// ��������� ����������������� ������� �� ������� �������
void MultTransposedMatrix(const std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs,
	std::vector<std::vector<std::complex<double>>> & result);

// ��������� ������� ������� �� ����������������� �������
void MultMatrixTransposed(const std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs,
	std::vector<std::vector<std::complex<double>>> & result);
