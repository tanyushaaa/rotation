#pragma once
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <utility>
///////////////////////////////  � � � � � � � � � � � �  //////////////////////////////

double det_(double **T, int N);
///////////////////////////////  � � � � �  � � � � � �  ///////////////////////////////

class Matrix {
public:
	Matrix();
	Matrix (std::istream& input, int n_);
	Matrix (const Matrix& m);
	Matrix& operator= (Matrix&& m);
	Matrix(std::pair<double**, int> m);
	bool main_diagonal() const;
	std::pair<double**, int> GetA() const;
	int GetN() const;
	double GetD() const;
	bool joint() const;
	~Matrix();
	friend std::ostream& operator<< (std::ostream &out, const Matrix& matrix);
private:
	int n; //���-�� �����������
	double** A = nullptr;
	double det;
};

Matrix operator*(const Matrix& m1, const Matrix& m2);
