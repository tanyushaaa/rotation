#pragma once
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <utility>
///////////////////////////////  О п р е д е л и т е л ь  //////////////////////////////

double det_(double **T, int N);
///////////////////////////////  К л а с с  м а т р и ц  ///////////////////////////////

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
	int n; //кол-во неизвестных
	double** A = nullptr;
	double det;
};

Matrix operator*(const Matrix& m1, const Matrix& m2);
