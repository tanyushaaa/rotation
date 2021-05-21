#include "matrix.h"

///////////////////////////////  О п р е д е л и т е л ь  //////////////////////////////

double det_(double **T, int N) {
	switch (N){
	case 1:
		return T[0][0];
	case 2:
		return T[0][0] * T[1][1] - T[0][1] * T[1][0];
	default:
		if (N < 1) {
			throw std::invalid_argument("There is no elements in this matrix.");
		}
		double **minor = new double*[N-1];		//массив указателей на столбцы исходной матрицы
		double det = 0;							//исключая заданные стр. и стл. (минор)
		int sign = 1;							//знак минора
		for (int i = 0; i < N; i++) {			//разложение по первому столбцу
			int sub_j = 0;
			for (int j = 0; j < N; j++) {		//заполнение "минорной" матрицы ссылками на исходные столбцы
				if (i != j)						//исключить i строку
					minor[sub_j++] = T[j] + 1;	//исключить 1й(0й) столбец
			}
			det += sign * T[i][0] * det_(minor, N-1);
			sign = -sign;
		}
		delete[] minor;
		return det;
	}
}

///////////////////////////////  К л а с с  м а т р и ц  ///////////////////////////////

	Matrix::Matrix() : n(3) {
		A = new double* [n];
		for (int i = 0; i < n; i++) {
			A[i] = new double [n];
		}

		A[0][0] = 68;
		A[0][1] = 4;
		A[0][2] = 3;

		A[1][0] = 4;
		A[1][1] = 32;
		A[1][2] = 4;

		A[2][0] = 3;
		A[2][1] = 4;
		A[2][2] = 84;

		det = det_(A, n);
	}

	Matrix::Matrix (std::istream& input, int n_) : n(n_) {
		A = new double* [n];
		for (int i = 0; i < n; i++) {
			A[i] = new double [n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				input >> A[i][j];
			}
		}

		det = det_(A, n);
	}

	Matrix::Matrix (const Matrix& m)
	: n(m.n),
	  det(m.det)
	{
		A = new double* [n];
		for(int i = 0; i < n; i++){
			A[i] = new double [n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = m.A[i][j];
			}
		}
	}

	Matrix& Matrix::operator= (Matrix&& m) {
		if (&m == this) {
			return *this;
		}

		for(int i = 0; i < n; i++){
			delete[] A[i];
		}

		n = m.n;
		det = m.det;

		for (int i = 0; i < n; i++) {
			A[i] = m.A[i];
		}

		for(int i = 0; i < n; i++){
			m.A[i] = nullptr;
		}

		return *this;
	}

	Matrix::Matrix(std::pair<double**, int> m){
		n = m.second;
		A = new double* [n];
		for(int i = 0; i < n; i++){
			A[i] = new double [n];
		}
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				A[i][j] = m.first[i][j];
			}
		}
		det = det_(A, n);
	}

	bool Matrix::main_diagonal() const { //проверка на то, являются ли элементы главной диагонали нулевыми (можно ли на них делить)
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				if (i == j){
					if (A[i][j] == 0)
						return false;
				}
			}
		}
		return true;
	}

	std::pair<double**, int> Matrix::GetA() const {
		return std::make_pair(A, n);
	}

	int Matrix::GetN() const { return n; }

	double Matrix::GetD() const { return det; }

	bool Matrix::joint() const { return (det != 0); }

	Matrix::~Matrix(){
		for (int i = 0; i < n; i++) {
			delete [] A[i];
		}
	}

std::ostream& operator<< (std::ostream &out, const Matrix& matrix){
	bool first = true;
	for (int i = 0; i < matrix.n; i++){
		if (!first)
			out << std::endl;

		for (int j = 0; j < matrix.n; j++){
			out << std::setw(10);
			out << std::fixed << std::setprecision(6);
			out << matrix.A[i][j] << " ";
		}

		first = false;
	}
	return out;
}

Matrix operator*(const Matrix& m1, const Matrix& m2){
	int n = m1.GetN();
	std::stringstream stream;

	double buff;
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++){
			buff = 0;
			for (int j = 0; j < n; j++) {
				buff += m1.GetA().first[i][j] * m2.GetA().first[j][k];
			}
			stream << buff << " ";
		}
	}
	Matrix m(stream, n);
	return m;
}
