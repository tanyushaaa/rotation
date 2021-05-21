#include "rotate.h"

double epsilon = 0.01;
Matrix E;

void E_mtrx(const Matrix& m){
	double** e = new double* [m.GetN()];
	for(int i = 0; i < m.GetN(); i++){
		e[i] = new double [m.GetN()];
	}

	for (int i = 0; i < m.GetN(); i++){
		for (int j = 0; j < m.GetN(); j++){
			if (i == j){
				e[i][j] = 1;
			} else {
				e[i][j] = 0;
			}
		}
	}
	E = Matrix(std::make_pair(e, m.GetN()));
}

std::pair<int, int> max_elem_up_diag(const Matrix& m){
	double **A_copy = new double* [m.GetN()];
	for (int i = 0; i < m.GetN(); i++) {
		A_copy[i] = new double [m.GetN()];
	}

	double max = 0;
	std::pair<int, int> max_indexes;
	for(int i = 0; i < m.GetN(); i++) {
		for(int j = 1; j < m.GetN(); j++) {
			if (j > i) {
				if (fabs(m.GetA().first[i][j]) >= max){
					max = fabs(m.GetA().first[i][j]);
					max_indexes = std::make_pair(i, j);
				}
			}
		}
	}
	return max_indexes;
}

double phi(std::pair<int, int> ind, const Matrix& m){
	return (0.5 * atan(2 * m.GetA().first[ind.first][ind.second] /
					  (m.GetA().first[ind.first][ind.first] - m.GetA().first[ind.second][ind.second])));
}

void transition_matrix(std::pair<int, int> ind, double phi, std::pair<double**, int> unit){
	for (int i = 0; i < unit.second; i++){
		for(int j = 0; j < unit.second; j++){
			if((i == ind.first && j == ind.first) || (i == ind.second && j == ind.second)){
				unit.first[i][j] = cos(phi);
			} else if ((i == ind.first && j == ind.second) || (i == ind.second && j == ind.first)){
				if (j > i){
					unit.first[i][j] = - sin(phi);
				} else {
					unit.first[i][j] = sin(phi);
				}
			} else if (i == j) {
				unit.first[i][j] = 1;
			} else if (i != j){
				unit.first[i][j] = 0;
			}
		}
	}
}

void transpose(std::pair<double**, int> unit){
	double tmp = 0;
	for (int i = 0; i < unit.second; i++){
		for (int j = 0; j < unit.second; j++){
			if (i > j){
				tmp = unit.first[i][j];
				unit.first[i][j] = unit.first[j][i];
				unit.first[j][i] = tmp;
			}
		}
	}
}

Matrix rotate_m(const Matrix& m){
	std::pair<int, int> indexes = max_elem_up_diag(m);
	double angle = phi(indexes, m);

	std::pair<double**, int> unit_m;
	unit_m.second = m.GetN();
	unit_m.first = new double* [unit_m.second];//выделили память под матрицу поворота
	for (int i = 0; i < unit_m.second; i++){
		unit_m.first[i] = new double [unit_m.second];
	}

	transition_matrix(indexes, angle, unit_m);//заполнили матрицу поворота
	Matrix transit(unit_m);
	E = E * transit;

	transpose(unit_m);
	Matrix transp(unit_m);

	Matrix result = transp * m * transit;
	for(int i = 0; i < result.GetN(); i++){
		for(int j = 0; j < result.GetN(); j++){
			if (i != j){
				if (result.GetA().first[i][j] > epsilon){
					return rotate_m(result);
				}
			}
		}
	}
	return result;
}

std::ostream& operator<<(std::ostream& out, const std::map<double, std::vector<double>>& mapa){
	int n = mapa.size() + 1;
	bool first = true;
	for(const auto& i : mapa){
		out << "lambda" << n - mapa.size() << " = " << std:: fixed
			<< std::setprecision(2) << std::setw(8) << std::left << i.first << "v" << n - mapa.size() << " = (";
		first = true;
		for (const auto& j : i.second){
			if(!first)
				out << ", ";
			out << std::setw(6) << std::right << j;
			first = false;
		}
		out << ")^T" << std::endl;
		n++;
	}
	return out;
}

std::map<double, std::vector<double>> eigenv(const Matrix& diag_m){
	std::map<double, std::vector<double>> eigen;//словарь пар собственное значение-собственный вектор
	int n = diag_m.GetN();
	double mult = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if(i == j){
				mult = 1 / E.GetA().first[i][j];//нормируем матрицу E - матрицу собственных векторов
				for (int k = 0; k < n; k++){
					//заполняем словарь: каждое с.значение - ключ словаря
					//в соответствие ставим вектор
					//(столбец матрицы полученной путем умножения всех матриц поворота)
					//(каждый столбец нормирован)
					eigen[diag_m.GetA().first[i][j]].push_back(E.GetA().first[k][j] * mult);
				}
				break;
			}
		}
	}
	std::cout << eigen;
	return eigen;
}
