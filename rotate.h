#pragma once
#include <cmath>
#include <utility>
#include <vector>
#include <map>
#include "matrix.h"

extern double epsilon;
extern Matrix E;

void E_mtrx(const Matrix& m);
std::pair<int, int> max_elem_up_diag(const Matrix& m);
double phi(std::pair<int, int> ind, const Matrix& m);
void transition_matrix(std::pair<int, int> ind, double phi, std::pair<double**, int> unit);
void transpose(std::pair<double**, int> unit);
Matrix rotate_m(const Matrix& m);
std::ostream& operator<<(std::ostream& out, const std::map<double, std::vector<double>>& mapa);
std::map<double, std::vector<double>> eigenv(const Matrix& diag_m);
