#include <iostream>
#include <fstream>
#include "rotate.h"

int main() {
	Matrix m1;
//	rotate_m(m1);
	E_mtrx(m1); //заполняем глобально объявленную единичную матрицу: нужна для последующего домножения
	// на матрицы поворота и нахождения собственных векторов
	eigenv(rotate_m(m1));
//	std::ifstream stream ("var1.txt");
//	Matrix m2(stream, 3);
//	Matrix x = m1 * m2;
//	std::cout << std::endl << x << std::endl;
	return 0;
}
