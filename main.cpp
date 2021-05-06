#include "Matrix.h"
#include "Complex.h"
#include "Extra.h"
#include "tests.h"


int main() {

	Matrix<double> A(3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
	cout << A << endl;

	return 0;
	
}
