#include "tests.h"
#include "Matrix.h"
#include <chrono>


int genRand(int low, int high, int s) {
    srand(time(0) + s*1000);
    return low + (rand() % high);
}

vector<pair<vector<vector<double>>, vector<double>>> testCross(int runs, int dim) {
    auto start = chrono::steady_clock::now();
    vector<pair<vector<vector<double>>, vector<double>>> fails {};
    auto getTime = [start](){return chrono::duration_cast<chrono::microseconds>(start - chrono::steady_clock::now()).count();};
    for(int i {0}; i < runs; ++i) {
        vector<vector<double>> vs;
        for(int k {0}; k < dim-1; ++k) {
            vs.push_back(vector<double>(dim));
            for(int j {0}; j < dim; ++j)
                vs[k][j] = genRand(0, 100, getTime());
        }

        vector<double> c = cross(vs);
        for(vector<double> v : vs) {
            if(dot(c, v) >= 1e-200) {
                fails.push_back(make_pair(vs, c));
            }
        }
    }
    cout << "Amount of fails: " << fails.size() << endl;
    return fails;
}


void demo() {
    
    	try {
		// Matrix construction
		Matrix <double> A(4, {1, 2, 3, 4, 3, 3, 4, 5, 3, 4, 5, 6});
		Matrix <double> B(3, {5, 5, 6, 7, 2, 4, 5, 7, 8, 1, 8, 2});
		vector <double> b {3, 6, 9};

		cout << "A = \n" << A << endl;
		cout << "B = \n" << B << endl;
		cout << "b = ";
		printVector(b);
		cout << "A_01 = " << A[0][1] << endl;	 // Indexing

		// Matrix operations
		Matrix <double> C = A*B;
		vector <double> prod = C*b;
		vector <double> solution = C.solve(b);

		cout << "\nC = A * B = \n" << C << endl; // Matrix multiplication
		cout << "A * b = ";						 // Matrix-vector multiplication
		printVector(prod);
		cout << "Solving Ax = b, x =";			 // Matrix solving (using cramers rule)
		printVector(solution);
		cout << "det(C) = " << C.det() << endl;  // Matrix determinant
		Matrix <double> AT = A.transpose();		 // Transpose
		cout << "\nA^T =\n" << AT << endl;

		// Other functionality
		vector <double> k {1, 1, 1, 1};
		vector <double> y = B.leastSquare(k);	 // Method of least squares
		cout << "y = ";
		printVector(y);
		vector <vector <double>> vecs = {b, y};
		cout << "b x y = ";
		printVector(cross(vecs));  // Cross product of vectors

		vector<double> Xdata {1, 2, 4, 5, 7, 8};
		vector<double> Ydata {2, 2, 5, 5, 8, 20};
		vector<double> result = regression(Xdata, Ydata, 3);
		cout << "Regression coefficients:";
		printVector(result);
	}
	catch (string error) {
		cout << error << endl;
		}

}