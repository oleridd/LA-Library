#pragma once
#include "Matrix.h"


double getVectorLength(vector<double> v);
vector<double> unify(vector<double> v);

template <class T>
T dot(vector<T> v1, vector<T> v2) {
    assert(v1.size() == v2.size());

    T prod;
    for(int i {0}; i < v1.size(); ++i) {
        prod += v1[i]*v2[i];
    }
    return prod;
}

template<class T>
vector<T> regression(vector<T> Xdata, vector<T> Ydata, int n) {
// Preforms regression or interpolation (depending on degree and amount of variables)
// on the data points Xdata and Ydata
// n is the amount of coefficients to solve for. Thit is, the degree - 1
    if(Xdata.size() != Ydata.size()) throw runtime_error("Datasets do not match");
    else if(n > Xdata.size()) 
        throw runtime_error("The polynomial can not have a degree greater than " + to_string(Xdata.size()));

    // Setting up the system of equations:
    vector<T> initVec {};
    for(int i = 0; i < Xdata.size(); i++) {
        for(int j = 0; j < n; j++) {
            initVec.push_back(pow(Xdata[i], n-(j+1)));
        }
    }
    Matrix<T> A(n, initVec);

    // Preforming regression/interpolation:
    if(A.isSquare()) return A.solve(Ydata); // If the matrix is square we don't need regression.
    else return A.leastSquare(Ydata);
}


template <class T>
pair<Matrix<T>, Matrix<T>> meshgrid(vector<T> xx, vector<T> yy) {
// Returns pair of matrices making a meshgrid
    int M = xx.size();
    int N = yy.size();
    Matrix<T> xxMatrix(N, M);
    Matrix<T> yyMatrix(N, M);
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            xxMatrix[i][j] = xx[j];
            yyMatrix[i][j] = yy[i];
        }
    }
    return make_pair(xxMatrix, yyMatrix);
}

template <class T>
Matrix<T> fillMatrix(vector<vector<T>> vs) {
// Fills matrix with the components of vs
    int dim = vs.size();
    Matrix<T> A(dim);
    for(int i {0}; i < dim; ++i) {
        for(int j {0}; j < dim; ++j) {
            A[i][j] = vs[i][j];
        }
    }
    return A;
}

template <class T>
vector<T> cross(vector<vector<T>> vs) {
// Experimental. Takes the cross product between n vectors
// All vectors must be one dimension higher than the amount of vectors.
    int dim = vs[0].size();
    assert(dim == vs.size()+1);
    for(vector<T> i : vs)
        assert(i.size() == dim);

    vector<T> newVec {};
    for(int i {0}; i < dim; ++i) {
        vector<vector<T>> vsExcl;
        for(vector<T> v : vs) {
            vector<T> vTemp = v;
            vTemp.erase(vTemp.begin()+i);
            vsExcl.push_back(vTemp);
        }
        Matrix<T> A = fillMatrix(vsExcl); // This will be transposed, but the determinant is unaffected.
        newVec.push_back(pow(-1, i) * A.det());
    }
    return newVec;
}


template <class T>
void printVector(vector<T> vec) {
    cout << "[";
    for (T el : vec) {
        cout << el << " ";
    }
    cout << "]" << endl;
}