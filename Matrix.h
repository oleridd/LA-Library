#pragma once
#include "libs.h"


template<class T>
class Matrix {
    private:
        int N; // # Rows
        int M; // # Columns
        T** mat;

        Matrix<T> Mult(const Matrix<T>& A, const Matrix<T>& B) const;
        unsigned int CountZeros(bool columns = false) const;
        Matrix<T> excludeRC(int n, int m) const;
        void switch_rows(int ind_1, int ind_2);

    public:
        Matrix();
        Matrix(int N, int M);
        explicit Matrix(int N, vector<T> init);
        explicit Matrix(int N);
        Matrix(const Matrix<T>& A);
        ~Matrix();

        T get(int n, int m) const;
        void set(int n, int m, T val);
        int getRows() const { return N; }
        int getColumns() const { return M; }
        bool isValid() const { return mat != nullptr; }
        bool isLinIndependent() { return det() != 0; }
        bool isSquare() { return N == M; }
        Matrix<T> transpose() const;
        vector<T> to_vector() const;
        vector<T> LU_factorize();
        // Following are only valid for square matrices!
        T det() const;
        vector<T> solve(vector<T> b) const;
        vector<T> leastSquare(vector<T> b) const;
    

        T*& operator[](const int index);
        Matrix<T> &operator=(Matrix<T> A);

        Matrix<T> operator+=(const Matrix<T>& A);
        Matrix<T> operator+(const Matrix<T>& A);
        Matrix<T> operator-=(const Matrix<T>& A);
        Matrix<T> operator-(const Matrix<T>& A);
        Matrix<T> operator*(const Matrix<T>& A);
        vector<T> operator*(const vector<T>& b);
        // Only for square matrices!
        Matrix<T> operator*=(const Matrix<T>& A);
        Matrix<T> operator^(const int a);

        friend ostream& operator<<(ostream& os, Matrix<T>& A) {
            if(!A.isValid()) {
                os << "Matrix is not valid" << endl;
            } else {
                for(int i = 0; i < A.N; i++) {
                    for(int j = 0; j < A.M; j++) {
                        os << A.mat[i][j] << " ";
                    }
                    os << endl;
                }
            }
            return os;
        }
};