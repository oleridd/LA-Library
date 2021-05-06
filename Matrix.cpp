#include "Matrix.h"


template<class T>
Matrix<T>::Matrix() : N(0), M(0), mat(nullptr) {}
// Constructor for empty matrix


template<class T>
Matrix<T>::Matrix(int N, int M) : N(N), M(M), mat(new T*[N]) {
// Constructor for custom dimensions (identity matrix)
    for(int i = 0; i < N; i++) {
        mat[i] = new T[M];
        for(int j = 0; j < M; j++) {
            if(i == j) mat[i][j] = 1;
            else mat[i][j] = 0;
        }
    }
}


template<class T>
Matrix<T>::Matrix(int M, vector<T> init) : Matrix<T>(init.size() / M, M) {
// Constructor for custom variables.
    if(init.size()%M) {
        mat = nullptr;
        throw range_error("Vector is not divisible by M");
    } else {
        int k = 0;
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                mat[i][j] = init[k++];
            }
        }
    }
}


template<class T>
Matrix<T>::Matrix(int N) :
// Constructor for square matrix (identity matrix)
    Matrix<T>(N, N) {}


template<class T>
Matrix<T>::Matrix(const Matrix<T>& A) : Matrix<T>(A.N, A.M) {
// Copy constructor
    if(!A.isValid()) mat = nullptr;
    else {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                mat[i][j] = A.mat[i][j];
            }
        }
    }
}


template<class T>
Matrix<T>::~Matrix() {
// Destructor
    if (mat != nullptr) {
        for(int i = 0; i < N; i++) delete[] mat[i];
        delete[] mat;
    }
}


template<class T>
T Matrix<T>::get(int n, int m) const {
// Get variables at index m, n
    if(n < N && m < M) return mat[n][m];
    else throw range_error("Index is out of range");
}


template<class T>
void Matrix<T>::set(int n, int m, T val) {
// Set variables at index m, n
    if(n < N && m < M) mat[n][m] = val;
    else throw range_error("Index is out of range");
}


template<class T>
void Matrix<T>::switch_rows(int ind_1, int ind_2) {
    T* placeholder = mat[ind_1];
    mat[ind_1] = mat[ind_2];
    mat[ind_2] = placeholder;
}


template<class T>
Matrix<T> Matrix<T>::Mult(const Matrix<T>& A, const Matrix<T>& B) const {
// Matrix multiplication function
    if(M != B.N) { 
        throw runtime_error("Invalid matrix dimensions");
        return Matrix<T>();
    } else {
        Matrix<T>newMat(A.N, B.M);
        for(int i = 0; i < A.N; i++) { // Through all rows in the left matrix.
            for(int j = 0; j < B.M; j++) { // Through all columns in the right matrix.
            newMat.mat[i][j] = 0;
                for(int k = 0; k < A.M; k++) { // Dot product through first row left matrix.
                    newMat.mat[i][j] += A.mat[i][k] * B.mat[k][j];
                }
            }
        }
        return newMat;
    }
}


template<class T>
unsigned int Matrix<T>::CountZeros(bool columns) const {
// Private function. Returns index to the row (or if columns = true, the column) that has the most zeros.
    Matrix<T> A;
    if(columns) A = this->transpose(); // Counts zeros for each column.
    else A = *this; // Counts zeros for each row.
    vector<int> zeros;
    for(int i = 0; i < N; i++) {
        zeros.push_back(0);
        for(int j = 0; j < M; j++) {
            if(A.mat[i][j] == 0) zeros[i]++;
        }
    }
    vector<int>::iterator it = max_element(zeros.begin(), zeros.end());
    return it - zeros.begin();
}


template<class T>
Matrix<T> Matrix<T>::excludeRC(int n, int m) const {
// Private function. Excludes the row and column containing the index at n, m.
    Matrix<T> B(N-1, N-1);
    int iScalar {0};
    int jScalar {0};
    for(int i = 0; i < N; i++) {
        if(i == n-1) { iScalar = 1; continue; }
        for(int j = 0; j < M; j++) {
            if(j == m-1) { jScalar = 1; continue; }
            B[i-iScalar][j-jScalar] = mat[i][j];
        }
        jScalar = 0;
    }
    return B;
}


template<class T>
Matrix<T> Matrix<T>::transpose() const {
// Returns transposed matrix
    Matrix<T> B(M, N);
    for(int i = 0; i < N; i++) { // Iterating through rows.
        for(int j = 0; j < M; j++) { // Iterating through columns.
            B[j][i] = mat[i][j];
        }
    }
    return B;
}


template<class T>
vector<T> Matrix<T>::to_vector() const {
// Converts 1xn og mx1 matrix to vector
    int m {0};
    int n {0};
    vector<T> res {};
    if(N > 1 && M > 1) throw runtime_error("Matrix must have only one row or column");
    else if(N == 1) m = 1;
    else if(M == 1) n = 1;
    for(int i = 0; i < (M*m + N*n); i++) {
        res.push_back(mat[n*i][m*i]);
    }
    return res;
}


template<class T>
vector<T> Matrix<T>::LU_factorize() {
// Factorizes the matrix to the LU factorization if possible.
// Returns a vector: indices for pivots. Also changes *this to be factorized.
    assert(N == M);

    pair <int, T> maxim(0, mat[0][0]);
    for (int n {0}; n < N; ++n) {
        if (mat[n][0] > maxim.second) {
            maxim = pair <int, T> (n, mat[n][0]);
        }
    }
    switch_rows(0, maxim.first);

    return {1, 2, 3};
}


template<class T>
T Matrix<T>::det() const {
// Returns the determinant of the matrix
    T det {0};
    int rowInd = CountZeros();
    for(int i = 0; i < N; i++) { // Iterates through columns in the selected row
        if(N != 2) {
            if(mat[rowInd][i] == 0) continue;
            det += pow(-1, i+rowInd) * mat[rowInd][i] * excludeRC(rowInd+1, i+1).det();
        } else det += pow(-1, i) * mat[i][0] * mat[1-i][1];
    }
    return det;
}


template<class T>
vector<T> Matrix<T>::solve(vector<T> b) const {
// Solves the matrix using cramers rule
    if(N != M) throw runtime_error("Matrix must be square to solve");
    if(det() == 0) throw runtime_error("Division by zero in solve()");
    vector<T> res {};
    for(int i = 0; i < M; i++) {
        Matrix<T> B = *this;
        for(int j = 0; j < N; j++) { // Need to replace column vector number i+1
            B[j][i] = b[j];
        }
        res.push_back(B.det()/det());
    }
    return res;
}


template<class T>
vector<T> Matrix<T>::leastSquare(vector<T> b) const {
// Least square method for overdetermined matrices
// b must be a vector that is not in the column space of the matrix
    if(M >= N) throw runtime_error("Matrix must be overdetermined in leastSquare()");
    Matrix<T> B = transpose() * *this;
    vector<T> bnew {transpose() * b};
    return B.solve(bnew);
}


template<class T>
T*& Matrix<T>::operator[](const int index) {
    if(index < N) return mat[index];
    else throw range_error("Index is out of range");
}


template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>A) {
    swap(mat, A.mat);
    return *this;
}


template<class T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& A) {
    if(!(M == A.M && N == A.N)) { 
        mat = nullptr;
        throw runtime_error("Matrix dimensions do not match");
    } else {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                mat[i][j] += A.mat[i][j];
            }
        }
    }
    return *this;
}


template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) {
    Matrix<T> B = *this;
    return B += A;
}


template<class T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T>& A) {
    if(!(M == A.M && N == A.N)) {
        mat = nullptr;
        throw runtime_error("Matrix dimensions do not match");
    } else {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                mat[i][j] -= A.mat[i][j];
            }
        }
    }
    return *this;
}


template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) {
    Matrix<T> B = *this;
    return B -= A;
}


template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& A) {
    return Mult(*this, A);
}


template<class T>
vector<T> Matrix<T>::operator*(const vector<T>& b) {
    if(M != size(b)) throw runtime_error("b must have " + to_string(M) + " elements");
    // if(M == Matrix(M.N)) return b;
    Matrix<T>B(1, b);
    Matrix<T>C = *this * B;
    return C.to_vector();
}


template<class T>
Matrix<T> Matrix<T>::operator*=(const Matrix<T>& A) {
    if(N != M) { 
        mat = nullptr;
        throw runtime_error("Matrix must be square");
    } else {
        *this = *this * A;
    }
    return *this;
}



template class Matrix<double long>;
template class Matrix<double>;
template class Matrix<float>;
template class Matrix<int>;
template class Matrix<long>;