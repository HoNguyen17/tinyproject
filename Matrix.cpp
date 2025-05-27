#include "Matrix.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

// Default constructor - creates an empty matrix
Matrix::Matrix() : mNumRows(0), mNumCols(0), mData(nullptr) {
}

// Constructor with dimensions
Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    assert(numRows > 0 && numCols > 0);
    
    // Allocate memory for rows (array of pointers)
    mData = new double*[numRows];
    
    // Allocate memory for each row and initialize to zero
    for (int i = 0; i < numRows; i++) {
        mData[i] = new double[numCols];
        for (int j = 0; j < numCols; j++) {
            mData[i][j] = 0.0;
        }
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& m) : mNumRows(m.mNumRows), mNumCols(m.mNumCols) {
    // Allocate memory for the new matrix
    if (mNumRows > 0 && mNumCols > 0) {
        mData = new double*[mNumRows];
        
        // Allocate memory for each row and copy data
        for (int i = 0; i < mNumRows; i++) {
            mData[i] = new double[mNumCols];
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = m.mData[i][j];
            }
        }
    } else {
        mData = nullptr;
    }
}

// Destructor
Matrix::~Matrix() {
    // Free memory for each row
    if (mData) {
        for (int i = 0; i < mNumRows; i++) {
            delete[] mData[i];
        }
        // Free memory for the array of pointers
        delete[] mData;
        mData = nullptr;
    }
}

// Assignment operator
Matrix& Matrix::operator=(const Matrix& m) {
    // Check for self-assignment
    if (this == &m) {
        return *this;
    }
    
    // Free existing memory
    if (mData) {
        for (int i = 0; i < mNumRows; i++) {
            delete[] mData[i];
        }
        delete[] mData;
    }
    
    // Copy dimensions
    mNumRows = m.mNumRows;
    mNumCols = m.mNumCols;
    
    // Allocate new memory and copy data
    if (mNumRows > 0 && mNumCols > 0) {
        mData = new double*[mNumRows];
        
        for (int i = 0; i < mNumRows; i++) {
            mData[i] = new double[mNumCols];
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = m.mData[i][j];
            }
        }
    } else {
        mData = nullptr;
    }
    
    return *this;
}

// Round bracket operator for one-based indexing (non-const version)
double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

// Round bracket operator for one-based indexing (const version)
double Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

// Negation (unary minus)
Matrix Matrix::operator-() const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = -mData[i][j];
        }
    }
    return result;
}

// Unary plus (returns copy)
Matrix Matrix::operator+() const {
    return Matrix(*this);
}

// Matrix addition
Matrix Matrix::operator+(const Matrix& m) const {
    assert(mNumRows == m.mNumRows && mNumCols == m.mNumCols);
    
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] + m.mData[i][j];
        }
    }
    return result;
}

// Matrix subtraction
Matrix Matrix::operator-(const Matrix& m) const {
    assert(mNumRows == m.mNumRows && mNumCols == m.mNumCols);
    
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] - m.mData[i][j];
        }
    }
    return result;
}

// Matrix multiplication
Matrix Matrix::operator*(const Matrix& m) const {
    assert(mNumCols == m.mNumRows);
    
    Matrix result(mNumRows, m.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < m.mNumCols; j++) {
            double sum = 0.0;
            for (int k = 0; k < mNumCols; k++) {
                sum += mData[i][k] * m.mData[k][j];
            }
            result.mData[i][j] = sum;
        }
    }
    return result;
}

// Matrix-vector multiplication
Vector Matrix::operator*(const Vector& v) const {
    assert(mNumCols == v.get_size());
    
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        double sum = 0.0;
        for (int j = 0; j < mNumCols; j++) {
            sum += mData[i][j] * v[j];
        }
        result[i] = sum;
    }
    return result;
}

// Matrix-scalar multiplication
Matrix Matrix::operator*(double scalar) const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return result;
}

// Matrix-scalar division
Matrix Matrix::operator/(double scalar) const {
    assert(scalar != 0.0);
    
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] / scalar;
        }
    }
    return result;
}

// Scalar-matrix multiplication (friend function)
Matrix operator*(double scalar, const Matrix& m) {
    return m * scalar;
}

// Vector-matrix multiplication (friend function)
Vector operator*(const Vector& v, const Matrix& m) {
    assert(v.get_size() == m.mNumRows);
    
    Vector result(m.mNumCols);
    for (int j = 0; j < m.mNumCols; j++) {
        double sum = 0.0;
        for (int i = 0; i < m.mNumRows; i++) {
            sum += v[i] * m.mData[i][j];
        }
        result[j] = sum;
    }
    return result;
}

// Compute determinant
double Matrix::Determinant() const {
    assert(mNumRows == mNumCols);
    
    // Base cases for small matrices
    if (mNumRows == 1) {
        return mData[0][0];
    }
    if (mNumRows == 2) {
        return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    }
    
    // For larger matrices, use cofactor expansion along first row
    double det = 0.0;
    for (int j = 0; j < mNumCols; j++) {
        // Create submatrix by excluding current row and column
        Matrix subMatrix(mNumRows - 1, mNumCols - 1);
        for (int row = 1; row < mNumRows; row++) {
            for (int col = 0, subCol = 0; col < mNumCols; col++) {
                if (col != j) {
                    subMatrix.mData[row-1][subCol++] = mData[row][col];
                }
            }
        }
        
        // Alternate sign based on column index
        double sign = (j % 2 == 0) ? 1.0 : -1.0;
        det += sign * mData[0][j] * subMatrix.Determinant();
    }
    
    return det;
}

// Compute inverse using adjoint method
Matrix Matrix::Inverse() const {
    assert(mNumRows == mNumCols);
    
    double det = Determinant();
    assert(fabs(det) > 1e-10);  // Check if matrix is invertible
    
    Matrix adjoint(mNumRows, mNumCols);
    
    // Special case for 1x1 matrix
    if (mNumRows == 1) {
        adjoint.mData[0][0] = 1.0;
        return adjoint / det;
    }
    
    // Special case for 2x2 matrix
    if (mNumRows == 2) {
        adjoint.mData[0][0] = mData[1][1];
        adjoint.mData[0][1] = -mData[0][1];
        adjoint.mData[1][0] = -mData[1][0];
        adjoint.mData[1][1] = mData[0][0];
        return adjoint / det;
    }
    
    // For larger matrices, compute cofactor matrix
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            // Create submatrix by excluding current row and column
            Matrix subMatrix(mNumRows - 1, mNumCols - 1);
            
            for (int row = 0, subRow = 0; row < mNumRows; row++) {
                if (row != i) {
                    for (int col = 0, subCol = 0; col < mNumCols; col++) {
                        if (col != j) {
                            subMatrix.mData[subRow][subCol++] = mData[row][col];
                        }
                    }
                    subRow++;
                }
            }
            
            // Compute cofactor and place in adjoint (with transposition)
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            adjoint.mData[j][i] = sign * subMatrix.Determinant();  // Note the j,i indexing for transpose
        }
    }
    
    return adjoint / det;
}

// Compute pseudo-inverse (Moore-Penrose)
Matrix Matrix::PseudoInverse() const {
    // For full rank square matrices, regular inverse
    if (mNumRows == mNumCols) {
        try {
            return Inverse();
        } catch (...) {
            // Continue with SVD approach if matrix is singular
        }
    }
    
    // For non-square matrices, use the formula: A⁺ = (Aᵀ·A)⁻¹·Aᵀ for full column rank
    // or A⁺ = Aᵀ·(A·Aᵀ)⁻¹ for full row rank
    
    if (mNumRows >= mNumCols) {  // More rows than columns or equal
        // Compute Aᵀ
        Matrix AT(mNumCols, mNumRows);
        for (int i = 0; i < mNumRows; i++) {
            for (int j = 0; j < mNumCols; j++) {
                AT.mData[j][i] = mData[i][j];
            }
        }
        
        // Compute Aᵀ·A
        Matrix ATA = AT * (*this);
        
        // Compute (Aᵀ·A)⁻¹
        Matrix ATAInv = ATA.Inverse();
        
        // Compute (Aᵀ·A)⁻¹·Aᵀ
        return ATAInv * AT;
    } else {  // More columns than rows
        // Compute Aᵀ
        Matrix AT(mNumCols, mNumRows);
        for (int i = 0; i < mNumRows; i++) {
            for (int j = 0; j < mNumCols; j++) {
                AT.mData[j][i] = mData[i][j];
            }
        }
        
        // Compute A·Aᵀ
        Matrix AAT = (*this) * AT;
        
        // Compute (A·Aᵀ)⁻¹
        Matrix AATInv = AAT.Inverse();
        
        // Compute Aᵀ·(A·Aᵀ)⁻¹
        return AT * AATInv;
    }
}

// Print matrix
void Matrix::Print() const {
    cout << "[" << endl;
    for (int i = 0; i < mNumRows; i++) {
        cout << "  ";
        for (int j = 0; j < mNumCols; j++) {
            cout << setw(10) << setprecision(4) << mData[i][j];
        }
        cout << endl;
    }
    cout << "]" << endl;
}
