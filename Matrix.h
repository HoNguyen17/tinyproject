#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#include <cassert>

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    double** mData;  // Pointer to pointer to a data element

public:
    // Constructors
    Matrix();  // Default constructor
    Matrix(int numRows, int numCols);  // Constructor with dimensions
    Matrix(const Matrix& m);  // Copy constructor
    ~Matrix();  // Destructor

    // Access methods
    int GetNumRows() const { return mNumRows; }
    int GetNumCols() const { return mNumCols; }

    // Overloaded operators
    Matrix& operator=(const Matrix& m);  // Assignment operator
    
    // Round bracket operator for one-based indexing
    double& operator()(int i, int j);  // Non-const version
    double operator()(int i, int j) const;  // Const version

    // Unary operators
    Matrix operator-() const;  // Negation
    Matrix operator+() const;  // Unary plus (returns copy)

    // Binary operations with matrices
    Matrix operator+(const Matrix& m) const;  // Matrix addition
    Matrix operator-(const Matrix& m) const;  // Matrix subtraction
    Matrix operator*(const Matrix& m) const;  // Matrix multiplication

    // Operations with vectors
    Vector operator*(const Vector& v) const;  // Matrix-vector multiplication

    // Operations with scalars
    Matrix operator*(double scalar) const;  // Matrix-scalar multiplication
    Matrix operator/(double scalar) const;  // Matrix-scalar division

    // Matrix operations
    double Determinant() const;  // Compute determinant
    Matrix Inverse() const;  // Compute inverse
    Matrix PseudoInverse() const;  // Compute pseudo-inverse (Moore-Penrose)

    // Helper methods
    void Print() const;  // Print matrix

    // Friend functions
    friend Matrix operator*(double scalar, const Matrix& m);  // Scalar-matrix multiplication
    friend Vector operator*(const Vector& v, const Matrix& m);  // Vector-matrix multiplication
};

#endif // MATRIX_H
