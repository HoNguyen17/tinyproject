#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Matrix.h"
#include "Vector.h"

class LinearSystem {
protected:
    int mSize;      // Size of the linear system
    Matrix* mpA;    // Pointer to coefficient matrix
    Vector* mpb;    // Pointer to right-hand side vector

public:
    // Constructor with matrix and vector
    LinearSystem(const Matrix& A, const Vector& b);
    
    // Virtual destructor
    virtual ~LinearSystem();
    
    // Delete the default constructor to prevent its use
    LinearSystem() = delete;
    
    // Delete the copy constructor to prevent its use
    LinearSystem(const LinearSystem& other) = delete;
    
    // Delete the assignment operator
    LinearSystem& operator=(const LinearSystem& other) = delete;
    
    // Get size of the system
    int GetSize() const { return mSize; }
    
    // Get coefficient matrix
    const Matrix& GetMatrix() const { return *mpA; }
    
    // Get right-hand side vector
    const Vector& GetRHS() const { return *mpb; }
    
    // Solve the linear system (using Gaussian elimination with pivoting)
    virtual Vector Solve() const;
};

// Derived class for positive definite symmetric linear systems
class PosSymLinSystem : public LinearSystem {
public:
    // Constructor with matrix and vector
    PosSymLinSystem(const Matrix& A, const Vector& b);
    
    // Virtual destructor
    virtual ~PosSymLinSystem();
    
    // Override the Solve method to use conjugate gradient method
    virtual Vector Solve() const override;
    
private:
    // Helper function to verify if a matrix is symmetric
    bool IsSymmetric(const Matrix& A) const;
};

// Class for non-square linear systems (under-determined or over-determined)
class NonSquareLinSystem : public LinearSystem {
public:
    // Constructor with matrix, vector, and regularization parameter
    NonSquareLinSystem(const Matrix& A, const Vector& b, double lambda = 0.0);
    
    // Virtual destructor
    virtual ~NonSquareLinSystem();
    
    // Override the Solve method to use pseudo-inverse or Tikhonov regularization
    virtual Vector Solve() const override;
    
    // Method to solve using Moore-Penrose pseudo-inverse
    Vector SolvePseudoInverse() const;
    
    // Method to solve using Tikhonov regularization
    Vector SolveTikhonov() const;
    
private:
    double mLambda; // Regularization parameter for Tikhonov regularization
};

#endif // LINEARSYSTEM_H
