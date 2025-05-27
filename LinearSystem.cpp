#include "LinearSystem.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

//==============================
// LinearSystem Implementation
//==============================

// Constructor with matrix and vector
LinearSystem::LinearSystem(const Matrix& A, const Vector& b) {
    // Check if matrix and vector dimensions are compatible
    if (A.GetNumRows() != b.get_size()) {
        throw std::invalid_argument("Matrix and vector dimensions are incompatible");
    }
    
    mSize = A.GetNumRows();
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

// Destructor
LinearSystem::~LinearSystem() {
    delete mpA;
    delete mpb;
}

// Solve the linear system using Gaussian elimination with pivoting
Vector LinearSystem::Solve() const {
    // Check if the matrix is square
    if (mpA->GetNumRows() != mpA->GetNumCols()) {
        throw std::invalid_argument("Base LinearSystem::Solve requires a square matrix");
    }
    // Create copies of the matrix and vector to avoid modifying the originals
    Matrix A = *mpA;
    Vector b = *mpb;
    Vector x(mSize);
    
    // Gaussian elimination with partial pivoting
    for (int k = 1; k <= mSize - 1; k++) {
        // Find pivot
        int pivot_row = k;
        double max_val = std::abs(A(k, k));
        
        for (int i = k + 1; i <= mSize; i++) {
            if (std::abs(A(i, k)) > max_val) {
                max_val = std::abs(A(i, k));
                pivot_row = i;
            }
        }
        
        // Swap rows if necessary
        if (pivot_row != k) {
            for (int j = k; j <= mSize; j++) {
                double temp = A(k, j);
                A(k, j) = A(pivot_row, j);
                A(pivot_row, j) = temp;
            }
            
            double temp = b(k);
            b(k) = b(pivot_row);
            b(pivot_row) = temp;
        }
        
        // Check for singularity
        if (std::abs(A(k, k)) < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }
        
        // Elimination
        for (int i = k + 1; i <= mSize; i++) {
            double factor = A(i, k) / A(k, k);
            
            for (int j = k; j <= mSize; j++) {
                A(i, j) -= factor * A(k, j);
            }
            
            b(i) -= factor * b(k);
        }
    }
    
    // Back substitution
    for (int i = mSize; i >= 1; i--) {
        double sum = 0.0;
        
        for (int j = i + 1; j <= mSize; j++) {
            sum += A(i, j) * x(j);
        }
        
        x(i) = (b(i) - sum) / A(i, i);
    }
    
    return x;
}

//==============================
// PosSymLinSystem Implementation
//==============================

// Constructor
PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b)
    : LinearSystem(A, b) {
    // Check if the matrix is symmetric
    if (!IsSymmetric(A)) {
        throw std::invalid_argument("Matrix must be symmetric for PosSymLinSystem");
    }
}

// Destructor
PosSymLinSystem::~PosSymLinSystem() {
    // Base class destructor will be called automatically
}

// Check if matrix is symmetric
bool PosSymLinSystem::IsSymmetric(const Matrix& A) const {
    if (A.GetNumRows() != A.GetNumCols()) {
        return false;
    }
    
    int n = A.GetNumRows();
    const double tolerance = 1e-10;
    
    for (int i = 1; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
            if (std::abs(A(i, j) - A(j, i)) > tolerance) {
                return false;
            }
        }
    }
    
    return true;
}

// Solve using conjugate gradient method
Vector PosSymLinSystem::Solve() const {
    const double tolerance = 1e-10;
    const int max_iterations = mSize * 10;
    
    // Create initial guess x = 0
    Vector x(mSize);
    
    // r = b - A*x
    Vector r = *mpb - (*mpA) * x;
    
    // p = r
    Vector p = r;
    
    // Initial residual norm squared
    double rsold = r * r;
    
    if (rsold < tolerance) {
        return x; // Already at solution
    }
    
    for (int iter = 0; iter < max_iterations; iter++) {
        // A*p
        Vector Ap = (*mpA) * p;
        
        // alpha = rsold / (p' * A * p)
        double alpha = rsold / (p * Ap);
        
        // x = x + alpha * p
        x = x + (p * alpha);
        
        // r = r - alpha * A * p
        r = r - (Ap * alpha);
        
        // Check convergence
        double rsnew = r * r;
        if (std::sqrt(rsnew) < tolerance) {
            break;
        }
        
        // beta = rsnew / rsold
        double beta = rsnew / rsold;
        
        // p = r + beta * p
        p = r + (p * beta);
        
        // Update rsold
        rsold = rsnew;
    }
    
    return x;
}

//==============================
// NonSquareLinSystem Implementation
//==============================

// Constructor with matrix, vector, and regularization parameter
NonSquareLinSystem::NonSquareLinSystem(const Matrix& A, const Vector& b, double lambda)
    : LinearSystem(A, b), mLambda(lambda) {
    // Note: LinearSystem constructor already checks if A and b are compatible
    // But it does not require A to be square anymore
}

// Destructor
NonSquareLinSystem::~NonSquareLinSystem() {
    // Base class destructor will be called automatically
}

// Solve method - choose between pseudo-inverse and Tikhonov based on lambda value
Vector NonSquareLinSystem::Solve() const {
    if (mLambda == 0.0) {
        return SolvePseudoInverse();
    } else {
        return SolveTikhonov();
    }
}

// Solve using Moore-Penrose pseudo-inverse
Vector NonSquareLinSystem::SolvePseudoInverse() const {
    // Get pseudo-inverse of A
    Matrix Aplus = mpA->PseudoInverse();
    
    // Return x = A⁺ * b
    return Aplus * (*mpb);
}

// Solve using Tikhonov regularization
Vector NonSquareLinSystem::SolveTikhonov() const {
    int n = mpA->GetNumCols(); // Number of columns in A
    
    // Create transpose of A
    Matrix AT(mpA->GetNumCols(), mpA->GetNumRows());
    for (int i = 1; i <= mpA->GetNumRows(); i++) {
        for (int j = 1; j <= mpA->GetNumCols(); j++) {
            AT(j, i) = (*mpA)(i, j);
        }
    }
    
    // Create identity matrix
    Matrix I(n, n);
    for (int i = 1; i <= n; i++) {
        I(i, i) = 1.0;
    }
    
    // Compute A^T * A + lambda * I
    Matrix ATA = AT * (*mpA);
    Matrix regularized = ATA + (I * mLambda);
    
    // Compute (A^T * A + lambda * I)^(-1)
    Matrix inv = regularized.Inverse();
    
    // Compute (A^T * A + lambda * I)^(-1) * A^T * b
    Vector ATb = AT * (*mpb);
    Vector x = inv * ATb;
    
    return x;
}
