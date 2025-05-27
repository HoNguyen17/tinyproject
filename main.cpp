#include "Vector.h"
#include "Matrix.h"
#include "LinearSystem.h"
#include <iostream>
#include <iomanip>
using namespace std;

// Compile with: cd "d:\CSE\nam 2\program2\tinyproject\part_A"; g++ -std=c++17 -Wall 
//                       main.cpp Vector.cpp Matrix.cpp LinearSystem.cpp -o linear_system_test
// Run with: cd "d:\CSE\nam 2\program2\tinyproject\part_A"; .\linear_system_test.exe
// Delete object file with: rm main.o
int main() {
    cout << "===== VECTOR TESTING =====" << endl;
    // Create vectors
    Vector v1(3);
    Vector v2(3);
    
    // Set values using zero-indexed operator
    v1[0] = 1.0;
    v1[1] = 2.0;
    v1[2] = 3.0;
    
    // Set values using one-indexed operator
    v2(1) = 4.0;
    v2(2) = 5.0;
    v2(3) = 6.0;
    
    // Print vectors
    cout << "v1 = ";
    v1.print();
    
    cout << "v2 = ";
    v2.print();
    
    // Vector addition
    Vector v3 = v1 + v2;
    cout << "v1 + v2 = ";
    v3.print();
    
    // Vector subtraction
    Vector v4 = v1 - v2;
    cout << "v1 - v2 = ";
    v4.print();
    
    // Dot product
    double dot = v1 * v2;
    cout << "v1 * v2 (dot product) = " << dot << std::endl;
    
    // Scalar multiplication
    Vector v5 = v1 * 2.0;
    cout << "v1 * 2.0 = ";
    v5.print();
    
    Vector v6 = 3.0 * v2;
    cout << "3.0 * v2 = ";
    v6.print();
    
    // Negation
    Vector v7 = -v1;
    cout << "-v1 = ";
    v7.print();
    
    // Division by scalar
    Vector v8 = v1 / 2.0;
    cout << "v1 / 2.0 = ";
    v8.print();
    
    // Copy constructor and assignment
    Vector v9 = v1;  // Copy constructor
    cout << "v9 (copy of v1) = ";
    v9.print();
    
    Vector v10(2);
    v10 = v2;  // Assignment operator
    cout << "v10 (assigned v2) = ";
    v10.print();
    
    // Access with bounds checking
    try {
        cout << "v1[1] = " << v1[1] << std::endl;  // Zero-indexed
        cout << "v1(2) = " << v1(2) << std::endl;  // One-indexed (same element)
        
        // This will throw an exception
        cout << "v1[3] = " << v1[3] << std::endl;
    }
    catch (const std::exception& e) {
        cout << "Exception caught: " << e.what() << std::endl;
    }
    
    cout << "\n===== MATRIX TESTING =====" << endl;
    
    // Create matrices
    Matrix A(2, 3);
    Matrix B(3, 2);
    
    // Set values using one-indexed operator
    A(1, 1) = 1.0; A(1, 2) = 2.0; A(1, 3) = 3.0;
    A(2, 1) = 4.0; A(2, 2) = 5.0; A(2, 3) = 6.0;
    
    B(1, 1) = 7.0; B(1, 2) = 8.0;
    B(2, 1) = 9.0; B(2, 2) = 10.0;
    B(3, 1) = 11.0; B(3, 2) = 12.0;
    
    // Print matrices
    cout << "Matrix A:" << endl;
    A.Print();
    
    cout << "Matrix B:" << endl;
    B.Print();
    
    // Matrix multiplication
    cout << "A * B:" << endl;
    Matrix C = A * B;
    C.Print();
    
    // Matrix-vector multiplication
    cout << "A * v1:" << endl;
    Vector v11 = A * v1;
    v11.print();
    
    // Create a square matrix for more operations
    Matrix D(3, 3);
    D(1, 1) = 2.0; D(1, 2) = 0.0; D(1, 3) = 1.0;
    D(2, 1) = 0.0; D(2, 2) = 3.0; D(2, 3) = 0.0;
    D(3, 1) = -1.0; D(3, 2) = 0.0; D(3, 3) = 4.0;
    
    cout << "Matrix D:" << endl;
    D.Print();
    
    // Matrix scalar operations
    cout << "D * 2.0:" << endl;
    Matrix E = D * 2.0;
    E.Print();
    
    cout << "3.0 * D:" << endl;
    Matrix F = 3.0 * D;
    F.Print();
    
    cout << "D / 2.0:" << endl;
    Matrix G = D / 2.0;
    G.Print();
    
    // Matrix determinant
    cout << "Determinant of D: " << D.Determinant() << endl;
    
    // Matrix inverse
    cout << "Inverse of D:" << endl;
    Matrix Dinv = D.Inverse();
    Dinv.Print();
    
    // Verify inverse: D * D^-1 should be identity
    cout << "D * D^-1 (should be identity):" << endl;
    Matrix I = D * Dinv;
    I.Print();
    
    // ===== LINEAR SYSTEM TESTING =====
    cout << "\n===== LINEAR SYSTEM TESTING =====" << endl;
    
    // Create a 3x3 square matrix
    Matrix M(3, 3);
    M(1, 1) = 4.0; M(1, 2) = 3.0; M(1, 3) = 2.0;
    M(2, 1) = 1.0; M(2, 2) = 4.0; M(2, 3) = 3.0;
    M(3, 1) = 2.0; M(3, 2) = 1.0; M(3, 3) = 4.0;
    
    // Right-hand side vector
    Vector b(3);
    b(1) = 10.0;
    b(2) = 13.0;
    b(3) = 14.0;
    
    cout << "Linear system matrix:" << endl;
    M.Print();
    
    cout << "Right-hand side vector:" << endl;
    b.print();
    
    // Solve using Gaussian elimination
    cout << "\nSolving using LinearSystem (Gaussian elimination):" << endl;
    try {
        LinearSystem ls(M, b);
        Vector x = ls.Solve();
        
        cout << "Solution vector:" << endl;
        x.print();
        
        // Verify solution
        Vector residual = b - M * x;
        cout << "Residual (should be close to zero):" << endl;
        residual.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in LinearSystem: " << e.what() << endl;
    }
    
    // Create a symmetric positive definite matrix
    Matrix S(3, 3);
    S(1, 1) = 4.0; S(1, 2) = 1.0; S(1, 3) = 0.0;
    S(2, 1) = 1.0; S(2, 2) = 5.0; S(2, 3) = 2.0;
    S(3, 1) = 0.0; S(3, 2) = 2.0; S(3, 3) = 3.0;
    
    // Right-hand side vector for symmetric system
    Vector c(3);
    c(1) = 6.0;
    c(2) = 15.0;
    c(3) = 11.0;
    
    cout << "\nSymmetric positive definite matrix:" << endl;
    S.Print();
    
    cout << "Right-hand side vector:" << endl;
    c.print();
    
    // Solve using conjugate gradient
    cout << "\nSolving using PosSymLinSystem (Conjugate Gradient):" << endl;
    try {
        PosSymLinSystem pls(S, c);
        Vector y = pls.Solve();
        
        cout << "Solution vector:" << endl;
        y.print();
        
        // Verify solution
        Vector residual = c - S * y;
        cout << "Residual (should be close to zero):" << endl;
        residual.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in PosSymLinSystem: " << e.what() << endl;
    }
    
    // Test with non-symmetric matrix
    cout << "\nTesting non-symmetric matrix with PosSymLinSystem:" << endl;
    try {
        PosSymLinSystem pls2(M, b);  // M is not symmetric
        Vector z = pls2.Solve();
        z.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in PosSymLinSystem: " << e.what() << endl;
    }
    
    // ===== NON-SQUARE LINEAR SYSTEM TESTING =====
    cout << "\n===== NON-SQUARE LINEAR SYSTEM TESTING =====" << endl;
    
    // Test over-determined system (more equations than unknowns)
    cout << "\nOver-determined system (more equations than unknowns):" << endl;
    Matrix O(4, 2);  // 4 equations, 2 unknowns
    O(1, 1) = 1.0; O(1, 2) = 1.0;
    O(2, 1) = 2.0; O(2, 2) = 1.0;
    O(3, 1) = 1.0; O(3, 2) = 2.0;
    O(4, 1) = 1.0; O(4, 2) = 1.0;
    
    Vector d(4);
    d(1) = 2.0;
    d(2) = 3.0;
    d(3) = 3.0;
    d(4) = 2.1;  // Added noise to make it inconsistent
    
    cout << "Matrix O:" << endl;
    O.Print();
    
    cout << "Right-hand side vector d:" << endl;
    d.print();
    
    // Solve using Moore-Penrose pseudo-inverse (least squares solution)
    cout << "\nSolving using Moore-Penrose pseudo-inverse:" << endl;
    try {
        NonSquareLinSystem nsls1(O, d);
        Vector x1 = nsls1.Solve();
        
        cout << "Solution vector:" << endl;
        x1.print();
        
        // Verify solution - compute residual
        Vector residual1 = d - O * x1;
        cout << "Residual (should be small):" << endl;
        residual1.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in NonSquareLinSystem: " << e.what() << endl;
    }
    
    // Solve using Tikhonov regularization
    cout << "\nSolving using Tikhonov regularization (lambda=0.1):" << endl;
    try {
        NonSquareLinSystem nsls2(O, d, 0.1);
        Vector x2 = nsls2.Solve();
        
        cout << "Solution vector:" << endl;
        x2.print();
        
        // Verify solution - compute residual
        Vector residual2 = d - O * x2;
        cout << "Residual:" << endl;
        residual2.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in NonSquareLinSystem: " << e.what() << endl;
    }
    
    // Test under-determined system (more unknowns than equations)
    cout << "\nUnder-determined system (more unknowns than equations):" << endl;
    Matrix U(2, 4);  // 2 equations, 4 unknowns
    U(1, 1) = 1.0; U(1, 2) = 2.0; U(1, 3) = 3.0; U(1, 4) = 4.0;
    U(2, 1) = 5.0; U(2, 2) = 6.0; U(2, 3) = 7.0; U(2, 4) = 8.0;
    
    Vector e(2);
    e(1) = 5.0;
    e(2) = 10.0;
    
    cout << "Matrix U:" << endl;
    U.Print();
    
    cout << "Right-hand side vector e:" << endl;
    e.print();
    
    // Solve using Moore-Penrose pseudo-inverse (minimum norm solution)
    cout << "\nSolving using Moore-Penrose pseudo-inverse:" << endl;
    try {
        NonSquareLinSystem nsls3(U, e);
        Vector x3 = nsls3.Solve();
        
        cout << "Solution vector (minimum norm):" << endl;
        x3.print();
        
        // Verify solution - should satisfy Ux = e
        Vector result3 = U * x3;
        cout << "U*x (should equal e):" << endl;
        result3.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in NonSquareLinSystem: " << e.what() << endl;
    }
    
    // Solve using Tikhonov regularization
    cout << "\nSolving using Tikhonov regularization (lambda=0.1):" << endl;
    try {
        NonSquareLinSystem nsls4(U, e, 0.1);
        Vector x4 = nsls4.Solve();
        
        cout << "Solution vector:" << endl;
        x4.print();
        
        // Verify solution
        Vector result4 = U * x4;
        cout << "U*x (should be close to e):" << endl;
        result4.print();
    }
    catch (const std::exception& e) {
        cout << "Exception in NonSquareLinSystem: " << e.what() << endl;
    }
    
    return 0;
}
