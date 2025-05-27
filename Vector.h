#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>

class Vector {
private:
    int mSize;    // Size of the vector
    double* mData; // Pointer to the data array

public:
    // Constructors and destructor
    Vector();                     // Default constructor
    Vector(int size);             // Constructor with size
    Vector(const Vector& v);      // Copy constructor
    ~Vector();                    // Destructor

    // Assignment operator
    Vector& operator=(const Vector& v);

    // Unary operators
    Vector operator-() const;   // Negation
    Vector operator+() const;   // Unary plus (returns copy)
    Vector operator++();        // Postfix increment
    Vector operator--();        // Postfix decrement

    // Binary operators with another vector
    Vector operator+(const Vector& v) const;
    Vector operator-(const Vector& v) const;
    double operator*(const Vector& v) const;  // Dot product

    // Binary operators with a scalar
    Vector operator*(double scalar) const;  // Vector * scalar
    Vector operator/(double scalar) const;  // Vector / scalar

    // Access operators
    double& operator[](int index);          // Zero-indexed access with bounds checking
    double operator[](int index) const;     // Const version
    double& operator()(int index);          // One-indexed access with bounds checking
    double operator()(int index) const;     // Const version    // Other useful methods
    int get_size() const {
        return mSize;  // Get size of vector
    }
    
    // Print the vector
    void print() const;
    // Friend functions for scalar operations
    friend Vector operator*(double scalar, const Vector& v);  // scalar * Vector
};

#endif // VECTOR_H
