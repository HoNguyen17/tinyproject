#include "Vector.h"
#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

// Default constructor
Vector::Vector() : mSize(0), mData(nullptr) {
}

// Constructor with size
Vector::Vector(int size) : mSize(size) {
    if (size < 0) {
        throw invalid_argument("Vector size must be non-negative");
    }
    
    if (size > 0) {
        mData = new double[size];
        // Initialize all elements to zero
        for (int i = 0; i < size; i++) {
            mData[i] = 0.0;
        }
    } else {
        mData = nullptr;
    }
}

// Copy constructor
Vector::Vector(const Vector& v) : mSize(v.mSize) {
    if (mSize > 0) {
        mData = new double[mSize];
        // Copy elements
        for (int i = 0; i < mSize; i++) {
            mData[i] = v.mData[i];
        }
    } else {
        mData = nullptr;
    }
}

// Destructor
Vector::~Vector() {
    delete[] mData;
    mData = nullptr;  
}

// Assignment operator
Vector& Vector::operator=(const Vector& v) {
    // Check for self-assignment
    if (this == &v) {
        return *this;
    }
    
    // Delete old data
    delete[] mData;
    
    // Copy size
    mSize = v.mSize;
    
    // Allocate new memory if needed
    if (mSize > 0) {
        mData = new double[mSize];
        // Copy elements
        for (int i = 0; i < mSize; i++) {
            mData[i] = v.mData[i];
        }
    } else {
        mData = nullptr;
    }
    
    return *this;
}

// Unary operator - (negation)
Vector Vector::operator-() const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = -mData[i];
    }
    return result;
}

// Unary operator + (returns a copy)
Vector Vector::operator+() const {
    return Vector(*this);
}

// Postfix increment
Vector Vector::operator++() {
    for (int i = 0; i < mSize; i++) {
        mData[i]++;
    }
    return *this;
}

// Postfix decrement
Vector Vector::operator--() {
    for (int i = 0; i < mSize; i++) {
        mData[i]--;
    }
    return *this;
}

// Vector addition
Vector Vector::operator+(const Vector& v) const {
    if (mSize != v.mSize) {
        throw invalid_argument("Vectors must be of the same size for addition");
    }
    
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] + v.mData[i];
    }
    return result;
}

// Vector subtraction
Vector Vector::operator-(const Vector& v) const {
    if (mSize != v.mSize) {
        throw invalid_argument("Vectors must be of the same size for subtraction");
    }
    
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] - v.mData[i];
    }
    return result;
}

// Vector dot product
double Vector::operator*(const Vector& v) const {
    if (mSize != v.mSize) {
        throw invalid_argument("Vectors must be of the same size for dot product");
    }
    
    double result = 0.0;
    for (int i = 0; i < mSize; i++) {
        result += mData[i] * v.mData[i];
    }
    return result;
}

// Vector multiplication by a scalar
Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] * scalar;
    }
    return result;
}

// Vector division by a scalar
Vector Vector::operator/(double scalar) const {
    if (scalar == 0.0) {
        throw invalid_argument("Division by zero is not allowed");
    }
    
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] / scalar;
    }
    return result;
}

// Zero-indexed access with bounds checking (non-const version)
double& Vector::operator[](int index) {
    if (index < 0 || index >= mSize) {
        throw out_of_range("Vector index out of range");
    }
    return mData[index];
}

// Zero-indexed access with bounds checking (const version)
double Vector::operator[](int index) const {
    if (index < 0 || index >= mSize) {
        throw out_of_range("Vector index out of range");
    }
    return mData[index];
}

// One-indexed access with bounds checking (non-const version)
double& Vector::operator()(int index) {
    if (index < 1 || index > mSize) {
        throw out_of_range("Vector index out of range");
    }
    return mData[index - 1];
}

// One-indexed access with bounds checking (const version)
double Vector::operator()(int index) const {
    if (index < 1 || index > mSize) {
        throw out_of_range("Vector index out of range");
    }
    return mData[index - 1];
}

// Print the vector
void Vector::print() const {
    cout << "[ ";
    for (int i = 0; i < mSize; i++) {
        cout << mData[i];
        if (i < mSize - 1) {
            cout << ", ";
        }
    }
    cout << " ]" << endl;
}

// Friend function for scalar * Vector
Vector operator*(double scalar, const Vector& v) {
    return v * scalar;  // Reuse the vector * scalar operation
}
