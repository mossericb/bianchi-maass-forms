//
// Created by Eric Moss on 6/26/23.
//

#include "Index.h"
#include <complex>
#include <cmath>
#include <algorithm>
#include "EricNT.h"

using namespace std;

Index::Index(int a, int b, int d) {
    this->a = a;
    this->b = b;
    this->d = d;

    std::complex<double> theta;
    //Initialize std::complex<double> theta
    if (EricNT::mod(-d,4) == 1) {
        theta = {1.0/2, sqrt(d)/2};
    } else if (EricNT:: mod(-d, 4) == 0) {
        throw(std::invalid_argument("d should be squarefree"));
    } else {
        theta = {0, sqrt(d)};
    }

    //Initialize std::complex<double> complex
    complex = std::complex<double> {(double)b,0};
    complex *= theta;
    complex += std::complex<double> {(double)a, 0};

    //Initialize double abs
    abs = std::abs(complex);

    //Initialize double angle
    angle = atan2(complex.imag(), complex.real());
    if (angle < 0) {
        angle += 2 * 3.141592653589793238462643383279502884197;
    }
}

std::ostream &operator<<(std::ostream &strm, const Index &index) {
    std::stringstream out;
    out << to_string(index.a) << " + " << to_string(index.b) << "*theta";
    return strm << std::move(out).str();
}

Index Index::rotate() const {
    if (d == 1) {
        return {-b, a, d};
    } else if (d == 3) {
        return {-b, a + b, d};
    } else {
        return {-a, -b, d};
    }
}

Index Index::conj() const {
    if (EricNT::mod(-d,4) == 1) {
        return {a + b, -b, d};
    } else {
        return {a, -b, d};
    }
}


bool operator<(const Index &index1, const Index &index2) {
    if (index1.abs < index2.abs) {
        return true;
    } else if (index1.abs == index2.abs && index1.angle < index2.angle) {
        return true;
    } else {
        return false;
    }
}

bool operator==(const Index &index1, const Index &index2) {
    return (index1.a == index2.a) && (index1.b == index2.b) && (index1.d == index2.d);
}

bool operator!=(const Index &index1, const Index &index2) {
    return !(index1 == index2);
}

Index operator-(const Index& lhs, const Index& rhs) {
    return Index(lhs.a - rhs.a, lhs.b - rhs.b, lhs.d);
}

string Index::to_string() const {
    string message = std::to_string(a) + " + " + std::to_string(b) + "*theta";
    return message;
}

const complex<double> &Index::getComplex() const {
    return complex;
}

double Index::getAbs() const {
    return abs;
}

double Index::getAngle() const {
    return angle;
}

int Index::getA() const {
    return a;
}

int Index::getB() const {
    return b;
}

int Index::getD() const {
    return d;
}


