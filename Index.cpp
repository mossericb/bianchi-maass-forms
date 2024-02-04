//
// Created by Eric Moss on 6/26/23.
//

#include "Index.h"
#include <complex>
#include <cmath>
#include <algorithm>
#include "Auxilliary.h"

using namespace std;

Index::Index(int a, int b) {
    this->a = a;
    this->b = b;
}

std::ostream &operator<<(std::ostream &strm, const Index &index) {
    std::stringstream out;
    out << to_string(index.a) << " + " << to_string(index.b) << "*theta";
    return strm << std::move(out).str();
}

Index Index::rotate(int d) const {
    if (d == 1) {
        return {-b, a};
    } else if (d == 3) {
        return {-b, a + b};
    } else {
        return {-a, -b};
    }
}

Index Index::conj(int d) const {
    if (Auxilliary::mod(-d, 4) == 1) {
        return {a + b, -b};
    } else {
        return {a, -b};
    }
}

bool operator==(const Index &index1, const Index &index2) {
    return (index1.a == index2.a) && (index1.b == index2.b);
}

bool operator!=(const Index &index1, const Index &index2) {
    return !(index1 == index2);
}

Index operator-(const Index& lhs, const Index& rhs) {
    return Index(lhs.a - rhs.a, lhs.b - rhs.b);
}

bool operator<(const Index& index1, const Index& index2) {
    if (index1.a < index2.a) {
        return true;
    } else if (index1.a == index2.a && index1.b < index2.b) {
        return true;
    } else {
        return false;
    }
}

string Index::to_string() const {
    string message = std::to_string(a) + " + " + std::to_string(b) + "*theta";
    return message;
}

std::complex<double> Index::getComplex(int d) const {
    std::complex<double> theta;
    //Initialize std::complex<double> theta
    if (Auxilliary::mod(-d, 4) == 1) {
        theta = {1.0/2, sqrt(d)/2};
    } else if (Auxilliary:: mod(-d, 4) == 0) {
        throw(std::invalid_argument("d should be squarefree"));
    } else {
        theta = {0, sqrt(d)};
    }
    return (double)a + ((double)b)*theta;
}

double Index::getAbs(int d) const {
    return abs(this->getComplex(d));
}

double Index::getAngle(int d) const {
    std::complex<double> complex = this->getComplex(d);
    double angle = atan2(complex.imag(), complex.real());
    if (angle < 0) {
        angle += 6.28318530717958647692528676655900576839;
    }
    return angle;
}

int Index::getA() const {
    return a;
}

int Index::getB() const {
    return b;
}


