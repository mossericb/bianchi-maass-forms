#include "../include/SL2C.h"

void SL2C::set(double &x, std::complex<double> &entry) {
    entry = {x, 0.0};
}

void SL2C::set(std::complex<double> &z, std::complex<double> &entry) {
    entry = z;
}

std::ostream &operator<<(std::ostream &strm, const SL2C &gamma) {
    std::stringstream out;
    out << gamma.a << ", " << gamma.b << ", " << gamma.c << ", " << gamma.d;
    return strm << std::move(out).str();
}