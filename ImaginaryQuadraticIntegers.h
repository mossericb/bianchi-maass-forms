//
// Created by Eric Moss on 2/3/24.
//

#ifndef BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
#define BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
#include <complex>

class ImaginaryQuadraticIntegers {
public:
    explicit ImaginaryQuadraticIntegers(int d);
    ImaginaryQuadraticIntegers();

    [[nodiscard]] int getd() const { return d; };
    [[nodiscard]] std::complex<double> getTheta() const { return theta; };
    [[nodiscard]] double getA() const { return A; };
    [[nodiscard]] double getY0() const { return Y0; };

private:
    int d;
    std::complex<double> theta;
    double A;
    double Y0;
};


#endif //BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
