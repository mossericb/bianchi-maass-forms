//
// Created by Eric Moss on 2/3/24.
//

#include "ImaginaryQuadraticIntegers.h"
#include "Auxilliary.h"

ImaginaryQuadraticIntegers::ImaginaryQuadraticIntegers(int d) {
    this->d = d;

    std::complex<double> theta = {0,0};
    if (Auxilliary::mod(-d, 4) == 1) {
        //theta = 1/2 + I*sqrt(d)/2
        theta = std::complex<double> {1.0/2, sqrt(d)/2};
    } else if (Auxilliary::mod(-d, 4) == 0) {
        throw std::invalid_argument("d is incorrect");
    } else {
        //theta = I*sqrt(d)
        theta = std::complex<double> {0, sqrt(d)};
    }

    A = theta.imag();

    //These results are hardcoded output from John Cremona's bianchi-progs
    if (d == 3) {
        Y0 = sqrt(2.0/3);
    } else if (d == 19) {
        Y0 = sqrt(2.0/19);
    } else if (d == 43) {
        Y0 = sqrt(2.0/43);
    } else if (d == 67) {
        Y0 = sqrt(2.0/67);
    } else if (d == 163) {
        Y0 = sqrt(2.0/163);
    } else {
        //Y0 = sqrt(1- (1/2)^2 - (theta.imag()/2)^2)
        //=sqrt(3/4 - theta.imag()^2/4)
        Y0 = .75 - pow(this->getTheta().imag(),2)/4;
        Y0 = sqrt(Y0);
    }
}

ImaginaryQuadraticIntegers::ImaginaryQuadraticIntegers() {
    d = 0;
    A = 0;
    theta = {0,0};
    Y0 = 0;
}


