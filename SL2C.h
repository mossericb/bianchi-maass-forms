//
// Created by Eric Moss on 9/27/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_SL2C_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_SL2C_H

#include <complex>
#include <iostream>

class SL2C {
public:
    template<typename T1, typename T2, typename T3, typename T4>
    SL2C(T1 a, T2 b, T3 c, T4 d);

    std::complex<double> a, b, c, d;
private:
    template<typename T>
    void set(T &n, std::complex<double>& entry);
    void set(double& x, std::complex<double>& entry);
    void set(std::complex<double>& z, std::complex<double>& entry);
    friend std::ostream& operator<<(std::ostream&, const SL2C&);
};

template<typename T1, typename T2, typename T3, typename T4>
SL2C::SL2C(T1 a, T2 b, T3 c, T4 d) {
    set(a, this->a);
    set(b, this->b);
    set(c, this->c);
    set(d, this->d);
}

template<typename T>
void SL2C::set(T &n, std::complex<double> &entry) {
    entry = {(double)n, 0.0};
}

#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_SL2C_H
