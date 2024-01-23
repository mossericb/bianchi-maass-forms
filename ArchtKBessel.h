#ifndef ARCHTKBESSEL_LIBRARY_H
#define ARCHTKBESSEL_LIBRARY_H

#include "archt/kbessel.c"

class ArchtKBessel {
public:
    ArchtKBessel();
    ~ArchtKBessel();

    void updateR(const double r);
    double evaluate(const double &x);
private:
    static int kbes_pr;
    mpfi_ptr f;
    mpfi_t mpfi_r, mpfi_x;
    mpfr_t mid;
};

#endif //ARCHTKBESSEL_LIBRARY_H
