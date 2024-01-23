#include "ArchtKBessel.h"

#include <iostream>

int ArchtKBessel::kbes_pr = 53;

ArchtKBessel::ArchtKBessel() {
    f = init_kbessel(kbes_pr);

    mpfi_init2(mpfi_r, mpfi_get_prec(f));
    mpfi_init2(mpfi_x, mpfi_get_prec(f));
    mpfr_init(mid);
}

ArchtKBessel::~ArchtKBessel() {
    mpfi_clear(mpfi_r);
    mpfi_clear(mpfi_x);
    mpfr_clear(mid);
    clear_kbessel(f);
}

void ArchtKBessel::updateR(const double r) {
    mpfi_set_d(mpfi_r, r);
}

double ArchtKBessel::evaluate(const double &x) {
    mpfi_set_d(mpfi_x, x);
    kbessel(f, mpfi_r, mpfi_x);
    mpfi_mid(mid, f);
    double ans = mpfr_get_d(mid, MPFR_RNDN);
    return ans;
}