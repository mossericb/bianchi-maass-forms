#include "ArchtKBessel.h"
#include "archt/kbessel.c"

#include <iostream>
#include <numbers>

ArchtKBessel::ArchtKBessel() {
    mpfi_ptr bess = init_kbessel(BEGINNING_BITS);
    vec_f.push_back(bess);
    r = 0;

    mpfi_t mpfi_r;
    mpfi_init2(mpfi_r, BEGINNING_BITS);
    vec_mpfi_r.push_back(*mpfi_r);

    mpfi_t mpfi_x;
    mpfi_init2(mpfi_x, BEGINNING_BITS);
    vec_mpfi_x.push_back(*mpfi_x);

    mpfr_t mid;
    mpfr_init2(mid, BEGINNING_BITS);
    vec_mid.push_back(*mid);
}

ArchtKBessel::~ArchtKBessel() {
    for (auto bess : vec_f) {
        clear_kbessel(bess);
    }
    vec_f.clear();

    for (auto mp_r : vec_mpfi_r) {
        mpfi_clear(&mp_r);
    }
    vec_mpfi_r.clear();

    for (auto mp_x : vec_mpfi_x) {
        mpfi_clear(&mp_x);
    }
    vec_mpfi_x.clear();

    for (auto mid : vec_mid) {
        mpfr_clear(&mid);
    }
    vec_mid.clear();
}

void ArchtKBessel::setR(const double &r) {
    this->r = r;
    for (int i = 0; i < vec_mpfi_r.size(); i++) {
        mpfr_set_d(&(vec_mpfi_r[i].left), r, MPFR_RNDN);
        mpfr_set_d(&(vec_mpfi_r[i].right), r, MPFR_RNDN);
    }
    zeroCutoff = (1136 + PI*r/2.0*log2(E) - 0.5*log2(E) + 0.5*log2(PI/2))/log2(E);
}

double ArchtKBessel::evaluate(const double &x) {
    if (x < 2.0 * r) {
        mpfr_set_d(&(vec_mpfi_x[0].left), x, MPFR_RNDN);
        mpfr_set_d(&(vec_mpfi_x[0].right), x, MPFR_RNDN);

        kbessel(vec_f[0], &vec_mpfi_r[0], &vec_mpfi_x[0]);
        mpfi_mid(&vec_mid[0], vec_f[0]);
        double ans = mpfr_get_d(&vec_mid[0], MPFR_RNDN);
        return ans;
    } else if (x > zeroCutoff) {
        return 0;
    }

    //Using the asymptotic formula for the K-Bessel function, we can expect to need
    //this many absolute bits to express the value to relative accuracy with 53 bits
    //The formula is int bits = std::ceil(-log2(sqrt(PI / (2.0 * x)) * exp(PI * r / 2.0 - x))) + 53;
    //but terms in the formula go beyond double precision pretty quickly, so we have to break it up
    int bits = std::ceil( -0.5*(log2(PI/2.0) - log2(x)) - (PI*r/2.0 - x)*log2(E) + 60);
    int bitsIndex = std::max(0, (int)std::ceil(log2(((double)bits)/BEGINNING_BITS)));

    for (int i = vec_f.size(); i <= bitsIndex; i++) {
        mpfi_ptr bess = init_kbessel(pow(2,i) * BEGINNING_BITS);
        vec_f.push_back(bess);

        mpfi_t mpfi_r;
        mpfi_init2(mpfi_r, pow(2,i) * BEGINNING_BITS);
        mpfr_set_d(&(mpfi_r->left), r, MPFR_RNDN);
        mpfr_set_d(&(mpfi_r->right), r, MPFR_RNDN);
        vec_mpfi_r.push_back(*mpfi_r);

        mpfi_t mpfi_x;
        mpfi_init2(mpfi_x, pow(2,i) * BEGINNING_BITS);
        vec_mpfi_x.push_back(*mpfi_x);

        mpfr_t mid;
        mpfr_init2(mid, pow(2,i) * BEGINNING_BITS);
        vec_mid.push_back(*mid);
    }

    mpfr_set_d(&(vec_mpfi_x[bitsIndex].left), x, MPFR_RNDN);
    mpfr_set_d(&(vec_mpfi_x[bitsIndex].right), x, MPFR_RNDN);

    kbessel(vec_f[bitsIndex], &vec_mpfi_r[bitsIndex], &vec_mpfi_x[bitsIndex]);
    //archt_out(stdout, vec_f[bitsIndex]);
    mpfi_mid(&vec_mid[bitsIndex], vec_f[bitsIndex]);
    double ans = mpfr_get_d(&vec_mid[bitsIndex], MPFR_RNDN);

    return ans;
}