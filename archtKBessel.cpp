#include "archtKBessel.h"
#include "archt/kbessel.c"

#include <iostream>

archtKBessel::archtKBessel(double r) {
    mpfr_init_set_d(acc, 1.0e-16, MPFR_RNDN);

    mpfi_ptr bess = init_kbessel(BEGINNING_BITS);
    vec_f.push_back(bess);

    mpfi_t mpfi_r;
    mpfi_init2(mpfi_r, BEGINNING_BITS);
    vec_mpfi_r.push_back(*mpfi_r);

    mpfi_t mpfi_x;
    mpfi_init2(mpfi_x, BEGINNING_BITS);
    vec_mpfi_x.push_back(*mpfi_x);

    mpfr_t mid;
    mpfr_init2(mid, BEGINNING_BITS);
    vec_mid.push_back(*mid);

    setR(r);
}

archtKBessel::~archtKBessel() {
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

    mpfr_free_cache();
}

void archtKBessel::setR(double r) {
    this->r = r;
    for (int i = 0; i < vec_mpfi_r.size(); i++) {
        mpfi_init_set_d(&(vec_mpfi_r[i]), r);
    }
    zeroCutoff = (1136 + PI*r/2.0*log2(E) - 0.5*log2(E) + 0.5*log2(PI/2))/log2(E);
}

double archtKBessel::evaluate(double x) {
    if (x > zeroCutoff) {
        return 0.0;
    } else {
        for (int i = 0; i < vec_f.size(); i++) {
            mpfi_set_d(&(vec_mpfi_x[i]), x);

            kbessel(vec_f[i], &vec_mpfi_r[i], &vec_mpfi_x[i]);
            /*
            mpfi_diam(&vec_mid[i], vec_f[i]);
            if (mpfr_lessequal_p(&vec_mid[i], acc)) {
                mpfi_mid(&vec_mid[i], vec_f[i]);
                double ans = mpfr_get_d(&vec_mid[i], MPFR_RNDN);
                return ans;
            }*/

            double ans1 = mpfr_get_d(&(vec_f[i]->left), MPFR_RNDN);
            double ans2 = mpfr_get_d(&(vec_f[i]->right), MPFR_RNDN);
            if (ans1 == ans2) {
                return ans1;
            }
        }

        /*At this point, none of the kbessel structs have enough precision
         * Make more
         */
        int i = vec_f.size();
        while (true) {
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

            mpfi_set_d(&(vec_mpfi_x[i]), x);

            kbessel(vec_f[i], &vec_mpfi_r[i], &vec_mpfi_x[i]);
            mpfi_diam(&vec_mid[i], vec_f[i]);
            if (mpfr_lessequal_p(&vec_mid[i], acc)) {
                mpfi_mid(&vec_mid[i], vec_f[i]);
                double ans = mpfr_get_d(&vec_mid[i], MPFR_RNDN);
                return ans;
            }
            ++i;
        }
    }
}