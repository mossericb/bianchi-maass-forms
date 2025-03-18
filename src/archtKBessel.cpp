#include "../include/archtKBessel.h"
#include "archt/kbessel.c"


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
        double bitsGuess = std::ceil( -0.5*(log2(PI/2.0) - log2(x)) - (PI*r/2.0 - x)*log2(E) + 60);
        int indexGuess = std::max(0, (int)std::ceil(log2(((double)bitsGuess)/BEGINNING_BITS)));
        indexGuess = std::min(indexGuess, (int)vec_f.size() - 1);
        for (int i = indexGuess; i < vec_f.size(); i++) {
            mpfi_set_d(&(vec_mpfi_x[i]), x);

            kbessel(vec_f[i], &vec_mpfi_r[i], &vec_mpfi_x[i]);

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
            vec_mpfi_r.push_back(*mpfi_r);
            mpfi_set_d(&vec_mpfi_r[i], r);

            mpfi_t mpfi_x;
            mpfi_init2(mpfi_x, pow(2,i) * BEGINNING_BITS);
            vec_mpfi_x.push_back(*mpfi_x);
            mpfi_set_d(&vec_mpfi_x[i], x);

            kbessel(vec_f[i], &vec_mpfi_r[i], &vec_mpfi_x[i]);

            double ans1 = mpfr_get_d(&(vec_f[i]->left), MPFR_RNDN);
            double ans2 = mpfr_get_d(&(vec_f[i]->right), MPFR_RNDN);
            if (ans1 == ans2) {
                return ans1;
            }

            ++i;
        }
    }
}