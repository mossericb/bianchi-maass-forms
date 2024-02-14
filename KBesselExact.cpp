//
// Created by Eric Moss on 2/13/24.
//

#include "KBesselExact.h"
#include "omp.h"
#include <acb_hypgeom.h>

KBesselExact::KBesselExact(double r, int bitsOfPrecision) {
    prec = bitsOfPrecision;
    this->r = r;

#pragma omp parallel default(none)
    {
        #pragma omp single
        threads = omp_get_num_threads();
    }

    realTemp.resize(threads);
    arbR.resize(threads);
    besselScale.resize(threads);
    preComputedBesselScale.resize(threads);
    complexTemp.resize(threads);
    acbNu.resize(threads);

    for (int i = 0; i < threads; i++) {
        arb_init(&realTemp[i]);
        arb_init(&arbR[i]);
        arb_init(&besselScale[i]);
        arb_init(&preComputedBesselScale[i]);
        acb_init(&complexTemp[i]);
        acb_init(&acbNu[i]);
    }

    for (int i = 0; i < threads; i++) {
        arb_set_d(&arbR[i],r);
        acb_set_d_d(&acbNu[i], 0, r);

        arb_const_pi(&preComputedBesselScale[i], 2*prec);
        arb_mul(&preComputedBesselScale[i], &preComputedBesselScale[i], &arbR[i], 2*prec);
        arb_div_ui(&preComputedBesselScale[i], &preComputedBesselScale[i], 2, 2*prec);
        arb_exp(&preComputedBesselScale[i], &preComputedBesselScale[i], 2*prec);
    }

}

KBesselExact::~KBesselExact() {
    for (int i = 0; i < threads; i++) {
        arb_clear(&realTemp[i]);
        arb_clear(&arbR[i]);
        arb_clear(&besselScale[i]);
        arb_clear(&preComputedBesselScale[i]);
        acb_clear(&complexTemp[i]);
        acb_clear(&acbNu[i]);
    }

    flint_cleanup_master();
}

double KBesselExact::exactKBessel(const double x) {
    if (x == 0) {
        return 0;
    }

    int threadNum = omp_get_thread_num();

    acb_set_d_d(&complexTemp[threadNum], x, 0);
    //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    //probably this precision will be enough, if not...
    int tryOnceBits = 2 * prec;

    acb_hypgeom_bessel_k(&complexTemp[threadNum],
                         &acbNu[threadNum],
                         &complexTemp[threadNum],
                         tryOnceBits);
    //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    acb_get_real(&realTemp[threadNum],&complexTemp[threadNum]);
    //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    arb_mul(&realTemp[threadNum],
            &realTemp[threadNum],
            &preComputedBesselScale[threadNum],
            tryOnceBits);
    //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    long accuracy = arb_rel_accuracy_bits(&realTemp[threadNum]);
    if (accuracy >= prec || accuracy >= 53) {
        double answer;
        try {
            answer = std::stod(arb_get_str(&realTemp[threadNum], 50, ARB_STR_NO_RADIUS));
        } catch (...) {
            answer = 0;
        }
        return answer;
    } else {
        //...then it enters this loop which will run until it reaches the desired precision
        for (int bits = 4 * prec; ; bits *= 2) {
            acb_set_d_d(&complexTemp[threadNum], x, 0);

            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_const_pi(&besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_mul(&besselScale[threadNum], &arbR[threadNum], &besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_div_ui(&besselScale[threadNum], &besselScale[threadNum], 2, bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_exp(&besselScale[threadNum], &besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            acb_hypgeom_bessel_k(&complexTemp[threadNum],
                                 &acbNu[threadNum],
                                 &complexTemp[threadNum],
                                 bits);
            //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            acb_get_real(&realTemp[threadNum],&complexTemp[threadNum]);
            //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_mul(&realTemp[threadNum],
                    &realTemp[threadNum],
                    &besselScale[threadNum],
                    bits);
            //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            accuracy = arb_rel_accuracy_bits(&realTemp[threadNum]);
            if (accuracy >= prec || accuracy >= 53) {
                break;
            }
        }
    }

    double answer;
    try {
        answer = std::stod(arb_get_str(&realTemp[threadNum], 50, ARB_STR_NO_RADIUS));
    } catch (...) {
        answer = 0;
    }

    return answer;
}


