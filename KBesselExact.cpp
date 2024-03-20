//
// Created by Eric Moss on 2/13/24.
//

#include "KBesselExact.h"
#include "omp.h"
#include <acb_hypgeom.h>
#include <boost/math/differentiation/finite_difference.hpp>

KBesselExact::KBesselExact(double r, int bitsOfPrecision) {
    prec = bitsOfPrecision;
    this->r = r;

#pragma omp parallel default(none) shared(threads)
    {
#pragma omp single
        threads = omp_get_max_threads();
    };

    K.reserve(threads);
    for (int i = 0; i < threads; i++) {
        K.push_back(new ArchtKBessel());
        K[i]->setR(r);
    }
}

double KBesselExact::exactKBessel(const double x) {
    int threadNum = omp_get_thread_num();

    return K[threadNum]->evaluate(x);
}

KBesselExact::~KBesselExact() {
    for (int i = 0; i < threads; i++) {
        delete K[i];
    }
    K.clear();
}

double KBesselExact::estimateDerivativeKBessel(const double x) {
    auto f = [this](double x) { return this->exactKBessel(x); };
    double dfdx = boost::math::differentiation::finite_difference_derivative(f, x);

    return dfdx;
}


