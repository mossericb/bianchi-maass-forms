//
// Created by Eric Moss on 2/13/24.
//

#include "KBesselExact.h"
#include "omp.h"
#include <boost/math/differentiation/finite_difference.hpp>

KBesselExact::KBesselExact(double r) {

#pragma omp parallel default(none) shared(threads)
    {
#pragma omp single
        threads = omp_get_max_threads();
    }

    K.reserve(threads);
    for (int i = 0; i < threads; i++) {
        K.emplace_back();
        K[i].setR(r);
    }
}

KBesselExact::KBesselExact() : KBesselExact(0) {
}

double KBesselExact::exactKBessel(const double x) {
    int threadNum = omp_get_thread_num();

    return K[threadNum].evaluate(x);
}

double KBesselExact::estimateDerivativeKBessel(const double x) {
    auto f = [this](double x) { return this->exactKBessel(x); };
    double dfdx = boost::math::differentiation::finite_difference_derivative(f, x);

    return dfdx;
}

void KBesselExact::updateR(const double r) {
    for (auto a : K) {
        a.setR(r);
    }
}


