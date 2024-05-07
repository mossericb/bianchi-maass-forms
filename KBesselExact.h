//
// Created by Eric Moss on 2/13/24.
//

#ifndef BIANCHI_MAASS_FORMS_KBESSELEXACT_H
#define BIANCHI_MAASS_FORMS_KBESSELEXACT_H

#include <ArchtKBessel.h>
#include <vector>

class KBesselExact {
public:
    explicit KBesselExact(double r);
    KBesselExact();

    double exactKBessel(const double x);
    double estimateDerivativeKBessel(const double x);
    void updateR(const double r);

private:
    int threads;

    std::vector<ArchtKBessel> K;
};


#endif //BIANCHI_MAASS_FORMS_KBESSELEXACT_H
