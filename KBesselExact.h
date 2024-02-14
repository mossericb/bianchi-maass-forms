//
// Created by Eric Moss on 2/13/24.
//

#ifndef BIANCHI_MAASS_FORMS_KBESSELEXACT_H
#define BIANCHI_MAASS_FORMS_KBESSELEXACT_H

#include <arb.h>
#include <acb.h>
#include <vector>

class KBesselExact {
public:
    KBesselExact(double r, int bitsOfPrecision = 53);
    ~KBesselExact();

    double exactKBessel(const double x);

private:
    int threads;
    int prec;

    double r;

    std::vector<arb_struct> realTemp;
    std::vector<arb_struct> arbR;
    std::vector<arb_struct> besselScale;
    std::vector<arb_struct> preComputedBesselScale;
    std::vector<acb_struct> complexTemp;
    std::vector<acb_struct> acbNu;
};


#endif //BIANCHI_MAASS_FORMS_KBESSELEXACT_H
