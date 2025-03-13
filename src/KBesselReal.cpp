#include "../include/KBesselReal.h"
#include "arb_hypgeom.h"

KBesselReal::KBesselReal(double r) {
    this->r = r;

    arb_init(&argument);
    arb_init(&order);
    arb_init(&answer);

    arb_set_d(&order, r);
}

KBesselReal::~KBesselReal() {
    arb_clear(&argument);
    arb_clear(&order);
    arb_clear(&answer);

    flint_cleanup();
}

double KBesselReal::evaluate(double x) {
    arb_set_d(&argument,x);

    int i = 1;

    while (true) {
        arb_hypgeom_bessel_k(&answer, &order, &argument, BEGINNING_BITS * pow(2,i));
        if (arb_rel_accuracy_bits(&answer) >= RELATIVE_ACC_BITS) {
            try {
                auto str= arb_get_str(&answer, 20, ARB_STR_NO_RADIUS);
                double result = atof(str);
                return result;
            } catch (...) {
                return 0;
            }
        }
        i++;
    }
}

void KBesselReal::setR(double r) {
    this->r = r;
    arb_set_d(&order, r);
}