//
// Created by Eric Moss on 6/22/24.
//

#ifndef BIANCHI_MAASS_FORMS_KBESSELREAL_H
#define BIANCHI_MAASS_FORMS_KBESSELREAL_H
#include "arb.h"

class KBesselReal {
public:
    explicit KBesselReal(double lambda);
    ~KBesselReal();

    void setLambda(double lambda);
    double evaluate(double x);

private:
    static constexpr double PI = 3.141592653589793238462643383279;
    static constexpr double E = 2.7182818284590452353602874713;
    static constexpr int BEGINNING_BITS = 54;
    static constexpr int RELATIVE_ACC_BITS = 53;

    double lambda;
    double r;

    arb_struct argument;
    arb_struct order;
    arb_struct answer;
};


#endif //BIANCHI_MAASS_FORMS_KBESSELREAL_H
