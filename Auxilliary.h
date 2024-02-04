//
// Created by Eric Moss on 6/27/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H

#include <iostream>
#include <vector>
#include "Index.h"

class Auxilliary {
public:
    static int mod(int x, int modulus);
    static int next(long double x);
    /*static void computeTheta(acb_t theta, int d, int bits);
    static void intervalReduce(arb_struct *answer, const arb_struct *x, const arb_struct *scale, int bits);
    static void vectorReduce(acb_struct *answer, const acb_struct *v, const acb_t toReduce, int bits);
    static void acbPlusArb(acb_t answer, const acb_t z, const arb_t x, int bits);
    static void acbTimesArb(acb_t answer, const acb_t z, const arb_t x, int bits);
    static void acbDivArb(acb_t answer, const acb_t z, const arb_t x, int bits);
    static void acbDivArb(acb_t answer, const acb_t z, int x, int bits);
    static void traceProduct(acb_struct *answer, const acb_struct *z, const acb_struct *w, int bits);*/
    //static double computeVolumeOfFD(int d);

    static double imagTheta(int d);
    static double multiPrecisionSummation(const std::vector<double>& numbers);
    //static double pi;
};


#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H
