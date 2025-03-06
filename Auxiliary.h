//
// Created by Eric Moss on 6/27/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H

#include <iostream>
#include <vector>
#include "Index.h"
#include "mpfr.h"
#include <boost/align/aligned_allocator.hpp>
#include <string>

using std::string;

class Auxiliary {
public:
    Auxiliary();
    ~Auxiliary();

    static long mod(long x, long modulus);
    static int next(double x);
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
    double multiPrecisionSummation(const std::vector<double>& numbers);
    static double kahanSummation(std::vector<double> &numbers);
    static double kahanSummationUnvectorized(std::vector<double> &numbers);
    static double nextWithinOrderOfMag(double x);
    static char legendreSymbol(long n, long p);
    static double rand(double a, double b);
    //static double pi;

    static int SIMD_DOUBLE_WIDTH;
    static string commonRound(string & s1, string & s2);
    static double roundN(const double r, int n);

private:
    std::vector<__mpfr_struct> sums;
    int threads = 0;
};


#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_AUXILLIARYMETHODS_H
