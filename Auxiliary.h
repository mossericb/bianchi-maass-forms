#pragma once

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

    static double imagTheta(int d);
    double multiPrecisionSummation(const std::vector<double>& numbers);
    static double kahanSummation(std::vector<double> &numbers);
    static double kahanSummationUnvectorized(std::vector<double> &numbers);
    static double nextWithinOrderOfMag(double x);
    static char legendreSymbol(long n, long p);
    static double rand(double a, double b);

    static int SIMD_DOUBLE_WIDTH;
    static string commonRound(string & s1, string & s2);
    static double roundN(const double r, int n);

private:
    std::vector<__mpfr_struct> sums;
    int threads = 0;
};