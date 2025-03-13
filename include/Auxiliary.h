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

    static long mod(long x, long modulus);

    static double kahanSummation(std::vector<double> &numbers);
    static double kahanSummationUnvectorized(std::vector<double> &numbers);
    static char legendreSymbol(long n, long p);
    static double rand(double a, double b);

    static string commonRound(string & s1, string & s2);
    static double roundN(const double r, int n);

private:
    int threads = 0;
};