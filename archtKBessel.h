#pragma once

#include <mpfr.h>
#include <mpfi.h>
#include <vector>

class archtKBessel {
public:
    explicit archtKBessel(double r);
    ~archtKBessel();

    void setR(double r);
    double evaluate(double x);
private:
    static constexpr double PI = 3.141592653589793238462643383279;
    static constexpr double E = 2.7182818284590452353602874713;
    static constexpr int BEGINNING_BITS = 54;
    double zeroCutoff;

    mpfr_t acc;

    double r;

    std::vector<mpfi_ptr> vec_f;
    std::vector<__mpfi_struct> vec_mpfi_r;
    std::vector<__mpfi_struct> vec_mpfi_x;
};