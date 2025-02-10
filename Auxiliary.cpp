//
// Created by Eric Moss on 6/27/23.
//

#include <stdexcept>

#include "Auxiliary.h"
//#include <math.h>
//#include "pari/pari.h"
//#include "flint/flint.h"
//#include "BianchiMaassComputer.h"
#include "omp.h"

#include <algorithm>
#include <chrono>
#include <random>

int Auxiliary::SIMD_DOUBLE_WIDTH = 2;

Auxiliary::Auxiliary() {

#pragma omp parallel default(none)
    {
        #pragma omp single
        threads = omp_get_num_threads();
    }

    sums.resize(threads);

    for (int i = 0; i < threads; i++) {
        mpfr_init2(&sums[i], 1076);
    }
}

Auxiliary::~Auxiliary() {
    for (int i = 0; i < threads; i++) {
        mpfr_clear(&sums[i]);
    }
    mpfr_free_cache();
    mpfr_mp_memory_cleanup();
}



long Auxiliary::mod(long x, long modulus) {
    long answer = x % modulus;
    if (answer < 0) {
        return answer + abs(modulus);
    } else {
        return answer;
    }
}

int Auxiliary::next(double x) {
    int answer = floor(x) + 1;
    return answer;
}

/*void Auxiliary::computeTheta(acb_t theta, int d, int bits) {
    acb_init(theta);
    if (Auxiliary::mod(-d,4) == 1) {
        //theta = 1/2 + I*sqrt(d)/2
        arb_t real, imag;
        arb_init(real);
        arb_init(imag);

        arb_set_si(real,1);
        arb_div_si(real,real,2, bits);

        arb_sqrt_ui(imag, d, bits);
        arb_div_si(imag, imag, 2, bits);

        acb_set_arb_arb(theta, real, imag);

        arb_clear(real);
        arb_clear(imag);
    } else if (Auxiliary::mod(-d,4) == 2 || Auxiliary::mod(-d,4) == 3) {
        //theta = I*sqrt(d)
        arb_t real, imag;
        arb_init(real);
        arb_init(imag);

        arb_zero(real);

        arb_sqrt_ui(imag, d, bits);

        acb_set_arb_arb(theta, real, imag);

        arb_clear(real);
        arb_clear(imag);
    } else {
        throw std::invalid_argument("d is incorrect");
    }
}

*//*
 * If the interval [-scale, scale] centered at 0 is translated to the rest of the real line,
 * where in the interval does x fall?
 * The answer is set to the variable answer.
 *
 * If x is positive, the right endpoint is favored.
 * If x is negative, the left endpoint is favored.
 * If x overlaps with an endpoint, we try to deal with it in a natural way, but there is no guarantee
 * that answer is a ball totally within [-scale, scale].
 *//*
void Auxiliary::intervalReduce(arb_struct *answer, const arb_struct *x, const arb_struct *scale, const int bits) {
    arb_t localCopyOfX;
    arb_init(localCopyOfX);
    arb_set(localCopyOfX,x);

    arb_t copyOfScale;
    arb_init(copyOfScale);
    arb_set(copyOfScale, scale);

    arb_t negScale;
    arb_init(negScale);
    arb_neg(negScale, scale);

    arb_t temp;
    arb_init(temp);

    if (arb_gt(localCopyOfX, copyOfScale)) {
        //use g = x - 2*scale*ceil(x/(2*scale) - 1/2)
        arb_sub(temp, localCopyOfX, copyOfScale, bits);
        arb_div_si(temp, temp, 2, bits);
        arb_div(temp, temp, copyOfScale, bits);
        arb_ceil(temp, temp, bits);
        arb_mul(temp, temp, copyOfScale, bits);
        arb_mul_si(temp,temp, 2, bits);
        arb_sub(answer, localCopyOfX, temp, bits);
    } else if (arb_lt(localCopyOfX, negScale)) {
        //use f = =x - 2*scale*floor(x/(2*scale) + 1/2)
        arb_add(temp, localCopyOfX, copyOfScale, bits);
        arb_div_si(temp, temp, 2, bits);
        arb_div(temp, temp, copyOfScale, bits);
        arb_floor(temp, temp, bits);
        arb_mul(temp, temp, copyOfScale, bits);
        arb_mul_si(temp,temp, 2, bits);
        arb_sub(answer, localCopyOfX, temp, bits);
    } else {
        //do something??
    }
    arb_clear(negScale);
    arb_clear(temp);
    arb_clear(localCopyOfX);
    arb_clear(copyOfScale);
}

*//*
 * Like interval reduce, but for the region defined by [-1,1]*v
 *
 * Project toReduce onto v to get the scalar projection arb_t s
 * final answer is toReduce - 2ceil((s-1)/2)*v
 *//*
void Auxiliary::vectorReduce(acb_struct *answer, const acb_struct *v, const acb_t toReduce, int bits) {
    acb_t copyOfToReduce;
    acb_init(copyOfToReduce);
    acb_set(copyOfToReduce, toReduce);

    acb_t copyOfV;
    acb_init(copyOfV);
    acb_set(copyOfV, v);

    arb_t vx, vy, rx, ry, s, denominator, temp;
    arb_init(vx);
    arb_init(vy);
    arb_init(rx);
    arb_init(ry);
    arb_init(s);
    arb_init(denominator);
    arb_init(temp);

    acb_get_real(vx,v);
    acb_get_imag(vy, v);
    acb_get_real(rx, copyOfToReduce);
    acb_get_imag(ry, copyOfToReduce);

    //Scalar projection of toReduce on to v
    arb_mul(s, vx, rx, bits);
    arb_addmul(s, vy, ry, bits);
    arb_sqr(denominator, vx, bits);
    arb_addmul(denominator, vy, vy, bits);
    arb_div(s, s, denominator, bits);

    //Do the rounding thing
    //temp = -2ceil((s-1)/2)
    arb_sub_si(temp, s, 1, bits);
    arb_div_si(temp, temp, 2, bits);
    arb_ceil(temp, temp, bits);
    arb_mul_si(temp, temp, 2, bits);
    arb_neg(temp, temp);

    //std::cout << "vector scaling is " << arb_get_str(temp, bits, 0) << std::endl;

    //Subtract the proper multiple of v from toReduce
    acb_mul_arb(answer, copyOfV, temp, bits);
    acb_add(answer, answer, copyOfToReduce, bits);

    acb_clear(copyOfToReduce);
    arb_clear(vx);
    arb_clear(vy);
    arb_clear(rx);
    arb_clear(ry);
    arb_clear(s);
    arb_clear(denominator);
    arb_clear(temp);
    acb_clear(copyOfV);
}

void Auxiliary::acbPlusArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_add(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxiliary::acbTimesArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_mul(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxiliary::acbDivArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_div(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxiliary::acbDivArb(acb_struct *answer, const acb_struct *z, int x, int bits) {
    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_set_si_si(complexEmbedding, x, 0);
    acb_div(answer, copyOfZ, complexEmbedding, bits);

    acb_clear(complexEmbedding);
    acb_clear(copyOfZ);
}

void Auxiliary::traceProduct(acb_struct *answer, const acb_struct *z, const acb_struct *w, int bits) {
    acb_t copyOfZ, copyOfW;
    acb_init(copyOfZ);
    acb_init(copyOfW);

    acb_set(copyOfZ, z);
    acb_set(copyOfW, w);

    acb_mul(answer, copyOfZ, copyOfW, bits);
    acb_conj(answer, answer);
    acb_addmul(answer, copyOfZ, copyOfW, bits);

    acb_clear(copyOfW);
    acb_clear(copyOfZ);
}*/

double Auxiliary::imagTheta(int d) {
    if (mod(-d, 4) == 1) {
        return sqrt(1.0*d)/2;
    } else {
        return sqrt(1.0*d);
    }
}

double Auxiliary::multiPrecisionSummation(const std::vector<double> &numbers) {
    int thread = omp_get_thread_num();

    mpfr_set_zero(&sums[thread], 1);
    for (const double number : numbers) {
        mpfr_add_d(&sums[thread], &sums[thread], number, MPFR_RNDN);
    }
    double answer = mpfr_get_d(&sums[thread], MPFR_RNDN);

    return answer;
}

/**
 * In scientific notation, if x = a * 10^b and a is a number in the interval [1,10),
 * returns next(a) * 10^b where next(a) is the closest integer strictly greater than a.
 * @param x = a * 10^b and a is a number in the interval [1,10)
 * @return next(a) * 10^b where next(a) is the closest integer strictly greater than a
 */
double Auxiliary::nextWithinOrderOfMag(double x) {
    double f = std::floor(std::log(x)/std::log(10.0));
    double g = floor(x/pow(10.0,f));
    double h = (g + 1)*pow(10.0,f);
    return h;
}

double Auxiliary::kahanSummation(std::vector<double> &numbers) {
    return kahanSummationUnvectorized(numbers);
}

double Auxiliary::kahanSummationUnvectorized(std::vector<double> &numbers) {
    std::sort(numbers.begin(), numbers.end(), [] (double left, double right) {
        return abs(left) < abs(right);
    });

    double sum = 0.0;
    double c = 0.0;
    for (const auto n : numbers) {
        double volatile y = n - c;
        double volatile t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

char Auxiliary::legendreSymbol(long n, long p) {
    n = n % p;
    if (n < 0) {
        n += p;
    }
    if (n == 0) {
        return 0;
    }
    if (p == 2) {
        return 1;
    }
    if (n == 2) {
        if (p % 8  == 1 || p % 8 == 7) {
            return 1;
        } else {
            return -1;
        }
    }
    if (mod(-1, p) == n) {
        if (p % 4 == 1) {
            return 1;
        } else {
            return -1;
        }
    }

    // 1 <= n < p
    for (int i = 1; i <= floor(p/2); i++) {
        if (mod(i*i, p) == n) {
            return 1;
        }
    }

    return -1;
}

double Auxiliary::rand(double a, double b) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dis(a, b);

    return dis(gen);
}




/*double Auxiliary::computeVolumeOfFD(int d) {
    pari_init(1000000,2);
    std::string arg = "L = lfuncreate(x^2+";
    arg.append(std::to_string(d));
    arg.append(");");
    gp_read_str(arg.c_str());
    arg = "lfun(L,2)";
    GEN res = gp_read_str(arg.c_str());
    char* out = GENtostr(res);
    pari_close();
    double zeta = std::stod(std::string(out));
    return pow(abs(d), 3.0/2)*zeta/(4*pow(Auxiliary::pi, 2));
}*/

//double Auxiliary::pi = std::numbers::pi_v<double>;


