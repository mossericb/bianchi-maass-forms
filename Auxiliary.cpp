#include "Auxiliary.h"
#include "omp.h"

#include <algorithm>
#include <random>
#include <sstream>
#include <iomanip>

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

string Auxiliary::commonRound(string & s1, string & s2) {

    if (s1 == s2) {
        return s1;
    }

    double r1 = std::stod(s1);
    double r2 = std::stod(s2);
    if (r1 == r2) {
        std::stringstream ss;
        ss << std::setprecision(16) << r1;
        return ss.str();
    }


    double last = roundN(r1,0);
    int n = 1;
    while (true) {
        double round1 = roundN(r1, n);
        double round2 = roundN(r2, n);
        if (round1 == round2) {
            last = round1;
        } else {
            std::stringstream ss;
            ss << std::setprecision(16) << last;
            return ss.str();
        }
        n++;
    }
}

double Auxiliary::roundN(const double r, int n) {
    double multiplier = std::pow(10, n);
    return std::round(r * multiplier) / multiplier;
}