/*TODO
     * For each m in indexTransversal, compute Q, then save a map that maps m to Q, and then have another map
     * that stores the test points for a particular Q.
     * There are lots of repeat values of Q, no sense in double-computing them.
     *
     * Compute the matrix one column at a time so that I can compute the test points once and use them all, then
     * throw them away and I haven't wasted anything.
     *
     * Optimize reduction for d=3
     *
     * Goal for this version
     *
     * Implement double-precision K-bessel evaluator with interpolation
     * Instead of secant method search, do a binary search. This should rely less on high precision.
     *
     * Stop double computing at beginning of a new interval.
     * Are my bounds of MY and M0 correct? Too large? I don't know.
     * Lagrange interpolation of the trisection step! To speed up interval bisection.
     * */

#include <iostream>
#include "CoefficientComputer.h"
#include "ArchtKBessel.h"
#include <chrono>

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl << std::flush

int main(int argc, char *argv[]) {
    if (argc != 3) {
        throw(std::invalid_argument("command line argument should be two numbers"));
    }
    double leftEndpoint = std::stod(argv[1]);
    double rightEndpoint = std::stod(argv[2]);

    /*auto start = std::chrono::high_resolution_clock::now();

    ArchtKBessel A = ArchtKBessel(53);
    A.setR(9.554);
    double y = A.evaluate(5.6);
    std::cout << std::setprecision(16);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = end - start;
    watch(y);
    watch(duration);

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < pow(10,4); i++) {
        double x = 0.01 + i * 50.0/pow(10,4);
        double y = A.evaluate(x);
    }
    end = std::chrono::high_resolution_clock::now();

    duration = end - start;
    watch(duration);

    A.setR(10);
    start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < pow(10,4); i++) {
        double x = 0.01 + i * 50.0/pow(10,4);
        double y = A.evaluate(x);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    watch(duration);

    KBesselApproximator K = KBesselApproximator(53);
    K.setRAndClearPrecompute(10);
    start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < pow(10,4); i++) {
        double x = 0.01 + i * 50.0/pow(10,4);
        double y = K.computeKBessel(x);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    watch(duration);

    double x = 7.51;
    watch(A.evaluate(x));
    watch(K.exactKBessel(x));*/

    int d = 19;
    char symClass = 'C';
    int D = 8;

    CoefficientComputer* c1 = new CoefficientComputer(d, D, symClass);
    c1->checkSingleEigenvalue2(6.0537);

    /*for (int D = 2; D <= 16; D++) {
        std::cout << "D = " << D << std::endl;
        CoefficientComputer* c1 = new CoefficientComputer(d, D, symClass, 6.05377197265625);
        c1->secantSearch(leftEndpoint, rightEndpoint);
        delete c1;
    }*/

    return 0;
}
