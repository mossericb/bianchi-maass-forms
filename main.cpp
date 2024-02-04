/*TODO
     * Add Weyl law spacing to search method
     *
     * Continue making pointwise checking work!
     *
     * Make reals-only arb code mpfr
     *
     * Implement "search for well conditioned value of Y" in pointwise methods
     * Don't need it right now, would be nice to have eventually.
     *
     * Implement a way to not store the M0 and MY indices separately (lots of duplicates). Sometimes I need only
     * indices up to M0, sometimes up to MY. Carefully figure these out and address over-precomputing.
     *
     * Get rid of unordered_map's for orbit retrieval
     *
     * Make include guards have the right project name.
     *
     * Optimize reduction for d=3
     *
     * Make QuadraticIntegers class that contains all the data inherent to the ring
     * -theta, A, Y0, indices
     * -methods that require d?
     * */

#include <iostream>
#include "ArchtKBessel.h"
#include "BianchiMaassPointwise.h"
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
    K.setRAndClear(10);
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
    int D = 12;

    BianchiMaassPointwise bmp = BianchiMaassPointwise(d, D, symClass);

    bmp.checkSingleEigenvalue(6.0537);

    /*for (int D = 2; D <= 16; D++) {
        std::cout << "D = " << D << std::endl;
        BianchiMaassComputer* c1 = new BianchiMaassComputer(d, D, symClass, 6.05377197265625);
        c1->secantSearch(leftEndpoint, rightEndpoint);
        delete c1;
    }*/

    return 0;
}
