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
#include "BianchiMaassSearch.h"
#include <chrono>

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl << std::flush

int main(int argc, char *argv[]) {
    if (argc != 3) {
        throw(std::invalid_argument("command line argument should be two numbers"));
    }
    double leftEndpoint = std::stod(argv[1]);
    double rightEndpoint = std::stod(argv[2]);

    int d = 19;
    char symClass = 'C';
    int D = 4;

    BianchiMaassSearch bms = BianchiMaassSearch(d, D, symClass);
    bms.searchForEigenvalues(6, 7);

    //BianchiMaassPointwise bmp = BianchiMaassPointwise(d, D, symClass);
    //bmp.checkSingleEigenvalue(6.011020660400391 ,0.065);

    /*for (int D = 2; D <= 16; D++) {
        std::cout << "D = " << D << std::endl;
        BianchiMaassComputer* c1 = new BianchiMaassComputer(d, D, symClass, 6.05377197265625);
        c1->secantSearch(leftEndpoint, rightEndpoint);
        delete c1;
    }*/

    return 0;
}
