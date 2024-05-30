/*TODO
     * Use the KBessel precomputation for the 4 symmetry classes, might as well! Define a special character for symClass
     * that triggers searching for all the classes at once.
     *
     * Add Weyl law spacing to search method
     *
     * Continue making pointwise checking work!
     *
     * Fully implement the rest of the fields
     *
     * Implement a way to not store the M0 and MY indices separately (lots of duplicates). Sometimes I need only
     * indices up to M0, sometimes up to MY. Carefully figure these out and address over-precomputing.
     *
     * Make include guards have the right project name.
     *
     * Optimize reduction for d=3
     * */

#include <iostream>
#include "BianchiMaassPointwise.h"
#include "BianchiMaassSearch.h"
#include <chrono>
#include <acb_hypgeom.h>

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl << std::flush
using namespace std::chrono;

int main(int argc, char *argv[]) {
    if (argc != 6) {
        throw (std::invalid_argument("command line arguments should be: d symClass D leftEndpoint rightEndpoint"));
    }
    int d = std::stoi(argv[1]);
    char symClass = argv[2][0];
    int D = std::stoi(argv[3]);
    double leftEndpoint = std::stod(argv[4]);
    double rightEndpoint = std::stod(argv[5]);

    std::cout << "d = " << d << '\n';
    std::cout << "symClass = " << symClass << '\n';
    std::cout << "D = " << D << '\n';
    std::cout << "leftEndpoint = " << leftEndpoint << '\n';
    std::cout << "rightEndpoint = " << rightEndpoint << '\n';

    BianchiMaassSearch bms = BianchiMaassSearch(d, D, symClass);
    bms.coarseSearchForEigenvalues(leftEndpoint, rightEndpoint);
    bms.refineEigenvalueIntervals();


    //BianchiMaassPointwise bmp = BianchiMaassPointwise(d, D, symClass);
    //bmp.checkSingleEigenvalue(6.011020660400391 ,0.065);

    /*for (int D = 2; D <= 16; D++) {
        std::cout << "D = " << D << std::endl;
        BianchiMaassComputer* c1 = new BianchiMaassComputer(d, D, symClass, 6.05377197265625);
        c1->secantMethod(leftEndpoint, rightEndpoint);
        delete c1;
    }*/

    return 0;
}
