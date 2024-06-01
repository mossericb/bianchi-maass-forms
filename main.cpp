/*TODO
     * Make include guards have the right project name.
     *
     * Optimize reduction for d=3
     *
     * Fix secant method where r1 = r2 leading to nan error
     * */

#include <iostream>
#include "BianchiMaassPointwise.h"
#include "BianchiMaassSearch.h"
#include <chrono>
#include <acb_hypgeom.h>

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl << std::flush
using namespace std::chrono;

int main(int argc, char *argv[]) {
    if (argc == 1) {
        throw std::invalid_argument("Please provide command line arguments.");
    }
    if (argc >= 2) {
        string mode = argv[1];
        if (mode == "0") { //this indicates coarse search
            if (argc != 7) {
                throw std::invalid_argument("Coarse search command line arguments should be: 0 d symClass D leftEndpoint rightEndpoint");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);
            double leftEndpoint = std::stod(argv[5]);
            double rightEndpoint = std::stod(argv[6]);

            std::cout << "mode = " << mode << ", coarse\n";
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';
            std::cout << "leftEndpoint = " << leftEndpoint << '\n';
            std::cout << "rightEndpoint = " << rightEndpoint << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode[0], d, D, symClass);
            bms.coarseSearchForEigenvalues(leftEndpoint, rightEndpoint);
        } else if (mode == "1") { //this indicates medium search
            if (argc != 5) {
                throw std::invalid_argument("Medium search command line arguments should be: 1 d symClass D");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);

            std::cout << "mode = " << mode << ", medium\n";
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode[0], d, D, symClass);
            bms.mediumSearchForEigenvalues();
        } else if (mode == "2") { //this indicates fine search
            if (argc != 5) {
                throw std::invalid_argument("Fine search command line arguments should be: 2 d symClass D where 10^-D is the final Hejhal truncation error for the secant method.");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);

            std::cout << "mode = " << mode << ", fine\n";
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode[0], d, D, symClass);
            bms.fineSearchForEigenvalues();
        } else {
            throw std::invalid_argument("First command line argument should be 0, 1, or 2 meaning coarse, medium, or fine search.");
        }
    }

    return 0;
}
