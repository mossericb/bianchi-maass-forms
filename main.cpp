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
        if (mode == "coarse") {
            if (argc != 7) {
                throw std::invalid_argument("Coarse search command line arguments should be: coarse d symClass D leftEndpoint rightEndpoint");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);
            double leftEndpoint = std::stod(argv[5]);
            double rightEndpoint = std::stod(argv[6]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';
            std::cout << "leftEndpoint = " << leftEndpoint << '\n';
            std::cout << "rightEndpoint = " << rightEndpoint << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.coarseSearchForEigenvalues(leftEndpoint, rightEndpoint);
        } else if (mode == "medium") { //this indicates medium search
            if (argc != 5) {
                throw std::invalid_argument("Medium search command line arguments should be: medium d symClass D");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.mediumSearchForEigenvalues();
        } else if (mode == "fine") { //this indicates fine search
            if (argc != 5) {
                throw std::invalid_argument("Fine search command line arguments should be: fine d symClass D where 10^-D is the final Hejhal truncation error for the secant method.");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.fineSearchForEigenvalues();
        } else if (mode == "extend") {
            if (argc != 4) {
                throw std::invalid_argument("Extend mode command line arguments should be: extend d symClass");
            }
            //Compute more coefficients. Compute a bunch. Up to M(Y) by default and can ask for more perhaps?
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, 0, symClass);
            bms.extendCoefficientComputation();
        } else if (mode == "test-modularity") {
            if (argc != 4) {
                throw std::invalid_argument("test-modularity mode command line arguments should be: test-modularity d symClass");
            }
            //Check modularity
            //Check Ramanujan
            //Produce Sato-Tate distribution
            //Check if the reflection-odd symmetry class is also modular
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, 0, symClass);
            bms.testForModularity();
        } else if (mode == "test-conjectures") {
            if (argc != 4) {
                throw std::invalid_argument("test-conjectures mode command line arguments should be: test-conjectures d symClass");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, 0, symClass);
            bms.testConjectures();
        } else if (mode == "Lfunction") {
            if (argc != 4) {
                throw std::invalid_argument("Lfunction mode command line arguments should be: Lfunction d");
            }
            //Produce the Dirichlet coefficients for all tested even forms
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, 0, symClass);
            bms.makeLFunctions();
        } else {
            throw std::invalid_argument(R"(First command line argument should be "coarse" "medium" "fine" "extend" "test-modularity" "test-conjectures" or "Lfunction".)");
        }
    }

    return 0;
}
