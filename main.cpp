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
        } else if (mode == "test-modularity-compute-reflections") {
            if (argc != 4) {
                throw std::invalid_argument("test-modularity mode command line arguments should be: test-modularity d symClass");
            }
            //Check modularity
            //Check Ramanujan
            //Produce Sato-Tate distribution
            //Check if the reflection-odd symmetry class is also modular
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];

            if (symClass != 'D' && symClass != 'C') {
                throw std::invalid_argument("symClass for this mode should be only D or C");
            }

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
        } else if (mode == "sandbox") {
            if (argc != 7) {
                throw std::invalid_argument("sandbox search command line arguments should be: sandbox d symClass D leftEndpoint rightEndpoint");
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
            bms.sandbox(leftEndpoint, rightEndpoint);
        } else if (mode == "sandbox2") {
            if (argc != 6) {
                throw std::invalid_argument("");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);
            double r = std::stod(argv[5]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';
            std::cout << std::setprecision(16) << "r = " << r << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.sandbox2(r);
        } else if (mode == "sandbox3") {
            if (argc != 6) {
                throw std::invalid_argument("");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);
            double r = std::stod(argv[5]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';
            std::cout << std::setprecision(16) << "r = " << r << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.sandbox3(r);
        } else if (mode == "sandbox4") {
            if (argc != 6) {
                throw std::invalid_argument("");
            }
            int d = std::stoi(argv[2]);
            char symClass = argv[3][0];
            double D = std::stod(argv[4]);
            double r = std::stod(argv[5]);

            std::cout << "mode = " << mode << '\n';
            std::cout << "d = " << d << '\n';
            std::cout << "symClass = " << symClass << '\n';
            std::cout << "D = " << D << '\n';
            std::cout << std::setprecision(16) << "r = " << r << '\n';

            BianchiMaassSearch bms = BianchiMaassSearch(mode, d, D, symClass);
            bms.sandbox4(r);
        } else {
            throw std::invalid_argument(R"(First command line argument should be "coarse" "medium" "fine" "test-modularity" "test-conjectures".)");
        }
    }

    return 0;
}
