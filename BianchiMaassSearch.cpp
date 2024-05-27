//
// Created by Eric Moss on 2/3/24.
//

#include "BianchiMaassSearch.h"

#include <cmath>
#include <stack>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <map>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <filesystem>
#include <eigen3/Eigen/SVD>
#include "archtKBessel.h"
//#include "Plotter.h"
//#include "PlotWindow.h"
//#include "FunctionToEvaluate.h"


//Containers
using std::find, std::reverse, std::stack;

//Math
using std::cos, std::sin, std::get, std::max, std::pow, std::abs, std::ceil, std::floor, std::min;

//IO
using std::cout, std::string, std::endl, std::flush, std::setprecision, std::to_string;

#define watch(x) cout << (#x) << " is " << (x) << endl << flush

BianchiMaassSearch::BianchiMaassSearch(int d, int D, char symClass) {

    //Check that d is valid
    vector<int> classNumberOne = {1,2,3,7,11,19,43,67,163};
    bool dIsClassNumberOne = find(classNumberOne.begin(), classNumberOne.end(), d) != classNumberOne.end();
    if (!dIsClassNumberOne) {
        throw(std::invalid_argument("d should be one of {1,2,3,7,11,19,43,67,163}"));
    }

    //Check that D is valid
    if (D <= 0) {
        throw(std::invalid_argument("D should be positive. 10^-D is the truncation relativeError."));
    }

    //Check that symclass is valid
    vector<char> symClasses = {'D','G','C','H'};
    bool symClassIsCorrect = find(symClasses.begin(), symClasses.end(), symClass) != symClasses.end();
    if (!symClassIsCorrect) {
        throw(std::invalid_argument("symClass should be one of D,G,C,H"));
    }

    this->d = d;
    this->D = D;
    this->symClass = symClass;

    Od = ImaginaryQuadraticIntegers(d);

    A = Od.getA();
    theta = Od.getTheta();
    Y0 = Od.getY0();

    truncation = pow(10,-D);
    tolerance = pow(10,-(D+6));

    //open file
    const std::string directory = "Output/"; // Change this to the desired directory
    const std::string prefix = "output_"
                               + to_string(d) + "_"
                               + symClass + "_";
    createOutputDirectory(directory);
    int maxNumber = findMaxFileNumber(directory, prefix);
    std::string outputFilename = directory
            + prefix
            + std::to_string(maxNumber + 1)
            + ".txt";

    outputFile.open(outputFilename);

    if (outputFile.is_open()) {
        // Write content to the file
        outputFile << "d = " << d << '\n';
        outputFile << "symClass = " << symClass << '\n';
        outputFile << "D = " << D << std::endl;
    } else {
        std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
    }
}



BianchiMaassSearch::~BianchiMaassSearch() {
    outputFile.close();
}



void BianchiMaassSearch::searchForEigenvalues(const double leftR, const double rightR) {
    if (leftR >= rightR) {
        throw(std::invalid_argument("leftR should be less than rightR"));
    }
    outputFile << "leftEndpoint = " << leftR << '\n';
    outputFile << "rightEndpoint = " << rightR << std::endl;

    //Magic number
    double numEigenValues = 0.5;

    vector<double> endpoints;
    double endpoint = leftR;
    while (endpoint < rightR) {
        endpoints.push_back(endpoint);
        endpoint = Od.eigenvalueIntervalRightEndpoint(endpoint, numEigenValues);
    }
    endpoints.push_back(rightR);

    for (int i = 0; i <= endpoints.size() - 2; i++) {
        computeAllConditionNumbers = true;
        stack<pair<double, double>> intervalsToSearch;
        pair<double, double> interval = {endpoints[i], endpoints[i + 1]};
        intervalsToSearch.push(interval);

        while (!intervalsToSearch.empty()) {
            pair<double,double> nextInterval = intervalsToSearch.top();
            intervalsToSearch.pop();
            vector<pair<double,double>> intervals = conditionedSearchForEigenvalues(nextInterval.first, nextInterval.second);
            for (auto itr = intervals.crbegin(); itr != intervals.crend(); itr++) {
                intervalsToSearch.emplace(itr->first, itr->second);
            }
        }
        outputFile << "Complete up to " << endpoints[i + 1] << endl;
    }
}



/***************************************************
 * These functions NARROW an already small interval containing a likely eigenvalue.
 ***************************************************/


/**
 * @brief Assumes there is one eigenvalue between the given parameters. Narrows interval down
 * until it can't any more and saves the result.
 *
 * @param leftR Left endpoint of interval containing one eigenvalue.
 * @param rightR Right endpoint of interval containing one eigenvalue.
 */
void BianchiMaassSearch::narrowLikelyInterval(const double leftR, const double rightR) {
    //do something
}

double BianchiMaassSearch::computeM0General(const double r) {
    archtKBessel bess(r);

    //Method: use binary search to find M0 such that
    //K(2*pi/A * M0 * Y0) = 10^-D * K(max(r,1))
    //with M0 > max(r,1)/(2*pi/A*Y0)
    //Reference: JAY JORGENSON, LEJLA SMAJLOVI  ́ C, AND HOLGER THEN
    //ON THE DISTRIBUTION OF EIGENVALUES OF MAASS FORMS ON CERTAIN MOONSHINE GROUPS

    //Step 1: Make sure the equality point happens between my left and right endpoints
    double minM0 = max(r, 1.0)/(2*pi/A*Y0);
    double maxM0 = 100;

    double peak = bess.evaluate(max(r,1.0));
    double evalLeft = truncation * peak - bess.evaluate(2*pi/A*minM0*Y0);
    double evalRight = truncation * peak - bess.evaluate(2*pi/A*maxM0*Y0);
    while (!areDifferentSign(evalLeft, evalRight)) {
        maxM0 *= 2;
        evalRight = truncation * peak - bess.evaluate(2*pi/A*maxM0*Y0);
    }


    //Step 2: Do binary search
    double left = minM0;
    double right = maxM0;
    evalLeft = truncation * peak - bess.evaluate(2*pi/A*left*Y0);

    //Super accuracy here doesn't really matter, but it's fast, so we go this far because we can
    while (right - left > 0.000001) {
        double center = (left + right)/2.0;
        double evalCenter = truncation * peak - bess.evaluate(2*pi/A*center*Y0);

        if (areDifferentSign(evalLeft, evalCenter)) {
            right = center;
            evalRight = evalCenter;
        } else {
            left = center;
            evalLeft = evalCenter;
        }
    }

    double M0 = (left + right)/2.0;
    M0++;
    return M0;
}

/**
 * @brief Compute M(Y) based on M(Y_0) and Y. Only assumes Y0 is computed. Does not modify anything. If D changes,
 * then M0 must be computed first, then this.
 *
 * @param M0 The given M0 bound.
 * @param Y Given horocycle height.
 * @return The bound MY = Y0/Y*M0.
 */
double BianchiMaassSearch::computeMYGeneral(const double M0, const double Y) {
    return Y0/Y * M0;
}

vector<TestPointOrbitData>
BianchiMaassSearch::getPointPullbackOrbits(const Index &m, const double Y, const double MY) {
    vector<TestPointOrbitData> answer;

    //These formulas provide bounds to guarantee that the DFT result is valid
    double absRealM = abs(m.getComplex(d).real());
    double absImagM = abs(m.getComplex(d).imag());
    double Q0double = 0;
    double Q1double = (MY + absRealM)/2.0;
    if (Auxiliary::mod(-d, 4) == 1) {
        double Q0doubleOption1 = MY + absRealM;
        double Q0doubleOption2 = (MY + absImagM)/(2.0*A);
        Q0double = max(Q0doubleOption1, Q0doubleOption2);
    } else {
        Q0double = (MY + absImagM)/(2.0*A);
    }

    int Q0 = ceil(Q0double);
    int Q1 = ceil(Q1double);

    //Tweak Q0 and Q1 to be exactly what we need for exploiting symmetry
    if (Auxiliary::mod(-d, 4) == 1 && d != 3) {
        //If -d = 1 mod 4 then we need Q0/Q1 to be an even integer for the map x -> -bar(x) to be defined on test points
        if (Auxiliary::mod(Q0, Q1) == 0 && Auxiliary::mod(Q0 / Q1, 2) == 0) {
            // Q0/Q1 is an even integer
            // there is nothing to do!
        } else {
            // Q0/Q1 is not an even integer
            int k = 0;
            double quotient = ((double)Q0)/Q1;

            while (!(2 * k < quotient && quotient < 2 * (k+1))) {
                k++;
            }

            if (k == 0) {
                Q0 = 2*Q1;
            } else {
                //It is true that 2k < Q0/Q1 < 2(k+1) and k > 0
                int firstOptionQ0 = ceil(Q1 * 2 * (k+1));
                while (!Auxiliary::mod(firstOptionQ0, 2 * (k + 1) == 0)) {
                    firstOptionQ0++;
                }
                int firstOptionQ1 = firstOptionQ0/(2*(k+1));

                int secondOptionQ0 = ceil(Q1 * 2 * k);
                while (!Auxiliary::mod(secondOptionQ0, 2 * k == 0)) {
                    secondOptionQ0++;
                }
                int secondOptionQ1 = secondOptionQ0/(2 * k);

                //Now both give a valid number of points, choose the one that will have fewest points overall (4*Q0*Q1 total)
                if (firstOptionQ0 * firstOptionQ1 < secondOptionQ0 * secondOptionQ1) {
                    Q0 = firstOptionQ0;
                    Q1 = firstOptionQ1;
                } else {
                    Q0 = secondOptionQ0;
                    Q1 = secondOptionQ1;
                }
            }
        }
    } else if (d == 3 || d == 1) {
        //So that we can do rotations by something other than -1, we need Q0=Q1
        //Add a little numerical wiggle room
        Q0 = max(Q0, Q1) + 1;
        Q1 = Q0;
    } else {
        //We don't need to fiddle with the divisibility of Q0 and Q1
        //Add some numerical wiggle room
        Q0 += 1;
        Q1 += 1;
    }


    short nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    short reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

    //Now generate the representatives
    //then generate the orbits
    if (Auxiliary::mod(-d, 4) == 1 && d != 3) {
        //iterate to make orbit representatives
        for (int l1 = 1; l1 <= Q1; l1++) {
            int lowerBound = ceil(0.5 - Q0*(l1-0.5)/(2.0*Q1));
            int upperBound = floor(0.5 + Q0 - Q0*(l1 - 0.5)/(2.0*Q1));
            for (int l0 = lowerBound; l0 <= upperBound; l0++) {
                //build the orbit
                TestPointOrbitData orbit;
                complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                //This magic quantity tells you about the real part of x
                int discriminant = 4*Q1*l0 - 2*Q1 + 2*Q0*l1 - Q0;

                //First case detects when real(x) = 0
                //Second case detects when real(x) has absolute value = 1/2
                if (discriminant == 0 || abs(discriminant) - 4*Q0*Q1 == 0) {
                    //The orbit in this case is {x, -x}
                    //So orbit/{+/-1} = {[x]}
                    //Therefore there are no proper translates mod +/-1
                } else {
                    //The orbit in this case is {x, -x, -bar(x), bar(x)}
                    //So orbit/{+/-1} = {[x], [-bar(x)]}
                    //So the only proper translate is -bar(x)
                    pair<complex<double>, short> tup (-conj(x), reflectionSign);
                    orbit.properTranslatesModSign.push_back(tup);
                }
                answer.push_back(orbit);
            }
        }
    } else if (d != 1) {
        //So -d = 2,3 mod 4 and d != 1

        //iterate to make orbit representatives
        for (int l0 = 1; l0 <= Q0; l0++) {
            for (int l1 = 1; l1 <= Q1; l1++) {
                //The orbit is {x, -x, -bar(x), bar(x)}
                //So orbit/{+/-1} is {[x], [-bar(x)]}
                //So the only proper translate is -bar(x)
                TestPointOrbitData orbit;
                complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                pair<complex<double>, short> tup (-conj(x), reflectionSign);
                orbit.properTranslatesModSign.push_back(tup);

                answer.push_back(orbit);
            }
        }
    } else if (d == 1) {
        //Iterate through orbit representatives
        for (int l0 = 1; l0 <= Q0; l0++) {
            for (int l1 = 1; l1 <= Q1; l1++) {
                TestPointOrbitData orbit;

                complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                if (abs(l0) == abs(l1)) {
                    //The orbit is {x, ix, -x, -ix}
                    //So orbit/{+/-1} is {[x], [ix]}
                    //So the only proper translate is ix

                    pair<complex<double>, short> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);
                } else {
                    //The orbit is {x, ix, -x, -ix, bar(x), ibar(x), -bar(x), -ibar(x)}
                    //So orbit/{+/-1} is {[x], [ix], [-bar(x)], [-ibar(x)]}
                    //So the proper translates are ix, -bar(x), and -ibar(x)

                    pair<complex<double>, short> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    get<0>(tup) = -conj(x);
                    get<1>(tup) = reflectionSign;
                    orbit.properTranslatesModSign.push_back(tup);

                    get<0>(tup) = -I*conj(x);
                    get<1>(tup) = reflectionSign * nonsquareUnitSign;
                    orbit.properTranslatesModSign.push_back(tup);
                }

                answer.push_back(orbit);
            }
        }
    } else {
        //d == 3
        //
        //In order for rotation to be defined on the set of testpoints, we do not use the shifted version
        //WWWWWWWAAAAARRRNNNIIINNNGGG
        //THIS CODE IS A MESS, IT'S WRONG AND WRONG AGAIN, NEEDS SERIOUS MATHEMATICAL REWORK

        int l1LowerBound = 0;
        int l1UpperBound = floor(2.0*Q1/3.0);
        for (int l1 = l1LowerBound; l1 <= l1UpperBound; l1++) {
            int l0LowerBound = l1; //should really be ceil(Q0*l1/(double)Q1), but Q0=Q1
            int l0UpperBound = floor(Q0 - l1/2.0); //should really be floor(Q0 - l1*Q0/(2*Q1)), but Q0=Q1
            for (int l0 = l0LowerBound; l0 <= l0UpperBound; l0++) {
                TestPointOrbitData orbit;

                complex<double> x = l0 / (2.0 * Q0) + theta * (double) l1 / (2.0 * Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                if (l0 == 0 && l1 == 0) {
                    //orbit is degenerate, just {0}
                    //is this okay?

                } else {

                    //For this comment let zeta=e^(2pi i /6) be represented by w
                    //The orbit is {x, wx, w^2x, -x, -wx, -w^2x}
                    //So orbit/{+/-1} is {[x], [wx], [w^2x]}
                    //So the proper translates are wx, w^2x

                    pair<complex<double>, short> tup (theta*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    get<0>(tup) = theta*theta*x;
                    get<1>(tup) = 1; //theta*theta is a square unit
                    orbit.properTranslatesModSign.push_back(tup);

                    //magic quantity that tells me about the real part of x
                    int discriminant = 2*Q1*l0 + Q1*l1 - 2*Q0*Q1;
                    bool orbitIsHalfSize = (l0 == 0 || l1 == 0 || l0 + l1 == 0 || discriminant == 0);

                    if (!orbitIsHalfSize) {
                        //For this comment let zeta=e^(2pi i /6) be represented by w
                        //The orbit is {x, wx, w^2x, -x, -wx, -w^2x, -bar(x), -wbar(x), -w^2bar(x), bar(x), wbar(x), w^2bar(x)}
                        //So orbit/{+/-1} is {[x], [wx], [w^2x], [-bar(x)], [-wbar(x)], [-w^2bar(x)]}
                        //So the proper translates are wx, w^2x, -bar(x), -wbar(x), -w^2bar(x)
                        get<0>(tup) = -conj(x);
                        get<1>(tup) = reflectionSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        get<0>(tup) = -theta*conj(x);
                        get<1>(tup) = reflectionSign * nonsquareUnitSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        get<0>(tup) = -theta*theta*conj(x);
                        get<1>(tup) = reflectionSign; //theta*theta is a square unit
                        orbit.properTranslatesModSign.push_back(tup);
                    }
                }
            }
        }
    }

    return answer;
}


double BianchiMaassSearch::traceProduct(const complex<double> &z, const complex<double> &w) {
    auto answer = z*w + conj(z)*conj(w);
    return answer.real();
}



bool BianchiMaassSearch::areDifferentSign(const double &a, const double &b) {
    return (a > 0 && b < 0) || (a < 0 && b > 0);
}

/*
void BianchiMaassSearch::recursiveSearchForEigenvalues(const double leftR, const double rightR,
                                                         vector<double>& leftG,
                                                         vector<double>& rightG) {

    //Replace this with something that selects Y1MATRIX and Y2MATRIX more carefully, dependent on r and Y0, etc.

    if (leftG.empty() || rightG.empty()) {
        if (leftG.empty() && rightG.empty()) {
            K.setRAndClear(rightR);
            //compute all the indices and points etc
            M0 = computeM0General(rightR);

            Y2 = std::max(1.0, rightR)/(2*pi/A*M0);
            Y1 = 0.95 * Y2;

            MY1 = computeMYGeneral(M0, Y1); //Y1MATRIX is the smaller one, produces larger MY
            MY2 = computeMYGeneral(M0, Y2);
            computeIndexData();
            searchPossibleSignChanges = indexTransversal.size() * 2;
            computeTestPointData();
        }

        if (leftG.empty()) {
            K.setRAndClear(leftR);
            //K precomputation is called in the getGVector() function
            leftG = getGVector();
        }

        if (rightG.empty()) {
            K.setRAndClear(rightR);
            //K precomputation is called in the getGVector() function
            rightG = getGVector();
        }
    }

    int signChanges = countSignChanges(leftG, rightG);

    if (sufficientSignChanges(leftG, rightG)) {
        double stop = 0.0000001;
        if (rightR - leftR < stop) {
            std::cout << "[" << leftR << ", " << rightR << "], final precision reached, eigenvalue possible" << std::endl;
            return;
        }

        double heckeThreshold = 0.00001;
        if (rightR - leftR < heckeThreshold) {
            double hecke = heckeCheck();
            if (abs(hecke) > 1) {
                std::cout << "[" << leftR << ", " << rightR << "], stopping, failed Hecke check" << std::endl;
                return;
            }
        }

        //Call on the left and right intervals
        std::cout << "[" << leftR << ", " << rightR << "] " << signChanges << "/" << searchPossibleSignChanges << std::endl;
        double centerR = (leftR + rightR)/2.0;
        vector<double> centerG = vector<double>();
        //In this function call, centerG gets computed and set to this variable
        recursiveSearchForEigenvalues(leftR, centerR, leftG, centerG);

        //In this function call, centerG is computed so it won't recompute
        recursiveSearchForEigenvalues(centerR, rightR, centerG, rightG);
        return;
    } else {
        std::cout << "[" << leftR << ", " << rightR << "] " << signChanges << "/" << searchPossibleSignChanges;
        std::cout << ", stopping, failed sign change check" << std::endl;
        return;
    }
}*/


MatrixXd BianchiMaassSearch::produceMatrix(const vector<Index> &indexTransversal,
                                           map<Index, vector<TestPointOrbitData>> &mToTestPointData,
                                           map<Index, vector<pair<Index, int>>> &ntoIndexOrbitData,
                                           double Y,
                                           KBessel &K) {

    int size = indexTransversal.size();
    MatrixXd answer;
    answer.resize(size, size);

#pragma omp parallel for schedule(dynamic) collapse(2) default(none) shared(indexTransversal, size, answer, mToTestPointData, ntoIndexOrbitData,Y, K)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Index m = indexTransversal[i];
            Index n = indexTransversal[j];
            double entry = computeEntry(m, n, K, mToTestPointData[m], Y, ntoIndexOrbitData[n]);

            answer(i,j) = entry;
        }
    }

    return answer;
}

pair<MatrixXd, double> BianchiMaassSearch::solveMatrix(const MatrixXd &matrix, const vector<Index> &indexTransversal,
                                         const int indexOfNormalization) {
    int size = indexTransversal.size();

    MatrixXd x;
    x.resize(size - 1, 1);
    MatrixXd b;
    b.resize(size - 1, 1);

    MatrixXd B;
    B.resize(size - 1, size - 1);

    //Finding kernel of M is the original problem.
    //Solving Bx=b is the system after moving a column over and deleting a row.
    //set b to negative of the indexOfNormalization column
    //set B to be M deleting row and column indexOfNormalization
    for (int row = 0; row < size; row++) {
        for (int column = 0; column < size; column++) {
            if (column == indexOfNormalization && row != indexOfNormalization) {
                //add it to b
                if (row > indexOfNormalization) {
                    b(row - 1) = matrix(row, column);
                } else {
                    b(row) = matrix(row, column);
                }
            } else if (column != indexOfNormalization && row != indexOfNormalization) {
                //add it to tempA somehow
                int modifiedRow = (row > indexOfNormalization) ? row - 1 : row;
                int modifiedColumn = (column > indexOfNormalization) ? column - 1 : column;
                B(modifiedRow, modifiedColumn) = matrix(row, column);
            }
        }
    }


    //make b negative
    b = -b;

    //solve the equation Ax=b for x
    x = B.colPivHouseholderQr().solve(b);

    //Now we need to turn x into xPrime by adding the normalization back in.
    Eigen::Matrix<double, Eigen::Dynamic, 1> xPrime;
    xPrime.resize(size);
    for (int i = 0; i < size; i++) {
        if (i < indexOfNormalization) {
            xPrime(i) = x(i);
        } else if (i == indexOfNormalization) {
            xPrime(i) = 1.0;
        } else {
            xPrime(i) = x(i - 1);
        }
    }
    Eigen::JacobiSVD<MatrixXd> svd(B);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);

    return {xPrime, cond};
}

double
BianchiMaassSearch::computeEntry(const Index &m, const Index &n, KBessel &K,
                                 const vector<TestPointOrbitData> &mTestPointOrbits, const double Y,
                                 const vector<pair<Index, int>> &nIndexOrbitDataModSign) {
    double answer = 0;

    long pointCount = 0;
    for (const auto& orbit : mTestPointOrbits) {
        pointCount += (orbit.properTranslatesModSign.size() + 1) * 2;
    }

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        vector<double> terms;
        for (const auto& testPointOrbit : mTestPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;


            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * pi / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : nIndexOrbitDataModSign) {
                Index l = tup.first;
                complex<double> xStar = pullback.getComplex();
                double ellTerm = tup.second; //This is a_l/a_n
                ellTerm *= cos(pi/A * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = cos(pi/A * traceProduct(-I* m.getComplex(d), x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                short sign = tup.second;
                testPointTerm += sign * cos(pi/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Aux.multiPrecisionSummation(terms);
        answer *= -4.0 / pointCount;
    } else {
        vector<double> terms;
        for (const auto& testPointOrbit : mTestPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;

            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * pi / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : nIndexOrbitDataModSign) {
                Index l = get<0>(tup);
                complex<double> xStar = pullback.getComplex();
                double ellTerm = get<1>(tup); //This is a_l/a_n
                ellTerm *= sin(pi/A * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = sin(pi/A * traceProduct(-I* m.getComplex(d), x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = get<0>(tup);
                short sign = get<1>(tup);
                testPointTerm += sign * sin(pi/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Aux.multiPrecisionSummation(terms);
        answer *= +4.0 / pointCount;
    }

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm = Y * K.exactKBessel(2 * pi / A * m.getAbs(d) * Y);

        answer = Aux.multiPrecisionSummation({answer, deltaTerm});
    }
    return answer;
}


int BianchiMaassSearch::countSignChanges(const vector<double> &v1, const vector<double> &v2) {
    int answer = 0;
    assert(v1.size() == v2.size());
    for (int i = 0; i < v1.size(); i++) {
        if (areDifferentSign(v1[i], v2[i])) {
            answer++;
        }
    }
    return answer;
}


bool BianchiMaassSearch::signChangeVectorIsIncreasing(vector<int> &v) {
    if (v.size() < 3) {
        return true;
    }

    int x3 = v.size() - 1;
    int x2 = x3 - 1;
    int x1 = x2 - 1;
    double y1 = v[x1];
    double y2 = v[x2];
    double y3 = v[x3];

    double xbar = x2;
    double ybar = (y1+y2+y3)/3.0;
    double beta = ((double)x1 - xbar)*(y1 - ybar) + ((double)x2 - xbar)*(y2 - ybar) + ((double)x3 - xbar)*(y3 - ybar);
    if (beta > 0) {
        return true;
    } else {
        return false;
    }
}
/**
 * If (x0, y0) and (x1, y1) are two distinct points then there is a unique line passing through them.
 * This function returns the unique zero of this line. Of course, if y0=y1 then there is no zero, and
 * the formula will return infinity if you do that.
 *
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 * @return The (x-coordinate of the) zero of the line passing through (x0, y0) and (x1, y1).
 */
double BianchiMaassSearch::findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1) {
    return -y0 * (x1 - x0)/(y1 - y0) + x0;
}

/*
 * TODO: make hecke check use two hecke relations
 * implement for all fields
 */
double BianchiMaassSearch::heckeCheck(map<Index, double>& coeffMap) {
    if (d == 1) {
        double avgSize = abs(coeffMap[Index(1,1)]);
        avgSize += abs(coeffMap[Index(3,0)]);
        avgSize += abs(coeffMap[Index(3,3)]);
        avgSize += abs(coeffMap[Index(1,2)]);
        avgSize += abs(coeffMap[Index(-1,3)]);
        avgSize /= 5;
        double hecke1 = coeffMap[Index(1,1)] * coeffMap[Index(3,0)] - coeffMap[Index(3,3)];
        double hecke2 = coeffMap[Index(1,1)] * coeffMap[Index(1,2)] - coeffMap[Index(-1,3)];
        hecke1 = abs(hecke1);
        hecke2 = abs(hecke2);
        double hecke = (hecke1 + hecke2)/avgSize;
        return hecke;
    } else if (d == 2) {
        double avgSize = abs(coeffMap[Index(0, 1)]);
        avgSize += abs(coeffMap[Index(-1, 1)]);
        avgSize += abs(coeffMap[Index(-2, -1)]);
        avgSize += abs(coeffMap[Index(5, 0)]);
        avgSize += abs(coeffMap[Index(0, 5)]);
        avgSize /= 5;
        double hecke1 = coeffMap[Index(0, 1)] * coeffMap[Index(-1, 1)] - coeffMap[Index(-2, -1)];
        double hecke2 = coeffMap[Index(0, 1)] * coeffMap[Index(5, 0)] - coeffMap[Index(0, 5)];
        hecke1 = abs(hecke1);
        hecke2 = abs(hecke2);
        double hecke = (hecke1 + hecke2)/avgSize;
        return hecke;
    } else if (d == 19) {
        double hecke1 = coeffMap[Index(2,0)] * coeffMap[Index(3,0)] - coeffMap[Index(6,0)];
        double hecke2 = coeffMap[Index(2, 0)] * coeffMap[Index(0, 1)] - coeffMap[Index(0, 2)];
        hecke1 = abs(hecke1);
        hecke2 = abs(hecke2);
        double hecke = hecke1 + hecke2;
        return hecke;
    } else {
        throw(std::invalid_argument("Hecke check not implemented."));
    }
}

bool BianchiMaassSearch::sufficientSignChanges(const vector<double> &v1, const vector<double> &v2) {
    int answer1 = 0;
    int answer2 = 0;

    for (int i = 0; i < v1.size()/2; i++) {
        if (areDifferentSign(v1[i], v2[i])) {
            answer1++;
        }
    }
    for (int i = v1.size()/2; i < v1.size(); i++) {
        if (areDifferentSign(v1[i], v2[i])) {
            answer2++;
        }
    }
    return (answer1 > v1.size()/6.0) && (answer2 > v1.size()/6.0);
}

double BianchiMaassSearch::minBess(KBessel &K, const vector<Index> &indexTransversal, double Y) {
    double min = +INFINITY;

    for(const auto index : indexTransversal) {
        double arg = 2 * pi / A * index.getAbs(d) * Y;
        if (arg > K.getR()) {
            break;
        }
        double bess = Y * K.exactKBessel(arg);
        if (abs(bess) < min) {
            min = abs(bess);
        }
    }
    double arg = 2 * pi / A * indexTransversal[indexTransversal.size() - 1].getAbs(d) * Y;
    double bess = Y * K.exactKBessel(arg);
    if (abs(bess) < min) {
        min = abs(bess);
    }

    return min;
}

vector<std::pair<double, double>> BianchiMaassSearch::conditionedSearchForEigenvalues(const double leftR, const double rightR) {
    //Assume that I will have to recompute everything multiple times here, so don't bother checking otherwise

    //First time through, compute the condition number everywhere to get an idea of how low you can go
    //In subsequent calls, start with Y1 and Y2 at the top end of the range and adjust a few times until
    //the condition number is close to the computed ballpark

    std::cout << std::setprecision(16)
              << "[" << leftR << ", " << rightR << "]" << std::endl;

    double M0 = computeM0General(rightR);

    KBessel K = KBessel(2 * pi / A /*usual factor*/
                        * 1 /*smallest magnitude of an index*/
                        * Y0 /*smallest height of a pullback*/
            , rightR);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

    int indexOfNormalization = 0;
    for (int i = 0; i < indexTransversal.size(); i++) {
        if (indexTransversal[i] == Index(1,0)) {
            indexOfNormalization = i;
            break;
        }
    }

    double Y1;
    double Y2;

    double MY1;
    double MY2;

    double maxYStar = 0;

    map<Index, vector<TestPointOrbitData>> mToY1TestPointOrbits;
    map<Index, vector<TestPointOrbitData>> mToY2TestPointOrbits;

    MatrixXd matrixY1;
    MatrixXd matrixY2;

    MatrixXd bY1;
    MatrixXd bY2;

    map<Index, double> coeffMap;

    //double upperY = 2.0 * max(1.0, rightR)/(2*pi/A*M0);
    //double lowerY = 1.7 * max(1.0, rightR)/(2*pi/A*M0);
    double c = max(1.0, rightR)/(2*pi/A*M0);
    double upperY = Y0;
    double lowerY = 0.5 * max(1.0, rightR)/(2*pi/A*M0);

    vector<double> heightsToTry;

    double tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.99;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    vector<double> r1ConditionNumbers;
    vector<double> r2ConditionNumbers;

    vector<double> r1MinDiag(heightsToTry.size(), 0);
    vector<double> r2MinDiag(heightsToTry.size(), 0);

    //compute Y1R2 matrix for all Y
    K.setRAndClear(rightR);
#pragma omp parallel for default(none) shared(heightsToTry, indexTransversal, r2MinDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r2MinDiag[i] = minBess(K, indexTransversal, Y);
    }

    K.setRAndClear(leftR);
#pragma omp parallel for default(none) shared(heightsToTry, indexTransversal, r1MinDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r1MinDiag[i] = minBess(K, indexTransversal, Y);
    }


    double ansY1 = 0;
    double ansY2 = 0;
    double max = 0;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i <= heightsToTry.size() - 2; i++) {//end i at second to last because Y2 has to be greater
        for (int j = i + 1; j < heightsToTry.size(); j++) {//j is always less than i, so j indexes a larger Y
            double maxMinBess = std::min({r1MinDiag[i],
                                  r1MinDiag[j],
                                  r2MinDiag[i],
                                  r2MinDiag[j]});
            if (max < maxMinBess) {
                max = maxMinBess;
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
            }
        }
    }

    argminY1 = ansY1;
    argminY2 = ansY2;

    /*Do it all over again but finer
     *and in a narrower range
     */

    upperY = ansY2 + 0.5 * c;
    lowerY = ansY1 - 0.5 * c;

    heightsToTry.clear();

    tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.9995;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    r1MinDiag.clear();
    r2MinDiag.clear();

    r1MinDiag.resize(heightsToTry.size(), 0);
    r2MinDiag.resize(heightsToTry.size(), 0);

    //compute Y1R2 matrix for all Y
    K.setRAndClear(rightR);
#pragma omp parallel for default(none) shared(heightsToTry, indexTransversal, r2MinDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r2MinDiag[i] = minBess(K, indexTransversal, Y);
    }

    K.setRAndClear(leftR);
#pragma omp parallel for default(none) shared(heightsToTry, indexTransversal, r1MinDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r1MinDiag[i] = minBess(K, indexTransversal, Y);
    }


    ansY1 = 0;
    ansY2 = 0;
    max = 0;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i <= heightsToTry.size() - 2; i++) {//end i at second to last because Y2 has to be greater
        for (int j = i + 1; j < heightsToTry.size(); j++) {//j is always less than i, so j indexes a larger Y
            double maxMinBess = std::min({r1MinDiag[i],
                                          r1MinDiag[j],
                                          r2MinDiag[i],
                                          r2MinDiag[j]});
            if (max < maxMinBess) {
                max = maxMinBess;
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
            }
        }
    }

    Y1 = ansY1;
    Y2 = ansY2;

    computeAllConditionNumbers = false;
    conditionGoal = 0;

    //compute the matrices again :)

    std::cout << Y1/c << " " << Y2/c << std::endl;

    MY1 = computeMYGeneral(M0, Y1);
    MY2 = computeMYGeneral(M0, Y2);

    mToY1TestPointOrbits.clear();
    mToY2TestPointOrbits.clear();
    maxYStar = 0;


#pragma omp parallel for default(none) shared(indexTransversal, mToY1TestPointOrbits, mToY2TestPointOrbits, Y1, MY1, Y2, MY2)
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index m = indexTransversal[i];
        auto orbitsY1 = getPointPullbackOrbits(m, Y1, MY1);
        auto orbitsY2 = getPointPullbackOrbits(m, Y2, MY2);
#pragma omp critical
        {
            mToY1TestPointOrbits[m] = orbitsY1;
            mToY2TestPointOrbits[m] = orbitsY2;
        }
    }

    for (const auto& m : indexTransversal) {
        for (const auto& itr : mToY1TestPointOrbits[m]) {
            if (maxYStar < itr.representativePullback.getJ()) {
                maxYStar = itr.representativePullback.getJ();
            }
        }
    }

    K.setRAndClear(leftR);
    K.extendPrecomputedRange(2 * pi / A * M0 * maxYStar);

    matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, K);
    matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, K);

    auto data1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto data2 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    std::cout << data1.second << " " << data2.second << " ";

    K.setRAndPrecompute(rightR, 2 * pi / A * M0 * maxYStar);

    matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, K);
    matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, K);

    data1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    data2 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    std::cout << data1.second << " " << data2.second << std::endl;

    /*
     * When you have to start over, for whatever reason, this means that condition numbers
     * blow up in that interval. The method resets the condition number search goal, but I
     * don't want this goal to stay artificially high. So whenever the algorithm has to start over,
     * it sets this flag to signal that next time it should re-restart over to make the condition
     * search goal lower again.
     */

    int searchPossibleSignChanges = indexTransversal.size() * 2;

    int signChanges = 0;

    for (int i = 0; i < indexTransversal.size(); i++) {
        Index n = indexTransversal[i];
        double coeff = data1.first(i);
        coeffMap[n] = coeff;
    }

    double hecke = heckeCheck(coeffMap);
    std::cout << "hecke " << hecke << std::endl;

    /* We expect most eigenvalues to appear like this
     */
    if (signChanges >= searchPossibleSignChanges * 0.975 && rightR - leftR < 0.001) {
        for (int i = 0; i < indexTransversal.size(); i++) {
            Index n = indexTransversal[i];
            double coeff = data1.first(i);
            coeffMap[n] = coeff;
        }

        double hecke = heckeCheck(coeffMap);

        std::cout << "--very high sign change count, eigenvalue probable\n";
        std::cout << "--sign changes " << signChanges << "/" << searchPossibleSignChanges << '\n';
        outputFile << setprecision(16) << "[" << leftR << ", " << rightR << "]" << std::endl;
        return {};
    }


    /* This is to ensure that we can zoom in far enough to find the good intervals.
     * Also provides a way to catch some eigenvalues that might not appear as clearly.
     *
     */
    if (signChanges > searchPossibleSignChanges/4.0) {
        double heckeThreshold = 0.0001;
        if (rightR - leftR < heckeThreshold) {
            //Compute the coeffMap
            for (int i = 0; i < indexTransversal.size(); i++) {
                Index n = indexTransversal[i];
                double coeff = data1.first(i);
                coeffMap[n] = coeff;
            }

            double hecke = heckeCheck(coeffMap);
            std::cout << "--hecke " << hecke << '\n';

            if (abs(hecke) > 1) {
                std::cout << "--sign changes " << signChanges << "/" << searchPossibleSignChanges << '\n';
                std::cout << "--stopping, failed Hecke check" << endl;
                return {};
            }
        }

        /*
         * Stopping condition determined by empirical observation, but is also
         * sensible! If you forge on past this point you start detecting fewer sign changes
         * than are really there, and the distinction between an interval having the eigenvalue
         * or not becomes blurry.
         */
        double stop = pow(10, -D + 2);
        bool finalPrecisionReached = (rightR - leftR) < stop;
        if (finalPrecisionReached && signChanges > searchPossibleSignChanges/2.0) {
            std::cout << "--final precision reached, eigenvalue possible\n";
            std::cout << "--sign changes " << signChanges << "/" << searchPossibleSignChanges << std::endl;
            outputFile << setprecision(16) << "[" << leftR << ", " << rightR << "]" << std::endl;
            return {};
        } else if (finalPrecisionReached) {
            std::cout << "--final precision reached, eigenvalue unlikely\n";
            std::cout << "--sign changes " << signChanges << "/" << searchPossibleSignChanges << std::endl;
            return {};
        }



        //Call on the left and right intervals
        std::cout << "--passed sign check\n";
        std::cout << "--sign changes " << signChanges << "/" << searchPossibleSignChanges << endl;
        double centerR = (leftR + rightR)/2.0;

        vector<pair<double,double>> intervals = {pair<double,double>{leftR, centerR},
                                                 pair<double,double>{centerR, rightR}};
        return intervals;
    } else {
        std::cout << "--stopping, failed sign change check, " << signChanges
        << "/" << searchPossibleSignChanges << std::endl;
        return {};
    }

}


int BianchiMaassSearch::findMaxFileNumber(const string &directory, const string &prefix) {
    int maxNumber = 0;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().stem().string();

            if (filename.find(prefix) == 0) {
                // Extract the number after prefix
                std::string numberStr = filename.substr(prefix.length());
                try {
                    int number = std::stoi(numberStr);
                    maxNumber = std::max(maxNumber, number);
                } catch (const std::invalid_argument& e) {
                    // Ignore files with invalid numbers in the filename
                }
            }
        }
    }

    return maxNumber;
}

void BianchiMaassSearch::createOutputDirectory(const std::string& directory) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directory(directory);
    }
}

vector<double> BianchiMaassSearch::mergeToVector(const MatrixXd& v1, const MatrixXd& v2) {
    if (v1.cols() != 1
        || v2.cols() != 1
        || v1.rows() != v2.rows()) {
        throw std::runtime_error("Solutions vectors do not have the right dimensions.");
    }

    vector<double> answer;
    answer.reserve(v1.rows() * 2);

    for (int i = 0; i < v1.rows(); i++) {
        answer.push_back(v1(i));
    }
    for (int i = 0; i < v2.rows(); i++) {
        answer.push_back(v1(i));
    }
    return answer;
}

/**
 * Given a relatively narrow interval [r1,r2] containing an eigenvalue,
 * this zooms in on the eigenvalue using the secant method (uses the first
 * Hecke relation as the cost function). r2-r1 should be at most 0.001.
 *
 * @param r1 Left endpoint of interval containing eigenvalue.
 * @param r2 Right endpoint of interval containing eigenvalue.
 */
void BianchiMaassSearch::secantMethod(double r1, double r2) {
    //find the conditioning goal the first time this is called
    //compute r1 matrix
    //compute hecke check for r1
    //compute hecke check for r2
    //secant
    //loop until
}

double BianchiMaassSearch::phi(map<Index, double> &coeffMap1, map<Index, double> &coeffMap2, vector<int> &signs) {
    if (d == 19) {
        double check1 = coeffMap1[Index(2,0)] - coeffMap2[Index(2,0)];
        double check2 = coeffMap1[Index(3,0)] - coeffMap2[Index(3,0)];
        double check3 = coeffMap1[Index(0,1)] - coeffMap2[Index(0,1)];
        check1 = abs(check1);
        check2 = abs(check2);
        check3 = abs(check3);
        return signs[0]*check1 + signs[1]*check2 + signs[2]*check3;
    } else {
        throw std::invalid_argument("not implemented for this d");
    }
}

