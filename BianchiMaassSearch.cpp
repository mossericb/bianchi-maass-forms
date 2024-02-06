//
// Created by Eric Moss on 2/3/24.
//

#include "BianchiMaassSearch.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <map>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <eigen3/Eigen/SVD>
//#include "Plotter.h"
//#include "PlotWindow.h"
//#include "FunctionToEvaluate.h"


//Containers
using std::find, std::reverse;

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
        throw(std::invalid_argument("D should be positive. 10^-D is the truncation error."));
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

    //outputFile.open("d" + to_string(d) + "_D" + to_string(D) + "_" + symClass + "_" + to_string(Y1Parameter) + ".txt");

    truncation = pow(10,-D);
    tolerance = pow(10,-(D+6));

    cout << setprecision(16);
    cout << "A: " << A << endl;
    cout << "theta: " << theta << endl;
    cout << "Y0: " << Y0 << endl;
    cout << "truncation: " << truncation << endl;
}



BianchiMaassSearch::~BianchiMaassSearch() {
    outputFile.close();
}



void BianchiMaassSearch::searchForEigenvalues(const double leftR, const double rightR) {
    //TODO: Initialize search private members.
    if (leftR >= rightR) {
        throw(std::invalid_argument("leftR should be less than rightR"));
    }
    vector<double> leftG = {};
    vector<double> rightG = {};
    recursiveSearchForEigenvalues(leftR, rightR, leftG, rightG);
}



void BianchiMaassSearch::clearSearchData() {
    //delete vectors and maps and stuff
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
    KBesselApproximator tempK = KBesselApproximator(53);
    tempK.setRAndClear(r);

    //Method: use binary search to find M0 such that
    //K(2*pi/A * M0 * Y0) = 10^-D * K(max(r,1))
    //with M0 > max(r,1)/(2*pi/A*Y0)
    //Reference: JAY JORGENSON, LEJLA SMAJLOVI  ́ C, AND HOLGER THEN
    //ON THE DISTRIBUTION OF EIGENVALUES OF MAASS FORMS ON CERTAIN MOONSHINE GROUPS

    //Step 1: Make sure the equality point happens between my left and right endpoints
    double minM0 = max(r, 1.0)/(2*pi/A*Y0);
    double maxM0 = 100;

    double peak = tempK.exactKBessel(max(r,1.0));
    double evalLeft = truncation * peak - tempK.exactKBessel(2*pi/A*minM0*Y0);
    double evalRight = truncation * peak - tempK.exactKBessel(2*pi/A*maxM0*Y0);
    while (!areDifferentSign(evalLeft, evalRight)) {
        maxM0 *= 2;
        evalRight = truncation * peak - tempK.exactKBessel(2*pi/A*maxM0*Y0);
    }


    //Step 2: Do binary search
    double left = minM0;
    double right = maxM0;
    evalLeft = truncation * peak - tempK.exactKBessel(2*pi/A*left*Y0);

    //Super accuracy here doesn't really matter, but it's fast, so we go this far because we can
    while (right - left > 0.000001) {
        double center = (left + right)/2.0;
        double evalCenter = truncation * peak - tempK.exactKBessel(2*pi/A*center*Y0);

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
    double maxYStar = 0;

    //These formulas provide bounds to guarantee that the DFT result is valid
    double absRealM = abs(m.getComplex(d).real());
    double absImagM = abs(m.getComplex(d).imag());
    double Q0double = 0;
    double Q1double = (MY + absRealM)/2.0;
    if (Auxilliary::mod(-d, 4) == 1) {
        double Q0doubleOption1 = MY + absRealM;
        double Q0doubleOption2 = (MY + absImagM)/(2.0*A);
        Q0double = max(Q0doubleOption1, Q0doubleOption2);
    } else {
        Q0double = (MY + absImagM)/(2.0*A);
    }

    int Q0 = ceil(Q0double);
    int Q1 = ceil(Q1double);

    //Tweak Q0 and Q1 to be exactly what we need for exploiting symmetry
    if (Auxilliary::mod(-d, 4) == 1 && d != 3) {
        //If -d = 1 mod 4 then we need Q0/Q1 to be an even integer for the map x -> -bar(x) to be defined on test points
        if (Auxilliary::mod(Q0, Q1) == 0 && Auxilliary::mod(Q0 / Q1, 2) == 0) {
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
                while (!Auxilliary::mod(firstOptionQ0, 2 * (k + 1) == 0)) {
                    firstOptionQ0++;
                }
                int firstOptionQ1 = firstOptionQ0/(2*(k+1));

                int secondOptionQ0 = ceil(Q1 * 2 * k);
                while (!Auxilliary::mod(secondOptionQ0, 2 * k == 0)) {
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


    char nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    char reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

    //Now generate the representatives
    //then generate the orbits
    if (Auxilliary::mod(-d, 4) == 1 && d != 3) {
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
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

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
                    tuple<complex<double>, char> tup (-conj(x), reflectionSign);
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
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                tuple<complex<double>, char> tup (-conj(x), reflectionSign);
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
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                if (abs(l0) == abs(l1)) {
                    //The orbit is {x, ix, -x, -ix}
                    //So orbit/{+/-1} is {[x], [ix]}
                    //So the only proper translate is ix

                    tuple<complex<double>, char> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);
                } else {
                    //The orbit is {x, ix, -x, -ix, bar(x), ibar(x), -bar(x), -ibar(x)}
                    //So orbit/{+/-1} is {[x], [ix], [-bar(x)], [-ibar(x)]}
                    //So the proper translates are ix, -bar(x), and -ibar(x)

                    tuple<complex<double>, char> tup (I*x, nonsquareUnitSign);
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
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

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

                    tuple<complex<double>, char> tup (theta*x, nonsquareUnitSign);
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

    if (Y == MY1) {
        testPointCountY1 = 2 * Q0 * 2 * Q1;
        maxY1Star = maxYStar;
    } else {
        testPointCountY2 = 2 * Q0 * 2 * Q1;
        maxY2Star = maxYStar;
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


void BianchiMaassSearch::recursiveSearchForEigenvalues(const double leftR, const double rightR,
                                                         vector<double>& leftG,
                                                         vector<double>& rightG) {

    //Batch out the interval. Replace this spacing variable with something from Weyl law later.
    double SPACING = 1.0 / 8;
    if (rightR - leftR < SPACING) {
        vector<double> endpoints;
        double endpoint = leftR;
        while (endpoint < rightR) {
            endpoints.push_back(endpoint);
            endpoint += SPACING;
        }
        for (int i = 0; i < endpoints.size() - 1; i++) {
            vector<double> newLeftG = vector<double>();
            vector<double> newRightG = vector<double>();
            //Since I'm passing two empty G vectors, the if-statement below will get triggered each time
            //this recursion is called for the "first" time
            recursiveSearchForEigenvalues(endpoints[i], endpoints[i + 1], newLeftG, newRightG);
        }
        return;
    }

    //Replace this with something that selects Y1Matrix and Y2Matrix more carefully, dependent on r and Y0, etc.
    Y1 = 0.065;
    Y2 = 0.07;

    if (leftG.empty() || rightG.empty()) {
        if (leftG.empty() && rightG.empty()) {
            K.setRAndClear(rightR);
            //compute all the indices and points etc
            M0 = computeM0General(rightR);
            MY1 = computeMYGeneral(M0, Y1); //Y1Matrix is the smaller one, produces larger MY
            MY2 = computeMYGeneral(M0, Y2);
            computeIndexData();
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

    if (signChanges >= searchPossibleSignChanges/3) {
        //Call on the left and right intervals

        double centerR = (leftR + rightR)/2.0;
        cout << "testing left and right subintervals" << endl;
        vector<double> centerG = vector<double>();
        //In this function call, centerG gets computed and set to this variable
        recursiveSearchForEigenvalues(leftR, centerR, leftG, centerG);

        //In this function call, centerG is computed so it won't recompute
        recursiveSearchForEigenvalues(centerR, rightR, centerG, rightG);
        return;
    } else {
        cout << "Stopping. Sign change count too low." << endl;
        cout << "-----------------" << endl;
        return;
    }



}

void BianchiMaassSearch::computeIndexData() {
    indicesM0.clear();
    indexTransversal.clear();
    indexOrbitData.clear();
    indexOrbitDataModMinusOne.clear();

    // Computed bounds using Lagrange multipliers
    // aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int aUpperBound = ceil(M0 * abs(theta)/ sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int aLowerBound = -aUpperBound;

    //bub = ceil(maxN / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int bUpperBound = ceil(M0 /sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int bLowerBound = -bUpperBound;

    for (int a = aLowerBound; a <= aUpperBound; a++) {
        for (int b = bLowerBound; b <= bUpperBound; b++) {
            if (a == 0 && b == 0) {
                continue;
            }
            Index index = Index(a, b);
            if (index.getAbs(d) <= M0) {
                indicesM0.push_back(index);
            }
        }
    }

    //Sort using lambda so that I don't have to store d in every Index object... waste of memory!
    auto indexComparator = [this](const Index& index1, const Index& index2) -> bool {
        if (index1.getAbs(d) < index2.getAbs(d)) {
            return true;
        } else if (index1.getAbs(d) == index2.getAbs(d) && index1.getAngle(d) < index2.getAngle(d)) {
            return true;
        } else {
            return false;
        }
    };

    /*Sort the indices by absolute value then by angle with positive real axis.*/
    std::sort(indicesM0.begin(), indicesM0.end(), indexComparator);

    vector<Index> tempIndices;
    copy(indicesM0.begin(), indicesM0.end(), back_inserter(tempIndices));

    int rotationCoeff = symClass == 'D' || symClass == 'G' ? 1 : -1;
    int conjCoeff = symClass == 'D' || symClass == 'C' ? 1 : -1;

    while (!tempIndices.empty()) {
        Index index = tempIndices.front();
        vector<tuple<Index, int>> orbit = {tuple<Index, int>(index, 1)};

        Index tempIndex = index.rotate(d);
        int tempCoeff = 1*rotationCoeff;
        tuple<Index, int> tempTuple = {tempIndex, tempCoeff};
        while (tempIndex != index) {
            orbit.push_back(tempTuple);

            tempIndex = tempIndex.rotate(d);
            tempCoeff = tempCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
        }

        Index conjIndex = index.conj(d);
        bool conjIsInRotations = false;
        for (const auto& tup : orbit) {
            if (conjIndex == get<0>(tup)) {
                conjIsInRotations = true;
                break;
            }
        }

        if (!conjIsInRotations) {
            /*add in the rotations of the conjugate*/
            tempTuple = {conjIndex, conjCoeff};
            orbit.push_back(tempTuple);

            tempIndex = conjIndex.rotate(d);
            tempCoeff = conjCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
            while (tempIndex != conjIndex) {
                orbit.push_back(tempTuple);

                tempIndex = tempIndex.rotate(d);
                tempCoeff = tempCoeff*rotationCoeff;
                tempTuple = {tempIndex, tempCoeff};
            }

        }
        /*Save our starting index to the transversal.*/
        indexTransversal.push_back(index);

        /*Save the orbit data.*/
        indexOrbitData[index] = orbit;

        /*Compute the orbit mod +-1*/
        vector<tuple<Index,int>> orbitModMinusOne;
        vector<Index> alreadyGotten;
        for (const auto& tup : orbit) {
            Index l = get<0>(tup);
            int pmOne = get<1>(tup);
            if (std::find(alreadyGotten.begin(), alreadyGotten.end(),l) == alreadyGotten.end()) {
                //add it
                tuple<Index,int> classModMinusOne = {l,pmOne};
                alreadyGotten.push_back(l);
                alreadyGotten.push_back(Index(-l.getA(), -l.getB()));
                orbitModMinusOne.push_back(classModMinusOne);
            }
        }
        indexOrbitDataModMinusOne[index] = orbitModMinusOne;

        /*Delete the indicesM0 in the orbit from the copied list of indicesM0.*/
        for (auto itr1 : orbit) {
            Index toDelete = get<0>(itr1);
            auto itr2 = tempIndices.begin();
            while (itr2 != tempIndices.end()) {
                if (*itr2 == toDelete) {
                    tempIndices.erase(itr2);
                    break;
                }
                ++itr2;
            }
        }
    }

    /*Save the index of n = 1 + 0*I. Should be 0, but just to check.*/
    for (int itr = 0; itr < indexTransversal.size(); itr++) {
        if (indexTransversal[itr].getComplex(d) == complex<double> {1, 0}) {
            searchIndexOfNormalization = itr;
            break;
        }
    }
}



void BianchiMaassSearch::populateMatrix(MatrixID matrixId) {

    int size = indexTransversal.size();
    searchMatrices[matrixId].resize(size, size);

    for (int i = 0; i < size; i++) {
        Index m = indexTransversal[i];

#pragma omp parallel for default(none) shared(i, m, size, matrixId)
        for (int j = 0; j < size; j++) {

            Index n = indexTransversal[j];
            double entry = computeEntry(m, n, matrixId);

            searchMatrices[matrixId](i,j) = entry;
        }
    }
}

void BianchiMaassSearch::solveMatrix(MatrixID matrixId) {
    int size = indexTransversal.size();

    searchMatrixG[matrixId].resize(size);
    searchMatrixSolutions[matrixId].resize(size - 1);
    searchMatrixB[matrixId].resize(size - 1);

    MatrixXd B;
    B.resize(size - 1, size - 1);

    //Finding kernel of M is the original problem.
    //Solving Bx=b is the system after moving a column over and deleting a row.
    auto M = searchMatrices[matrixId];
    auto b = searchMatrixB[matrixId];
    auto x = searchMatrixSolutions[matrixId];
    //set b to negative of the indexOfNormalization column
    //set B to be M deleting row and column indexOfNormalization
    for (int row = 0; row < size; row++) {
        for (int column = 0; column < size; column++) {
            if (column == searchIndexOfNormalization && row != searchIndexOfNormalization) {
                //add it to b
                if (row > searchIndexOfNormalization) {
                    b(row - 1) = M(row, column);
                } else {
                    b(row) = M(row, column);
                }
            } else if (column != searchIndexOfNormalization && row != searchIndexOfNormalization) {
                //add it to tempA somehow
                int modifiedRow = (row > searchIndexOfNormalization) ? row - 1 : row;
                int modifiedColumn = (column > searchIndexOfNormalization) ? column - 1 : column;
                B(modifiedRow, modifiedColumn) = M(row, column);
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
        if (i < searchIndexOfNormalization) {
            xPrime(i) = x(i);
        } else if (i == searchIndexOfNormalization) {
            xPrime(i) = 1.0;
        } else {
            xPrime(i) = x(i - 1);
        }
    }

    //What I really want is M*xPrime = g
    searchMatrixG[matrixId] = M * xPrime;
}

vector<double> BianchiMaassSearch::getGVector() {
    int size = indexTransversal.size();
    vector<double> gVector;
    gVector.resize(2*size);

    K.extendPrecomputedRange(2 * pi / A * 1 * Y1, 2 * pi / A * M0 * maxY1Star);

    populateMatrix(MatrixID::Y1Matrix);
    populateMatrix(MatrixID::Y2Matrix);

    //solve both matrices
    //apply each matrix to the solution of the other
    solveMatrix(MatrixID::Y1Matrix);
    solveMatrix(MatrixID::Y2Matrix);

    //concatenate the output vectors (coherently!)
    for (int i = 0; i < size; i++) {
        gVector.push_back(searchMatrixG[MatrixID::Y1Matrix](i));
    }
    for (int i = 0; i < size; i++) {
        gVector.push_back(searchMatrixG[MatrixID::Y2Matrix](i));
    }

    return gVector;
}

double
BianchiMaassSearch::computeEntry(const Index &m, const Index &n, MatrixID matrixId) {
    double answer = 0;

    int pointCount = 0;
    vector<TestPointOrbitData> testPointOrbits;
    if (matrixId == MatrixID::Y1Matrix) {
        testPointOrbits = mToY1TestPointOrbits[m];
        pointCount = testPointCountY1;
    } else {
        testPointOrbits = mToY2TestPointOrbits[m];
        pointCount = testPointCountY2;
    }

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        vector<double> terms;
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;


            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * pi / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : indexOrbitDataModMinusOne[n]) {
                Index l = get<0>(tup);
                complex<double> xStar = pullback.getComplex();
                double ellTerm = get<1>(tup); //This is a_l/a_n
                ellTerm *= cos(pi/A * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = cos(pi/A * traceProduct(-I* m.getComplex(d), x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = get<0>(tup);
                char sign = get<1>(tup);
                testPointTerm += sign * cos(pi/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Auxilliary::multiPrecisionSummation(terms);
        answer *= -4.0 / pointCount;
    } else {
        vector<double> terms;
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;

            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * pi / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : indexOrbitDataModMinusOne[n]) {
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
                char sign = get<1>(tup);
                testPointTerm += sign * sin(pi/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Auxilliary::multiPrecisionSummation(terms);
        answer *= +4.0 / pointCount;
    }

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm;
        if (matrixId == MatrixID::Y1Matrix) {
            deltaTerm = Y1 * K.exactKBessel(2 * pi / A * m.getAbs(d) * Y1);
        } else {
            deltaTerm = Y2 * K.exactKBessel(2 * pi / A * m.getAbs(d) * Y2);
        }

        answer += deltaTerm;
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

double BianchiMaassSearch::findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1) {
    return -y0 * (x1 - x0)/(y1 - y0) + x0;
}

void BianchiMaassSearch::computeTestPointData() {
    mToY1TestPointOrbits.clear();
    mToY2TestPointOrbits.clear();
    for (const auto& m : indexTransversal) {
        mToY1TestPointOrbits[m] = getPointPullbackOrbits(m, Y1, MY1);
        mToY2TestPointOrbits[m] = getPointPullbackOrbits(m, Y2, MY2);
    }
}
