#include "BianchiMaassPointwise.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iterator>
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


/**
 * @brief Initializes a BianchiMaassPointwise object. Can be used for pointwise calculations.
 *
 * @param d Q(sqrt(-d)) is the field that the forms are defined for. d should be in {1,2,3,7,11,19,43,67,163}.
 * @param D 10^(-D) is the magnitude of the relativeError term of the truncation of Bianchi-Maass form Fourier-Bessel series.
 * @param symClass 'D', 'G', 'C', or 'H' are symmetry classes as defined in the literature.
 */
BianchiMaassPointwise::BianchiMaassPointwise(int d, int D, char symClass)
    : K(computeK(d, D, symClass), 0) {

    auto start = std::chrono::high_resolution_clock::now();

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

    K = KBessel(2 * PI / A * 1 * Y0, 0);

    //outputFile.open("d" + to_string(d) + "_D" + to_string(D) + "_" + symClass + "_" + to_string(Y1Parameter) + ".txt");

    truncation = pow(10,-D);
    tolerance = pow(10,-(D+6));

    cout << setprecision(16);
    cout << "A: " << A << endl;
    cout << "theta: " << theta << endl;
    cout << "Y0: " << Y0 << endl;
    cout << "truncation: " << truncation << endl;
}

BianchiMaassPointwise::~BianchiMaassPointwise() {
    outputFile.close();
}

/***************************************************
 * These functions work POINTWISE, as in, they are for seeing
 * behavior at a single point, not for finding eigenvalues.
 ***************************************************/

/**
 * @brief Conducts basic tests to see if r is a genuine eigenvalue. If Y is not specified, the method will
 * search for a well-conditioned value of Y, but it will take much longer than if you specify Y.
 *
 * @param r Eigenvalue/spectral parameter to check.
 * @param Y Height of horosphere. If Y = 0 or Y is not specified, will search for a good value of Y. Will take much longer.
 */
void BianchiMaassPointwise::checkSingleEigenvalue(const double r, const double Y) {

    //If any of these are zero or empty, intialize/reinitialize everything
    if (M0 == 0 || MY == 0 || indexTransversal.empty()) {
        K.setRAndClear(r);
        this->Y = Y;
        M0 = computeM0General(r);
        MY = computeMYGeneral(M0, Y);
        computeIndexData();
    } else {
        if (K.getR() == r && this->Y == Y) {
            //do nothing, we're all set
        } else if (K.getR() != r && this->Y != Y) {
            this->Y = Y;
            K.setRAndClear(r);
            M0 = computeM0General(r);
            MY = computeMYGeneral(M0, Y);
            computeIndexData();
        } else if (K.getR() != r && this->Y == Y) {
            K.setRAndClear(r);
            M0 = computeM0General(r);
            MY = computeMYGeneral(M0, Y);
            computeIndexData();
        } else {
            //K.getR() == r && Y != Y
            this->Y = Y;
            MY = computeMYGeneral(M0, Y);
            computeIndexData();
        }
    }

    populateMatrix();
    solveMatrix();

    int numPoints = 20;
    vector<Quaternion> samplePoints;
    vector<Quaternion> samplePointPullbacks;
    while (samplePoints.size() < numPoints) {
        double height = Y0 + (1-Y0) * rand01();
        double x0 = rand01()*10.0 - 5.0;
        double x1 = rand01()*10.0 - 5.0;
        Quaternion point = Quaternion(x0, x1, height, 0);
        Quaternion pullback = Quaternion(point);
        pullback.reduce(d);
        if ((point - pullback).abs() > 0.000001) {
            samplePoints.push_back(point);
            samplePointPullbacks.push_back(pullback);
        }
    }

    double maxDifference = 0;
    for (int i = 0; i < samplePoints.size(); i++) {
        complex<double> testPointEvaluation = evaluate(r, samplePoints[i], Y);
        complex<double> pullbackEvaluation = evaluate(r, samplePointPullbacks[i], Y);
        double difference = abs(testPointEvaluation - pullbackEvaluation);
        if (difference > maxDifference) {
            maxDifference = difference;
        }
    }

    cout << "Max of |f(z) - f(z^*)| over " << numPoints << " points: " << maxDifference << setprecision(16) << endl;
    //TODO: put hecke relation check in here
}


/**
 * @brief Returns a map of indices and Fourier-Bessel coefficients.
 *
 * The matrix from Hejhal's identity provides coefficients up to bound M0, and then we do a post-computation
 * that computes coefficients up to bound MY.
 * @param r Eigenvalue/spectral parameter to use.
 * @param Y Height of horosphere. If Y = 0 or Y is not specified, will search for a good value of Y. Will take much longer.
 * @return A map of all coefficients up to bound MY and the corresponding Fourier-Bessel coefficient.
 */
unordered_map<Index, double> BianchiMaassPointwise::computeFourierCoefficients(const double r, const double Y) {
    return unordered_map<Index, double>();
}

double BianchiMaassPointwise::evaluate(const double r, const Quaternion& z, const double Y) {

    if (z.getJ() < Y0) {
        throw(std::invalid_argument("Coefficients only good for j-part above Y0"));
    }

    if (coefficientMap.empty()) {
        K.setRAndClear(r);
        this->Y = Y;
        M0 = computeM0General(r);
        MY = computeMYGeneral(M0, Y);
        computeIndexData();
        populateMatrix();
        solveMatrix();
    }

    ArchtKBessel bess(r);

    vector<double> terms;
    for (auto itr : indexOrbitDataModMinusOne) {
        auto n = itr.first;
        auto a_n = coefficientMap[n];

        double term = a_n*z.getJ();

        //kappa(2*pi/A*|n|*y)
        term *= bess.evaluate(2 * PI / A * n.getAbs(d) * z.getJ());

        double cs = 0;
        if (symClass == 'D' || symClass == 'G' || d == 1) {
            for (auto lClass : itr.second) {
                cs += lClass.second * 2 * cos(PI/A * traceProduct(I * (lClass.first).getComplex(d),
                                                                   z.getComplex()));
            }
        } else {
            for (auto lClass : itr.second) {
                cs += lClass.second * 2 * sin(PI/A * traceProduct(I * (lClass.first).getComplex(d),
                                                                   z.getComplex()));
            }
        }

        term *= cs;

        terms.push_back(term);
    }
    double answer = Aux.multiPrecisionSummation(terms);
    return answer;
}


/**
 * @brief Computes the Sato-Tate histogram with coefficients up to MY.
 *
 * @param r Eigenvalue/spectral parameter to use.
 * @param Y Height of horosphere. If Y = 0 or Y is not specified, will search for a good value of Y. Will take much longer.
 */
void BianchiMaassPointwise::checkSatoTate(const double r, const double Y) {
    //check that coefficients for r have been computed
    //check sato tate
    //output a text file, probably?
}


/**
 * @brief Checks the Ramanujan-Petersson conjecture for the specified r and Y.
 *
 * Checks coefficients up to bound MY.
 * @param r Eigenvalue/spectral parameter to use.
 * @param Y Height of horosphere. If Y = 0 or Y is not specified, will search for a good value of Y. Will take much longer.
 * @return Returns true if the conjecture holds, false if there is a counterexample.
 */
bool BianchiMaassPointwise::checkRamanujanPetersson(const double r, const double Y) {
    //check that coefficients for r have been computed
    //verify that |a(p)| \leq 2*N(p)^(-1/2) <--ideal norm

    //make a list of ideals which are within the bounds of computed coefficients
    //for indexM0
    //if prime
    //check |a(p)| \leq 2 * N(p)^(-1/2)

    //we will compute prime elements up to the bound M0
    vector<int> rationalPrimes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};
    return false;
}

/**
 * @brief Produces the L-function of the Bianchi-Maass form for the specified r and Y.
 *
 * @param r Eigenvalue/spectral parameter to use.
 * @param Y Height of horosphere. If Y = 0 or Y is not specified, will search for a good value of Y. Will take much longer.
 */
void BianchiMaassPointwise::produceLFunction(const double r, const double Y) {
    //check that coefficients for r have been computed
    //idk
    int a = 0;
}

void BianchiMaassPointwise::clearData() {
    //clear all the vectors
    int a = 0;
}

/**
 * @brief Computes M(Y_0) for the spectral parameter r. Assumes Y0, D, A are already computed.
 * Does not interfere with other objects for computing the K-Bessel function.
 *
 * @param r Spectral parameter
 */
double BianchiMaassPointwise::computeM0General(const double r) {
    ArchtKBessel bess(r);

    //Method: use binary search to find M0 such that
    //K(2*pi/A * M0 * Y0) = 10^-D * K(max(r,1))
    //with M0 > max(r,1)/(2*pi/A*Y0)
    //Reference: JAY JORGENSON, LEJLA SMAJLOVI  ́ C, AND HOLGER THEN
    //ON THE DISTRIBUTION OF EIGENVALUES OF MAASS FORMS ON CERTAIN MOONSHINE GROUPS

    //Step 1: Make sure the equality point happens between my left and right endpoints
    double minM0 = max(r, 1.0)/(2*PI/A*Y0);
    double maxM0 = 100;

    double peak = bess.evaluate(max(r,1.0));
    double evalLeft = truncation * peak - bess.evaluate(2*PI/A*minM0*Y0);
    double evalRight = truncation * peak - bess.evaluate(2*PI/A*maxM0*Y0);
    while (!areDifferentSign(evalLeft, evalRight)) {
        maxM0 *= 2;
        evalRight = truncation * peak - bess.evaluate(2*PI/A*maxM0*Y0);
    }


    //Step 2: Do binary search
    double left = minM0;
    double right = maxM0;
    evalLeft = truncation * peak - bess.evaluate(2*PI/A*left*Y0);

    //Super accuracy here doesn't really matter, but it's fast, so we go this far because we can
    while (right - left > 0.000001) {
        double center = (left + right)/2.0;
        double evalCenter = truncation * peak - bess.evaluate(2*PI/A*center*Y0);

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
double BianchiMaassPointwise::computeMYGeneral(const double M0, const double Y) {
    return Y0/Y * M0;
}

vector<TestPointOrbitData>
BianchiMaassPointwise::getPointPullbackOrbits(const Index &m, const double Y, const double MY) {
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
    testPointCount = 2 * Q0 * 2 * Q1;

    char nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    char reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

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
                    pair<complex<double>, char> tup (-conj(x), reflectionSign);
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

                pair<complex<double>, char> tup (-conj(x), reflectionSign);
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

                    pair<complex<double>, char> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);
                } else {
                    //The orbit is {x, ix, -x, -ix, bar(x), ibar(x), -bar(x), -ibar(x)}
                    //So orbit/{+/-1} is {[x], [ix], [-bar(x)], [-ibar(x)]}
                    //So the proper translates are ix, -bar(x), and -ibar(x)

                    pair<complex<double>, char> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    tup.first = -conj(x);
                    tup.second = reflectionSign;
                    orbit.properTranslatesModSign.push_back(tup);

                    tup.first = -I*conj(x);
                    tup.second = reflectionSign * nonsquareUnitSign;
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

                    pair<complex<double>, char> tup (theta*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    tup.first = theta*theta*x;
                    tup.second = 1; //theta*theta is a square unit
                    orbit.properTranslatesModSign.push_back(tup);

                    //magic quantity that tells me about the real part of x
                    int discriminant = 2*Q1*l0 + Q1*l1 - 2*Q0*Q1;
                    bool orbitIsHalfSize = (l0 == 0 || l1 == 0 || l0 + l1 == 0 || discriminant == 0);

                    if (!orbitIsHalfSize) {
                        //For this comment let zeta=e^(2pi i /6) be represented by w
                        //The orbit is {x, wx, w^2x, -x, -wx, -w^2x, -bar(x), -wbar(x), -w^2bar(x), bar(x), wbar(x), w^2bar(x)}
                        //So orbit/{+/-1} is {[x], [wx], [w^2x], [-bar(x)], [-wbar(x)], [-w^2bar(x)]}
                        //So the proper translates are wx, w^2x, -bar(x), -wbar(x), -w^2bar(x)
                        tup.first = -conj(x);
                        tup.second = reflectionSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        tup.first = -theta*conj(x);
                        tup.second = reflectionSign * nonsquareUnitSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        tup.first = -theta*theta*conj(x);
                        tup.second = reflectionSign; //theta*theta is a square unit
                        orbit.properTranslatesModSign.push_back(tup);
                    }
                }
            }
        }
    }

    return answer;
}


/***************************************************
 * Private methods used in POINTWISE calculations.
 ***************************************************/

void BianchiMaassPointwise::populateMatrix() {

    int size = indexTransversal.size();
    matrix.resize(size, size);

    for (int i = 0; i < size; i++) {
        Index m = indexTransversal[i];
        watch(m);

        testPointOrbits = getPointPullbackOrbits(m, Y, MY);
        double maxYStar = 0;
        for (const auto& orbit : testPointOrbits) {
            double Y = orbit.representativePullback.getJ();
            if (Y > maxYStar) {
                maxYStar = Y;
            }
        }

        K.extendPrecomputedRange(2 * PI / A * M0 * maxYStar);
#pragma omp parallel for default(none) shared(matrix, indexTransversal, i, m, size)
        for (int j = 0; j < size; j++) {

            Index n = indexTransversal[j];
            double entry = computeEntry(m, n);

            matrix(i, j) = entry;
        }
    }
}


double BianchiMaassPointwise::computeEntry(const Index &m, const Index &n) {
    double answer = 0;

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        vector<double> terms;
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;


            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * PI / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : indexOrbitDataModMinusOne[n]) {
                Index l = tup.first;
                complex<double> xStar = pullback.getComplex();
                double ellTerm = tup.second; //This is a_l/a_n
                ellTerm *= cos(PI/A * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = cos(PI/A * traceProduct(-I* m.getComplex(d), x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                char sign = tup.second;
                testPointTerm += sign * cos(PI/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Aux.multiPrecisionSummation(terms);
        answer *= -4.0 / testPointCount;
    } else {
        vector<double> terms;
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;

            double yStar = pullback.getJ();
            double term = yStar * K.approxKBessel(2 * PI / A * n.getAbs(d) * yStar);

            double indexTerm = 0;
            for (const auto& tup : indexOrbitDataModMinusOne[n]) {
                Index l = tup.first;
                complex<double> xStar = pullback.getComplex();
                double ellTerm = tup.second; //This is a_l/a_n
                ellTerm *= sin(PI/A * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = sin(PI/A * traceProduct(-I* m.getComplex(d), x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                char sign = tup.second;
                testPointTerm += sign * sin(PI/A * traceProduct(-I* m.getComplex(d), etaX));
            }
            term *= testPointTerm;

            terms.push_back(term);
        }
        answer = Aux.multiPrecisionSummation(terms);
        answer *= +4.0 / testPointCount;
    }

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm = Y * K.approxKBessel(2 * PI / A * m.getAbs(d) * Y);

        answer += deltaTerm;
    }
    return answer;
}


void BianchiMaassPointwise::solveMatrix() {
    Eigen::Matrix<double, Eigen::Dynamic, 1> x;
    Eigen::Matrix<double, Eigen::Dynamic, 1> b;
    MatrixXd subMatrix;

    coefficientMap.clear();

    long size = indexTransversal.size() - 1;
    x.resize(size);
    b.resize(size);
    subMatrix.resize(size, size);

    //set b to negative of the indexOfNormalization column
    //set tempA to be matrixA deleting row and column indexOfNormalization
    for (int row  = 0; row < indexTransversal.size(); row++) {
        for (int column = 0; column < indexTransversal.size(); column++) {
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
                subMatrix(modifiedRow, modifiedColumn) = matrix(row, column);
            }
        }
    }


    //make b negative
    b = -b;

    //solve the equation tempA * x = b for x
    x = subMatrix.colPivHouseholderQr().solve(b);

    //put the coefficients in the coefficient map, remembering to add back in the coefficients for indexOfNormalization
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index index = indexTransversal[i];
        if (i == indexOfNormalization) {
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = tup.second;
                coefficientMap[tup.first] = coeff;
            }
        } else {
            //look through the orbit and assign the same value to all of them using the +/- factor
            int modifiedi = (i > indexOfNormalization) ? i - 1 : i;
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = x(modifiedi);
                coeff *= tup.second;
                coefficientMap[tup.first] = coeff;
            }
        }
    }
}


/**
 * @brief Computes indices up to bounds M0 and MY. Also computes orbit data from symmetry relations.
 * Assumes that M0 and MY are already computed.
 */
void BianchiMaassPointwise::computeIndexData() {
    indicesM0.clear();
    indicesMY.clear();
    indexTransversalMY.clear();
    indexTransversal.clear();
    indexOrbitData.clear();
    indexOrbitDataModMinusOne.clear();

    // Computed bounds using Lagrange multipliers
    // aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int aUpperBound = ceil(MY * abs(theta) / sqrt(pow(abs(theta), 2) - pow(theta.real(), 2)));
    int aLowerBound = -aUpperBound;

    //bub = ceil(maxN / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int bUpperBound = ceil(MY / sqrt(pow(abs(theta), 2) - pow(theta.real(), 2)));
    int bLowerBound = -bUpperBound;

    for (int a = aLowerBound; a <= aUpperBound; a++) {
        for (int b = bLowerBound; b <= bUpperBound; b++) {
            if (a == 0 && b == 0) {
                continue;
            }
            Index index = Index(a, b);
            if (index.getAbs(d) <= MY) {
                indicesMY.push_back(index);
                if (index.getAbs(d) <= M0) {
                    indicesM0.push_back(index);
                }
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
    std::sort(indicesMY.begin(), indicesMY.end(), indexComparator);

    vector<Index> tempIndices;
    copy(indicesMY.begin(), indicesMY.end(), back_inserter(tempIndices));

    int rotationCoeff = symClass == 'D' || symClass == 'G' ? 1 : -1;
    int conjCoeff = symClass == 'D' || symClass == 'C' ? 1 : -1;

    while (!tempIndices.empty()) {
        Index index = tempIndices.front();
        vector<pair<Index, int>> orbit = {pair<Index, int>(index, 1)};

        Index tempIndex = index.rotate(d);
        int tempCoeff = 1*rotationCoeff;
        pair<Index, int> tempTuple = {tempIndex, tempCoeff};
        while (tempIndex != index) {
            orbit.push_back(tempTuple);

            tempIndex = tempIndex.rotate(d);
            tempCoeff = tempCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
        }

        Index conjIndex = index.conj(d);
        bool conjIsInRotations = false;
        for (const auto& tup : orbit) {
            if (conjIndex == tup.first) {
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
        if (index.getAbs(d) <= M0) {
            indexTransversal.push_back(index);
            indexTransversalMY.push_back(index);
        } else {
            indexTransversalMY.push_back(index);
        }


        /*Save the orbit data.*/
        indexOrbitData[index] = orbit;

        /*Compute the orbit mod +-1*/
        vector<pair<Index,int>> orbitModMinusOne;
        vector<Index> alreadyGotten;
        for (const auto& tup : orbit) {
            Index l = tup.first;
            int pmOne = get<1>(tup);
            if (std::find(alreadyGotten.begin(), alreadyGotten.end(),l) == alreadyGotten.end()) {
                //add it
                pair<Index,int> classModMinusOne = {l,pmOne};
                alreadyGotten.push_back(l);
                alreadyGotten.push_back(Index(-l.getA(), -l.getB()));
                orbitModMinusOne.push_back(classModMinusOne);
            }
        }
        indexOrbitDataModMinusOne[index] = orbitModMinusOne;

        /*Delete the indicesM0 in the orbit from the copied list of indicesM0.*/
        for (auto itr1 : orbit) {
            Index toDelete = itr1.first;
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
            indexOfNormalization = itr;
            break;
        }
    }
}


double BianchiMaassPointwise::rand01() {
    return std::rand()/ (RAND_MAX + 1.0);
}

double BianchiMaassPointwise::traceProduct(const complex<double> &z, const complex<double> &w) {
    auto answer = z*w + conj(z)*conj(w);
    return answer.real();
}

bool BianchiMaassPointwise::areDifferentSign(const double &a, const double &b) {
    return (a > 0 && b < 0) || (a < 0 && b > 0);
}

double BianchiMaassPointwise::computeK(int d, int D, char symClass) {
    Od = ImaginaryQuadraticIntegers(d);

    A = Od.getA();
    theta = Od.getTheta();
    Y0 = Od.getY0();

    return 2 * PI/A * 1 * Y0;
}




