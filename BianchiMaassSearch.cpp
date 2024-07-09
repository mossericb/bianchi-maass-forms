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
#include <fstream>
#include "archtKBessel.h"
#include <boost/align/aligned_allocator.hpp>
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

BianchiMaassSearch::BianchiMaassSearch(string mode, int d, double D, char symClass) {

    //Check that mode is valid
    if (mode != "coarse" && mode != "medium" && mode != "fine") {
        throw std::invalid_argument("mode should be coarse, medium, or fine");
    }
    this->mode = mode;

    //Check that d is valid
    vector<int> classNumberOne = {1,2,3,7,11,19,43,67,163};
    bool dIsClassNumberOne = find(classNumberOne.begin(), classNumberOne.end(), d) != classNumberOne.end();
    if (!dIsClassNumberOne) {
        throw(std::invalid_argument("d should be one of {1,2,3,7,11,19,43,67,163}"));
    }
    this->d = d;

    //Check that symclass is valid
    vector<char> symClasses = {'D','G','C','H'};
    bool symClassIsCorrect = find(symClasses.begin(), symClasses.end(), symClass) != symClasses.end();
    if (!symClassIsCorrect) {
        throw(std::invalid_argument("symClass should be one of D,G,C,H"));
    }
    this->symClass = symClass;

    /*
     * These are values that determine what precision the calculation runs at.
     * Determined experimentally to achieve as good results as possible within
     * reasonable computation time spans.
     */
    dToDMap[1] = {6,6,8,16};
    dToDMap[2] = {6,6,8,16};
    dToDMap[3] = {6,6,8,16};
    dToDMap[7] = {6,6,8,16};
    dToDMap[11] = {6,6,8,16};
    dToDMap[19] = {5,6,8,16};
    dToDMap[43] = {4,4,4,6};
    dToDMap[67] = {3,3,3,4};
    dToDMap[163] = {2,2,2,3};

    this->D = D;
    if (D == 0) {
        autoPrecision = true;
    } else {
        autoPrecision = false;
        truncation = pow(10.0, -(this->D));
    }

    setUpOutputLogFiles();

    Od = ImaginaryQuadraticIntegers(d);
    A = Od.getA();
    theta = Od.getTheta();
    Y0 = Od.getY0();
    twoPiOverA = 2 * pi / A;
    piOverA = pi/A;

    dToPrimes[1] = {Index(1,1), Index(2,1), Index(3,0)};
    dToPrimes[2] = {Index(0,1), Index(1,1), Index(3,1)};
    dToPrimes[3] = {Index(1,1), Index(2,0), Index(2,1)};
    dToPrimes[7] = {Index(0,1), Index(-1,2), Index(3,0)};
    dToPrimes[11] = {Index(0,1), Index(2,0), Index(1,1)};
    dToPrimes[19] = {Index(2,0), Index(0,1), Index(1,1)};
    dToPrimes[43] = {Index(2,0), Index(3,0), Index(0,1)};
    dToPrimes[67] = {Index(2,0), Index(3,0), Index(0,1)};
    dToPrimes[163] = {Index(2,0), Index(3,0), Index(5,0)};
}


void BianchiMaassSearch::coarseSearchForEigenvalues(const double leftLambda, const double rightLambda) {
    if (leftLambda >= rightLambda) {
        throw(std::invalid_argument("leftR should be less than rightR"));
    }

    auto intervalsFromFile = getIntervalsForCoarseSearch(leftLambda, rightLambda);

    stack<pair<double,double>> intervals;

    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    double lastPrecisionRecompute = intervals.top().first;
    if (autoPrecision) {
        computeMaximumD(lastPrecisionRecompute, 180);
        std::cout << "D = " << D << std::endl;
    }

    KBessel* leftLambdaK = new KBessel(twoPiOverA * 1 * Y0, intervals.top().first);
    KBessel* rightLambdaK = new KBessel(twoPiOverA * 1 * Y0, intervals.top().second);
    while (!intervals.empty()) {
        double left = intervals.top().first;
        double right = intervals.top().second;
        intervals.pop();

        if (autoPrecision) {
            if (lastPrecisionRecompute + 1 < left) {
                computeMaximumD(left, 180);
                lastPrecisionRecompute = left;
                std::cout << "D = " << D << std::endl;
            }
        }

        rightLambdaK = new KBessel(twoPiOverA * 1 * Y0, right);
        if (possiblyContainsEigenvalue(left, right, leftLambdaK, rightLambdaK)) {
            coarseOutputFile << std::setprecision(16) << "[" << left << ", " << right << "]" << std::endl;
        }
        delete leftLambdaK;
        leftLambdaK = rightLambdaK;

        coarseOutputFile << std::setprecision(16) << "Complete up to " << right << std::endl;
    }
    delete rightLambdaK;
}

void BianchiMaassSearch::mediumSearchForEigenvalues() {
    /*
     * The coarse log looks like
     *
     * Complete up to 6.000000
     * [6.5, 6.75]
     * Complete up to 7.0000
     *
     * So the medium log should also look like
     *
     * Complete up to 6.75
     * [6.845454545, 6.84545455]
     * Complete up to 7.0000
     *
     * Then the fine log will comb through this and grab the [inter, vals]
     */


    auto intervalsFromFile = getIntervalsForMediumSearch();

    stack<pair<double, double>> intervals;
    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    /*
     *
     *
     * refine likely intervals down to a specified narrowness
     *
     *
     */

    //Magic number!
    double weylEigenvalueCount = 0.0005;

    map<double, KBessel*> KBessMap;

    while (!intervals.empty()) {
        stack<pair<double, double>> spawnedIntervals;
        double upper = intervals.top().second;
        spawnedIntervals.push(intervals.top());
        intervals.pop();

        if (autoPrecision) {
            computeMaximumD(spawnedIntervals.top().first, 180);
            std::cout << "D = " << D << std::endl;
        }

        while (!spawnedIntervals.empty()) {
            auto interval = spawnedIntervals.top();
            spawnedIntervals.pop();

            bool zoomIn = true;
            double weylRightEndpoint = Od.eigenvalueIntervalRightEndpoint(interval.first, weylEigenvalueCount);

            //Cutoff interval should be not wider than 0.001
            weylRightEndpoint = min(weylRightEndpoint, interval.first + 0.001);

            //Cutoff interval should be not narrower than 0.00000001
            weylRightEndpoint = max(weylRightEndpoint, interval.first + 0.00000001);
            if (interval.second <= weylRightEndpoint) {
                mediumOutputFile << std::setprecision(16) << "[" << interval.first << ", " << interval.second << "]" << std::endl;
                cout << "Likely eigenvalue " << std::setprecision(16) << "[" << interval.first << ", " << interval.second << "]" << std::endl;
            } else {
                double left = interval.first;
                double right = interval.second;
                double center = (left + right) / 2;

                if (KBessMap.find(left) == KBessMap.end()) {
                    KBessMap[left] = new KBessel(twoPiOverA * 1 * Y0, left);
                }
                if (KBessMap.find(center) == KBessMap.end()) {
                    KBessMap[center] = new KBessel(twoPiOverA * 1 * Y0, center);
                }
                if (KBessMap.find(right) == KBessMap.end()) {
                    KBessMap[right] = new KBessel(twoPiOverA * 1 * Y0, right);
                }

                bool leftInterval = possiblyContainsEigenvalue(left, center, KBessMap[left], KBessMap[center]);
                bool rightInterval = possiblyContainsEigenvalue(center, right, KBessMap[center], KBessMap[right]);
                if (!leftInterval && !rightInterval) {
                    delete KBessMap[center];
                    KBessMap.erase(center);
                    continue;
                } else if (leftInterval && !rightInterval) {
                    spawnedIntervals.emplace(left, center);
                    continue;
                } else if (!leftInterval && rightInterval) {
                    spawnedIntervals.emplace(center, right);
                    continue;
                } else {
                    spawnedIntervals.emplace(left, center);
                    spawnedIntervals.emplace(center, right);
                    continue;
                }
            }
        }
        for (auto itr = KBessMap.begin(); itr != KBessMap.end(); ) {
            if (itr->first < upper) {
                delete itr->second;
                KBessMap.erase(itr++);
            } else {
                ++ itr;
            }
        }
        mediumOutputFile << "Complete up to " << std::setprecision(16) << upper << std::endl;
    }

    for (auto itr = KBessMap.begin(); itr != KBessMap.end(); ) {
        delete itr->second;
        KBessMap.erase(itr++);
    }

}

void BianchiMaassSearch::fineSearchForEigenvalues() {

    auto intervalsFromFile = getIntervalsForFineSearch();

    stack<pair<double, double>> intervals;
    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    while(!intervals.empty()) {
        auto interval = intervals.top();
        intervals.pop();

        if (autoPrecision) {
            computeMaximumD(intervals.top().first, 180);
            D *= 3.0/4;
            truncation = pow(10.0, -D);
            std::cout << "D = " << D << std::endl;
        }

        auto firstPassOutput = fineSecantMethod(interval.first, interval.second);
        auto firstPassHecke = get<0>(firstPassOutput);
        double lastHecke = firstPassHecke[firstPassHecke.size() - 1].second;
        if (lastHecke > 1) { //not modular
            fineOutputFile << "Complete up to " << std::setprecision(16) << interval.second << std::endl;
            continue;
        }

        double firstLambda = get<1>(firstPassOutput);
        double secondLambda = get<2>(firstPassOutput);

        if (autoPrecision) {
            computeMaximumD(intervals.top().first, 180);
            std::cout << "D = " << D << std::endl;
        }

        auto secondPassOutput = fineSecantMethod(firstLambda, secondLambda);
        auto secondPassHecke = get<0>(secondPassOutput);
        lastHecke = secondPassHecke[secondPassHecke.size() - 1].second;
        if (lastHecke > 1) {
            fineOutputFile << "Complete up to " << std::setprecision(16) << interval.second << std::endl;
            continue;
        }

        const std::string directory = "Output/Final/"; // Change this to the desired directory
        const std::string prefix = to_string(d) + "_"
                                   + symClass + "_";
        createOutputDirectory(directory);
        int number = findMaxFileNumber(directory, prefix);
        std::string outputFilename = directory
                                     + prefix
                                     + to_string(number + 1)
                                     + ".txt";

        ofstream out;
        out.open(outputFilename, std::ofstream::out | std::ofstream::app);

        out << "First pass\n";

        for (auto pair : firstPassHecke) {
            out << std::setprecision(16) << pair.first << ", " << pair.second << '\n';
        }

        out << "Second pass\n";

        for (auto pair : secondPassHecke) {
            out << std::setprecision(16) << pair.first << ", " << pair.second << '\n';
        }
        out << std::flush;

        vector<pair<double, double>> combined;
        for (auto pair : firstPassHecke) {
            combined.push_back(pair);
        }
        for (auto pair : secondPassHecke) {
            combined.push_back(pair);
        }

        std::sort(combined.begin(), combined.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });

        string r0;
        string r1;
        string r2;
        out << "Least hecke:\n";
        for (int i = 0; i < min((int)combined.size(), 3); i++) {
            out << std::setprecision(16) << combined[i].first << ", " << combined[i].second << std::endl;
            if (i == 0) {
                r0 = to_string(combined[i].first);
            } else if (i == 1) {
                r1 = to_string(combined[i].first);
            } else if (i == 2) {
                r2 = to_string(combined[i].first);
            }
        }
        out << "These digits are correct:\n";
        int size = r0.length();
        size = min(size, (int)r1.length());
        size = min(size, (int)r2.length());
        for (int n = 0; n < size; n++) {
            if (r0[n] == r1[n] && r0[n] == r2[n]) {
                out << r0[n];
                continue;
            }
        }


        fineOutputFile << "Complete up to " << std::setprecision(16) << interval.second << std::endl;
    }
}


bool BianchiMaassSearch::possiblyContainsEigenvalue(const double leftLambda, const double rightLambda, KBessel *leftLambdaK, KBessel *rightLambdaK) {
    //Assume that I will have to recompute everything multiple times here, so don't bother checking otherwise

    //First time through, compute the condition number everywhere to get an idea of how low you can go
    //In subsequent calls, start with Y1 and Y2 at the top end of the range and adjust a few times until
    //the condition number is close to the computed ballpark

    std::cout << std::setprecision(16)
              << "[" << leftLambda << ", " << rightLambda << "]" << std::endl;

    double M0 = computeM0General(rightLambda);

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

    auto pairOfYs = computeTwoWellConditionedY(leftLambdaK, rightLambdaK, leftLambda, rightLambda, M0, indexTransversal);

    double Y1 = pairOfYs.first;
    double Y2 = pairOfYs.second;

    /*
     * Y1 and Y2 have been computed. Proceed with computing the actual matrices
     */
    double MY1 = computeMYGeneral(M0, Y1);
    double MY2 = computeMYGeneral(M0, Y2);

    double maxYStar = 0;

    map<Index, vector<TestPointOrbitData>> mToY1TestPointOrbits;
    map<Index, vector<TestPointOrbitData>> mToY2TestPointOrbits;

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

    leftLambdaK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    MatrixXd matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *leftLambdaK);
    MatrixXd matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *leftLambdaK);

    auto solAndCondY1Lambda1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2Lambda1 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    rightLambdaK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *rightLambdaK);
    matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *rightLambdaK);

    auto solAndCondY1Lambda2 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2Lambda2 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    std::cout << solAndCondY1Lambda1.second << " " << solAndCondY2Lambda1.second << " ";
    std::cout << solAndCondY1Lambda2.second << " " << solAndCondY2Lambda2.second << std::endl;

    /*
     * When you have to start over, for whatever reason, this means that condition numbers
     * blow up in that interval. The method resets the condition number search goal, but I
     * don't want this goal to stay artificially high. So whenever the algorithm has to start over,
     * it sets this flag to signal that next time it should re-restart over to make the condition
     * search goal lower again.
     */

    vector<double> coeffDiffLambda1(indexTransversal.size(), NAN);
    vector<double> coeffDiffLambda2(indexTransversal.size(), NAN);

    for (int i = 0; i < indexTransversal.size(); i++) {
        coeffDiffLambda1[i] = solAndCondY1Lambda1.first(i) - solAndCondY2Lambda1.first(i);
        coeffDiffLambda2[i] = solAndCondY1Lambda2.first(i) - solAndCondY2Lambda2.first(i);
    }

    int signChanges = countSignChanges(coeffDiffLambda1, coeffDiffLambda2);
    std::cout << signChanges << "/" << coeffDiffLambda1.size() << std::endl;

    int possibleSignChanges = coeffDiffLambda1.size() - 1;
    double proportion = ((double)signChanges)/possibleSignChanges;
    if (proportion >=  0.5) {
        return true;
    } else {
        return false;
    }
}

double BianchiMaassSearch::computeM0General(const double lambda) {
    KBessel bess(1, lambda);

    //Method: use binary search to find M0 such that
    //K(2*pi/A * M0 * Y0) = 10^-D * K(max(r,1))
    //with M0 > max(r,1)/(2*pi/A*Y0)
    //Reference: JAY JORGENSON, LEJLA SMAJLOVI  ́ C, AND HOLGER THEN
    //ON THE DISTRIBUTION OF EIGENVALUES OF MAASS FORMS ON CERTAIN MOONSHINE GROUPS

    //Step 1: Make sure the equality point happens between my left and right endpoints
    double minM0 = 0;
    if (lambda <= 1) {
        minM0 = 1.0/(2 * pi / A * Y0);
    } else {
        double r = sqrt(lambda - 1);
        minM0 = max(r, 1.0) / (2 * pi / A * Y0);
    }
    double maxM0 = 700.0 / (2 * pi / A * Y0);

    double peak = 0;
    if (lambda <= 1) {
        peak = bess.exactKBessel(2*pi/A*1*Y0);
    } else {
        double r = sqrt(lambda - 1);
        peak = bess.exactKBessel(max(r, 1.0));
    }

    double evalLeft = truncation * peak - bess.exactKBessel(2*pi/A*minM0*Y0);
    double evalRight = truncation * peak - bess.exactKBessel(2*pi/A*maxM0*Y0);
    while (!areDifferentSign(evalLeft, evalRight)) {
        maxM0 *= 2;
        evalRight = truncation * peak - bess.exactKBessel(2*pi/A*maxM0*Y0);
    }


    //Step 2: Do binary search
    double left = minM0;
    double right = maxM0;
    evalLeft = truncation * peak - bess.exactKBessel(2*pi/A*left*Y0);

    //Super accuracy here doesn't really matter, but it's fast, so we go this far because we can
    while (right - left > 0.000001) {
        double center = (left + right)/2.0;
        double evalCenter = truncation * peak - bess.exactKBessel(2*pi/A*center*Y0);

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

unsigned long long int BianchiMaassSearch::countPointPullbackOrbits(const Index &m, double Y, double MY) {
    unsigned long long int answer = 0;

    auto Q0Q1 = computeQ0Q1(m, MY);
    int Q0 = Q0Q1.first;
    int Q1 = Q0Q1.second;

    //Now generate the representatives
    //then generate the orbits
    if (Auxiliary::mod(-d, 4) == 1 && d != 3) {
        //iterate to make orbit representatives
        answer = Q1 * (Q0 + 1);

    } else if (d != 1) {
        //So -d = 2,3 mod 4 and d != 1

        //iterate to make orbit representatives
        answer = Q0 * Q1;

    } else if (d == 1) {

        answer = Q0 * Q1;

    } else {
        //d == 3
        //
        //In order for rotation to be defined on the set of testpoints, we do not use the shifted version
        //WWWWWWWAAAAARRRNNNIIINNNGGG
        //THIS CODE IS A MESS, IT'S WRONG AND WRONG AGAIN, NEEDS SERIOUS MATHEMATICAL REWORK

        throw std::invalid_argument("Code not ready yet for d = 3");
    }

    return answer;
}

vector<TestPointOrbitData>
BianchiMaassSearch::getPointPullbackOrbits(const Index &m, const double Y, const double MY) {
    vector<TestPointOrbitData> answer;

    //These formulas provide bounds to guarantee that the DFT result is valid
    double absRealM = abs(m.getComplex(d).real());
    double absImagM = abs(m.getComplex(d).imag());
    double Q0double = (MY + absImagM)/A;
    double Q1double = MY + absRealM;

    int Q0 = floor(Q0double) + 1;
    int Q1 = floor(Q1double) + 1;

    //Method 1, only quotient by negation
    long method1 = 0;
    if (Aux.mod(Q0, 2) == 1 and Aux.mod(Q1,2) == 1) {
        method1 = (Q0*Q1 - 1)/2 + 1;
    } else {
        method1 = Q0*Q1/2;
    }

    //Method 2, quotient by reflection
    long method2 = 0;
    int method2AdjustedQ0 = Q0;
    int method2AdjustedQ1 = Q1;

    if (method2AdjustedQ0 % 2 == 1) {
        method2AdjustedQ0++;
    }
    if (method2AdjustedQ1 % 2 == 1) {
        method2AdjustedQ1++;
    }

    if (method2AdjustedQ0 % method2AdjustedQ1 == 0 && (method2AdjustedQ0/method2AdjustedQ1) % 2 == 0) {
        //do nothing
    } else {
        if ((1.0*method2AdjustedQ0)/method2AdjustedQ1 <= 2) {
            method2AdjustedQ0 = 2 * method2AdjustedQ1;
        } else {
            int k = 2;
            while ( !(k <= ((1.0*Q0)/Q1) && ((1.0*Q0)/Q1) <= k + 2) ) {
                k += 2;
            }

            //either make the ratio go up to k+2 by pumping up Q0
            int goUpQ0 = (k + 2) * method2AdjustedQ1;
            int goUpQ1 = method2AdjustedQ1;

            //or make the ratio go down to k by pumping up Q1
            int goDownQ0 = method2AdjustedQ0;
            while (goDownQ0 % k != 0) {
                goDownQ0 += 2;
            }

            int goDownQ1 = goDownQ0/k;

            if (goUpQ0 * goUpQ1 > goDownQ0 * goDownQ1) {
                method2AdjustedQ0 = goDownQ0;
                method2AdjustedQ1 = goDownQ1;
            } else {
                method2AdjustedQ0 = goUpQ0;
                method2AdjustedQ1 = goUpQ1;
            }
        }
    }

    method2 = method2AdjustedQ0 * method2AdjustedQ1 / 4;

    //Method 3, quotient by negation and then reflection as much as possible
    long method3 = 0;

    int method3AdjustedQ0 = Q0;
    int method3AdjustedQ1 = Q1;

    if (method2AdjustedQ1 % 2 == 0) {
        method3AdjustedQ1++;
    }

    if (method3AdjustedQ0 % method3AdjustedQ1 == 0) {
        // do nothing
    } else {
        if ((1.0*method3AdjustedQ0)/method3AdjustedQ1 <= 1) {
            method3AdjustedQ0 = method3AdjustedQ1;
        } else {
            int k = 1;
            while ( !(k <= (1.0*Q0)/Q1 && (1.0*Q0)/Q1 <= k + 1)) {
                k++;
            }

            //either make the ratio go up to k+1 by pumping up Q0
            int goUpQ0 = (k + 1) * method3AdjustedQ1;
            int goUpQ1 = method3AdjustedQ1;

            //or make the ratio go down to k by pumping up Q1
            int goDownQ0 = Q0;
            while (goDownQ0 % k != 0 || (goDownQ0/k) % 2 != 1) {
                goDownQ0++;
            }

            int goDownQ1 = goDownQ0/k;

            if (goUpQ0*goUpQ1 + goUpQ0 + goUpQ1 > goDownQ0*goDownQ1 + goDownQ0 + goDownQ1) {
                method3AdjustedQ0 = goDownQ0;
                method3AdjustedQ1 = goDownQ1;
            } else {
                method3AdjustedQ0 = goUpQ0;
                method3AdjustedQ1 = goUpQ1;
            }
        }
    }
    method3 = method3AdjustedQ0*method3AdjustedQ1 + method3AdjustedQ0 + method3AdjustedQ1;
    method3 /= 4;

    char method = 1;
    if (method2 < method1 && method2 < method3) {
        method = 2;
        Q0 = method2AdjustedQ0;
        Q1 = method2AdjustedQ1;
    } else if (method3 < method1 && method3 < method2) {
        method = 3;
        Q0 = method3AdjustedQ0;
        Q1 = method3AdjustedQ1;
    } else {
        //Q0 and Q1 stay the same
    }
    double xi0 = Q0 % 2 == 0 ? 0.5 : 0;
    double xi1 = Q1 % 2 == 0 ? 0.5 : 0;
    int L0 = Q0 % 2 == 0 ? 1 - Q0/2 : (1 - Q0)/2;
    int U0 = Q0 % 2 == 0 ? Q0/2 : (Q0 - 1)/2;
    int L1 = Q1 % 2 == 0 ? 1 - Q1/2 : (1 - Q1)/2;
    int U1 = Q1 % 2 == 0 ? Q1/2 : (Q1 - 1)/2;

    short nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    short reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

    if (Aux.mod(-d, 4) == 1) {
        for (int l0 = L0; l0 <= U0; l0++) {
            for (int l1 = L1; l1 <= U1; l1++) {
                complex<double> x = (l0 - xi0)/(1.0*Q0) + theta * (l1 - xi1)/(1.0*Q1);
            }
        }


    } else {

    }

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


double BianchiMaassSearch::traceProduct(complex<double> z, complex<double> w) {
    double answer = 2*(z*w).real(); //this is the same as z*w + bar(z*w)
    return answer;
}



bool BianchiMaassSearch::areDifferentSign(double a, double b) {
    return (a > 0 && b < 0) || (a < 0 && b > 0);
}


MatrixXd BianchiMaassSearch::produceMatrix(const vector<Index> &indexTransversal,
                                           map<Index, vector<TestPointOrbitData>> &mToTestPointData,
                                           map<Index, vector<pair<Index, int>>> &ntoIndexOrbitData,
                                           double Y,
                                           KBessel &K) {

    int size = indexTransversal.size();
    MatrixXd answer;
    answer.resize(size, size);

    unsigned long long int terms = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Index m = indexTransversal[i];
            terms += mToTestPointData[m].size();
        }
    }

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic) default(none) shared(indexTransversal, size, answer, mToTestPointData, ntoIndexOrbitData,Y, K)
    for (int i = 0; i < size; i++) {
        Index m = indexTransversal[i];
        for (int j = 0; j < size; j++) {
            Index n = indexTransversal[j];
            double entry = computeEntry(m, n, K, mToTestPointData[m], Y, ntoIndexOrbitData[n]);

            answer(i,j) = entry;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);

    nanosecondsPerTerm = ((double)duration.count())/((double)terms);

    watch(nanosecondsPerTerm);

    return answer;
}

pair<MatrixXd, double> BianchiMaassSearch::solveMatrix(const MatrixXd &matrix, const vector<Index> &indexTransversal,
                                         int indexOfNormalization) {
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
    auto nComplex = n.getComplex(d);
    auto mComplex = m.getComplex(d);
    auto nAbs = abs(nComplex);
    auto mAbs = abs(mComplex);
    auto twoPiOverATimesNAbs = twoPiOverA * nAbs;
    size_t size = mTestPointOrbits.size();

    unsigned long long int pointCount = 0;
    for (const auto& orbit : mTestPointOrbits) {
        pointCount += (orbit.properTranslatesModSign.size() + 1);
    }
    pointCount *= 2;

    vector<double> yStars;
    vector<double> bessels;
    vector<double> xStarTerms;
    vector<double> xTerms;
    yStars.reserve(size);
    bessels.reserve(size);
    xStarTerms.reserve(size);
    xTerms.reserve(size);

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        for (const auto& testPointOrbit : mTestPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;
            auto xStar = pullback.getComplex();

            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            yStars.push_back(yStar);
            double bess = K.approxKBessel(twoPiOverATimesNAbs * yStar);
            bessels.push_back(bess);

            double indexTerm = 0;
            for (const auto& tup : nIndexOrbitDataModSign) {
                Index l = tup.first;
                double ellTerm = tup.second; //This is a_l/a_n
                ellTerm *= cos(piOverA * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            xStarTerms.push_back(indexTerm);

            double testPointTerm = cos(piOverA * traceProduct(-I* mComplex, x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                short sign = tup.second;
                testPointTerm += sign * cos(piOverA * traceProduct(-I* mComplex, etaX));
            }
            xTerms.push_back(testPointTerm);
        }
    } else {
        for (const auto& testPointOrbit : mTestPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;
            auto xStar = pullback.getComplex();

            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            yStars.push_back(yStar);
            double bess = K.approxKBessel(twoPiOverATimesNAbs * yStar);
            bessels.push_back(bess);

            double indexTerm = 0;
            for (const auto& tup : nIndexOrbitDataModSign) {
                Index l = tup.first;
                double ellTerm = tup.second; //This is a_l/a_n
                ellTerm *= sin(piOverA * traceProduct(I* l.getComplex(d), xStar));
                indexTerm += ellTerm;
            }
            xStarTerms.push_back(indexTerm);

            double testPointTerm = sin(piOverA * traceProduct(-I* mComplex, x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                short sign = tup.second;
                testPointTerm += sign * sin(piOverA * traceProduct(-I* mComplex, etaX));
            }
            xTerms.push_back(testPointTerm);
        }
    }

    vector<double> terms;
    terms.resize(size, 0);

    //Do it this way so the multiplications are vectorizable
    for (size_t i = 0; i < size; ++i) {
        terms[i] = yStars[i] * bessels[i] * xStarTerms[i] * xTerms[i];
    }

    double answer = Aux.kahanSummation(terms);

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        answer *= -4.0 / pointCount;
    } else {
        answer *= 4.0 / pointCount;
    }

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm = Y * K.exactKBessel(twoPiOverA * mAbs * Y);

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
    return (x0 * y1 - y0 * x1)/(y1 - y0);
}

double BianchiMaassSearch::heckeCheck(map<Index, double>& coeffMap) {
    double coeff0 = coeffMap.at(dToPrimes.at(d)[0]);
    double coeff1 = coeffMap.at(dToPrimes.at(d)[1]);
    double coeff2 = coeffMap.at(dToPrimes.at(d)[2]);

    Index prod01 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[1], d);
    Index prod02 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[2], d);

    double hecke1 = coeff0 * coeff1 - coeffMap.at(prod01);
    double hecke2 = coeff0 * coeff2 - coeffMap.at(prod02);

    return abs(hecke1) + abs(hecke2);
}

double BianchiMaassSearch::approxCondition(KBessel *K, const vector<Index> &indexTransversal, double Y) {

    if (K->getLambda() <= 1) {
        double min = +INFINITY;
        double max = -INFINITY;


        for(const auto index : indexTransversal) {
            double arg = 2 * pi / A * index.getAbs(d) * Y;
            double bess = Y * K->exactKBessel(arg);
            if (abs(bess) < min) {
                min = abs(bess);
            }
            if (abs(bess) > max) {
                max = abs(bess);
            }
        }

        return max/min;
    } else {
        double min = +INFINITY;

        for(const auto index : indexTransversal) {
            double arg = 2 * pi / A * index.getAbs(d) * Y;
            if (arg > K->getLambda()) {
                break;
            }
            double bess = Y * K->exactKBessel(arg);
            if (abs(bess) < min) {
                min = abs(bess);
            }
        }
        double arg = 2 * pi / A * indexTransversal[indexTransversal.size() - 1].getAbs(d) * Y;
        double bess = Y * K->exactKBessel(arg);
        if (abs(bess) < min) {
            min = abs(bess);
        }

        return 1.0/min;
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
        std::filesystem::create_directories(directory);
    }
}

bool BianchiMaassSearch::isFileEmpty(const string &filename) {
    return std::filesystem::is_empty(filename);
}


tuple<vector<pair<double,double>>, double, double> BianchiMaassSearch::fineSecantMethod(double leftLambda, double rightLambda) {
    std::cout << "Starting final refinement." << std::endl;
    std::cout << std::setprecision(16) << "[" << leftLambda << ", " << rightLambda << "]" << std::endl;

    double M0 = computeM0General(rightLambda);

    KBessel* leftLambdaK = new KBessel(2 * pi / A /*usual factor*/
                        * 1 /*smallest magnitude of an index*/
                        * Y0 /*smallest height of a pullback*/
            , leftLambda);

    KBessel* rightLambdaK = new KBessel(twoPiOverA * 1 * Y0, rightLambda);

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

    double Y = computeWellConditionedY(leftLambdaK, rightLambdaK, leftLambda, rightLambda, M0, indexTransversal);

    double MY = computeMYGeneral(M0, Y);

    double maxYStar = 0;

    map<Index, vector<TestPointOrbitData>> mToYTestPointOrbits;

#pragma omp parallel for schedule(dynamic) default(none) shared(indexTransversal, mToYTestPointOrbits, Y, MY)
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index m = indexTransversal[i];
        auto orbitsY = getPointPullbackOrbits(m, Y, MY);
#pragma omp critical
        {
            mToYTestPointOrbits[m] = orbitsY;
        }
    }

    for (const auto& m : indexTransversal) {
        for (const auto& itr : mToYTestPointOrbits[m]) {
            if (maxYStar < itr.representativePullback.getJ()) {
                maxYStar = itr.representativePullback.getJ();
            }
        }
    }

    leftLambdaK->setLambdaAndPrecompute(leftLambda, 2 * pi / A * M0 * maxYStar);
    MatrixXd matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, (*leftLambdaK));
    auto solAndCondYLambda1 = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

    rightLambdaK->setLambdaAndPrecompute(rightLambda, 2 * pi / A * M0 * maxYStar);
    matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, (*rightLambdaK));
    auto solAndCondYLambda2 = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

    std::cout << solAndCondYLambda1.second << " " << solAndCondYLambda2.second << std::endl;


    map<Index, double> coeffMapLambda1, coeffMapLambda2;

    for (int i = 0; i < indexTransversal.size(); i++) {
        Index index = indexTransversal[i];
        coeffMapLambda1[index] = solAndCondYLambda1.first(i);
        coeffMapLambda2[index] = solAndCondYLambda2.first(i);
    }

    //We have the solutions from everything, now we cook up the phi function

    double coeff0Lambda1 = coeffMapLambda1.at(dToPrimes.at(d)[0]);
    double coeff1Lambda1 = coeffMapLambda1.at(dToPrimes.at(d)[1]);
    double coeff2Lambda1 = coeffMapLambda1.at(dToPrimes.at(d)[2]);

    double coeff0Lambda2 = coeffMapLambda2.at(dToPrimes.at(d)[0]);
    double coeff1Lambda2 = coeffMapLambda2.at(dToPrimes.at(d)[1]);
    double coeff2Lambda2 = coeffMapLambda2.at(dToPrimes.at(d)[2]);

    Index prod01 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[1], d);
    Index prod02 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[2], d);

    double hecke1Lambda1 = coeff0Lambda1 * coeff1Lambda1 - coeffMapLambda1.at(prod01);
    double hecke2Lambda1 = coeff0Lambda1 * coeff2Lambda1 - coeffMapLambda1.at(prod02);
    hecke1Lambda1 = abs(hecke1Lambda1);
    hecke2Lambda1 = abs(hecke2Lambda1);

    double hecke1Lambda2 = coeff0Lambda2 * coeff1Lambda2 - coeffMapLambda2.at(prod01);
    double hecke2Lambda2 = coeff0Lambda2 * coeff2Lambda2 - coeffMapLambda2.at(prod02);
    hecke1Lambda2 = abs(hecke1Lambda2);
    hecke2Lambda2 = abs(hecke2Lambda2);


    double phi1;
    double phi2;
    double eps1;
    double eps2;
    bool signChange = false;
    int n = 1;
    while (!signChange) {
        vector<double> epsilon;
        for (int i = -n; i <= n; i++) {
            if (i == 0) {
                continue;
            }
            epsilon.push_back(i);
        }

        for (auto epsItr1 : epsilon) {
            for (auto epsItr2 : epsilon) {
                phi1 = epsItr1 * hecke1Lambda1 + epsItr2 * hecke2Lambda1;
                phi2 = epsItr1 * hecke1Lambda2 + epsItr2 * hecke2Lambda2;
                if (areDifferentSign(phi1, phi2)) {
                    signChange = true;
                    eps1 = epsItr1;
                    eps2 = epsItr2;
                    break;
                }
            }
            if (signChange) {
                break;
            }
        }
        n++;
    }

    double lambdaN = leftLambda;
    double phiN = phi1;
    double lambdaNPlus1 = rightLambda;
    double phiNPlus1 = phi2;

    std::cout << std::setprecision(16) << "r = " << leftLambda << ", " << heckeCheck(coeffMapLambda1) << std::endl;
    std::cout << std::setprecision(16) << "r = " << rightLambda << ", " << heckeCheck(coeffMapLambda2) << std::endl;

    vector<pair<double, double>> heckeValues;
    heckeValues.emplace_back(lambdaN, heckeCheck(coeffMapLambda1));
    heckeValues.emplace_back(lambdaNPlus1, heckeCheck(coeffMapLambda2));

    int minIterations = 4;
    int maxIterations = 30;
    int iterations = 0;

    while (true) {
        iterations++;

        double nextLambda = findZeroOfLinearInterpolation(lambdaN, phiN, lambdaNPlus1, phiNPlus1);

        KBessel K = KBessel(twoPiOverA * 1 * Y0, nextLambda);
        K.extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, K);
        auto solAndCondYNextLambda = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

        std::cout << solAndCondYNextLambda.second << " " << solAndCondYNextLambda.second << std::endl;

        map<Index, double> coeffMapNextLambda;
        for (int i = 0; i < indexTransversal.size(); i++) {
            Index index = indexTransversal[i];
            coeffMapNextLambda[index] = solAndCondYNextLambda.first(i);
        }

        double hecke = heckeCheck(coeffMapNextLambda);
        heckeValues.emplace_back(nextLambda, hecke);
        std::cout << std::setprecision(16) << "r = " << nextLambda << ", " << hecke << std::endl;

        //We have the solutions from everything, now we cook up the phi function

        double coeff0NextLambda = coeffMapNextLambda.at(dToPrimes.at(d)[0]);
        double coeff1NextLambda = coeffMapNextLambda.at(dToPrimes.at(d)[1]);
        double coeff2NextLambda = coeffMapNextLambda.at(dToPrimes.at(d)[2]);

        double hecke1NextLambda = coeff0NextLambda * coeff1NextLambda - coeffMapNextLambda.at(prod01);
        double hecke2NextLambda = coeff0NextLambda * coeff2NextLambda - coeffMapNextLambda.at(prod02);
        hecke1NextLambda = abs(hecke1NextLambda);
        hecke2NextLambda = abs(hecke2NextLambda);

        double nextPhi = eps1 * hecke1NextLambda + eps2 * hecke2NextLambda;

        /*
         * This is the Anderson-Bjorck bracketed root finding algorithm
         */
        if (areDifferentSign(phiNPlus1, nextPhi)) {
            lambdaN = lambdaNPlus1;
            phiN = phiNPlus1;

            lambdaNPlus1 = nextLambda;
            phiNPlus1 = nextPhi;
        } else {
            double m = 1 - nextPhi/phiNPlus1;
            if (m <= 0) {
                m = 0.5;
            }
            phiN *= m;

            lambdaNPlus1 = nextLambda;
            phiNPlus1 = nextPhi;
        }

        if (iterations < minIterations) {
            continue;
        }
        if (hecke > 1 //not modular
            || heckeHasConverged(heckeValues) //computation succeeded and is complete
            || iterations > maxIterations //computation not converging as expected
            || phiN == phiNPlus1 //numerical error incoming, could just be very good convergence
            || lambdaN == lambdaNPlus1) {

            delete leftLambdaK;
            delete rightLambdaK;
            return {heckeValues, lambdaN, lambdaNPlus1};

        }
    }
}

bool BianchiMaassSearch::heckeHasConverged(const vector<pair<double, double>> &heckeValues) {
    int size = heckeValues.size();
    if (size < 4) {
        return false;
    }

    vector<double> values;
    for (int i = size - 1; i >= size - 4; i--) {
        values.push_back(heckeValues[i].second);
    }

    double max = *std::max_element(values.begin(), values.end());
    double min = *std::min_element(values.begin(), values.end());

    double ratio = max/min;

    if (ratio < 10) {
        return true;
    } else {
        return false;
    }
}

double BianchiMaassSearch::computeWellConditionedY(KBessel *K, double M0, vector<Index> &indexTransversal) {
    double c = 0;
    if (K->getLambda() <= 1) {
        c = 1 / (2*pi/A*M0);
    } else {
        double r = sqrt(K->getLambda() - 1);
        c = r / (2*pi/A*M0);
    }
    double upperY = Y0;
    double lowerY = 0.5 * c;

    vector<double> heightsToTry;

    double tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.99;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    vector<double> approxConditionNumbers(heightsToTry.size(), 0);


#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, approxConditionNumbers, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        approxConditionNumbers[i] = approxCondition(K, indexTransversal, Y);
    }


    for (int i = heightsToTry.size() - 1; i >= 0; i--)  {
        if (approxConditionNumbers[i] <= 100) {
            return heightsToTry[i];
        }
    }

    double Y = 0;
    double minCondSoFar = +INFINITY;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        if (minCondSoFar > approxConditionNumbers[i]) {
            minCondSoFar = approxConditionNumbers[i];
            Y = heightsToTry[i];
        }
    }

    /*Do it all over again but finer
     *and in a narrower range
     */

    upperY = Y + 0.5 * c;
    lowerY = Y - 0.5 * c;
    if (lowerY <= 0) {
        lowerY = Y * 0.5;
    }

    heightsToTry.clear();

    tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.9995;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    approxConditionNumbers.clear();

    approxConditionNumbers.resize(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, approxConditionNumbers, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        approxConditionNumbers[i] = approxCondition(K, indexTransversal, testY);
    }

    for (int i = heightsToTry.size() - 1; i >= 0; i--)  {
        if (approxConditionNumbers[i] <= 100) {
            return heightsToTry[i];
        }
    }

    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        if (minCondSoFar > approxConditionNumbers[i]) {
            minCondSoFar = approxConditionNumbers[i];
            Y = heightsToTry[i];
        }
    }
    return Y;
}

double
BianchiMaassSearch::computeWellConditionedY(KBessel *leftLambdaK, KBessel *rightLambdaK, double leftLambda, double rightLambda, double M0,
                                            const vector<Index> &indexTransversal) {
    //double upperY = 2.0 * max(1.0, rightR)/(2*pi/A*M0);
    //double lowerY = 1.7 * max(1.0, rightR)/(2*pi/A*M0);

    double c = 0;
    if (rightLambdaK->getLambda() <= 1) {
        c = 1 / (2*pi/A*M0);
    } else {
        double r = sqrt(rightLambdaK->getLambda() - 1);
        c = r / (2*pi/A*M0);
    }
    double upperY = Y0;
    double lowerY = 0.5 * c;

    vector<double> heightsToTry;

    double tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.99;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    vector<double> lambda1ApproxCond(heightsToTry.size(), 0);
    vector<double> lambda2ApproxCond(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda2ApproxCond, rightLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda2ApproxCond[i] = approxCondition(rightLambdaK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda1ApproxCond, leftLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda1ApproxCond[i] = approxCondition(leftLambdaK, indexTransversal, Y);
    }

    for (int i = heightsToTry.size() - 1; i >= 0; i--) {
        double condHere = std::max({lambda1ApproxCond[i],
                                    lambda2ApproxCond[i]});
        if (condHere <= 100) {
            return heightsToTry[i];
        }
    }

    double Y = 0;
    double minCondSoFar = +INFINITY;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double condHere = std::max({lambda1ApproxCond[i],
                                      lambda2ApproxCond[i]});
        if (condHere < minCondSoFar) {
            minCondSoFar = condHere;
            Y = heightsToTry[i];
        }
    }

    /*Do it all over again but finer
     *and in a narrower range
     */

    upperY = Y + 0.5 * c;
    lowerY = Y - 0.5 * c;
    if (lowerY <= 0) {
        lowerY = 0.5 * Y;
    }

    heightsToTry.clear();

    tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.9995;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    lambda1ApproxCond.clear();
    lambda2ApproxCond.clear();

    lambda1ApproxCond.resize(heightsToTry.size(), 0);
    lambda2ApproxCond.resize(heightsToTry.size(), 0);


#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda2ApproxCond, rightLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        lambda2ApproxCond[i] = approxCondition(rightLambdaK, indexTransversal, testY);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda1ApproxCond, leftLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        lambda1ApproxCond[i] = approxCondition(leftLambdaK, indexTransversal, testY);
    }

    for (int i = heightsToTry.size() - 1; i >= 0; i--) {
        double condHere = std::max({lambda1ApproxCond[i],
                                    lambda2ApproxCond[i]});
        if (condHere <= 100) {
            return heightsToTry[i];
        }
    }

    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double condHere = std::max({lambda1ApproxCond[i],
                                      lambda2ApproxCond[i]});
        if (condHere < minCondSoFar) {
            minCondSoFar = condHere;
            Y = heightsToTry[i];
        }
    }
    return Y;
}

pair<double, double>
BianchiMaassSearch::computeTwoWellConditionedY(KBessel *leftLambdaK, KBessel *rightLambdaK, double leftLambda, double rightLambda,
                                               double M0, const vector<Index> &indexTransversal) {

    double c = 0;
    if (rightLambdaK->getLambda() <= 1) {
        c = 1 / (2*pi/A*M0);
    } else {
        double r = sqrt(rightLambdaK->getLambda() - 1);
        c = r / (2*pi/A*M0);
    }
    double upperY = Y0;
    double lowerY = 0.5 * c;

    vector<double> heightsToTry;

    double tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.99;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    vector<double> lambda1ApproxCond(heightsToTry.size(), 0);
    vector<double> lambda2ApproxCond(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda2ApproxCond, rightLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda2ApproxCond[i] = approxCondition(rightLambdaK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda1ApproxCond, leftLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda1ApproxCond[i] = approxCondition(leftLambdaK, indexTransversal, Y);
    }


    double ansY1 = 0;
    double ansY2 = 0;
    double minCondSoFar = +INFINITY;
    //heightsToTry is sorted smallest to largest
    for (int i = heightsToTry.size() - 1; i >= 2; i--) {//end i at second to last because Y2 has to be greater
        for (int j = i - 2; j >= 0; j--) {//j is always less than i, so j indexes a larger Y
            double condHere = std::max({lambda1ApproxCond[i],
                                          lambda1ApproxCond[j],
                                          lambda2ApproxCond[i],
                                          lambda2ApproxCond[j]});

            if (condHere <= 100) {
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
                return {min(ansY1, ansY2), max(ansY1, ansY2)};
            }
            if (condHere < minCondSoFar) {
                minCondSoFar = condHere;
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
            }
        }
    }

    /*Do it all over again but finer
     *and in a narrower range
     */

    upperY = ansY2 + 0.5 * c;
    lowerY = ansY1 - 0.5 * c;
    if (lowerY <= 0) {
        lowerY = 0.5 * ansY1;
    }

    heightsToTry.clear();

    tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.9995;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    lambda1ApproxCond.clear();
    lambda2ApproxCond.clear();

    lambda1ApproxCond.resize(heightsToTry.size(), 0);
    lambda2ApproxCond.resize(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda2ApproxCond, rightLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda2ApproxCond[i] = approxCondition(rightLambdaK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, lambda1ApproxCond, leftLambdaK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        lambda1ApproxCond[i] = approxCondition(leftLambdaK, indexTransversal, Y);
    }

    //heightsToTry is sorted smallest to largest
    for (int i = heightsToTry.size() - 1; i >= 2; i--) {//end i at second to last because Y2 has to be greater
        for (int j = i - 2; j >= 0; j--) {//j is always less than i, so j indexes a larger Y
            double condHere = std::max({lambda1ApproxCond[i],
                                        lambda1ApproxCond[j],
                                        lambda2ApproxCond[i],
                                        lambda2ApproxCond[j]});

            if (condHere <= 100) {
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
                return {min(ansY1, ansY2), max(ansY1, ansY2)};
            }
            if (condHere < minCondSoFar) {
                minCondSoFar = condHere;
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
            }
        }
    }
    return {min(ansY1, ansY2), max(ansY1, ansY2)};
}

void BianchiMaassSearch::setUpOutputLogFiles() {
    if (mode == "coarse") {
        const std::string directory = "Output/Coarse/";
        const std::string prefix = "coarse_"
                                   + to_string(d) + "_"
                                   + symClass;
        createOutputDirectory(directory);

        std::string outputFilename = directory
                                     + prefix
                                     + ".txt";

        coarseOutputFile.open(outputFilename, std::ofstream::out | std::ofstream::app);

        if (coarseOutputFile.is_open()) {
            if (isFileEmpty(outputFilename)) {
                coarseOutputFile << "d = " << d << '\n';
                coarseOutputFile << "symClass = " << symClass << '\n';
                coarseOutputFile << "D = " << D << std::endl;
            }
        } else {
            std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
        }
    } else if (mode == "medium") {

        const std::string directoryOut = "Output/Medium/";
        const std::string prefixOut = "medium_"
                                      + to_string(d) + "_"
                                      + symClass;
        createOutputDirectory(directoryOut);

        std::string outputFilename = directoryOut
                                     + prefixOut
                                     + ".txt";

        mediumOutputFile.open(outputFilename, std::ofstream::out | std::ofstream::app);

        if (mediumOutputFile.is_open()) {
            if (isFileEmpty(outputFilename)) {
                mediumOutputFile << "d = " << d << '\n';
                mediumOutputFile << "symClass = " << symClass << '\n';
                mediumOutputFile << "D = " << D << std::endl;
            }
        } else {
            std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
        }

    } else if (mode == "fine") {
        const std::string directoryOut = "Output/Fine/";
        const std::string prefixOut = "fine_"
                                      + to_string(d) + "_"
                                      + symClass;
        createOutputDirectory(directoryOut);

        std::string outputFilename = directoryOut
                                     + prefixOut
                                     + ".txt";

        fineOutputFile.open(outputFilename, std::ofstream::out | std::ofstream::app);

        if (fineOutputFile.is_open()) {
            if (isFileEmpty(outputFilename)) {
                fineOutputFile << "d = " << d << '\n';
                fineOutputFile << "symClass = " << symClass << std::endl;
                fineOutputFile << "D =" << D << std::endl;
            }
        } else {
            std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
        }
    }

}

vector<pair<double, double>> BianchiMaassSearch::getIntervalsForCoarseSearch(double startLambda, double endLambda) {

    const std::string directoryIn = "Output/Coarse/";
    const std::string prefixIn = "coarse_"
                                  + to_string(d) + "_"
                                  + symClass;

    std::string inputFilename = directoryIn
                                 + prefixIn
                                 + ".txt";

    //scan through and find the last value we're complete to
    std::ifstream coarseInFile;
    coarseInFile.open(inputFilename);
    string line;
    while (getline(coarseInFile, line)) {
        if (line.substr(0, 15) == "Complete up to ") {
            coarseComplete = std::stod(line.substr(15));
        }
    }

    //Magic number
    //Supposed to reflect that we are dispatching the searcher on intervals that we expect
    //to have an eigenvalue in them 2/10th of the time. Allows for easier conditioning and
    //guards against unusually close eigenvalues, but not conclusively.
    double numEigenValues = 0.2;

    vector<double> endpoints;
    double endpoint = max(startLambda, coarseComplete);
    while (endpoint < endLambda) {
        endpoints.push_back(endpoint);
        endpoint = Od.eigenvalueIntervalRightEndpoint(endpoint, numEigenValues);
    }
    endpoints.push_back(endLambda);

    if (endpoints.size() < 2) {
        return {};
    }

    vector<pair<double, double>> intervals;
    for (int i = 0; i <= endpoints.size() - 2; i++) {
        intervals.emplace_back(endpoints[i], endpoints[i + 1]);
    }

    return intervals;
}

vector<pair<double, double>> BianchiMaassSearch::getIntervalsForMediumSearch() {

    /*
     *
     * Open the output file and record how far the medium searching has gone
     *
     */

    const std::string directoryOut = "Output/Medium/";
    const std::string prefixOut = "medium_"
                                  + to_string(d) + "_"
                                  + symClass;
    createOutputDirectory(directoryOut);

    std::string outputFilename = directoryOut
                                 + prefixOut
                                 + ".txt";

    //scan through and find the last value we're complete to
    std::ifstream mediumInFile;
    mediumInFile.open(outputFilename);
    string line;
    while (getline(mediumInFile, line)) {
        if (line.substr(0, 15) == "Complete up to ") {
            mediumComplete = std::stod(line.substr(15));
        }
    }

    /*
     *
     * Open the coarse search file and gather the intervals that need to be medium searched still
     *
     */

    const std::string directoryIn = "Output/Coarse/"; // Change this to the desired directory
    const std::string prefixIn = "coarse_"
                                 + to_string(d) + "_"
                                 + symClass;

    std::string filenameIn = directoryIn
                             + prefixIn
                             + ".txt";

    std::ifstream inputFile;
    inputFile.open(filenameIn);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening file \"" << filenameIn << "\"" << std::endl;
    }

    vector<pair<double, double>> intervalsFromFile;

    while (std::getline(inputFile, line)) {
        if (line[0] == '[') {
            //remove brackets
            string numbers = line.substr(1, line.size() - 2);
            int comma = numbers.find(',');
            string leftStr = numbers.substr(0,comma);
            string rightStr = numbers.substr(comma + 2);
            double left = std::stod(leftStr);
            double right = std::stod(rightStr);
            if (mediumComplete < right) {
                intervalsFromFile.emplace_back(left,right);
            }
        }
    }

    return intervalsFromFile;
}

vector<pair<double,double>> BianchiMaassSearch::getIntervalsForFineSearch() {


    /*
     *
     * Open the output file and record how far the fine searching has gone
     *
     */
    const std::string directoryOut = "Output/Fine/";
    const std::string prefixOut = "fine_"
                                  + to_string(d) + "_"
                                  + symClass;
    createOutputDirectory(directoryOut);

    std::string outputFilename = directoryOut
                                 + prefixOut
                                 + ".txt";

    std::ifstream fineInFile;
    fineInFile.open(outputFilename);
    string line;
    while (getline(fineInFile, line)) {
        if (line.substr(0, 15) == "Complete up to ") {
            fineComplete = std::stod(line.substr(15));
        }
    }

    /*
     *
     * Open the medium search file and gather the intervals that need to be searched still
     *
     */

    const std::string directoryIn = "Output/Medium/";
    const std::string prefixIn = "medium_"
                                 + to_string(d) + "_"
                                 + symClass;

    std::string filenameIn = directoryIn
                             + prefixIn
                             + ".txt";

    std::ifstream inputFile;
    inputFile.open(filenameIn);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening file \"" << filenameIn << "\"" << std::endl;
    }

    vector<pair<double, double>> intervalsFromFile;
    while (std::getline(inputFile, line)) {
        if (line[0] == '[') {
            //remove brackets
            string numbers = line.substr(1, line.size() - 2);
            int comma = numbers.find(',');
            string leftStr = numbers.substr(0,comma);
            string rightStr = numbers.substr(comma + 2);
            double left = std::stod(leftStr);
            double right = std::stod(rightStr);
            if (fineComplete < right) {
                intervalsFromFile.emplace_back(left,right);
            }
        }
    }

    return intervalsFromFile;
}

void BianchiMaassSearch::computeMaximumD(double lambda, int timeLimitSeconds) {
    double maxNanoseconds = timeLimitSeconds * pow(10,9);

    double maxTerms = maxNanoseconds / nanosecondsPerTerm;

    double minD = 3;
    double maxD = 20;
    if (d == 163) {
        maxD = 12;
    }
    double leftD = minD; //I don't want to compute anything with any less precision than this
    double rightD = maxD; //This is the most precision I will ever work to


    unsigned long long int leftTerms = 0;
    unsigned long long int rightTerms = 0;

    truncation = pow(10.0, -leftD);
    double M0 = computeM0General(lambda);

    KBessel* K = new KBessel(twoPiOverA * 1 * Y0, lambda);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);

    double Y = computeWellConditionedY(K, M0, indexTransversal);
    delete K;
    double MY = computeMYGeneral(M0, Y);

    for (auto m : indexTransversal) {
        leftTerms += countPointPullbackOrbits(m, Y, MY);
    }
    leftTerms *= indexTransversal.size();


    /*
     *
     *
     *
     */

    truncation = pow(10.0, -rightD);
    M0 = computeM0General(lambda);

    indicesM0 = Od.indicesUpToM(M0);
    data = Od.indexOrbitQuotientData(indicesM0, symClass);
    indexTransversal = get<0>(data);

    MY = computeMYGeneral(M0, Y);

    for (auto m : indexTransversal) {
        rightTerms += countPointPullbackOrbits(m, Y, MY) * indexTransversal.size();
    }

    //In this case, both precisions exceed or meet the time limit, but we're not willing
    //to work in precision less than D = 2
    if (leftTerms >= maxTerms && rightTerms >= maxTerms) {
        D = minD;
        truncation = pow(10.0, -D);
        return;
    }

    //In this case, both precisions will finish on time, so we return the largest
    //precision that seems sensible. This may never happen.
    if (leftTerms <= maxTerms && rightTerms <= maxTerms) {
        double answer;
        if (mode == "coarse") {
            answer = rightD * 0.5;
        } else if (mode == "medium") {
            answer =  rightD * 0.75;
        } else {
            //mode == "fine"
            answer =  rightD;
        }
        D = max(minD, answer);
        truncation = pow(10.0, -D);
        return;
    }

    while (rightD - leftD > 0.1) {
        double centerD = (leftD + rightD)/2;
        unsigned long long int centerTerms = 0;
        truncation = pow(10.0, -centerD);
        M0 = computeM0General(lambda);

        indicesM0 = Od.indicesUpToM(M0);
        data = Od.indexOrbitQuotientData(indicesM0, symClass);
        indexTransversal = get<0>(data);

        MY = computeMYGeneral(M0, Y);

        for (auto m : indexTransversal) {
            centerTerms += countPointPullbackOrbits(m, Y, MY) * indexTransversal.size();
        }

        if (centerTerms > maxTerms) {
            rightD = centerD;
            rightTerms = centerTerms;
        } else {
            leftD = centerD;
            leftTerms = centerTerms;
        }
    }

    double answer;

    if (mode == "coarse") {
        answer = rightD * 0.5;
    } else if (mode == "medium") {
        answer =  rightD * 0.75;
    } else {
        //mode == "fine"
        answer =  rightD;
    }

    D = max(minD, answer);
    truncation = pow(10.0, -D);
    return;
}

pair<int, int> BianchiMaassSearch::computeQ0Q1(Index m, double MY) {
    if (d == 3 || d == 1) {
        throw std::invalid_argument("not implemented for d = " + to_string(d));
    }

    complex<double> mComplex = m.getComplex(d);
    int Q0 = ceil((MY + abs(mComplex.imag()))/(2 * A));
    int Q1 = ceil((MY + abs(mComplex.real()))/2);
    bool adjusted = true;

    while (adjusted) {

        adjusted = false;
        bool circleContainsNoPoints = true;
        while (circleContainsNoPoints) {
            int tempQ0 = Q0 - 1;

            vector<complex<double>> pointsToMiss;
            pointsToMiss.emplace_back(2.0 * Q1, 0);
            pointsToMiss.emplace_back(-2.0 * Q1, 0);
            pointsToMiss.push_back(2.0 * tempQ0 * conj(theta));
            pointsToMiss.push_back(-2.0 * tempQ0 * conj(theta));
            pointsToMiss.push_back(pointsToMiss[3] + pointsToMiss[0]);
            pointsToMiss.push_back(pointsToMiss[2] + pointsToMiss[1]);

            for (const auto& p : pointsToMiss) {
                bool circleContainsP = abs((p-mComplex)) <= MY;
                if (circleContainsP) {
                    circleContainsNoPoints = false;
                    break;
                }
            }

            if (circleContainsNoPoints) {
                Q0 = tempQ0;
                adjusted = true;
            }

        }

        circleContainsNoPoints = true;
        while (circleContainsNoPoints) {
            int tempQ1 = Q1 - 1;

            vector<complex<double>> pointsToMiss;
            pointsToMiss.emplace_back(2.0 * tempQ1, 0);
            pointsToMiss.emplace_back(-2.0 * tempQ1, 0);
            pointsToMiss.push_back(2.0 * Q0 * conj(theta));
            pointsToMiss.push_back(-2.0 * Q0 * conj(theta));
            pointsToMiss.push_back(pointsToMiss[3] + pointsToMiss[0]);
            pointsToMiss.push_back(pointsToMiss[2] + pointsToMiss[1]);

            for (const auto& p : pointsToMiss) {
                bool circleContainsP = abs((p-mComplex)) <= MY;
                if (circleContainsP) {
                    circleContainsNoPoints = false;
                    break;
                }
            }

            if (circleContainsNoPoints) {
                Q1 = tempQ1;
                adjusted = true;
            }
        }
    }

    if (Auxiliary::mod(-d,4) != 1) {
        return {Q0, Q1};
    }

    if (Q0 % Q1 == 0 && Q0/Q1 % 2 == 0) {
        return {Q0, Q1};
    }

    int k = 0;
    double quotient = ((double)Q0)/Q1;

    while (!(2 * k < quotient && quotient < 2 * (k+1))) {
        k++;
    }

    int testQ0;
    int testQ1;

    if (k == 0) {
        testQ1 = Q1;
        testQ0 = 2 * Q1;
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
            testQ0 = firstOptionQ0;
            testQ1 = firstOptionQ1;
        } else {
            testQ0 = secondOptionQ0;
            testQ1 = secondOptionQ1;
        }
    }

    long pointsIfAdjust = (2 * testQ0) * (2 * testQ1) / 4;
    long pointsIfNoAdjust = (2 * Q0) * (2 * Q1) / 2;

    if (pointsIfAdjust < pointsIfNoAdjust) {
        return {testQ0, testQ1};
    } else {
        return {Q0, Q1};
    }
}

