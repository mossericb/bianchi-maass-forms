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
#include <boost/align/aligned_allocator.hpp>
#include <sstream>
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


void BianchiMaassSearch::coarseSearchForEigenvalues(const double leftR, const double rightR) {
    if (leftR >= rightR) {
        throw(std::invalid_argument("leftR should be less than rightR"));
    }

    auto intervalsFromFile = getIntervalsForCoarseSearch(leftR, rightR);
    if (intervalsFromFile.empty()) {
        std::cout << "Range already complete." << std::endl;
        return;
    }
    stack<pair<double,double>> intervals;

    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    double lastPrecisionRecompute = intervals.top().first;
    if (autoPrecision) {
        computeMaximumD(lastPrecisionRecompute, 180);
        std::cout << "D = " << D << std::endl;
    }

    KBessel* leftRK = new KBessel(twoPiOverA * 1 * Y0, intervals.top().first);
    KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, intervals.top().second);
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

        rightRK = new KBessel(twoPiOverA * 1 * Y0, right);
        if (possiblyContainsEigenvalue(left, right, leftRK, rightRK)) {
            coarseOutputFile << std::setprecision(16) << "[" << left << ", " << right << "]" << std::endl;
        }
        delete leftRK;
        leftRK = rightRK;

        coarseOutputFile << std::setprecision(16) << "Complete up to " << right << std::endl;
    }
    delete rightRK;
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

    if (intervalsFromFile.empty()) {
        std::cout << "Range already complete." << std::endl;
        return;
    }

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
            //weylRightEndpoint = max(weylRightEndpoint, interval.first + 0.00000001);
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
                    spawnedIntervals.emplace(center, right);
                    spawnedIntervals.emplace(left, center);
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
    if (intervalsFromFile.empty()) {
        std::cout << "Range already complete." << std::endl;
        return;
    }


    stack<pair<double, double>> intervals;
    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    while (!intervals.empty()) {
        auto interval = intervals.top();
        intervals.pop();
        if (autoPrecision) {
            computeMaximumD(interval.second, 180);
            std::cout << "D = " << D << std::endl;
        }
        medianIllinoisSearch(interval.first, interval.second);
        fineOutputFile << "Complete up to " << std::setprecision(16) << interval.second << std::endl;
    }
}

void BianchiMaassSearch::medianIllinoisSearch(double leftR, double rightR) {

    std::cout << std::setprecision(16)
              << "[" << leftR << ", " << rightR << "]" << std::endl;

    double M0 = computeM0General(rightR);

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

    KBessel* leftRK = new KBessel(twoPiOverA * 1 * Y0, leftR);
    KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, rightR);

    auto pairOfYs = computeTwoWellConditionedY(leftRK, rightRK, leftR, rightR, M0, indexTransversal);

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

    leftRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    MatrixXd matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *leftRK);
    MatrixXd matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *leftRK);

    auto solAndCondY1R1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2R1 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    rightRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *rightRK);
    matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *rightRK);

    delete leftRK;
    delete rightRK;

    auto solAndCondY1R2 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2R2 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    std::cout << solAndCondY1R1.second << " " << solAndCondY2R1.second << " ";
    std::cout << solAndCondY1R2.second << " " << solAndCondY2R2.second << std::endl;

    vector<double> coeffDiffa(indexTransversal.size(), NAN);
    vector<double> coeffDiffb(indexTransversal.size(), NAN);

    for (int i = 0; i < indexTransversal.size(); i++) {
        coeffDiffa[i] = solAndCondY1R1.first(i) - solAndCondY2R1.first(i);
        coeffDiffb[i] = solAndCondY1R2.first(i) - solAndCondY2R2.first(i);
    }

    int signChangeIndex = 0;
    double maxDiff = 0;
    for (int i = 0; i < indexTransversal.size(); i++) {
        if (areDifferentSign(coeffDiffa[i], coeffDiffb[i])) {
            double diffHere = std::min(abs(coeffDiffa[i]), abs(coeffDiffb[i]));
            if (diffHere > maxDiff) {
                signChangeIndex = i;
                maxDiff = diffHere;
            }
        }
    }

    double a = leftR;
    double b = rightR;
    double c = 0;

    pair<Eigen::Matrix<double, -1, -1, 0>, double> solAndCondY1NextR;

    std::cout << std::setprecision(16) << a << "\n";
    std::cout << std::setprecision(16) << b << std::endl;

    int side = 0;

    double cutoff = pow(10.0, -15);
    while (abs((1 - b/a)) > cutoff) {
        vector<double> c_vec;
        for (int i = 0; i < indexTransversal.size(); i++) {
            if (areDifferentSign(coeffDiffa[i], coeffDiffb[i])) {
                double c_here = (coeffDiffa[i] * b - coeffDiffb[i] * a) / (coeffDiffa[i] - coeffDiffb[i]);
                c_vec.push_back(c_here);
            }
        }
        std::sort(c_vec.begin(), c_vec.end());
        c = c_vec.size() % 2 == 1 ? c_vec[c_vec.size()/2] : (c_vec[c_vec.size()/2 - 1] + c_vec[c_vec.size()/2])/2.0;

        KBessel* nextRK = new KBessel(twoPiOverA * 1 * Y0, c);
        nextRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *nextRK);
        matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *nextRK);

        delete nextRK;

        solAndCondY1NextR = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
        auto solAndCondY2NextR = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

        vector<double> coeffDiffc(indexTransversal.size(), NAN);

        for (int i = 0; i < indexTransversal.size(); i++) {
            coeffDiffc[i] = solAndCondY1NextR.first(i) - solAndCondY2NextR.first(i);
        }

        std::cout << std::setprecision(16) << c << std::endl;

        int leftSignChanges = countSignChanges(coeffDiffa, coeffDiffc);
        int rightSignChanges = countSignChanges(coeffDiffc, coeffDiffb);

        std::cout << leftSignChanges << "/" << indexTransversal.size() << ", ";
        std::cout << rightSignChanges << "/" << indexTransversal.size() << std::endl;

        double leftRatio = leftSignChanges*1.0/indexTransversal.size();
        double rightRatio = rightSignChanges*1.0/indexTransversal.size();

        if (leftRatio > 0.5 && rightRatio < 0.3) {
            b = c;
            coeffDiffb = coeffDiffc;
            if (side == -1) {
                for (auto & itr : coeffDiffa) {
                    itr *= 0.5;
                }
            }
            side = -1;
        } else if (rightRatio > 0.5 && leftRatio < 0.3) {
            a = c;
            coeffDiffa = coeffDiffc;
            if (side == +1) {
                for (auto & itr : coeffDiffb) {
                    itr *= 0.5;
                }
            }
            side = +1;
        } else if (rightRatio > 0.5 && leftRatio > 0.5) {
            std::cout << "!!!!!!!!!!!!!!Interval might contain two spec params!!!!!!!!!!!!!\n";
            std::cout << setprecision(16) << "[" << a << ", " << c << "] and [" << c << ", " << b << "]" << std::endl;
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            break;
        } else {
            std::cout << "Sign changes no longer useful, stopping." << std::endl;
            fineOutputFile << std::setprecision(16) << "[" << a << ", " << b << "], " << c << std::endl;

            vector<double> coeffs;
            for (int i = 0; i < indexTransversal.size(); i++) {
                double coeff = solAndCondY1NextR.first(i);
                coeffs.push_back(coeff);
            }

            saveSearchResult({a, b, c}, indexTransversal, coeffs);
            break;
        }
    }
    if (abs((1 - b/a)) <= cutoff) {
        std::cout << "Reached final precision, stopping." << std::endl;
        fineOutputFile << std::setprecision(16) << "[" << a << ", " << b << "], " << c << std::endl;

        vector<double> coeffs;
        for (int i = 0; i < indexTransversal.size(); i++) {
            double coeff = solAndCondY1NextR.first(i);
            coeffs.push_back(coeff);
        }

        saveSearchResult({a, b, c}, indexTransversal, coeffs);
    }
}


bool BianchiMaassSearch::possiblyContainsEigenvalue(const double leftR, const double rightR, KBessel *leftRK, KBessel *rightRK) {
    //Assume that I will have to recompute everything multiple times here, so don't bother checking otherwise

    //First time through, compute the condition number everywhere to get an idea of how low you can go
    //In subsequent calls, start with Y1 and Y2 at the top end of the range and adjust a few times until
    //the condition number is close to the computed ballpark

    std::cout << std::setprecision(16)
              << "[" << leftR << ", " << rightR << "]" << std::endl;

    double M0 = computeM0General(rightR);

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

    auto pairOfYs = computeTwoWellConditionedY(leftRK, rightRK, leftR, rightR, M0, indexTransversal);

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

    leftRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    MatrixXd matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *leftRK);
    MatrixXd matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *leftRK);

    auto solAndCondY1R1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2R1 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    rightRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y1, *rightRK);
    matrixY2 = produceMatrix(indexTransversal, mToY2TestPointOrbits, indexOrbitDataModSign, Y2, *rightRK);

    auto solAndCondY1R2 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);
    auto solAndCondY2R2 = solveMatrix(matrixY2, indexTransversal, indexOfNormalization);

    std::cout << solAndCondY1R1.second << " " << solAndCondY2R1.second << " ";
    std::cout << solAndCondY1R2.second << " " << solAndCondY2R2.second << std::endl;

    /*
     * When you have to start over, for whatever reason, this means that condition numbers
     * blow up in that interval. The method resets the condition number search goal, but I
     * don't want this goal to stay artificially high. So whenever the algorithm has to start over,
     * it sets this flag to signal that next time it should re-restart over to make the condition
     * search goal lower again.
     */

    vector<double> coeffDiffR1(indexTransversal.size(),NAN);
    vector<double> coeffDiffR2(indexTransversal.size(),NAN);

    for (int i = 0; i < indexTransversal.size(); i++) {
        coeffDiffR1[i] = solAndCondY1R1.first(i) - solAndCondY2R1.first(i);
        coeffDiffR2[i] = solAndCondY1R2.first(i) - solAndCondY2R2.first(i);
    }

    int signChanges = countSignChanges(coeffDiffR1, coeffDiffR2);
    std::cout << signChanges << "/" << coeffDiffR1.size() << std::endl;

    int possibleSignChanges = coeffDiffR1.size() - 1;
    double proportion = ((double)signChanges)/possibleSignChanges;
    if (proportion >=  0.3) {
        return true;
    } else {
        return false;
    }
}

double BianchiMaassSearch::computeM0General(const double r) {
    KBessel bess(1, r);

    //Method: use binary search to find M0 such that
    //K(2*pi/A * M0 * Y0) = 10^-D * K(max(r,1))
    //with M0 > max(r,1)/(2*pi/A*Y0)
    //Reference: JAY JORGENSON, LEJLA SMAJLOVI  ́ C, AND HOLGER THEN
    //ON THE DISTRIBUTION OF EIGENVALUES OF MAASS FORMS ON CERTAIN MOONSHINE GROUPS

    //Step 1: Make sure the equality point happens between my left and right endpoints
    double minM0 = max(r, 1.0)/(2*pi/A*Y0);
    double maxM0 = 100;

    double peak = bess.exactKBessel(max(r,1.0));
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
    M0 += 2;
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

    //These formulas provide bounds to guarantee that the DFT result is valid
    double absRealM = abs(m.getComplex(d).real());
    double absImagM = abs(m.getComplex(d).imag());
    double Q0double = (MY + absImagM)/A;
    double Q1double = MY + absRealM;

    int Q0 = floor(Q0double) + 1;
    int Q1 = floor(Q1double) + 1;
    Q0 += 2;
    Q1 += 2;

    if (Auxiliary::mod(-d, 4) == 2 || Auxiliary::mod(-d, 4) == 3) {
        if (Q0 % 2 == 1) {
            Q0++;
        }
        if (Q1 % 2 == 1) {
            Q1++;
        }
        return Q0 * Q1 / 4;
    }

    //testcode
    if (Q0 % 2 == 1) {
        Q0++;
    }
    if (Q1 % 2 == 1) {
        Q1++;
    }

    return Q0 * Q1 / 2;
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
    Q0 += 2;
    Q1 += 2;

    short nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    short reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

    if (Auxiliary::mod(-d, 4) == 2 || Auxiliary::mod(-d, 4) == 3) {
        if (Q0 % 2 == 1) {
            Q0++;
        }
        if (Q1 % 2 == 1) {
            Q1++;
        }

        double xi0 = 0.5;
        double xi1 = 0.5;
        int L0 = 1 - Q0/2;
        int U0 = Q0/2;
        int L1 = 1 - Q1/2;
        int U1 = Q1/2;

        for (int l0 = 1; l0 <= U0; l0++) {
            for (int l1 = 1; l1 <= U1; l1++) {
                TestPointOrbitData orbit;
                complex<double> x = (l0 - xi0)/(1.0*Q0) + theta * (l1 - xi1)/(1.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;
                orbit.properTranslatesModSign = vector<pair<complex<double>, short>>();

                orbit.properTranslatesModSign.emplace_back(-conj(x), reflectionSign);

                answer.push_back(orbit);
            }
        }
        return answer;
    }

    //testcode
    if (Q0 % 2 == 1) {
        Q0++;
    }
    if (Q1 % 2 == 1) {
        Q1++;
    }

    double xi0, xi1;
    int L0, U0, L1, U1;
    if (Q0 % 2 == 0) {
        xi0 = 0.5;
        L0 = 1 - Q0/2;
        U0 = Q0/2;
    } else {
        xi0 = 0;
        L0 = (1 - Q0)/2;
        U0 = (Q0 - 1)/2;
    }
    if (Q1 % 2 == 0) {
        xi1 = 0.5;
        L1 = 1 - Q1/2;
        U1 = Q1/2;
    } else {
        xi1 = 0;
        L1 = (1 - Q1)/2;
        U1 = (Q1 - 1)/2;
    }

    for (int l0 = L0; l0 <= U0; l0++) {
        for (int l1 = 1; l1 <= U1; l1++) {
            TestPointOrbitData orbit;
            complex<double> x = (l0 - xi0)/(1.0*Q0) + theta * (l1 - xi1)/(1.0*Q1);
            Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
            pullback.reduce(d);

            orbit.representativeComplex = x;
            orbit.representativePullback = pullback;
            orbit.properTranslatesModSign = vector<pair<complex<double>, short>>();

            answer.push_back(orbit);
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
        auto mTestPointData = mToTestPointData[m];
        for (int j = 0; j < size; j++) {
            Index n = indexTransversal[j];
            double entry = computeEntry(m, n, K, mTestPointData, Y, ntoIndexOrbitData[n]);

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
        if (orbit.representativeComplex.real() == 0 && orbit.representativeComplex.imag() == 0) {
            pointCount++;
        } else {
            pointCount += (orbit.properTranslatesModSign.size() + 1) * 2;
        }
    }

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

            if (x.real() == 0 && x.imag() == 0) {

                double indexTerm = 0;
                for (const auto& tup : nIndexOrbitDataModSign) {
                    indexTerm += tup.second * 2;
                }
                xStarTerms.push_back(indexTerm);

                xTerms.push_back(1.0/4);
                continue;
            }

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

            if (x.real() == 0 && x.imag() == 0) {
                xStarTerms.push_back(0);
                xTerms.push_back(0);
                continue;
            }
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

double BianchiMaassSearch::minBess(KBessel *K, const vector<Index> &indexTransversal, double Y) {
    double min = +INFINITY;

    for(const auto index : indexTransversal) {
        double arg = 2 * pi / A * index.getAbs(d) * Y;
        if (arg > K->getR()) {
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

    return min;
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


tuple<vector<pair<double,double>>, double, double> BianchiMaassSearch::fineSecantMethod(double leftR, double rightR) {
    std::cout << "Starting final refinement." << std::endl;
    std::cout << std::setprecision(16) << "[" << leftR << ", " << rightR << "]" << std::endl;

    double M0 = computeM0General(rightR);

    KBessel* leftRK = new KBessel(2 * pi / A /*usual factor*/
                        * 1 /*smallest magnitude of an index*/
                        * Y0 /*smallest height of a pullback*/
            , leftR);

    KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, rightR);

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

    double Y = computeWellConditionedY(leftRK, rightRK, leftR, rightR, M0, indexTransversal);

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

    leftRK->setRAndPrecompute(leftR, 2 * pi / A * M0 * maxYStar);
    MatrixXd matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, (*leftRK));
    auto solAndCondYR1 = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

    rightRK->setRAndPrecompute(rightR, 2 * pi / A * M0 * maxYStar);
    matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, (*rightRK));
    auto solAndCondYR2 = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

    std::cout << solAndCondYR1.second << " " << solAndCondYR2.second << std::endl;


    map<Index, double> coeffMapR1, coeffMapR2;

    for (int i = 0; i < indexTransversal.size(); i++) {
        Index index = indexTransversal[i];
        coeffMapR1[index] = solAndCondYR1.first(i);
        coeffMapR2[index] = solAndCondYR2.first(i);
    }

    //We have the solutions from everything, now we cook up the phi function

    double coeff0R1 = coeffMapR1.at(dToPrimes.at(d)[0]);
    double coeff1R1 = coeffMapR1.at(dToPrimes.at(d)[1]);
    double coeff2R1 = coeffMapR1.at(dToPrimes.at(d)[2]);

    double coeff0R2 = coeffMapR2.at(dToPrimes.at(d)[0]);
    double coeff1R2 = coeffMapR2.at(dToPrimes.at(d)[1]);
    double coeff2R2 = coeffMapR2.at(dToPrimes.at(d)[2]);

    Index prod01 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[1], d);
    Index prod02 = dToPrimes.at(d)[0].mul(dToPrimes.at(d)[2], d);

    double hecke1R1 = coeff0R1 * coeff1R1 - coeffMapR1.at(prod01);
    double hecke2R1 = coeff0R1 * coeff2R1 - coeffMapR1.at(prod02);
    hecke1R1 = abs(hecke1R1);
    hecke2R1 = abs(hecke2R1);

    double hecke1R2 = coeff0R2 * coeff1R2 - coeffMapR2.at(prod01);
    double hecke2R2 = coeff0R2 * coeff2R2 - coeffMapR2.at(prod02);
    hecke1R2 = abs(hecke1R2);
    hecke2R2 = abs(hecke2R2);


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
                phi1 = epsItr1*hecke1R1 + epsItr2*hecke2R1;
                phi2 = epsItr1*hecke1R2 + epsItr2*hecke2R2;
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

    double rN = leftR;
    double phiN = phi1;
    double rNPlus1 = rightR;
    double phiNPlus1 = phi2;

    std::cout << std::setprecision(16) << "r = " << leftR << ", " << heckeCheck(coeffMapR1) << std::endl;
    std::cout << std::setprecision(16) << "r = " << rightR << ", " << heckeCheck(coeffMapR2) << std::endl;

    vector<pair<double, double>> heckeValues;
    heckeValues.emplace_back(rN, heckeCheck(coeffMapR1));
    heckeValues.emplace_back(rNPlus1, heckeCheck(coeffMapR2));

    int minIterations = 4;
    int maxIterations = 30;
    int iterations = 0;

    while (true) {
        iterations++;

        double nextR = findZeroOfLinearInterpolation(rN, phiN, rNPlus1, phiNPlus1);

        KBessel K = KBessel(twoPiOverA * 1 * Y0, nextR);
        K.extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, K);
        auto solAndCondYNextR = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

        std::cout << solAndCondYNextR.second << " " << solAndCondYNextR.second << std::endl;

        map<Index, double> coeffMapNextR;
        for (int i = 0; i < indexTransversal.size(); i++) {
            Index index = indexTransversal[i];
            coeffMapNextR[index] = solAndCondYNextR.first(i);
        }

        double hecke = heckeCheck(coeffMapNextR);
        heckeValues.emplace_back(nextR, hecke);
        std::cout << std::setprecision(16) << "r = " << nextR << ", " << hecke << std::endl;

        //We have the solutions from everything, now we cook up the phi function

        double coeff0NextR = coeffMapNextR.at(dToPrimes.at(d)[0]);
        double coeff1NextR = coeffMapNextR.at(dToPrimes.at(d)[1]);
        double coeff2NextR = coeffMapNextR.at(dToPrimes.at(d)[2]);

        double hecke1NextR = coeff0NextR * coeff1NextR - coeffMapNextR.at(prod01);
        double hecke2NextR = coeff0NextR * coeff2NextR - coeffMapNextR.at(prod02);
        hecke1NextR = abs(hecke1NextR);
        hecke2NextR = abs(hecke2NextR);

        double nextPhi = eps1 * hecke1NextR + eps2 * hecke2NextR;

        /*
         * This is the Anderson-Bjorck bracketed root finding algorithm
         */
        if (areDifferentSign(phiNPlus1, nextPhi)) {
            rN = rNPlus1;
            phiN = phiNPlus1;

            rNPlus1 = nextR;
            phiNPlus1 = nextPhi;
        } else {
            double m = 1 - nextPhi/phiNPlus1;
            if (m <= 0) {
                m = 0.5;
            }
            phiN *= m;

            rNPlus1 = nextR;
            phiNPlus1 = nextPhi;
        }

        if (iterations < minIterations) {
            continue;
        }
        if (hecke > 1 //not modular
            || heckeHasConverged(heckeValues) //computation succeeded and is complete
            || iterations > maxIterations //computation not converging as expected
            || phiN == phiNPlus1 //numerical error incoming, could just be very good convergence
            || rN == rNPlus1) {
            watch(hecke);
            watch(heckeHasConverged(heckeValues));
            watch(iterations > maxIterations);
            watch(phiN == phiNPlus1);
            watch(rN == rNPlus1);

            for (auto const& itr : coeffMapNextR) {
                std::cout << std::setprecision(16) << itr.first << ", " << itr.second << "\n";
            }
            std::cout << std::flush;
            delete leftRK;
            delete rightRK;
            return {heckeValues, rN, rNPlus1};

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

double BianchiMaassSearch::computeWellConditionedY(KBessel *K, double r, double M0, vector<Index> &indexTransversal) {
    double c = max(1.0, r)/(2*pi/A*M0);
    double upperY = Y0;
    double lowerY = 0.5 * max(1.0, r)/(2*pi/A*M0);

    vector<double> heightsToTry;

    double tempY = upperY;
    while (tempY > lowerY) {
        heightsToTry.push_back(tempY);
        tempY *= 0.99;
    }
    std::sort(heightsToTry.begin(), heightsToTry.end());

    vector<double> minDiag(heightsToTry.size(), 0);


#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, minDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        minDiag[i] = minBess(K, indexTransversal, Y);
    }


    double Y = 0;
    double maxSoFar = 0;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double maxMinBess = minDiag[i];
        if (maxSoFar < maxMinBess) {
            maxSoFar = maxMinBess;
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

    minDiag.clear();

    minDiag.resize(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, minDiag, K)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        minDiag[i] = minBess(K, indexTransversal, testY);
    }

    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double maxMinBess = minDiag[i];
        if (maxSoFar < maxMinBess) {
            maxSoFar = maxMinBess;
            Y = heightsToTry[i];
        }
    }
    return Y;
}

double
BianchiMaassSearch::computeWellConditionedY(KBessel *leftRK, KBessel *rightRK, double leftR, double rightR, double M0,
                                            const vector<Index> &indexTransversal) {
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

    vector<double> r1MinDiag(heightsToTry.size(), 0);
    vector<double> r2MinDiag(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r2MinDiag, rightRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r2MinDiag[i] = minBess(rightRK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r1MinDiag, leftRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r1MinDiag[i] = minBess(leftRK, indexTransversal, Y);
    }


    double Y = 0;
    double maxSoFar = 0;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double maxMinBess = std::min({r1MinDiag[i],
                                      r2MinDiag[i]});
        if (maxSoFar < maxMinBess) {
            maxSoFar = maxMinBess;
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

    r1MinDiag.clear();
    r2MinDiag.clear();

    r1MinDiag.resize(heightsToTry.size(), 0);
    r2MinDiag.resize(heightsToTry.size(), 0);


#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r2MinDiag, rightRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        r2MinDiag[i] = minBess(rightRK, indexTransversal, testY);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r1MinDiag, leftRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double testY = heightsToTry[i];
        r1MinDiag[i] = minBess(leftRK, indexTransversal, testY);
    }

    //heightsToTry is sorted smallest to largest
    for (int i = 0; i < heightsToTry.size(); i++) {//end i at second to last because Y2 has to be greater
        double maxMinBess = std::min({r1MinDiag[i],
                                      r2MinDiag[i]});
        if (maxSoFar < maxMinBess) {
            maxSoFar = maxMinBess;
            Y = heightsToTry[i];
        }
    }
    return Y;
}

pair<double, double>
BianchiMaassSearch::computeTwoWellConditionedY(KBessel *leftRK, KBessel *rightRK, double leftR, double rightR,
                                               double M0, const vector<Index> &indexTransversal) {
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

    vector<double> r1MinDiag(heightsToTry.size(), 0);
    vector<double> r2MinDiag(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r2MinDiag, rightRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r2MinDiag[i] = minBess(rightRK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r1MinDiag, leftRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r1MinDiag[i] = minBess(leftRK, indexTransversal, Y);
    }


    double ansY1 = 0;
    double ansY2 = 0;
    double maxSoFar = 0;
    //heightsToTry is sorted smallest to largest
    for (int i = 0; i <= heightsToTry.size() - 2; i++) {//end i at second to last because Y2 has to be greater
        for (int j = i + 1; j < heightsToTry.size(); j++) {//j is always less than i, so j indexes a larger Y
            double maxMinBess = std::min({r1MinDiag[i],
                                          r1MinDiag[j],
                                          r2MinDiag[i],
                                          r2MinDiag[j]});
            if (maxSoFar < maxMinBess) {
                maxSoFar = maxMinBess;
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

    r1MinDiag.clear();
    r2MinDiag.clear();

    r1MinDiag.resize(heightsToTry.size(), 0);
    r2MinDiag.resize(heightsToTry.size(), 0);

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r2MinDiag, rightRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r2MinDiag[i] = minBess(rightRK, indexTransversal, Y);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, r1MinDiag, leftRK)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        r1MinDiag[i] = minBess(leftRK, indexTransversal, Y);
    }

    //heightsToTry is sorted smallest to largest
    for (int i = 0; i <= heightsToTry.size() - 2; i++) {//end i at second to last because Y2 has to be greater
        for (int j = i + 1; j < heightsToTry.size(); j++) {//j is always less than i, so j indexes a larger Y
            double maxMinBess = std::min({r1MinDiag[i],
                                          r1MinDiag[j],
                                          r2MinDiag[i],
                                          r2MinDiag[j]});
            if (maxSoFar < maxMinBess) {
                maxSoFar = maxMinBess;
                ansY1 = heightsToTry[i];
                ansY2 = heightsToTry[j];
            }
        }
    }
    return {min(ansY1, ansY2), max(ansY1, ansY2)};
}

void BianchiMaassSearch::saveSearchResult(vector<double> spectralParameters,
                                          vector<Index>& indexTransversal,
                                          vector<double>& coeffs) {
    string directory = "Output/SearchResults/" + to_string(d) + "/" + symClass + "/";
    std::stringstream ss;
    ss << std::setprecision(16) << spectralParameters[2];
    string prefix = ss.str() + ".txt";

    createOutputDirectory(directory);

    ofstream outputFile;
    outputFile.open(directory + prefix, std::ofstream::out | std::ofstream::app);

    outputFile << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << ", D = " << to_string(D) << "\n";
    outputFile << "[" << spectralParameters[0] << ", " << spectralParameters[1] << "], ";
    outputFile << spectralParameters[2];

    for (int i = 0; i < indexTransversal.size(); i++) {
        auto index = indexTransversal[i];
        double coeff = coeffs[i];
        outputFile << "\n" << index.getA() << ", " << index.getB() << ", " << coeff;
    }
    outputFile << std::flush;
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

vector<pair<double, double>> BianchiMaassSearch::getIntervalsForCoarseSearch(double startR, double endR) {

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
    double numEigenValues = 0.1;

    vector<double> endpoints;
    double endpoint = max(startR, coarseComplete);
    while (endpoint < endR) {
        endpoints.push_back(endpoint);
        endpoint = Od.eigenvalueIntervalRightEndpoint(endpoint, numEigenValues);
    }
    if (endR > 0 && endR <= 1e-11) {
        endpoints.push_back(2e-11);
    } else {
        endpoints.push_back(endR);
    }


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
        } else if (line[0] == '[') {
            //remove brackets
            string numbers = line.substr(1, line.size() - 2);
            int comma = numbers.find(',');
            string leftStr = numbers.substr(0, comma);
            string rightStr = numbers.substr(comma + 2);
            double left = std::stod(leftStr);
            double right = std::stod(rightStr);
            mediumComplete = right;
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
            if (left < mediumComplete && mediumComplete < right) {
                intervalsFromFile.emplace_back(mediumComplete, right);
            } else if (mediumComplete < right) {
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

void BianchiMaassSearch::computeMaximumD(double r, int timeLimitSeconds) {
    double maxNanoseconds = timeLimitSeconds * pow(10,9);

    double maxTerms = maxNanoseconds / nanosecondsPerTerm;

    double minD = 3;
    double maxD = 16;
    if (d == 163) {
        maxD = 12;
    }
    double leftD = minD; //I don't want to compute anything with any less precision than this
    double rightD = maxD; //This is the most precision I will ever work to


    unsigned long long int leftTerms = 0;
    unsigned long long int rightTerms = 0;

    truncation = pow(10.0, -leftD);
    double M0 = computeM0General(r);

    KBessel* K = new KBessel(twoPiOverA * 1 * Y0, r);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);

    double Y = computeWellConditionedY(K,r, M0, indexTransversal);
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
    M0 = computeM0General(r);

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
            answer = std::min(6.0, rightD * 0.3);
        } else if (mode == "medium") {
            answer =  std::min(6.0, rightD * 0.5);
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
        M0 = computeM0General(r);

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
        answer = std::min(6.0, rightD * 0.3);
    } else if (mode == "medium") {
        answer =  std::min(6.0, rightD * 0.5);
    } else {
        //mode == "fine"
        answer =  rightD;
    }

    D = max(minD, answer);
    truncation = pow(10.0, -D);
}

