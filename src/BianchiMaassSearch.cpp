#include "../include/BianchiMaassSearch.h"

#include <cmath>
#include <stack>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <map>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <eigen3/Eigen/SVD>
#include <fstream>
#include <boost/align/aligned_allocator.hpp>
#include <sstream>


//Containers
using std::find, std::reverse, std::stack;

//Math
using std::cos, std::sin, std::get, std::max, std::pow, std::abs, std::ceil, std::floor, std::min;

//IO
using std::cout, std::string, std::endl, std::flush, std::setprecision, std::to_string;

#define watch(x) cout << (#x) << " is " << (x) << endl << flush

BianchiMaassSearch::BianchiMaassSearch(string mode, int d, double D, char symClass) {

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


    this->D = D;
    if (D == 0) {
        autoPrecision = true;
    } else {
        autoPrecision = false;
        truncation = pow(10.0, -(this->D));
    }

    if (mode == "test-modularity-compute-reflections") {
        autoPrecision = false;
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

    int indexOfNormalization = getIndexOfNormalization(indexTransversal);

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

    pair<Eigen::Matrix<double, -1, -1, 0>, double> solAndCondY1NextR = solAndCondY1R1;

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

        //this is the median
        if (c_vec.size() == 0) {

            vector<double> coeffs;
            for (int i = 0; i < indexTransversal.size(); i++) {

                double coeff = solAndCondY1NextR.first(i);
                coeffs.push_back(coeff);
            }
            c = a;
            std::cout << "ZERO SIGN CHANGES DETECTED. CONDUCT FURTHER INVESTIGATION\n";
            std::cout << setprecision(16) << "[" << a << ", " << b << "]" << std::endl;
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            saveSearchResult({a, b, c}, indexTransversal, coeffs, NAN);
            break;
        } else if (c_vec.size() % 2 == 1) {

            c = c_vec[(c_vec.size() - 1)/2];
        } else {

            c = (c_vec[c_vec.size()/2 - 1] + c_vec[c_vec.size()/2])/2.0;
        }


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

            saveSearchResult({a, b, c}, indexTransversal, coeffs, Y1);
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

        saveSearchResult({a, b, c}, indexTransversal, coeffs, Y1);
    }
}


bool BianchiMaassSearch::possiblyContainsEigenvalue(const double leftR, const double rightR, KBessel *leftRK, KBessel *rightRK) {


    std::cout << std::setprecision(16)
              << "[" << leftR << ", " << rightR << "]" << std::endl;

    double M0 = computeM0General(rightR);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

    int indexOfNormalization = getIndexOfNormalization(indexTransversal);

    auto pairOfYs = computeTwoWellConditionedY(leftRK, rightRK, leftR, rightR, M0, indexTransversal);

    double Y1 = pairOfYs.first;
    double Y2 = pairOfYs.second;

    /*
     * Y1 and Y2 have been computed. Proceed with computing the actual matrices
     */
    double MY1 = computeMYGeneral(M0, Y1);
    double MY2 = computeMYGeneral(M0, Y2);

    watch(M0);
    watch(MY1);
    watch(MY2);
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
        //minM0 /= 2;
        //evalLeft = truncation * peak - bess.exactKBessel(2*pi/A*minM0*Y0);
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

int BianchiMaassSearch::getIndexOfNormalization(const vector<Index> & indexTransversal) {
    if (d < 7) {
        throw std::invalid_argument("Not implemented for d = 1,2,3");
    }

    if (symClass == 'D' || symClass == 'C') {
        for (int i = 0; i < indexTransversal.size(); i++) {
            auto n = indexTransversal[i];
            if (Index(1,0) == n) {
                return i;
            }
        }
        throw std::invalid_argument("Index for n = 1 not in the list of indices provided.");
    }

    for (int i = 0; i < indexTransversal.size(); i++) {
        auto n = indexTransversal[i];
        auto nComplex = n.getComplex(d);
        if (nComplex != conj(nComplex) && nComplex != -conj(nComplex)) {
            return i;
        }
    }
    throw std::invalid_argument("Indices provided must be incorrect");
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

            if (m == n) {
                answer(i,j) = Y * K.exactKBessel(twoPiOverA * m.getAbs(d) * Y) - entry;
            } else {
                answer(i,j) = -entry;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);

    nanosecondsPerTerm = ((double)duration.count())/((double)terms);

    //watch(nanosecondsPerTerm);
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


            double indexTerm = 0;
            for (const auto& tup : nIndexOrbitDataModSign) {
                Index l = tup.first;
                double ellTerm = tup.second; //This is a_l/a_n
                indexTerm += ellTerm * cos(piOverA * traceProduct(l.getComplex(d), I*xStar));
            }
            xStarTerms.push_back(indexTerm);

            double testPointTerm = cos(piOverA * traceProduct(-mComplex, I*x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                short sign = tup.second;
                testPointTerm += sign * cos(piOverA * traceProduct(-mComplex, I*etaX));
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
                indexTerm += ellTerm * sin(piOverA * traceProduct(l.getComplex(d), I*xStar));
            }
            xStarTerms.push_back(indexTerm);

            double testPointTerm = sin(piOverA * traceProduct(-mComplex, I*x));
            for (const auto& tup : testPointOrbit.properTranslatesModSign) {
                complex<double> etaX = tup.first;
                short sign = tup.second;
                testPointTerm += sign * sin(piOverA * traceProduct(-mComplex, I*etaX));
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
        answer *= 4.0 / pointCount;
    } else {
        answer *= -4.0 / pointCount;
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

void BianchiMaassSearch::createOutputDirectory(const std::string& directory) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }
}

bool BianchiMaassSearch::isFileEmpty(const string &filename) {
    return std::filesystem::is_empty(filename);
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
    upperY = min(Y0, upperY);
    lowerY = Y - 0.5 * c;
    lowerY = max(lowerY, 0.5 * Y);

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
    upperY = min(Y0, upperY);
    lowerY = Y - 0.5 * c;
    lowerY = max(lowerY, 0.5 * Y);

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
    upperY = min(Y0, upperY);
    lowerY = ansY1 - 0.5 * c;
    lowerY = max(lowerY, 0.5 * ansY1);

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

void BianchiMaassSearch::saveSearchResult(vector<double> spectralParameters, vector<Index> &indexTransversal,
                                          vector<double> &coeffs, double Y) {
    string directory = "Output/Search/Final/" + to_string(d) + "/" + symClass + "/";
    std::stringstream ss;
    ss << std::setprecision(16) << spectralParameters[2];
    string prefix = ss.str() + ".txt";

    createOutputDirectory(directory);

    ofstream outputFile;
    outputFile.open(directory + prefix, std::ofstream::out | std::ofstream::app);

    outputFile << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << ", D = " << to_string(D) << ", Y = " << Y << "\n";
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
        const std::string directory = "Output/Search/Coarse/";
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

        const std::string directoryOut = "Output/Search/Medium/";
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
        const std::string directoryOut = "Output/Search/Fine/";
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

    const std::string directoryIn = "Output/Search/Coarse/";
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

    const std::string directoryOut = "Output/Search/Medium/";
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

    const std::string directoryIn = "Output/Search/Coarse/"; // Change this to the desired directory
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
    const std::string directoryOut = "Output/Search/Fine/";
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

    const std::string directoryIn = "Output/Search/Medium/";
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

void BianchiMaassSearch::testForModularity() {

    //The symClass has been verified to by one of D or C
    //The first step is to test the forms in symClass for modularity
    //Gather which ones have already been verified modular
    string outDirectory = "Output/TestedData/TestedForms/" + to_string(d) + "/" + symClass + "/";
    createOutputDirectory(outDirectory);

    vector<string> completeSpecParams;
    vector<string> incompleteSpecParams;
    for (const auto& entry : std::filesystem::directory_iterator(outDirectory)) {
        if (entry.is_regular_file()) {
            auto filename = entry.path().filename().string();
            filename = filename.substr(0, filename.find(".txt"));
            completeSpecParams.push_back(filename);
        }
    }

    //Gather the ones that haven't been done yet (or which failed to be verified)
    string inDirectory = "Output/Search/Final/" + to_string(d) + "/" + symClass + "/";
    for (const auto& entry : std::filesystem::directory_iterator(inDirectory)) {
        if (entry.is_regular_file()) {
            auto filename = entry.path().filename().string();
            filename = filename.substr(0, filename.find(".txt"));
            if (std::find(completeSpecParams.begin(),
                          completeSpecParams.end(),
                          filename)
                != completeSpecParams.end()) {
                continue;
            }

            incompleteSpecParams.push_back(filename);
        }
    }



    //Read in the coefficients for each incompleteSpecParam and evaluate at points to check modularity
    //If it is modular.................
    for (const auto s : incompleteSpecParams) {
        std::stringstream ss;
        ss << "Output/Search/Final/" << to_string(d) << "/" << symClass << "/";
        ss << s << ".txt";
        string inFilename = ss.str();
        ss.clear();

        double r = stod(s);

        std::ifstream inFile(inDirectory + s + ".txt");
        string stringD;
        getline(inFile, stringD); //this line has D = .... in it
        int startD = stringD.find("D = ") + 4;
        int nextComma = stringD.find(",", startD);
        if (nextComma != string::npos) {
            stringD = stringD.substr(stringD.find("D = ") + 4, nextComma - startD);
        } else {
            stringD = stringD.substr(stringD.find("D = ") + 4, stringD.length() - startD);
        }

        double DForThisForm = std::stod(stringD);

        string line;
        getline(inFile, line); //this line has the interval in it
        int leftBracket = 0;
        int firstComma = line.find(",", 0);
        string leftEndpoint = line.substr(leftBracket + 1, firstComma - leftBracket - 1);
        int rightBracket = line.find("]", 0);
        string rightEndpoint = line.substr(firstComma + 2, rightBracket - (firstComma + 2));
        string confidentR = Auxiliary::commonRound(leftEndpoint, rightEndpoint);


        D = DForThisForm;
        truncation = pow(10,-(this->D));

        KBessel* Kb = new KBessel(twoPiOverA * 1 * Y0, r);

        double M0 = computeM0General(r);

        vector<Index> indicesM0 = Od.indicesUpToM(M0);
        auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
        vector<Index> indexTransversal = get<0>(data);
        map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

        int indexOfNormalization = getIndexOfNormalization(indexTransversal);

        double Y = computeWellConditionedY(Kb, r, M0, indexTransversal);

        double MY = computeMYGeneral(M0, Y);

        double maxYStar = 0;

        map<Index, vector<TestPointOrbitData>> mToYTestPointOrbits;

#pragma omp parallel for default(none) shared(indexTransversal, mToYTestPointOrbits, Y, MY)
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

        Kb->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        MatrixXd matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, *Kb);

        auto solAndCondY = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

        vector<double> coeffs;
        for (int i = 0; i < indexTransversal.size(); i++) {
            coeffs.push_back(solAndCondY.first(i));
        }

        vector<pair<Quaternion, Quaternion>> zAndZStar;
        int numPointsToCheck = 100;
        while (zAndZStar.size() < numPointsToCheck) {
            Quaternion z = Quaternion(Auxiliary::rand(-0.5, 0.5), Auxiliary::rand(-A/2.0, A/2.0), Y0 - Auxiliary::rand(0, Y0/5), 0);
            Quaternion zstar = Quaternion(z);
            zstar.reduce(d);

            if ((z - zstar).abs() > 1) {
                zAndZStar.emplace_back(z, zstar);
            }
        }

        KBessel K = KBessel(twoPiOverA * 1 * Y0, r);
        vector<pair<double, double>> fzAndfZStar;
        for (const auto& zZStar : zAndZStar) {
            double fz = evaluate(zZStar.first, K, indexTransversal, coeffs);
            double fZStar = evaluate(zZStar.second, K, indexTransversal, coeffs);
            pair<double,double> fzfZStar = {fz, fZStar};
            fzAndfZStar.push_back(fzfZStar);
        }

        double maxDiff = 0;
        for (const auto& fzFZStar : fzAndfZStar) {
            double diff = fzFZStar.first - fzFZStar.second;
            if (maxDiff < abs(diff)) {
                maxDiff = abs(diff);
            }
        }

        if (maxDiff < 0.01) {
            std::cout << std::setprecision(16) << s << " is modular for symClass " << symClass << std::endl;
            //It's modular!

            string destination = "Output/TestedData/TestedForms/" + to_string(d) + "/" + symClass + "/" + s + ".txt";
            std::ofstream dest(destination, std::ios::binary);

            dest << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << ", r = " << r << ", D = " << D << ", Y = " << Y;

            dest << "\nmax of |f(z) - f(z^*)| for " << to_string(numPointsToCheck) << " points = " << std::setprecision(16) << maxDiff << "; ";

            /*
             *
             * Hecke
             */

            auto index0 = std::find(indexTransversal.begin(), indexTransversal.end(),dToPrimes[d][0]) - indexTransversal.begin();
            auto index1 = std::find(indexTransversal.begin(), indexTransversal.end(),dToPrimes[d][1]) - indexTransversal.begin();
            auto index2 = std::find(indexTransversal.begin(), indexTransversal.end(),dToPrimes[d][2]) - indexTransversal.begin();

            auto indexProd01 = dToPrimes[d][0].mul(dToPrimes[d][1],d);
            auto indexProd12 = dToPrimes[d][1].mul(dToPrimes[d][2],d);
            auto indexProd02 = dToPrimes[d][0].mul(dToPrimes[d][2],d);

            auto indexItr = find(indexTransversal.begin(), indexTransversal.end(), indexProd01);
            if (indexItr != indexTransversal.end()) {
                auto index = indexItr - indexTransversal.begin();
                double heckeRelation = coeffs[index0] * coeffs[index1] - coeffs[index];
                dest << std::setprecision(16) << "a(" << dToPrimes[d][0] << ")*a(" << dToPrimes[d][1] <<") - a(" << indexProd01 << ") = " << heckeRelation << ", ";
            }

            indexItr = find(indexTransversal.begin(), indexTransversal.end(), indexProd02);
            if (indexItr != indexTransversal.end()) {
                auto index = indexItr - indexTransversal.begin();
                double heckeRelation = coeffs[index0] * coeffs[index2] - coeffs[index];
                dest << std::setprecision(16) << "a(" << dToPrimes[d][0] << ")*a(" << dToPrimes[d][2] <<") - a(" << indexProd02 << ") = " << heckeRelation << ", ";
            }

            indexItr = find(indexTransversal.begin(), indexTransversal.end(), indexProd12);
            if (indexItr != indexTransversal.end()) {
                auto index = indexItr - indexTransversal.begin();
                double heckeRelation = coeffs[index1] * coeffs[index2] - coeffs[index];
                dest << std::setprecision(16) << "a(" << dToPrimes[d][1] << ")*a(" << dToPrimes[d][2] <<") - a(" << indexProd12 << ") = " << heckeRelation << ", ";
            }

            /*
             *
             *
             * end Hecke
             */

            dest << "\nconfident_r = " << std::setprecision(16) << confidentR;

            for (int i = 0; i < indexTransversal.size(); i++) {
                dest << "\n" << std::setprecision(16) << indexTransversal[i].getA() << ", " << indexTransversal[i].getB() << ", " << coeffs[i];
            }

            dest.close();
        } else {
            std::cout << std::setprecision(16) << s << " is NOT modular for symClass " << symClass << std::endl;
        }

    }

    /**
     * Now, go through the tested forms and check if the reflection-odd counterpart is also modular
     * If it is, save it to the respective folder
     *
     *
     *
     *
     *
     */

    string testedDirectory = "Output/TestedData/TestedForms/" + to_string(d) + "/" + symClass + "/";

    vector<string> testedSpecParams;
    for (const auto& entry : std::filesystem::directory_iterator(testedDirectory)) {
        if (entry.is_regular_file()) {
            auto filename = entry.path().filename().string();
            auto itr = filename.find(".txt");
            if (itr != string::npos) {
                filename = filename.substr(0, itr);
                testedSpecParams.push_back(filename);
            }
        }
    }

    std::cout << "Going to test these spectral parameters for modularity in the corresponding symmetry class of reflection-odd forms" << std::endl;

    for (auto& s : testedSpecParams) {
        std::cout << s << std::endl;
    }
    char originalSymClass = symClass;

    if (originalSymClass == 'D') {
        symClass = 'G';
    } else {
        //symClass is 'C'
        symClass = 'H';
    }

    for (auto s : testedSpecParams) {
        double r = std::stod(s);

        std::ifstream inFile(inDirectory + s + ".txt");
        string stringD;
        getline(inFile, stringD); //this line has D = .... in it
        int startD = stringD.find("D = ") + 4;
        int nextComma = stringD.find(",", startD);
        if (nextComma != string::npos) {
            stringD = stringD.substr(stringD.find("D = ") + 4, nextComma - startD);
        } else {
            stringD = stringD.substr(stringD.find("D = ") + 4, stringD.length() - startD);
        }

        double DForThisForm = std::stod(stringD);

        string line;
        getline(inFile, line); //this line has the interval in it
        int leftBracket = 0;
        int firstComma = line.find(",", 0);
        string leftEndpoint = line.substr(leftBracket + 1, firstComma - leftBracket - 1);
        int rightBracket = line.find("]", 0);
        string rightEndpoint = line.substr(firstComma + 2, rightBracket - (firstComma + 2));
        string confidentR = Auxiliary::commonRound(leftEndpoint, rightEndpoint);


        D = DForThisForm;
        truncation = pow(10,-(this->D));

        KBessel* Kb = new KBessel(twoPiOverA * 1 * Y0, r);

        double M0 = computeM0General(r);

        vector<Index> indicesM0 = Od.indicesUpToM(M0);
        auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
        vector<Index> indexTransversal = get<0>(data);
        map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

        int indexOfNormalization = getIndexOfNormalization(indexTransversal);

        double Y = computeWellConditionedY(Kb, r, M0, indexTransversal);

        double MY = computeMYGeneral(M0, Y);

        double maxYStar = 0;

        map<Index, vector<TestPointOrbitData>> mToYTestPointOrbits;

#pragma omp parallel for default(none) shared(indexTransversal, mToYTestPointOrbits, Y, MY)
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

        Kb->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        MatrixXd matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, *Kb);

        auto solAndCondY = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

        vector<double> coeffs;
        for (int i = 0; i < indexTransversal.size(); i++) {
            coeffs.push_back(solAndCondY.first(i));
        }

        double coeffavg = 0;
        for (const auto & coeff : coeffs) {
            coeffavg += abs(coeff);
        }
        coeffavg /= (double)indexTransversal.size();
        if (coeffavg > 10) {
            std::cout << "warning ";
            watch(coeffavg);
        }

        vector<pair<Quaternion, Quaternion>> zAndZStar;
        int numPointsToCheck = 100;
        while (zAndZStar.size() < numPointsToCheck) {
            Quaternion z = Quaternion(Auxiliary::rand(-0.5, 0.5), Auxiliary::rand(-A/2.0, A/2.0), Y0 - Auxiliary::rand(0, Y0/5), 0);
            Quaternion zstar = Quaternion(z);
            zstar.reduce(d);

            if ((z - zstar).abs() > 0.1) {
                zAndZStar.emplace_back(z, zstar);
            }
        }

        vector<pair<double, double>> fzAndfZStar;
        for (const auto& zZStar : zAndZStar) {
            double fz = evaluate(zZStar.first, *Kb, indexTransversal, coeffs);
            double fZStar = evaluate(zZStar.second, *Kb, indexTransversal, coeffs);
            pair<double,double> fzfZStar = {fz, fZStar};
            fzAndfZStar.push_back(fzfZStar);
        }

        double maxDiff = 0;
        for (const auto& fzFZStar : fzAndfZStar) {
            double diff = fzFZStar.first - fzFZStar.second;
            if (maxDiff < abs(diff)) {
                maxDiff = abs(diff);
            }
        }

        watch(maxDiff);
        if (maxDiff < 0.01) {
            //It's modular!
            std::cout << std::setprecision(16) << s << " is modular for symClass " << symClass << std::endl;
            string directory = "Output/TestedData/TestedForms/" + to_string(d) + "/" + symClass + "/";
            createOutputDirectory(directory);

            string filename = directory + s + ".txt";
            ofstream outputFile;
            outputFile.open(filename, std::ofstream::out);

            outputFile << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << ", r = " << r << ", D = " << D << ", Y = " << Y << std::endl;
            outputFile << "max of |f(z) - f(z^*)| for " << to_string(numPointsToCheck) << " points = " << std::setprecision(16) << maxDiff;
            outputFile << "\nconfident_r = " << std::setprecision(16) << confidentR;

            for (int i = 0; i < indexTransversal.size(); i++) {
                auto index = indexTransversal[i];
                double coeff = coeffs[i];
                outputFile << std::setprecision(16) << "\n" << index.getA() << ", " << index.getB() << ", " << coeff;
            }
            outputFile << std::flush;
            outputFile.close();
        } else {
            std::cout << std::setprecision(16) << s << " is NOT modular for symClass " << symClass << std::endl;
        }
    }

}

void BianchiMaassSearch::testConjectures() {


    //Gather the ones that haven't been done yet
    string inDirectory = "Output/TestedData/TestedForms/" + to_string(d) + "/" + symClass + "/";
    vector<string> incompleteSpecParams;
    for (const auto& entry : std::filesystem::directory_iterator(inDirectory)) {
        if (entry.is_regular_file()) {
            auto filename = entry.path().filename().string();
            filename = filename.substr(0, filename.find(".txt"));
            double specParam = stod(filename);
            /*if (std::find(completeSpecParams.begin(),
                          completeSpecParams.end(),
                          specParam)
                != completeSpecParams.end()) {
                continue;
            }*/

            incompleteSpecParams.push_back(filename);
        }
    }

    for (const auto& s : incompleteSpecParams) {
        std::ifstream inFile(inDirectory + s + ".txt");


        //read in coefficients
        string line;
        std::getline(inFile, line); // read the modularity info
        std::getline(inFile, line); //read info like d = 19, etc.
        std::getline(inFile, line); //read the eigenvalue interval

        //each line now contains a, b, a_n where index = a + b*theta and a_n is Fourier coefficient
        vector<Index> indices;
        vector<double> coeffs;
        while (std::getline(inFile, line)) {
            int firstComma = line.find(',');
            int secondComma = line.find(',', firstComma + 1);
            string a = line.substr(0, firstComma);
            string b = line.substr(firstComma + 2, secondComma - firstComma - 2);
            string coeff = line.substr(secondComma + 2);

            Index n(stoi(a), stoi(b));
            double an = stod(coeff);

            indices.push_back(n);
            coeffs.push_back(an);
        }

        string satoOutDirectory = "Output/Conjectures/SatoTate/" + to_string(d) + "/" + symClass + "/";
        createOutputDirectory(satoOutDirectory);

        string ramanujanOutDirectory = "Output/Conjectures/Ramanujan/" + to_string(d) + "/" + symClass + "/";
        createOutputDirectory(ramanujanOutDirectory);

        string satoOutFilename = satoOutDirectory + s + ".txt";
        string ramanujanOutFilename = ramanujanOutDirectory + s + ".txt";

        std::ofstream satoOut(satoOutFilename);
        std::ofstream ramanujanOut(ramanujanOutFilename);

        satoOut << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << std::endl;
        satoOut << std::setprecision(16) << s << std::endl;
        satoOut << std::endl;

        ramanujanOut << std::setprecision(16) << "d = " << d << ", symClass = " << symClass << std::endl;
        ramanujanOut << std::setprecision(16) << s << std::endl;
        ramanujanOut << "Listed are indices, coefficient that fail the Ramanujan conjecture" << std::endl;

        for (int i = 0; i < indices.size(); i++) {
            Index index = indices[i];
            if (!Od.isPrime(index)) {
                continue;
            }

            double coeff = coeffs[i];

            if (coeff > 2.0) {
                ramanujanOut << std::setprecision(16) << index.getA() << ", " << index.getB() << ", " << coeff << std::endl;
            } else {
                double sato = acos(coeff/2.0);
                satoOut << std::setprecision(16) << sato << std::endl;
            }
        }
    }
}

double BianchiMaassSearch::evaluate(const Quaternion &z, KBessel &K, const vector<Index> &indexTransversal,
                                    const vector<double> &coeffs) {
    auto x = z.getComplex();
    auto y = z.getJ();
    vector<double> terms;
    terms.resize(indexTransversal.size(), 0);
    if (d == 1 || d == 3 || d == 2) {
        throw std::invalid_argument("d = 1, 3, 2 evaluation not implemented");
    }

#pragma omp parallel for default(none) shared(indexTransversal, coeffs, x, y, K, terms)
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index n = indexTransversal[i];
        double an = coeffs[i];
        double term = an * K.exactKBessel(twoPiOverA * n.getAbs(d) * y);
        auto nComplex = n.getComplex(d);

        double cs = 0;
        if (symClass == 'D') {
            cs = cos(pi/A * traceProduct(nComplex, x * I));
            if (n.getB() != 0 && n.getA() + 0.5 * n.getB() != 0.0) {
                cs += cos(pi / A * traceProduct(-conj(nComplex), x * I));
            }
        } else if (symClass == 'C') {
            cs = sin(pi/A * traceProduct(nComplex, x * I));
            if (n.getB() != 0 && n.getA() + 0.5 * n.getB() != 0.0) {
                cs -= sin(pi / A * traceProduct(-conj(nComplex), x * I));
            }
        } else if (symClass == 'G') {
            cs = cos(pi/A * traceProduct(nComplex, x * I));
            if (n.getB() != 0 && n.getA() + 0.5 * n.getB() != 0.0) {
                cs -= cos(pi / A * traceProduct(-conj(nComplex), x * I));
            }
        } else if (symClass == 'H') {
            cs = sin(pi/A * traceProduct(nComplex, x * I));
            if (n.getB() != 0 && n.getA() + 0.5 * n.getB() != 0.0) {
                cs += sin(pi / A * traceProduct(-conj(nComplex), x * I));
            }
        }
        term *= cs;
        terms[i] = term;
    }
    double answer = Auxiliary::kahanSummation(terms);
    answer *= y;
    return answer;
}

/***
 * Computes condition numbers of many matrices V_m,n.
 * For demonstration purposes only.
 * @param leftR
 * @param rightR
 */
void BianchiMaassSearch::sandbox(const double leftR, const double rightR) {

    KBessel* leftRK = new KBessel(twoPiOverA * 1 * Y0, leftR);
    KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, rightR);


    ofstream outFile;
    string filename = "Output/Misc";
    createOutputDirectory(filename);
    filename += "/cond_" + to_string(d) + ".txt";
    outFile.open(filename);

    double M0 = computeM0General(rightR);
    watch(M0);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

    int indexOfNormalization = getIndexOfNormalization(indexTransversal);


    auto bestY = computeWellConditionedY(rightRK, rightR, M0, indexTransversal);
    watch(Y0);
    std::cout << "estimated best Y is " << std::setprecision(16) << bestY << std::endl;
    outFile << Y0 << std::endl;
    outFile << "estimated best Y is " << std::setprecision(16) << bestY << std::endl;

    double val = 0.06790910071637833;
    double shift = 0.005;

    vector<double> bruteY = {};
    bruteY.push_back(bestY);
    for (double Y = val + shift; Y > val - shift; Y -= (2.0*shift)/2000.0) {
        bruteY.push_back(Y);
    }
    bruteY.push_back(0.005);

    for (auto Y : bruteY) {
        double MY = computeMYGeneral(M0, Y);

        double maxYStar = 0;

        map<Index, vector<TestPointOrbitData>> mToY1TestPointOrbits;

        maxYStar = 0;


#pragma omp parallel for default(none) shared(indexTransversal, mToY1TestPointOrbits, Y, MY)
        for (int i = 0; i < indexTransversal.size(); i++) {
            Index m = indexTransversal[i];
            auto orbitsY1 = getPointPullbackOrbits(m, Y, MY);
#pragma omp critical
            {
                mToY1TestPointOrbits[m] = orbitsY1;
            }
        }

        for (const auto& m : indexTransversal) {
            for (const auto& itr : mToY1TestPointOrbits[m]) {
                if (maxYStar < itr.representativePullback.getJ()) {
                    maxYStar = itr.representativePullback.getJ();
                }
            }
        }

        rightRK->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
        MatrixXd matrixY1 = produceMatrix(indexTransversal, mToY1TestPointOrbits, indexOrbitDataModSign, Y, *rightRK);

        auto solAndCondY1R1 = solveMatrix(matrixY1, indexTransversal, indexOfNormalization);

        std::cout << std::setprecision(16) << "(" << Y << ", " << solAndCondY1R1.second << ")" << std::endl;
        outFile << std::setprecision(16) << "(" << Y << ", " << solAndCondY1R1.second << ")" << std::endl;
    }
}

/***
 * Computes the size of some pieces of data in an example computation.
 * For demonstration purposes only.
 * @param r
 */
void BianchiMaassSearch::sandbox2(double r) {

    std::ofstream outFile;
    outFile.open("Output/Misc/size_trend.txt", std::ofstream::out | std::ofstream::app);
    KBessel* Kb = new KBessel(twoPiOverA * 1 * Y0, r);

    double M0 = computeM0General(r);
    watch(M0);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

    int indexOfNormalization = getIndexOfNormalization(indexTransversal);

    double Y = computeWellConditionedY(Kb, r, M0, indexTransversal);
    watch(Y);
    double MY = computeMYGeneral(M0, Y);
    watch(MY);
    double maxYStar = 0;

    std::cout << "Matrix ix " << indexTransversal.size() << "x" << indexTransversal.size() << std::endl;

    unsigned long long int totalTerms = 0;
    for (int i = 0; i < indexTransversal.size(); i++) {
        auto m = indexTransversal[i];
        totalTerms += countPointPullbackOrbits(m, Y, MY) * indexTransversal.size();
    }
    std::cout << "Total terms " << totalTerms << std::endl;

    double avgNumTerms = ((double)totalTerms)/(indexTransversal.size()*indexTransversal.size());

    std::cout << "Avg num terms per entry " << std::setprecision(16) << avgNumTerms << std::endl;

    outFile << "d = " << d << ", ";
    outFile << "r = " << std::setprecision(16) << r << ", ";
    outFile << "D = " << D << ", ";
    outFile << "M0 = " << std::setprecision(16) << M0 << ", ";
    outFile << "Y = " << std::setprecision(16) << Y << ", ";
    outFile << "MY = " << std::setprecision(16) << MY << ", ";
    outFile << "matrix = " << indexTransversal.size() << ", ";
    outFile << "total_terms = " << totalTerms << ", ";
    outFile << "avg_terms = " << std::setprecision(16) << avgNumTerms;
    outFile << '\n';

}

/***
 * This method times a single instance of computing the matrix V_m,n.
 * For demonstration purposes only.
 * @param r
 */
void BianchiMaassSearch::sandbox3(double r) {

    /**
     * Compute M0
     *
     */
    auto start = std::chrono::high_resolution_clock::now();
    double M0 = computeM0General(r);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute M0 " << duration.count() << " ms" << std::endl;




    /**
     *
     * Compute indices
     */

    start = std::chrono::high_resolution_clock::now();
    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute index transversal and index orbit data " << duration.count() << " ms" << std::endl;

    /**
     *
     * Compute well-conditioned Y
     */
    start = std::chrono::high_resolution_clock::now();
    KBessel* Kb = new KBessel(twoPiOverA * 1 * Y0, r);
    double Y = computeWellConditionedY(Kb, r, M0, indexTransversal);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute well-conditioned Y " << duration.count() << " ms" << std::endl;


    /**
     *
     * Compute test points and pullbacks
     */
    start = std::chrono::high_resolution_clock::now();
    double MY = computeMYGeneral(M0, Y);
    map<Index, vector<TestPointOrbitData>> mToYTestPointOrbits;

#pragma omp parallel for default(none) shared(indexTransversal, mToYTestPointOrbits, Y, MY)
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index m = indexTransversal[i];
        auto orbitsY = getPointPullbackOrbits(m, Y, MY);
#pragma omp critical
        {
            mToYTestPointOrbits[m] = orbitsY;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute test points X " << duration.count() << " ms" << std::endl;


    /**
     *
     * Do K-Bessel precomputation
     */
    start = std::chrono::high_resolution_clock::now();
    double maxYStar = 0;
    for (const auto& m : indexTransversal) {
        for (const auto& itr : mToYTestPointOrbits[m]) {
            if (maxYStar < itr.representativePullback.getJ()) {
                maxYStar = itr.representativePullback.getJ();
            }
        }
    }
    Kb->extendPrecomputedRange(2 * pi / A * M0 * maxYStar);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute KBessel " << duration.count() << " ms" << std::endl;


    /**
     *
     * Compute the matrix from Hejhal's identity
     */
    start = std::chrono::high_resolution_clock::now();
    MatrixXd matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, *Kb);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Compute matrix " << duration.count() << " ms" << std::endl;

    /**
     *
     * Solve the matrix equation
     */
    start = std::chrono::high_resolution_clock::now();
    int indexOfNormalization = getIndexOfNormalization(indexTransversal);
    auto solAndCondY = solveMatrix(matrixY, indexTransversal, indexOfNormalization);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Solve matrix " << duration.count() << " ms" << std::endl;
}

/***
 * This method computes many values of the function which we maximize for heuristic conditioning.
 * For demonstration purposes only.
 * @param r
 */
void BianchiMaassSearch::sandbox4(double r) {

    KBessel* Kb = new KBessel(twoPiOverA * 1 * Y0, r);

    double M0 = computeM0General(r);

    vector<Index> indicesM0 = Od.indicesUpToM(M0);
    auto data = Od.indexOrbitQuotientData(indicesM0, symClass);
    vector<Index> indexTransversal = get<0>(data);
    map<Index, vector<pair<Index, int>>> indexOrbitDataModSign = get<1>(data);

    int indexOfNormalization = getIndexOfNormalization(indexTransversal);
    double upperY = Y0;
    double lowerY = 0.0001;

    vector<double> heightsToTry;
    for (double Y = upperY; Y > lowerY; Y -= Y0/5000.0) {
        heightsToTry.push_back(Y);
    }

    vector<pair<double, double>> heuristics(heightsToTry.size());

#pragma omp parallel for schedule(dynamic) default(none) shared(heightsToTry, indexTransversal, heuristics, Kb)
    for (int i = 0; i < heightsToTry.size(); i++) {
        double Y = heightsToTry[i];
        pair<double,double> my_pair = {Y, minBess(Kb, indexTransversal, Y)};
        heuristics[i] = my_pair;
    }

    std::ofstream outFile;
    outFile.open("Output/Misc/heuristic_" + to_string(d) + ".txt", std::ofstream::out);
    for (const auto my_pair : heuristics) {
        outFile << std::setprecision(16) << "(" << my_pair.first << ", " << my_pair.second << ")\n";
    }
}

/***
 * This method computes many values of the function which we maximize for heuristic conditioning.
 * For demonstration purposes only.
 * @param r
 */
void BianchiMaassSearch::sandbox5() {

    for (int thisd : vector<int>{19, 43, 67, 163}) {
        ImaginaryQuadraticIntegers O = ImaginaryQuadraticIntegers(thisd);
        double thisY0 = sqrt(2.0/thisd);
        A = O.getA();
        for (double thisD : vector<double>{4.0, 8.0, 16.0, 20.0}) {
            truncation = pow(10.0,-thisD);
            for (double r : vector<double>{10.0, 2*10.0, pow(2,2)*10.0, pow(2,3)*10.0}) {
                double M0 = computeM0General(r);
                KBessel* Kb = new KBessel(2*pi/A * 1 * thisY0, r);

                vector<Index> indicesM0 = O.indicesUpToM(M0);
                auto data = Od.indexOrbitQuotientData(indicesM0, 'D');
                vector<Index> indexTransversal = get<0>(data);

                auto bestY = computeWellConditionedY(Kb, r, M0, indexTransversal);

                std::cout << std::setprecision(16) << thisd << ", " << thisD << ", " << r << ", " << M0 << ", " << bestY << std::endl;
            }
        }
    }
}

