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

BianchiMaassSearch::BianchiMaassSearch(char mode, int d, double D, char symClass) {

    //Check that d is valid
    vector<int> classNumberOne = {1,2,3,7,11,19,43,67,163};
    bool dIsClassNumberOne = find(classNumberOne.begin(), classNumberOne.end(), d) != classNumberOne.end();
    if (!dIsClassNumberOne) {
        throw(std::invalid_argument("d should be one of {1,2,3,7,11,19,43,67,163}"));
    }

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
    dToDMap[19] = {5,5,8,16};
    dToDMap[43] = {4,4,4,6};
    dToDMap[67] = {3,3,3,4};
    dToDMap[163] = {2,2,2,3};

    switch (mode) {
        case '0':
            this->D = dToDMap[d][0];
            break;
        case '1':
            this->D = dToDMap[d][1];
            break;
        case '2':
            this->D = dToDMap[d][2];
            break;
        default:
            throw std::invalid_argument("mode incorrect");
    }
    truncation = pow(10.0, -(this->D));

    //Check that symclass is valid
    vector<char> symClasses = {'D','G','C','H'};
    bool symClassIsCorrect = find(symClasses.begin(), symClasses.end(), symClass) != symClasses.end();
    if (!symClassIsCorrect) {
        throw(std::invalid_argument("symClass should be one of D,G,C,H"));
    }

    this->d = d;
    this->symClass = symClass;

    Od = ImaginaryQuadraticIntegers(d);

    A = Od.getA();
    theta = Od.getTheta();
    Y0 = Od.getY0();
    twoPiOverA = 2 * pi / A;

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
            coarseOutputFile << "D = " << D << '\n';
            coarseOutputFile << "Complete up to " << 0.75 << std::endl;
        }
    } else {
        std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
    }

    //Magic number
    //Supposed to reflect that we are dispatching the searcher on intervals that we expect
    //to have an eigenvalue in them 2/10th of the time. Allows for easier conditioning and
    //guards against unusually close eigenvalues, but not conclusively.
    double numEigenValues = 0.2;

    vector<double> endpoints;
    double endpoint = leftR;
    while (endpoint < rightR) {
        endpoints.push_back(endpoint);
        endpoint = Od.eigenvalueIntervalRightEndpoint(endpoint, numEigenValues);
    }
    endpoints.push_back(rightR);

    KBessel* leftRK = new KBessel(twoPiOverA * 1 * Y0, endpoints[0]);
    KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, endpoints[1]);
    for (size_t i = 0; i <= endpoints.size() - 2; i++) {
        double left = endpoints[i];
        double right = endpoints[i + 1];
        rightRK = new KBessel(twoPiOverA * 1 * Y0, right);
        possiblyContainsEigenvalue(left, right, leftRK, rightRK);
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


    const std::string directoryOut = "Output/Medium/";
    const std::string prefixOut = "medium_"
                               + to_string(d) + "_"
                               + symClass;
    createOutputDirectory(directoryOut);

    std::string outputFilename = directoryOut
                                 + prefixOut
                                 + ".txt";

    mediumOutputFile.open(outputFilename, std::ofstream::out | std::ofstream::app);

    double mediumComplete = 0.75;

    if (mediumOutputFile.is_open()) {
        if (isFileEmpty(outputFilename)) {
            mediumOutputFile << "d = " << d << '\n';
            mediumOutputFile << "symClass = " << symClass << '\n';
            mediumOutputFile << "D = " << D << '\n';
        } else {
            //scan through and find the last value we're complete to
            std::ifstream mediumInFile;
            mediumInFile.open(outputFilename, std::ofstream::out | std::ofstream::app);
            string line;
            while (getline(mediumInFile, line)) {
                if (line.substr(0, 15) == "Complete up to ") {
                    mediumComplete = std::stod(line.substr(15));
                }
            }
            mediumInFile.close();
        }
    } else {
        std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
    }

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
    string line;
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
    inputFile.close();

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

    while (!intervals.empty()) {
        stack<pair<double, double>> spawnedIntervals;
        double upper = intervals.top().second;
        spawnedIntervals.push(intervals.top());
        intervals.pop();

        while (!spawnedIntervals.empty()) {
            auto interval = spawnedIntervals.top();
            spawnedIntervals.pop();

            bool zoomIn = true;
            double weylRightEndpoint = Od.eigenvalueIntervalRightEndpoint(interval.first, weylEigenvalueCount);
            weylRightEndpoint = min(weylRightEndpoint, interval.first + 0.001);
            weylRightEndpoint = max(weylRightEndpoint, interval.first + 0.00000001);
            if (interval.second <= weylRightEndpoint) {
                mediumOutputFile << std::setprecision(16) << "[" << interval.first << ", " << interval.second << "]" << std::endl;
            } else {
                double left = interval.first;
                double right = interval.second;
                double center = (left + right) / 2;

                KBessel* leftRK = new KBessel(twoPiOverA * 1 * Y0, left);
                KBessel* centerRK = new KBessel(twoPiOverA * 1 * Y0, center);
                KBessel* rightRK = new KBessel(twoPiOverA * 1 * Y0, right);
                bool leftInterval = possiblyContainsEigenvalue(left, center, leftRK, centerRK);
                bool rightInterval = possiblyContainsEigenvalue(center, right, centerRK, rightRK);
                delete leftRK;
                delete centerRK;
                delete rightRK;
                if (!leftInterval && !rightInterval) {
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
        mediumOutputFile << "Complete up to " << std::setprecision(16) << upper << std::endl;
    }
    mediumOutputFile.close();
}

void BianchiMaassSearch::fineSearchForEigenvalues() {

    const std::string directoryOut = "Output/Fine/";
    const std::string prefixOut = "fine_"
                                  + to_string(d) + "_"
                                  + symClass;
    createOutputDirectory(directoryOut);

    std::string outputFilename = directoryOut
                                 + prefixOut
                                 + ".txt";

    fineOutputFile.open(outputFilename, std::ofstream::out | std::ofstream::app);

    double fineComplete = 0.75;

    if (fineOutputFile.is_open()) {
        if (isFileEmpty(outputFilename)) {
            fineOutputFile << "d = " << d << '\n';
            fineOutputFile << "symClass = " << symClass << '\n';
        } else {
            //scan through and find the last value we're complete to
            std::ifstream fineInFile;
            fineInFile.open(outputFilename);
            string line;
            while (getline(fineInFile, line)) {
                if (line[0] == '[') {
                    //remove brackets
                    string numbers = line.substr(1, line.size() - 2);
                    int comma = numbers.find(',');
                    string leftStr = numbers.substr(0,comma);
                    string rightStr = numbers.substr(comma + 2);
                    double left = std::stod(leftStr);
                    double right = std::stod(rightStr);
                    if (fineComplete < right) {
                        fineComplete = right;
                    }
                }
            }
            fineInFile.close();
        }
    } else {
        std::cerr << "Error creating file \"" << outputFilename << "\"" << std::endl;
    }


    /*
     *
     * Read in the medium file and collect the intervals there
     *
     *
     */

    const std::string directoryIn = "Output/Medium/"; // Change this to the desired directory
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
    string line;
    while (std::getline(inputFile, line)) {
        if (line[0] == '[') {
            //remove brackets
            string numbers = line.substr(1, line.size() - 2);
            int comma = numbers.find(',');
            string leftStr = numbers.substr(0,comma);
            string rightStr = numbers.substr(comma + 2);
            double left = std::stod(leftStr);
            double right = std::stod(rightStr);
            if (fineComplete < left) {
                intervalsFromFile.emplace_back(left,right);
            }
        }
    }
    inputFile.close();

    stack<pair<double, double>> intervals;
    for (auto itr = intervalsFromFile.rbegin(); itr != intervalsFromFile.rend(); itr++) {
        intervals.push(*itr);
    }

    while(!intervals.empty()) {
        auto interval = intervals.top();
        intervals.pop();

        truncation = pow(10.0, -dToDMap[d][2]);
        auto firstPassOutput = fineSecantMethod(interval.first, interval.second, D);
        if (firstPassOutput.first == firstPassOutput.second) {
            fineOutputFile << "[" << std::setprecision(16) << interval.first << ", " << interval.second << "]" << std::endl;
            return;
        }

        truncation = pow(10.0, -dToDMap[d][3]);
        auto secondPassOutput = fineSecantMethod(firstPassOutput.first, firstPassOutput.second, finalPrecision);
        fineOutputFile << "[" << std::setprecision(16) << interval.first << ", " << interval.second << "]" << std::endl;
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
    if (proportion >=  0.5) {
        coarseOutputFile << setprecision(16) << "[" << leftR << ", " << rightR << "]" << std::endl;
        return true;
    }
    return false;
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
        std::filesystem::create_directory(directory);
    }
}

bool BianchiMaassSearch::isFileEmpty(const string &filename) {
    return std::filesystem::is_empty(filename);
}

pair<double, double> BianchiMaassSearch::fineSecantMethod(double leftR, double rightR, int secantD) {
    std::cout << "Starting final refinement." << std::endl;
    std::cout << std::setprecision(16) << "[" << leftR << ", " << rightR << "]" << std::endl;

    D = secantD;
    truncation = pow(10.0, -D);

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

    /*
     * Y1 and Y2 have been computed. Proceed with computing the actual matrices
     */
    double MY = computeMYGeneral(M0, Y);

    double maxYStar = 0;

    map<Index, vector<TestPointOrbitData>> mToYTestPointOrbits;

    mToYTestPointOrbits.clear();


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

    double thisR = rightR;
    double thisPhi = phi2;
    double nextR = findZeroOfLinearInterpolation(leftR, phi1, rightR, phi2);
    double nextPhi;

    std::cout << std::setprecision(16) << "r = " << leftR << std::endl;
    std::cout << std::setprecision(16) << "r = " << rightR << std::endl;
    std::cout << std::setprecision(16) << "r = " << nextR << std::endl;

    vector<pair<double, double>> eigenvalueHeckePairs;
    bool recomputeTruncation = false;

    int maxIterations = 30;
    int iterations = 0;
    bool convergence = true;

    KBessel K = KBessel(twoPiOverA * 1 * Y0, nextR);

    while (true) {
        iterations++;

        K.setRAndPrecompute(nextR, 2 * pi / A * M0 * maxYStar);
        matrixY = produceMatrix(indexTransversal, mToYTestPointOrbits, indexOrbitDataModSign, Y, K);
        auto solAndCondYNextR = solveMatrix(matrixY, indexTransversal, indexOfNormalization);

        std::cout << solAndCondYNextR.second << " " << solAndCondYNextR.second << std::endl;

        map<Index, double> coeffMapNextR;

        for (int i = 0; i < indexTransversal.size(); i++) {
            Index index = indexTransversal[i];
            coeffMapNextR[index] = solAndCondYNextR.first(i);
        }

        //We have the solutions from everything, now we cook up the phi function

        double coeff0NextR = coeffMapNextR.at(dToPrimes.at(d)[0]);
        double coeff1NextR = coeffMapNextR.at(dToPrimes.at(d)[1]);
        double coeff2NextR = coeffMapNextR.at(dToPrimes.at(d)[2]);

        double hecke1NextR = coeff0NextR * coeff1NextR - coeffMapNextR.at(prod01);
        double hecke2NextR = coeff0NextR * coeff2NextR - coeffMapNextR.at(prod02);
        hecke1NextR = abs(hecke1NextR);
        hecke2NextR = abs(hecke2NextR);

        nextPhi = eps1 * hecke1NextR + eps2 * hecke2NextR;

        if (thisPhi == nextPhi) { //End right now so we don't get arithmetic errors
            const std::string directory = "Output/Final/"; // Change this to the desired directory
            const std::string prefix = to_string(d) + "_"
                                       + symClass + "_";
            createOutputDirectory(directory);

            std::string outputFilename = directory
                                         + prefix
                                         + ".txt";

            ofstream out;
            out.open(outputFilename, std::ofstream::out | std::ofstream::app);

            out << "[" << std::setprecision(16) << leftR << ", " << rightR << "]" << std::endl;
            for (auto p : eigenvalueHeckePairs) {
                out << std::setprecision(16) << p.first << ", " << p.second << '\n';
            }

            return {min(thisR, nextR), max(thisR, nextR)};
        }

        double tempR = nextR;
        double tempPhi = nextPhi;

        nextR = findZeroOfLinearInterpolation(thisR, thisPhi, nextR, nextPhi);
        thisR = tempR;
        thisPhi = nextPhi;

        double hecke = heckeCheck(coeffMapNextR);
        std::cout << std::setprecision(16) << "r = " << nextR << std::endl;
        std::cout << "hecke " << std::setprecision(16) << hecke << std::endl;

        eigenvalueHeckePairs.emplace_back(nextR, hecke);

        if (hecke > 1) {
            std::cout << "hecke check failed " << std::setprecision(16) << hecke << std::endl;
            return {0,0};
        } else if (heckeHasConverged(eigenvalueHeckePairs)) {

            const std::string directory = "Output/Final/"; // Change this to the desired directory
            const std::string prefix = to_string(d) + "_"
                                       + symClass + "_";
            createOutputDirectory(directory);

            std::string outputFilename = directory
                                         + prefix
                                         + ".txt";

            ofstream out;
            out.open(outputFilename, std::ofstream::out | std::ofstream::app);

            out << "[" << std::setprecision(16) << leftR << ", " << rightR << "]" << std::endl;
            for (auto p : eigenvalueHeckePairs) {
                out << std::setprecision(16) << p.first << ", " << p.second << '\n';
            }

            return {min(thisR, nextR), max(thisR, nextR)};
        } else if (iterations > maxIterations) {

            const std::string directory = "Output/Final/"; // Change this to the desired directory
            const std::string prefix = to_string(d) + "_"
                                       + symClass + "_";
            createOutputDirectory(directory);

            std::string outputFilename = directory
                                         + prefix
                                         + ".txt";

            ofstream out;
            out.open(outputFilename, std::ofstream::out | std::ofstream::app);

            out << "[" << std::setprecision(16) << leftR << ", " << rightR << "]" << std::endl;
            out << "Did not converge as expected" << std::endl;
            for (auto p : eigenvalueHeckePairs) {
                out << std::setprecision(16) << p.first << ", " << p.second << '\n';
            }

            return {min(thisR, nextR), max(thisR, nextR)};
        }
    }
}

bool BianchiMaassSearch::heckeHasConverged(const vector<pair<double, double>> &eigenvalueHeckePairs) {
    int size = eigenvalueHeckePairs.size();
    if (size < 4) {
        return false;
    }

    vector<double> values;
    for (int i = size - 1; i >= size - 4; i--) {
        values.push_back(eigenvalueHeckePairs[i].second);
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


    Y = 0;
    maxSoFar = 0;
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

    vector<double> r1ConditionNumbers;
    vector<double> r2ConditionNumbers;

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


    ansY1 = 0;
    ansY2 = 0;
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

