//
// Created by Eric Moss on 6/27/23.
//

#include <cmath>
#include <stdexcept>
#include "CoefficientComputer.h"
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


#define watch(x) std::cout << (#x) << " is " << (x) << std::endl << std::flush

using Eigen::MatrixXd;

double CoefficientComputer::SUBINTERVAL_SEARCH_WIDTH = 0.5;
double CoefficientComputer::LOWER_BOUND_FOR_EIGENVALUE = 0.75;
double CoefficientComputer::MINIMUM_SPACING = 0.01;
double CoefficientComputer::SCALE_WEYL_LAW = 0.1;
double CoefficientComputer::BINARY_SEARCH_GRID_SPACING = 0.1;
double CoefficientComputer::RATIO_OF_Y1_TO_Y0 = 0.3;
double CoefficientComputer::RATIO_OF_Y2_TO_Y0 = 0.8;
double CoefficientComputer::EIGENVALUE_INTERVAL_WIDTH_CUTOFF = pow(10, -7);
double CoefficientComputer::CONDITION_CUTOFF = 200.0;

CoefficientComputer::CoefficientComputer(int d, int D, char symClass) :
matrixR1Y1(0,0), matrixR1Y2(0,0), matrixR2Y1(0,0), matrixR2Y2(0,0) {
    //I want the accuracy to indeed be 10^(-D), the truncation quantity of the fourier series
    //i.e. the size of the tail
    //I first compute the tail until the terms are [ +/- ...]
    //I need enough bits so that when I lose precision during computing fourier terms
    //the error doesn't start inter
    auto start = std::chrono::high_resolution_clock::now();

    this->d = d;
    this->D = D;
    this->symClass = symClass;

    outputFile.open("d" + std::to_string(d) + "_D" + std::to_string(D) + "_" + symClass + "_" + std::to_string(Y1Parameter) + ".txt");

    /*Check that d and symClass are correct*/
    std::vector<int> classOne = {1,2,3,7,11,19};
    bool dIsClassOne = std::find(classOne.begin(), classOne.end(), d) != classOne.end();
    assert(dIsClassOne);

    std::vector<char> symClasses = {'D','G','C','H'};
    bool symClassIsCorrect = std::find(symClasses.begin(), symClasses.end(), symClass) != symClasses.end();
    assert(symClassIsCorrect);

    /*Initialize values that are derived from d, D, and symClass*/
    //volumeOfFD = EricNT::computeVolumeOfFD(d);
    computeA();
    computeTheta();
    computeY0(); //Y0 depends on theta
    truncation = pow(10,-D);
    tolerance = pow(10,-(D+6));

    /*double Y = 0.1;
    while (Y > 0.05) {
        YGrid.push_back(Y);
        Y *= 0.99;
    }
    std::reverse(YGrid.begin(), YGrid.end());*/
    YGrid.push_back(0.07);
    YGrid.push_back(0.065);
    YGrid.push_back(0.06);
    YGrid.push_back(0.055);
    YGrid.push_back(0.05);
    YGrid.push_back(0.045);
    YGrid.push_back(0.04);
    std::reverse(YGrid.begin(), YGrid.end());

    K = new KBesselApproximator(0);
    T = new TrigApproximator();

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "A: " << A << std::endl;
    std::cout << "theta: " << theta << std::endl;
    std::cout << "Y0: " << Y0 << std::endl;
    std::cout << "truncation: " << truncation << std::endl;

    matrices.resize(4);
    conditionNumbers.resize(4);
    coefficientMaps.resize(4);
}

/**
 * @brief Sets the private member A to teh volume of the fundamental parallelogram for the specified number field.
 */
void CoefficientComputer::computeA() {
    if (EricNT::mod(-d,4) == 1) {
        A = sqrt(d)/2;
    } else if (EricNT::mod(-d,4) == 0) {
        throw std::invalid_argument("d is incorrect");
    } else {
        A = sqrt(d);
    }
}

/**
 * @brief Sets the private member theta to the generator for the specified ring of integers.
 */
void CoefficientComputer::computeTheta() {
    if (EricNT::mod(-d,4) == 1) {
        //theta = 1/2 + I*sqrt(d)/2
        theta = std::complex<double> {1.0/2, sqrt(d)/2};
    } else if (EricNT::mod(-d,4) == 0) {
        throw std::invalid_argument("d is incorrect");
    } else {
        //theta = I*sqrt(d)
        theta = std::complex<double> {0, sqrt(d)};
    }
}

/**
 * @brief Sets the private member Y0 to the lowest point of the fundamental domain for the specified ring of integers.
 */
void CoefficientComputer::computeY0() {
    if (d == 3) {
        Y0 = sqrt(2.0/3);
    } else if (d == 19) {
        Y0 = sqrt(2.0/19);
    } else if (d == 43) {
        Y0 = sqrt(2.0/43);
    } else if (d == 67) {
        Y0 = sqrt(2.0/67);
    } else if (d == 163) {
        Y0 = sqrt(2.0/163);
    } else {
        //Y0 = sqrt(1- (1/2)^2 - (theta.imag()/2)^2)
        //=sqrt(3/4 - theta.imag()^2/4)
        Y0 = .75 - pow(theta.imag(),2)/4;
        Y0 = sqrt(Y0);
    }
}


/*void CoefficientComputer::computeM(const std::string& inY) {
    auto start = std::chrono::high_resolution_clock::now();
    double localY;
    if (inY == "Y") {
        localY = Y;
    } else if (inY == "Y0") {
        localY = Y0;
    } else {
        throw std::invalid_argument("computeM string should be Y or Y0");
    }

    double modestImpliedConstant = 1.0;

    *//*Binary search to find where the terms become zero to 16 digits.*//*
    //startN is where BesselK starts to be exponential
    //startN = r/(2*pi/A*localY) = r*A/(2*pi*localY)
    //maxN is a random guess, will adjust later
    //K.updateRAndPrecompute(r);

    double startN = r;
    double maxN = 100;

    double bess = kappaExact(2*pi/A*startN*localY);
    double temp = modestImpliedConstant * localY * bess;

    while (temp < truncation) {
        startN /= 2;
        bess = kappaExact(2*pi/A*startN*localY);
        temp = modestImpliedConstant * localY * bess;
    }

    //double temp = modestImpliedConstant * right * localY * kappaExact(2*pi/A*right*localY);
    bess = kappaExact(2*pi/A*maxN*localY);
    temp = modestImpliedConstant * localY * bess;

    while (temp > tolerance) {
        maxN *= 2;

        //temp = modestImpliedConstant * right * localY * kappaExact(2*pi/A*right*localY);
        bess = kappaExact(2*pi/A*maxN*localY);
        temp = modestImpliedConstant * localY * bess;
    }

    *//*std::cout << "binary search starting right point ended with term evaluation = ";
    arb_printn(temp, bits, 0);
    printf("\n");*//*

    //Now the terms evaluated at "right" are 0
    //Begin the binary search to find where the
    double left = startN;
    double right = maxN;
    double center;

    //while (right - left).abs() > 1/2
    while (right - left > 0.001) {
        center = (left + right)/2;
        //temp = modestImpliedConstant * center * localY * kappaExact(2*pi/A*center*localY);
        bess = kappaExact(2*pi/A*center*localY);
        temp = modestImpliedConstant * localY * bess;

        if (temp > tolerance) {
            left = center;
        } else {
            //temp < tolerance
            right = center;
        }
    }

    double toleranceCutoff = right;

    //////start over with truncation bound finding

    left = startN;
    right = maxN;

    //double temp = modestImpliedConstant * right * localY * kappaExact(2*pi/A*right*localY);
    bess = kappaExact(2*pi/A*right*localY);
    temp = modestImpliedConstant * localY * bess;

    *//*std::cout << "binary search starting right point ended with term evaluation = ";
    arb_printn(temp, bits, 0);
    printf("\n");*//*

    //Now the terms evaluated at "right" are 0
    //Begin the binary search to find where the

    //while (right - left).abs() > 1/2
    while (right - left > 0.001) {
        center = (left + right)/2;
        //temp = modestImpliedConstant * center * localY * kappaExact(2*pi/A*center*localY);
        bess = kappaExact(2*pi/A*center*localY);
        temp = modestImpliedConstant * localY * bess;

        if (temp > truncation) {
            left = center;
        } else {
            //temp < tolerance
            right = center;
        }
    }

    double truncationCutoff = right;
    *//*std::cout << "binary search ended with term evaluation =  ";
    arb_printn(temp, bits, 0);
    printf("\n");*//*
    *//*include some code to verify this answer is correct?*//*


    *//*Compute the actual tail of the series by computing all the indicesM0 up to our found bound.*//*
    *//*Start by generating the indicesM0 that go into the tail.*//*

    //aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));

    int aUpperBound = ceil(toleranceCutoff * abs(theta)/ sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int aLowerBound = -aUpperBound;

    //bub = ceil(maxN / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));
    int bUpperBound = ceil(toleranceCutoff /sqrt(pow(abs(theta),2) - pow(theta.real(),2)));
    int bLowerBound = -bUpperBound;

    std::vector<Index> tailIndices;

    for (int a = aLowerBound; a <= aUpperBound; a++) {
        for (int b = bLowerBound; b <= bUpperBound; b++) {
            Index index = Index(a,b,d);
            //std::cout << index << std::endl;
            //arb_printn(index.angle, bits, 0);
            //printf("\n");
            //arb_printn(index.abs, bits, 0);
            //printf("\n");
            if (index.getAbs() == 0 || index.getAbs() < truncationCutoff || index.getAbs() > toleranceCutoff) {
                continue;
            }
            tailIndices.push_back(index);
        }
    }

    //Sort indicesM0 by absolute value then positive angle w.r.t. pos real axis
    std::sort(tailIndices.begin(), tailIndices.end());

    *//*Start computing the tail.*//*
    double tail = 0;
    double term;

    for (const auto& itr : tailIndices) {
        bess = kappaExact(2*pi/A*itr.getAbs()*localY);
        term = modestImpliedConstant * localY * bess;
        tail += term;
    }

    double M;
    //Subtract off the early terms until the tail is below the accuracy bound
    for (const auto &itr : tailIndices) {
        //tail < accuracy
        if (tail < truncation) {
            M = itr.getAbs();
            break;
        } else {
            bess = kappaExact(2*pi/A*itr.getAbs()*localY);
            term = modestImpliedConstant * localY * bess;
            tail -= term;
        }
    }

    if (inY == "Y") {
        MY = M;
        //std::cout << "MY: " << MY << std::endl;
    } else if (inY == "Y0") {
        M0 = M;
        //std::cout << "M0: " << M0 << std::endl;
    } else {
        throw std::invalid_argument("computeM string should be Y or Y0");
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    computeMTime += duration.count();
}*/

void CoefficientComputer::computeY1() {
    //Computed so that kappa(2*pi/A*M0*Y1) is half of the max of kappa(...)
    //Then Y2 can be slightly larger and it's still okay... I Hope
    std::vector<double> point = K->maximize();
    double minY1 = point[0]/(2*pi/A*M0);
    double maxY1 = Y0;
    double bessMax = point[1];

    double left = minY1;
    double right = maxY1;
    double center = (left + right)/2.0;

    double leftValue = kappaExact(2*pi/A*M0*left) - bessMax*0.4;
    double centerValue = kappaExact(2*pi/A*M0*center) - bessMax*0.4;
    while (right - left > 0.00001) {
        if (areDifferentSign(leftValue, centerValue)) {
            //zero is in left interval
            right = center;
        } else {
            //zero is in right interval
            left = center;
            leftValue = centerValue;
        }
        center = (left + right)/2.0;
        centerValue = kappaExact(2*pi/A*M0*center) - bessMax*0.4;
    }
    Y1 = center;
    watch(Y1);
}

void CoefficientComputer::computeY2() {
    std::vector<double> point = K->maximize();
    double minY2 = point[0]/(2*pi/A*M0);
    double maxY2 = Y0;
    double bessMax = point[1];

    double left = minY2;
    double right = maxY2;
    double center = (left + right)/2.0;

    double leftValue = kappaExact(2*pi/A*M0*left) - bessMax*0.3;
    double centerValue = kappaExact(2*pi/A*M0*center) - bessMax*0.3;
    while (right - left > 0.00001) {
        if (areDifferentSign(leftValue, centerValue)) {
            //zero is in left interval
            right = center;
        } else {
            //zero is in right interval
            left = center;
            leftValue = centerValue;
        }
        center = (left + right)/2.0;
        centerValue = kappaExact(2*pi/A*M0*center) - bessMax*0.3;
    }
    Y2 = center;
    watch(Y2);
}

/**
 * @brief Computes indices and sets vectors up to bounds M0 and MY. Also computes symmetry relations between
 * coefficients.
 */
void CoefficientComputer::computeIndexData() {
    indicesM0.clear();
    indicesMY.clear();
    indexTransversalMY.clear();
    indexTransversal.clear();
    indexOrbitData.clear();
    indexOrbitDataModMinusOne.clear();

    //aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));

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
            Index index = Index(a, b, d);
            if (index.getAbs() <= MY) {
                indicesMY.push_back(index);
                if (index.getAbs() <= M0) {
                    indicesM0.push_back(index);
                }
            }
        }
    }

    /*Sort the indicesM0 by absolute value then by angle with positive real axis.*/
    std::sort(indicesMY.begin(), indicesMY.end());

    std::vector<Index> tempIndices;
    copy(indicesMY.begin(), indicesMY.end(), back_inserter(tempIndices));

    int rotationCoeff = symClass == 'D' || symClass == 'G' ? 1 : -1;
    int conjCoeff = symClass == 'D' || symClass == 'C' ? 1 : -1;

    while (!tempIndices.empty()) {
        Index index = tempIndices.front();
        std::vector<std::tuple<Index, int>> orbit = {std::tuple<Index, int>(index, 1)};

        Index tempIndex = index.rotate();
        int tempCoeff = 1*rotationCoeff;
        std::tuple<Index, int> tempTuple = {tempIndex, tempCoeff};
        while (tempIndex != index) {
            orbit.push_back(tempTuple);

            tempIndex = tempIndex.rotate();
            tempCoeff = tempCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
        }

        Index conjIndex = index.conj();
        bool conjIsInRotations = false;
        for (auto tup : orbit) {
            if (conjIndex == std::get<0>(tup)) {
                conjIsInRotations = true;
                break;
            }
        }

        if (!conjIsInRotations) {
            /*add in the rotations of the conjugate*/
            tempTuple = {conjIndex, conjCoeff};
            orbit.push_back(tempTuple);

            tempIndex = conjIndex.rotate();
            tempCoeff = conjCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
            while (tempIndex != conjIndex) {
                orbit.push_back(tempTuple);

                tempIndex = tempIndex.rotate();
                tempCoeff = tempCoeff*rotationCoeff;
                tempTuple = {tempIndex, tempCoeff};
            }

        }
        /*Save our starting index to the transversal.*/
        if (index.getAbs() <= M0) {
            indexTransversal.push_back(index);
            indexTransversalMY.push_back(index);
        } else {
            indexTransversalMY.push_back(index);
        }


        /*Save the orbit data.*/
        indexOrbitData[index] = orbit;

        /*Compute the orbit mod +-1*/
        std::vector<std::tuple<Index,int>> orbitModMinusOne;
        std::vector<Index> alreadyGotten;
        for (auto tup : orbit) {
            Index l = std::get<0>(tup);
            int pmOne = std::get<1>(tup);
            if (std::find(alreadyGotten.begin(), alreadyGotten.end(),l) == alreadyGotten.end()) {
                //add it
                std::tuple<Index,int> classModMinusOne = {l,pmOne};
                alreadyGotten.push_back(l);
                alreadyGotten.push_back(Index(-l.getA(), -l.getB(), l.getD()));
                orbitModMinusOne.push_back(classModMinusOne);
            }
        }
        indexOrbitDataModMinusOne[index] = orbitModMinusOne;

        /*Delete the indicesM0 in the orbit from the copied list of indicesM0.*/
        for (auto itr1 : orbit) {
            Index toDelete = std::get<0>(itr1);
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
        if (indexTransversal[itr].getComplex() == std::complex<double> {1,0}) {
            indexOfNormalization = itr;
            break;
        }
    }

    //watch(indexTransversal.size());

    /*PlotWindow window = PlotWindow(-MY, MY, -MY, MY, 4000);
    auto origin = window.getOrigin();
    auto domainPoints = window.getPoints();
    FunctionToEvaluate func = FunctionToEvaluate();

    auto evaluatedPoints = func.indexPoints(domainPoints, indexTransversal);
    std::string filename = "indexTransversal.png";
    Plotter plotter = Plotter(evaluatedPoints);
    plotter.drawAxes(origin[0], origin[1]);
    plotter.histogramSmoothing = false;
    plotter.makeImage(filename, "greyscale");

    evaluatedPoints = func.indexPoints(domainPoints, indicesM0);
    filename = "indicesM0.png";
    plotter = Plotter(evaluatedPoints);
    plotter.drawAxes(origin[0], origin[1]);
    plotter.histogramSmoothing = false;
    plotter.makeImage(filename, "greyscale");

    evaluatedPoints = func.indexPoints(domainPoints, indicesMY);
    filename = "indicesMY.png";
    plotter = Plotter(evaluatedPoints);
    plotter.drawAxes(origin[0], origin[1]);
    plotter.histogramSmoothing = false;
    plotter.makeImage(filename, "greyscale");*/
}

/**
 * @brief Computes test points on two horospheres specified by Y1 and Y2; saves the points and pullbacks to a private map.
 */
void CoefficientComputer::computeTestPointData() {
    auto start = std::chrono::high_resolution_clock::now();
    mToTestPointPairsY1.clear();
    mToTestPointPairsY2.clear();
    maxYStar = 0;

    bool postCompute = false;

    for (auto m : indexTransversalMY) {
        if (postCompute == false && m.getAbs() > M0) {
            continue;
        }

        int Q0 = 0;
        int Q1 = (MY + abs(m.getComplex().real()))/2.0;
        if (EricNT::mod(-d,4) == 1) {
            Q0 = std::max(MY + abs(m.getComplex().real()), (MY + abs(m.getComplex().imag()))/(2.0*A));
        } else {
            Q0 = (MY + abs(m.getComplex().imag()))/(2*A);
        }
        Q0 += 2;
        Q1 += 2;

        int totalNumberOfPoints = (2*Q0)*(2*Q1)/2;

         std::vector<Quaternion> pointsY1;
         std::vector<Quaternion> pointsY2;
         pointsY1.reserve(totalNumberOfPoints);
         pointsY2.reserve(totalNumberOfPoints);
         //Quaternion point = {-1.0L/2 + 1.0L/(2*P) + 1.0L*l0/P, -imagTheta / 2.0L + imagTheta / (2 * Q) + imagTheta * l1 / Q, Y, 0};
         //= {(-P + 1)/(2*P) + l0/P, (-imagTheta*Q + imagTheta)/(2*Q) + imagTheta*l1/Q, Y, 0}

         for (int l0 = 1-Q0; l0 <= Q0; l0++) {
             for (int l1 = 1; l1 <= Q1; l1++) {
                 std::complex<double> x = ((double)(l0-0.5))/(2.0*Q0) + theta*((double)(l1-0.5))/((double)2.0*Q1);
                 double x0 = x.real();
                 double x1 = x.imag();

                 Quaternion pointY1 = {x0, x1, Y1, 0};
                 pointsY1.push_back(pointY1);
                 Quaternion pointY2 = {x0, x1, Y2, 0};
                 pointsY2.push_back(pointY2);
             }
         }

         std::vector<std::vector<Quaternion>> pointAndPullbackPairsY1;
         std::vector<std::vector<Quaternion>> pointAndPullbackPairsY2;
         pointAndPullbackPairsY1.resize(totalNumberOfPoints);
         pointAndPullbackPairsY2.resize(totalNumberOfPoints);

#pragma omp parallel for
         for (int i = 0; i < pointsY1.size(); i++) {
             Quaternion pointY1 = pointsY1[i];
             Quaternion pointY2 = pointsY2[i];
             Quaternion pullbackY1 = Quaternion(pointY1);
             Quaternion pullbackY2 = Quaternion(pointY2);
             pullbackY1.reduce(d);
             pullbackY2.reduce(d);
             assert(pullbackY1.getJ() > Y0);
             assert(pullbackY2.getJ() > Y0);
             std::vector<Quaternion> pairY1 = {pointY1,pullbackY1};
             std::vector<Quaternion> pairY2 = {pointY2,pullbackY2};
             pointAndPullbackPairsY1[i] = pairY1;
             pointAndPullbackPairsY2[i] = pairY2;
         }

         for (int i = 0; i < pointsY1.size(); i++) {
             Quaternion pullbackY1 = pointAndPullbackPairsY1[i][1];
             Quaternion pullbackY2 = pointAndPullbackPairsY2[i][1];
             if (pullbackY1.getJ() > maxYStar) {
                 maxYStar = pullbackY1.getJ();
             }
             if (pullbackY2.getJ() > maxYStar) {
                 maxYStar = pullbackY2.getJ();
             }
         }

         mToTestPointPairsY1[m] = pointAndPullbackPairsY1;
         mToTestPointPairsY2[m] = pointAndPullbackPairsY2;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    computeTestPointDataTime += duration.count();
}

/**
 * @brief Computes the Hejhal matrix with Y and r specified via matrixID. The 4 possible matrices are independently stored.
 *
 * Also computes the condition number for each matrix computed.
 *
 * @param matrixID Y1R1, Y1R2, Y2R1, or Y2R2
 */
void CoefficientComputer::populateMatrix(short matrixID) {

    auto tempMatrix = matrices[matrixID];
    tempMatrix.resize(indexTransversal.size(), indexTransversal.size());

#pragma omp parallel for collapse(2)
    for (int i = 0; i < indexTransversal.size(); i++) {
        for (int j = 0; j < indexTransversal.size(); j++) {
            Index m = indexTransversal[i];
            Index n = indexTransversal[j];
            double entry = computeEntry(matrixID, m, n);
            //std::cout << entry << std::endl;

            tempMatrix(i,j) = entry;
        }
    }

    Eigen::MatrixXd subMatrix;
    subMatrix.resize(indexTransversal.size() - 1, indexTransversal.size() - 1);
    for (int i = 0; i < indexTransversal.size(); i++) {
        for (int j = 0; j < indexTransversal.size(); j++) {
            if (i == indexOfNormalization || j == indexOfNormalization) {
                continue;
            }
            int modifiedI;
            int modifiedJ;
            if (i < indexOfNormalization) {
                modifiedI = i;
            } else {
                modifiedI = i - 1;
            }

            if (j < indexOfNormalization) {
                modifiedJ = j;
            } else {
                modifiedJ = j - 1;
            }

            subMatrix(modifiedI, modifiedJ) = tempMatrix(i, j);
        }
    }
    std::cout << "-------------" << std::endl;
    std::cout << subMatrix << std::endl;
    std::cout << "-------------" << std::endl;
    Eigen::JacobiSVD<MatrixXd> svd(subMatrix);
    conditionNumbers[matrixID] = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);

    matrices[matrixID] = tempMatrix;


   /* std::cout << "-----------------------------------" << std::endl;
    std::cout << matrix << std::endl;
    std::cout << "-----------------------------------" << std::endl;*/


    /*if (cond * truncation > 0.001) {
        double offDiag = 0;
        for (int i = 0; i < matrix.rows(); i++) {
            for (int j = 0; j < matrix.cols(); j++) {
                if (i == j) {
                    std::cout << matrix(i,j) << " ";
                } else {
                    offDiag += abs(matrix(i,j));
                }
            }
        }
        std::cout << std::endl << "off diagonal average: " << offDiag/(matrix.rows()*matrix.cols() - matrix.rows()) << std::endl;
    }*/

}

/*
 * m is the row
 * n is the column
 */
/*std::complex<double> CoefficientComputer::computeEntry(const Index &m, const Index &n) {
    std::complex<double> answer = {0,0};

    std::vector<std::vector<Quaternion>> pointAndPullbackPairs;
    double Y;
    if (whichYToUse == 1) {
        pointAndPullbackPairs = testPointDataY1[m];
        Y = Y1;
    } else {
        pointAndPullbackPairs = testPointDataY2[m];
        Y = Y2;
    }

    auto numberOfPoints = pointAndPullbackPairs.size();

    for (auto pair : pointAndPullbackPairs) {
        Quaternion point = pair[0];
        Quaternion pullback = pair[1];

        std::complex<double> x = point.getComplex();
        std::complex<double> xStar = pullback.getComplex();

        double y = point.getJ();
        double yStar = pullback.getJ();

        //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
        std::complex<double> term = yStar * kappaApprox(2*pi/A*n.getAbs()*yStar);

        //e(symClass, n, x^*)
        //sum_{l \in S_n} +/- 1 * exp(pi I /A * <I*l, xstar>)
        std::complex<double> symExpFactor = {0,0};
        for (auto data : indexOrbitData[n]) {
            Index l = std::get<0>(data);
            int pmOne = std::get<1>(data);
            symExpFactor += (double)pmOne * std::exp(pi*I/A * traceProduct(I*l.getComplex(), xStar));
        }
       term *= symExpFactor;

        //exp(-pi*I/A * <I*m,x>)
        std::complex<double> expFactor = std::exp(-pi*I/A* traceProduct(I*m.getComplex(),x));
        term *= expFactor;

        answer += term;
    }
    //Right now answer = sum_(points) Y^* kappa(2*pi/A*|n|Y^*)e(symClass, n, x^*) exp(-pi *I/A *<I*m, x>)
    //multiply answer by -im(theta)/#pts
    answer *= -theta.imag()/numberOfPoints;

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm;
        if (whichYToUse == 1) {
            deltaTerm = Y1 * kappaApprox(2*pi/A *m.getAbs() * Y1);
        } else {
            deltaTerm = Y2 * kappaApprox(2*pi/A *m.getAbs() * Y2);
        }

        answer += deltaTerm;
    }
    return answer;
}*/

/**
 * @brief Computes matrix entry (m,n) with Y and r specified by matrixID. Requires points to be precomputed
 * @param matrixID
 * @param m
 * @param n
 * @return
 */
double CoefficientComputer::computeEntry(short matrixID, const Index &m, const Index &n) {
    double answer = 0;

    std::vector<std::vector<Quaternion>> pointAndPullbackPairs;
    double localY;
    if (matrixID == Y1R1 || matrixID == Y1R2) {
        pointAndPullbackPairs = mToTestPointPairsY1[m];
        localY = Y1;
    } else {
        pointAndPullbackPairs = mToTestPointPairsY2[m];
        localY = Y2;
    }

    for (auto pair : pointAndPullbackPairs) {
        auto point = pair[0];
        auto pullback = pair[1];

        std::complex<double> x = point.getComplex();
        std::complex<double> xStar = pullback.getComplex();

        double yStar = pullback.getJ();
        assert(yStar > Y0);

        //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
        double term = yStar * kappaApprox(2*pi/A*n.getAbs()*yStar);

        if (symClass== 'D' || symClass == 'G' || d == 1) {
            term *= std::cos(pi/A * traceProduct(-I*m.getComplex(),x));
        } else {
            term *= -std::sin(pi/A * traceProduct(-I*m.getComplex(),x));
        }


        //e(symClass, n, x^*)
        //sum_{l \in S_n} +/- 1 * exp(pi I /A * <I*l, xstar>)
        double cs = 0;
        if (symClass == 'D' || symClass == 'G' || d == 1) {
            for (auto tup : indexOrbitDataModMinusOne[n]) {
                Index l = std::get<0>(tup);
                double ellTerm = std::get<1>(tup);
                ellTerm *= std::cos(pi/A * traceProduct(I*l.getComplex(),xStar));
                cs += ellTerm;
            }
        } else {
            for (auto tup : indexOrbitDataModMinusOne[n]) {
                Index l = std::get<0>(tup);
                double ellTerm = std::get<1>(tup);
                ellTerm *= std::sin(pi/A * traceProduct(I*l.getComplex(),xStar));
                cs += ellTerm;
            }
        }
        term *= cs;

        answer += term;
    }

    //Right now answer = sum_(points) Y^* kappa(2*pi/A*|n|Y^*)e(symClass, n, x^*) exp(-pi *I/A *<I*m, x>)
    answer *= -4.0/(2.0*pointAndPullbackPairs.size());

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm = localY * kappaExact(2*pi/A *m.getAbs() * localY);

        answer += deltaTerm;
    }
    return answer;
}

void CoefficientComputer::solveMatrixAndComputeCoefficientMap(short matrixID) {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> x;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> b;
    Eigen::MatrixXcd subMatrix;

    auto tempCoefficientMap = coefficientMaps[matrixID];
    tempCoefficientMap.clear();
    long size = indexTransversal.size() - 1;
    x.resize(size);
    b.resize(size);
    subMatrix.resize(size, size);

    //set b to negative of the indexOfNormalization column
    //set tempA to be matrixA deleting row and column indexOfNormalization
    auto tempMatrix = matrices[matrixID];
    for (int row  = 0; row < indexTransversal.size(); row++) {
        for (int column = 0; column < indexTransversal.size(); column++) {
            if (column == indexOfNormalization && row != indexOfNormalization) {
                //add it to b
                if (row > indexOfNormalization) {
                    b(row - 1) = tempMatrix(row, column);
                } else {
                    b(row) = tempMatrix(row, column);
                }
            } else if (column != indexOfNormalization && row != indexOfNormalization) {
                //add it to tempA somehow
                int modifiedRow = (row > indexOfNormalization) ? row - 1 : row;
                int modifiedColumn = (column > indexOfNormalization) ? column - 1 : column;
                subMatrix(modifiedRow, modifiedColumn) = tempMatrix(row, column);
            }
        }
    }


    //make b negative
    b = -b;

    //solve the equation tempA * x = b for x
    x = subMatrix.colPivHouseholderQr().solve(b);

    /*std::cout << "computed coefficients: " << std::endl;
    std::cout << x << std::endl;
    std::cout << "--------------" << std::endl;*/

    //put the coefficients in the coefficient map, remembering to add back in the coefficients for indexOfNormalization
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index index = indexTransversal[i];
        if (i == indexOfNormalization) {
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = std::get<1>(tup);
                tempCoefficientMap[std::get<0>(tup)] = coeff;
            }
        } else {
            //look through the orbit and assign the same value to all of them using the +/- factor
            int modifiedi = (i > indexOfNormalization) ? i - 1 : i;
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = x(modifiedi).real();
                coeff *= std::get<1>(tup);
                tempCoefficientMap[std::get<0>(tup)] = coeff;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    double coeffAverage = 0;
    for (auto itr : tempCoefficientMap) {
        coeffAverage += std::abs(itr.second);
    }
    //std::cout << "coeff average " << matrixID << ": " << coeffAverage/tempCoefficientMap.size() << std::endl;
    coefficientMaps[matrixID] = tempCoefficientMap;
}

double CoefficientComputer::evaluate(const Quaternion &z) {
    //check that complex is above Y
    assert(!coefficientMaps[Y1R1].empty());
    auto coeffMap = coefficientMaps[Y1R1];

    double answer = 0;

    for (auto itr : coeffMap) {
        auto n = itr.first;
        auto a_n = itr.second;

        double term = a_n*z.getJ();

        //kappa(2*pi/A*|n|*y)
        term *= kappaApprox(2*pi/A * n.getAbs() * z.getJ());

        //exp(pi*I/A <I*n, x>)
        if (symClass == 'D' || symClass == 'G' || d == 1) {
            term *= 2.0 * std::cos(pi/A * traceProduct(I*n.getComplex(), z.getComplex()));
        } else {
            term *= 2.0 * std::sin(pi/A * traceProduct(I*n.getComplex(), z.getComplex()));
        }

        answer += term;
    }
    return answer;
}

/*void CoefficientComputer::searchForEigenvalues(double leftR, double rightR) {
    auto start = std::chrono::high_resolution_clock::now();



    //TODO: compute "optimal" Y1 and Y2 based on where the left and right bessel functions decay



    //translate these into Y-intervals


    //Initialize the rest of the data that depends on the eigenvalue
    //We use the greatest eigenvalue in the interval to compute M0, and consequently the indicesM0 of
    // our truncation, which must remain constant throughout.
    //M(Y) (and thus test points) is/are recomputed with each new eigenvalue in order to reduce total
    // number of K-bessel computations.
    this->leftR = leftR;
    this->rightR = rightR;

    computeM0();
    computeIndexTransversal();
    computeY1Y2();
    computeMY();
    computeIndexData();
    computeTestPointData();


    assert (leftR < rightR);

    //split up interval into chunks (how to decide? might have to guess and adjust.)
    auto left = (leftR < LOWER_BOUND_FOR_EIGENVALUE) ? LOWER_BOUND_FOR_EIGENVALUE : leftR;
    auto right = (rightR < LOWER_BOUND_FOR_EIGENVALUE) ? LOWER_BOUND_FOR_EIGENVALUE : rightR;
    assert (left < right);

    std::vector<double> rGrid;

    double itr = left;
    while (itr < right) {
        rGrid.push_back(itr);
        itr += BINARY_SEARCH_GRID_SPACING;
    }
    rGrid.push_back(rGrid[rGrid.size() - 1] + BINARY_SEARCH_GRID_SPACING);
    for (auto itr : rGrid) {
        std::cout << itr << " ";
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    std::cout << "searchForEigenvalues setup: " << duration.count() << " seconds." << std::endl;

    //iterate over chunks (subintervals), where for each chunk you recompute
    //--M(Y) and M0
    //--Indices
    //--Test points
    //Call function which executes a search in the subinterval and does not have to worry about recomputing anything.
    std::vector<double> leftZeros;
    std::vector<double> rightZeros = setNewRAndComputeZeroesVector(rGrid[0]);

    for (int i = 1; i < rGrid.size(); i++) {
        std::cout << "Starting [" << rGrid[i-1] << ", " << rGrid[i] << "]...";
        leftZeros = rightZeros;
        rightZeros = setNewRAndComputeZeroesVector(rGrid[i]);
        int signChanges = countSignChanges(leftZeros, rightZeros);
        std::cout << signChanges << "/" << leftZeros.size() << " sign changes." << std::endl << std::flush;
        if (signChanges >= leftZeros.size()/3) {
            //Eigenvalue likely!
            std::cout << "Above interval might have an eigenvalue!" << std::endl;
            auto leftEndpoint = rGrid[i - 1];
            auto rightEndpoint = rGrid[i];
            auto subLeftZeros = leftZeros;
            auto subRightZeros = rightZeros;

            bool exitedEarly = false;
            while (rightEndpoint - leftEndpoint > 0.000001) {

                *//*std::vector<double> centerZeros;
                for (int j = 0; j < subLeftZeros.size(); j++) {
                    double x0 = leftEndpoint;
                    double x1 = rightEndpoint;
                    double y0 = subLeftZeros[j];
                    double y1 = subRightZeros[j];

                    if (areDifferentSign(y0, y1)) {
                        centerZeros.push_back(findZeroOfLinearInterpolation(x0, y0, x1, y1));
                    }
                }
                std::sort(centerZeros.begin(), centerZeros.end());
                double centerEndpoint = centerZeros[floor(centerZeros.size()/2)];*//*
                double centerEndpoint = (leftEndpoint + rightEndpoint)/2.0;

                watch(leftEndpoint);
                watch(centerEndpoint);
                watch(rightEndpoint);

                auto subCenterZeros = setNewRAndComputeZeroesVector(centerEndpoint);
                int leftSignChanges = countSignChanges(subLeftZeros, subCenterZeros);
                int rightSignChanges = countSignChanges(subCenterZeros, subRightZeros);

                watch(leftSignChanges);
                watch(rightSignChanges);

                if (std::max(leftSignChanges, rightSignChanges) < leftZeros.size()/3) {
                    std::cout << "number of sign changes went too low, eigenvalue unlikely" << std::endl;
                    exitedEarly = true;
                    break;
                }

                if (leftSignChanges < rightSignChanges) {
                    subLeftZeros = subCenterZeros;
                    leftEndpoint = centerEndpoint;
                    *//*if (rightSignChanges < leftZeros.size()/2.0) {
                        std::cout << "Stopping because the number of sign changes got too small." << std::endl << std::flush;
                        break;
                    }*//*
                } else if (leftSignChanges > rightSignChanges) {
                    subRightZeros = subCenterZeros;
                    rightEndpoint = centerEndpoint;
                    *//*if (leftSignChanges < leftZeros.size()/2.0) {
                        std::cout << "Stopping because the number of sign changes got too small." << std::endl << std::flush;
                        break;
                    }*//*
                } else {
                    //equal number of sign changes
                    std::cout << "equal number of sign changes; inconclusive" << std::endl;
                    exitedEarly = true;
                    break;
                }
                watch(heckeCheck());

            }
            if (exitedEarly) {
                continue;
            }
            outputFile << std::setprecision(16) << leftEndpoint << " " << rightEndpoint;
            outputFile << " Hecke check: " << heckeCheck() << std::endl;
            std::cout << "--------------------------------\n";
            std::cout << std::flush;
        }
    }
}*/


void CoefficientComputer::recursiveSearchForEigenvalues(double localLeftR, double localRightR) {
    //set up all the parameters
    //Choose Y1 and Y2
    //If one of the condition numbers is too high,
    //recompute Y1 and Y2

    double spacing = 1.0/8;
    if (localRightR - localLeftR < spacing) {
        std::vector<double> endpoints;
        double endpoint = localLeftR;
        while (endpoint < localRightR) {
            endpoints.push_back(endpoint);
            endpoint += spacing;
        }
        for (int i = 0; i < endpoints.size() - 1; i++) {
            recursiveSearchForEigenvalues(endpoints[i], endpoints[i+1]);
        }
    }

    double stop = pow(10,-16);
    if (localRightR - localLeftR < stop) {
        std::cout << "Search finished to accuracy of " << stop << std::endl;
        std::cout << heckeCheck();
        std::cout << "-----------------" << std::endl;

        outputFile << "Search finished to accuracy of " << stop << std::endl;
        outputFile << std::setprecision(16) << "[" << localLeftR << ", " << localRightR << "]" << std::endl;
        outputFile << std::setprecision(16) <<  heckeCheck() << std::endl;
        outputFile << "--------------------------" << std::endl;
        return;
    }

    this->leftR = localLeftR;
    this->rightR = localRightR;
    watch(localLeftR);
    watch(localRightR);

    //general setup
    computeM0();
    computeIndexTransversal();
    computeY1Y2();


    //make a vector of {Y,max(condr1, condr2)}
    //Once the Y and cond are loosely correlated
    //Choose Y1 and Y2 to be the two smallest condition numbers

    if (auto search = rToConditionMap.find(localLeftR); search == rToConditionMap.end()) {
        std::vector<double> R1ConditionNumbers;
        K->updateR(localLeftR);

        Y1 = YGrid[0];
        watch(Y1);
        computeMY();
        computeIndexData();
        computeTestPointData();
        refreshKBesselPrecompute(Y1R1);
        populateMatrix(Y1R1);
        R1ConditionNumbers.push_back(conditionNumbers[Y1R1]);
        watch(conditionNumbers[Y1R1]);
        //std::cout << tempStruct.Y << " " << tempStruct.condition << std::endl;
        for (int i = 1; i < YGrid.size(); i ++) {
            Y1 = YGrid[i];
            watch(Y1);
            computeMY();
            computeIndexData();
            computeTestPointData();
            populateMatrix(Y1R1);
            R1ConditionNumbers.push_back(conditionNumbers[Y1R1]);
            watch(conditionNumbers[Y1R1]);
        }
        rToConditionMap[localLeftR] = R1ConditionNumbers;
    }

    if (auto search = rToConditionMap.find(localRightR); search == rToConditionMap.end()) {
        std::vector<double> R2ConditionNumbers;
        K->updateR(localRightR);

        Y1 = YGrid[0];
        watch(Y1);
        computeMY();
        computeIndexData();
        computeTestPointData();
        refreshKBesselPrecompute(Y1R2);
        populateMatrix(Y1R2);
        R2ConditionNumbers.push_back(conditionNumbers[Y1R2]);
        watch(conditionNumbers[Y1R2]);
        //std::cout << tempStruct.Y << " " << tempStruct.condition << std::endl;
        for (int i = 1; i < YGrid.size(); i ++) {
            Y1 = YGrid[i];
            watch(Y1);
            computeMY();
            computeIndexData();
            computeTestPointData();
            populateMatrix(Y1R2);
            R2ConditionNumbers.push_back(conditionNumbers[Y1R2]);
            watch(conditionNumbers[Y1R2]);
        }
        rToConditionMap[localRightR] = R2ConditionNumbers;
    }

    auto leftConditions = rToConditionMap[localLeftR];
    auto rightConditions = rToConditionMap[localRightR];
    std::vector<double> combinedConditionNumbers;
    for (int i = 0; i < leftConditions.size(); i++) {
        combinedConditionNumbers.push_back(std::max(leftConditions[i], rightConditions[i]));
    }

    std::vector<heightConditionPair> limInf;

    double minCondition = (double)pow(10,10);
    int minIndex = 0;
    double secondMinCondition = minCondition;
    int secondMinIndex = 0;

    for (int i = 0; i < combinedConditionNumbers.size(); i++) {
        auto Y = YGrid[i];
        auto condition = combinedConditionNumbers[i];

        watch(Y);
        watch(condition);
        if (condition < minCondition) {
            minCondition = condition;
            minIndex = i;
        }
    }
    for (int i = 0; i < combinedConditionNumbers.size(); i++) {
        auto condition = combinedConditionNumbers[i];

        if (condition < secondMinCondition && condition > minCondition) {
            secondMinCondition = condition;
            secondMinIndex = i;
        }
    }

    Y1 = std::min(YGrid[minIndex], YGrid[secondMinIndex]);
    Y2 = std::max(YGrid[minIndex], YGrid[secondMinIndex]);

    watch(Y1);
    watch(Y2);

    computeMY();
    computeIndexData();
    computeTestPointData();
    refreshKBesselPrecompute(Y1R1);
    populateMatrix(Y1R1);
    populateMatrix(Y2R1);
    refreshKBesselPrecompute(Y1R2);
    populateMatrix(Y1R2);
    populateMatrix(Y2R2);

    std::cout << "-------------SUMMARY-------------" << std::endl;
    std::cout << "final parameters" << std::endl;
    watch(localLeftR);
    watch(localRightR);
    watch(Y1);
    watch(conditionNumbers[Y1R1]);
    watch(conditionNumbers[Y1R2]);
    watch(Y2);
    watch(conditionNumbers[Y2R1]);
    watch(conditionNumbers[Y2R2]);

    //See if there are enough sign changes to proceed
    int signChanges = solveComputeZerosVectorsCountSignChanges();
    int possibleSignChanges = indexTransversal.size() * 2;
    std::cout << "sign changes " << signChanges << "/" << possibleSignChanges << std::endl;
    std::cout << heckeCheck();

    outputFile << "r1: " << localLeftR << " r2: " << localRightR << std::endl;
    outputFile << "Y1: " << Y1 << " Y2: " << Y2 << std::endl;
    outputFile << "sign changes: " << signChanges << "/" << possibleSignChanges << std::endl;
    outputFile << "---------------------------------" << std::endl << std::flush;

    if (signChanges >= possibleSignChanges/3) {
        //Call on the left and right intervals

        double centerR = (localLeftR + localRightR)/2;
        std::cout << "testing left and right subintervals" << std::endl;
        //TODO: passing these parameters like this then globablly changing them doesn't work when working "recursively"
        recursiveSearchForEigenvalues(localLeftR, centerR);
        recursiveSearchForEigenvalues(centerR, localRightR);
        return;
    } else {
        std::cout << "Stopping. Sign change count too low." << std::endl;
        std::cout << "-----------------" << std::endl;
        return;
    }
}



int CoefficientComputer::countSignChanges(const std::vector<double> &v1, const std::vector<double> &v2) {
    int answer = 0;
    assert(v1.size() == v2.size());
    for (int i = 0; i < v1.size(); i++) {
        if (areDifferentSign(v1[i], v2[i])) {
            answer++;
        }
    }
    return answer;
}

/*std::vector<double> CoefficientComputer::setNewRAndComputeZeroesVector(double newR) {
    //auto start = std::chrono::high_resolution_clock::now();

    coefficientMap.clear();
    conditionY1 = 0;
    conditionY2 = 0;

    K->updateR(newR);

    double maxKappaArg = 2*pi/A*M0*maxYStar;
    double minKappaArg = 2*pi/A*1*Y1;
    K->updateRAndPrecompute(newR, minKappaArg, maxKappaArg);

    Eigen::Matrix<double, Eigen::Dynamic, 1> coefficientsFromY1, coefficientsFromY2, zeros1, zeros2;
    coefficientsFromY1.resize(indexTransversal.size());
    coefficientsFromY2.resize(indexTransversal.size());
    zeros1.resize(indexTransversal.size());
    zeros2.resize(indexTransversal.size());

    whichYToUse = 1;
    populateMatrix();
    if (conditionY1 > CONDITION_CUTOFF) {
        return {};
    }
    solveMatrixAndComputeCoefficientMap();

    for (int i = 0; i < indexTransversal.size(); i++) {
        coefficientsFromY1(i) =  coefficientMap[indexTransversal[i]];
    }

    whichYToUse = 2;
    populateMatrix();
    if (conditionY2 > CONDITION_CUTOFF) {
        return {};
    }
    solveMatrixAndComputeCoefficientMap();

    for (int i = 0; i < indexTransversal.size(); i++) {
        coefficientsFromY2(i) =  coefficientMap[indexTransversal[i]];
    }

    zeros1 = matrix2*coefficientsFromY1;
    zeros2 = matrix1*coefficientsFromY2;

    std::vector<double> answer;
    for (int i = 0; i < indexTransversal.size(); i++) {
        answer.push_back(zeros1(i));
    }
    for (int i = 0; i < indexTransversal.size(); i++) {
        answer.push_back(zeros2(i));
    }

    //std::cout << "double zeros" << std::endl;
    //for (auto itr: answer) {
    //    std::cout << std::setprecision(53) << itr << std::endl;
    //}

    return answer;
}*/

bool CoefficientComputer::signChangeVectorIsIncreasing(std::vector<int> &v) {
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


double CoefficientComputer::traceProduct(const std::complex<double> &z, const std::complex<double> &w) {
    auto answer = z*w + conj(z)*conj(w);
    return answer.real();
}

void CoefficientComputer::printTime() {
    watch(computeMTime);
    watch(computeIndexDataTime);
    watch(computeTestPointDataTime);
    watch(populateMatrixTime);
    watch(solveMatrixAndComputeCoefficientMapTime);
}

double CoefficientComputer::kappaApprox(const double &x) {
    return K->approxKBessel(x);
}

double CoefficientComputer::kappaExact(const double &x) {
    return K->computeKBessel(x);
}

CoefficientComputer::~CoefficientComputer() {
    delete K;
    outputFile.close();
}

double CoefficientComputer::findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1) {
    return -y0 * (x1 - x0)/(y1 - y0) + x0;
}

bool CoefficientComputer::areDifferentSign(const double &a, const double &b) {
    return (a > 0 && b < 0) || (a < 0 && b > 0);
}

void CoefficientComputer::computeM0() {
    //The goal is to find the point M0 where
    // \epsilon*kappa(max(t,1)) - kappa(2*pi/A*M0*Y0) == 0
    K->updateR(rightR);
    double minM0 = std::max(rightR, 1.0)/(2*pi/A*Y0);
    double maxM0 = 100;

    double thing = std::min(1.0, kappaExact(rightR));
    double evalLeft = truncation * thing - kappaExact(2*pi/A*minM0*Y0);
    double evalRight = truncation * thing - kappaExact(2*pi/A*maxM0*Y0);
    while (!areDifferentSign(evalLeft, evalRight)) {
        maxM0 *= 2;
        evalRight = truncation * thing - kappaExact(2*pi/A*maxM0*Y0);
    }

    //double temp = modestImpliedConstant * right * localY * kappaExact(2*pi/A*right*localY);
    double left = minM0;
    double right = maxM0;
    evalLeft = truncation * thing - kappaExact(2*pi/A*left*Y0);
    //evalRight = truncation * kappaExact(std::max(r, 1.0)) - kappaExact(2*pi/A*right*Y0);

    while (right - left > 0.000001) {
        double center = (left + right)/2.0;
        double evalCenter = truncation * thing - kappaExact(2*pi/A*center*Y0);

        if (areDifferentSign(evalLeft, evalCenter)) {
            right = center;
            evalRight = evalCenter;
        } else {
            left = center;
            evalLeft = evalCenter;
        }
    }

    M0 = (left + right)/2.0;
    M0++;
    watch(M0);
}

void CoefficientComputer::computeMY() {
    MY = Y0/Y1 * M0;
    //watch(MY);
}

std::vector<std::vector<Quaternion>> CoefficientComputer::getPointAndPullbackPairsY1(const Index &m) {
    double P = 0;
    double Q = 0;
    if (EricNT::mod(-d, 4) == 1) {
        double temp1 = 2*MY/sqrt(4*pow(abs(theta),2) - 1) - m.getB();
        double temp2 = -2*MY/sqrt(4*pow(abs(theta),2) - 1) - m.getB();
        P = std::max(abs(temp1), abs(temp2));
        P = ceil(2*P) + 1;

        temp1 = (2*(MY - m.getA()) - m.getB())/sqrt(d);
        temp2 = (2*(-MY - m.getA()) - m.getB())/sqrt(d);
        Q = std::max(abs(temp1), abs(temp2));
        Q = ceil(2*Q) + 1;
    } else {
        double temp1 = MY/(sqrt(d) - m.getB());
        double temp2 = MY/(sqrt(d) - m.getB());
        P = std::max(abs(temp1), abs(temp2));
        P = ceil(2*P) + 1;

        temp1 = (MY - m.getA())/ sqrt(d);
        temp2 = (-MY - m.getA())/ sqrt(d);
        Q = std::max(abs(temp1), abs(temp2));
        Q = ceil(2*Q) + 1;
    }
    watch(P*Q);


    std::vector<Quaternion> pointsY1;
    pointsY1.reserve(P*Q);
    //Quaternion point = {-1.0L/2 + 1.0L/(2*P) + 1.0L*l0/P, -imagTheta / 2.0L + imagTheta / (2 * Q) + imagTheta * l1 / Q, Y, 0};
    //= {(-P + 1)/(2*P) + l0/P, (-imagTheta*Q + imagTheta)/(2*Q) + imagTheta*l1/Q, Y, 0}

    //baseX0 = (-P+1)/(2*P) = -1/2 + 1/2*P
    double baseX0 = -1.0/2 + 1.0/2*P;

    //baseX1 = (-Q+1)*imagTheta/(2*Q) = -imagTheta/2 + imagTheta/2*Q
    double baseX1 = -theta.imag()/2 + theta.imag()/(2*Q);

    for (int l0 = 0; l0 < P; l0++) {
        for (int l1 = 0; l1 < Q; l1++) {
            double x0 = baseX0 + (double)l0/P;
            double x1 = baseX1 + theta.imag()*l1/Q;

            Quaternion pointY1 = {x0, x1, Y1, 0};
            pointsY1.push_back(pointY1);
        }
    }

    std::vector<std::vector<Quaternion>> pointAndPullbackPairsY1;
    pointAndPullbackPairsY1.resize(P*Q);

    maxYStar = 0;
#pragma omp parallel for
    for (int i = 0; i < pointsY1.size(); i++) {
        Quaternion pointY1 = pointsY1[i];
        Quaternion pullbackY1 = Quaternion(pointY1);
        pullbackY1.reduce(d);
        if (pullbackY1.getJ() > maxYStar) {
            maxYStar = pullbackY1.getJ();
        }
        std::vector<Quaternion> pairY1 = {pointY1,pullbackY1};
        pointAndPullbackPairsY1[i] = pairY1;
    }

    return pointAndPullbackPairsY1;
}

std::vector<std::vector<Quaternion>> CoefficientComputer::getPointAndPullbackPairsY2(const Index &m) {
    double P = 0;
    double Q = 0;
    if (EricNT::mod(-d, 4) == 1) {
        double temp1 = 2*MY/sqrt(4*pow(abs(theta),2) - 1) - m.getB();
        double temp2 = -2*MY/sqrt(4*pow(abs(theta),2) - 1) - m.getB();
        P = std::max(abs(temp1), abs(temp2));
        P = ceil(2*P) + 1;

        temp1 = (2*(MY - m.getA()) - m.getB())/sqrt(d);
        temp2 = (2*(-MY - m.getA()) - m.getB())/sqrt(d);
        Q = std::max(abs(temp1), abs(temp2));
        Q = ceil(2*Q) + 1;
    } else {
        double temp1 = MY/(sqrt(d) - m.getB());
        double temp2 = MY/(sqrt(d) - m.getB());
        P = std::max(abs(temp1), abs(temp2));
        P = ceil(2*P) + 1;

        temp1 = (MY - m.getA())/ sqrt(d);
        temp2 = (-MY - m.getA())/ sqrt(d);
        Q = std::max(abs(temp1), abs(temp2));
        Q = ceil(2*Q) + 1;
    }
    watch(P*Q);


    std::vector<Quaternion> pointsY2;
    pointsY2.reserve(P*Q);
    //Quaternion point = {-1.0L/2 + 1.0L/(2*P) + 1.0L*l0/P, -imagTheta / 2.0L + imagTheta / (2 * Q) + imagTheta * l1 / Q, Y, 0};
    //= {(-P + 1)/(2*P) + l0/P, (-imagTheta*Q + imagTheta)/(2*Q) + imagTheta*l1/Q, Y, 0}

    //baseX0 = (-P+1)/(2*P) = -1/2 + 1/2*P
    double baseX0 = -1.0/2 + 1.0/2*P;

    //baseX1 = (-Q+1)*imagTheta/(2*Q) = -imagTheta/2 + imagTheta/2*Q
    double baseX1 = -theta.imag()/2 + theta.imag()/(2*Q);

    for (int l0 = 0; l0 < P; l0++) {
        for (int l1 = 0; l1 < Q; l1++) {
            double x0 = baseX0 + (double)l0/P;
            double x1 = baseX1 + theta.imag()*l1/Q;

            Quaternion pointY2 = {x0, x1, Y2, 0};
            pointsY2.push_back(pointY2);
        }
    }

    std::vector<std::vector<Quaternion>> pointAndPullbackPairsY2;
    pointAndPullbackPairsY2.resize(P*Q);

#pragma omp parallel for
    for (int i = 0; i < pointsY2.size(); i++) {
        Quaternion pointY2 = pointsY2[i];
        Quaternion pullbackY2 = Quaternion(pointY2);
        pullbackY2.reduce(d);
        std::vector<Quaternion> pairY2 = {pointY2,pullbackY2};
        pointAndPullbackPairsY2[i] = pairY2;
    }

    return pointAndPullbackPairsY2;
}

void CoefficientComputer::computeY1Y2() {

    //The largest I can afford is where it's 0.01 the max of the K-bessel
    //So Y2 is where K_left(--) is 0.01 of its max
    //Y1 is 0.99 of Y2?

    Y1 = Y0;
    Y2 = Y0;
    return;

    /*//Search for Y2 first
    double proportion = 0.01;

    K->updateR(localLeftR);
    auto maxLeft = K->maximize();
    std::cout << "left maximized at (" << maxLeft[0] << ", " << maxLeft[1] << ")" << std::endl;

    double bessMax = maxLeft[1];

    double left = maxLeft[0] / (2 * pi / A * M0);
    double right = Y0;
    double center = (left + right) / 2.0;

    double leftValue = kappaExact(2 * pi / A * M0 * left) - bessMax * proportion;
    double centerValue = kappaExact(2 * pi / A * M0 * center) - bessMax * proportion;
    while (right - left > 0.00001) {
        if (areDifferentSign(leftValue, centerValue)) {
            //zero is in left interval
            right = center;
        } else {
            //zero is in right interval
            left = center;
            leftValue = centerValue;
        }
        center = (left + right) / 2.0;
        centerValue = kappaExact(2 * pi / A * M0 * center) - bessMax * proportion;
    }
    Y2 = center;
    Y1 = 0.9 * Y2;*/
}

std::string CoefficientComputer::heckeCheck() {
    if (d == 1) {
        return std::string();
        //return coefficientMap[Index(1,1,d)]*coefficientMap[Index(3,0,d)] - coefficientMap[Index(3,3,d)];
    } else if (d == 2) {
        return std::string();
        //return coefficientMap[Index(0,1,d)]*coefficientMap[Index(1,1,d)] - coefficientMap[Index(-2,1,d)];
    }
    else if (d == 19) {
        std::stringstream ss;
        for (int i = 0; i < 4; i++) {
            auto temp = coefficientMaps[i];
            auto hecke =  temp[Index(2,0,d)]*temp[Index(3,0,d)] - temp[Index(6,0,d)];
            ss << std::setprecision(16) << "hecke--matrixID " << i << " " << hecke << std::endl;
        }
        return ss.str();
    } else {
        throw(std::invalid_argument("Hecke check not implemented."));
    }
}

void CoefficientComputer::computeIndexTransversal() {
    auto start = std::chrono::high_resolution_clock::now();
    indicesM0.clear();
    indexTransversal.clear();

    //aub = ceil(maxN * thetaModulus / sqrtl(pow(thetaModulus, 2) - pow(thetaReal, 2)));

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
            Index index = Index(a, b, d);
            if (index.getAbs() <= M0) {
                    indicesM0.push_back(index);
            }
        }
    }

    /*Sort the indicesM0 by absolute value then by angle with positive real axis.*/
    std::sort(indicesM0.begin(), indicesM0.end());

    std::vector<Index> tempIndices;
    copy(indicesM0.begin(), indicesM0.end(), back_inserter(tempIndices));

    int rotationCoeff = symClass == 'D' || symClass == 'G' ? 1 : -1;
    int conjCoeff = symClass == 'D' || symClass == 'C' ? 1 : -1;

    while (!tempIndices.empty()) {
        Index index = tempIndices.front();
        std::vector<std::tuple<Index, int>> orbit = {std::tuple<Index, int>(index, 1)};

        Index tempIndex = index.rotate();
        int tempCoeff = 1*rotationCoeff;
        std::tuple<Index, int> tempTuple = {tempIndex, tempCoeff};
        while (tempIndex != index) {
            orbit.push_back(tempTuple);

            tempIndex = tempIndex.rotate();
            tempCoeff = tempCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
        }

        Index conjIndex = index.conj();
        bool conjIsInRotations = false;
        for (auto tup : orbit) {
            if (conjIndex == std::get<0>(tup)) {
                conjIsInRotations = true;
                break;
            }
        }

        if (!conjIsInRotations) {
            /*add in the rotations of the conjugate*/
            tempTuple = {conjIndex, conjCoeff};
            orbit.push_back(tempTuple);

            tempIndex = conjIndex.rotate();
            tempCoeff = conjCoeff*rotationCoeff;
            tempTuple = {tempIndex, tempCoeff};
            while (tempIndex != conjIndex) {
                orbit.push_back(tempTuple);

                tempIndex = tempIndex.rotate();
                tempCoeff = tempCoeff*rotationCoeff;
                tempTuple = {tempIndex, tempCoeff};
            }

        }
        /*Save our starting index to the transversal.*/
        indexTransversal.push_back(index);

        /*Save the orbit data.*/
        indexOrbitData[index] = orbit;

        /*Delete the indicesM0 in the orbit from the copied list of indicesM0.*/
        for (auto itr1 : orbit) {
            Index toDelete = std::get<0>(itr1);
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

}

void CoefficientComputer::adjustY1() {
    Y1 *= 0.99;
}

void CoefficientComputer::adjustY2() {
    Y2  = Y1 + 0.99*(Y2-Y1);
}

void CoefficientComputer::refreshKBesselPrecompute(short matrixID) {
    double tempR;
    if (matrixID == Y1R1 || matrixID == Y2R1) {
        tempR = leftR;
    } else {
        tempR = rightR;
    }
    K->updateR(tempR);

    double maxKappaArg = 2*pi/A*M0*maxYStar;
    double minKappaArg = 2*pi/A*1*Y1;
    K->updateRAndPrecompute(tempR, minKappaArg, maxKappaArg);
}

int CoefficientComputer::solveComputeZerosVectorsCountSignChanges() {
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> coefficients;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> zeros;

    coefficients.resize(4);
    zeros.resize(4);
    for (auto &temp : coefficients) {
        temp.resize(indexTransversal.size());
    }

    for (auto &temp : zeros) {
        temp.resize(indexTransversal.size());
    }

    for (short i = 0; i < 4; i++) {
        solveMatrixAndComputeCoefficientMap(i);
    }

    for (short i = 0; i < 4; i++) {
        auto tempCoefficientMap = coefficientMaps[i];
        for (int j = 0; j < indexTransversal.size(); j++) {
            auto index = indexTransversal[j];
            auto coeff = tempCoefficientMap[index];
            coefficients[i](j) = coeff;
        }
    }

    zeros[Y1R1] = matrices[Y1R1] * coefficients[Y2R1];
    zeros[Y2R1] = matrices[Y2R1] * coefficients[Y1R1];

    zeros[Y1R2] = matrices[Y1R2] * coefficients[Y2R2];
    zeros[Y2R2] = matrices[Y2R2] * coefficients[Y1R2];

    std::vector<double> left;
    std::vector<double> right;
    for (int i = 0; i < zeros[Y1R1].size(); i++) {
        left.push_back(zeros[Y1R1][i]);
    }
    for (int i = 0; i < zeros[Y2R1].size(); i++) {
        left.push_back(zeros[Y2R1][i]);
    }

    for (int i = 0; i < zeros[Y1R2].size(); i++) {
        right.push_back(zeros[Y1R2][i]);
    }
    for (int i = 0; i < zeros[Y2R2].size(); i++) {
        right.push_back(zeros[Y2R2][i]);
    }

    return countSignChanges(left, right);
}

void CoefficientComputer::checkSingleEigenvalue(double r) {
    K->updateR(r);
    leftR = r;
    rightR = r;

    computeM0();

    /*std::vector<double> rConditionNumbers;

    Y1 = YGrid[0];
    Y2 = Y0;
    watch(Y1);
    computeMY();
    computeIndexData();
    computeTestPointData();
    refreshKBesselPrecompute(Y1R1);
    populateMatrix(Y1R1);
    watch(conditionNumbers[Y1R1]);
    rConditionNumbers.push_back(conditionNumbers[Y1R1]);
    for (int i = 1; i < YGrid.size(); i ++) {
        Y1 = YGrid[i];
        watch(Y1);
        computeMY();
        computeIndexData();
        computeTestPointData();
        populateMatrix(Y1R1);
        rConditionNumbers.push_back(conditionNumbers[Y1R1]);
        watch(conditionNumbers[Y1R1]);
    }
    rToConditionMap[r] = rConditionNumbers;

    double minCondition = pow(10,10);
    double minConditionY = Y0;
    for (int i = 0; i < YGrid.size(); i++) {
        if (rConditionNumbers[i] < minCondition) {
            minCondition = rConditionNumbers[i];
            minConditionY = YGrid[i];
        }
    }

    std::cout << std::setprecision(16) << "minCondition is " << minCondition << " for Y " << minConditionY << std::endl;*/

    Y1 = 0.065;
    Y2 = Y1;
    computeMY();
    computeIndexData();
    computeTestPointData();
    refreshKBesselPrecompute(Y1R1);
    populateMatrix(Y1R1);
    solveMatrixAndComputeCoefficientMap(Y1R1);

    //generate a bunch of random points above Y
    //pull back
    //evaluate at paired points
    //compute differences
    std::vector<Quaternion> testPoints;
    std::vector<Quaternion> testPointPullbacks;

    /*PlotWindow window = PlotWindow(16.0/9, 1.5, {0,Y1},1080);
    auto domainPoints = window.getPoints();
    FunctionToEvaluate func = FunctionToEvaluate();

    std::vector<std::vector<double>> evaluatedPoints;
    for (auto row : domainPoints) {
        std::vector<double> evaluatedRow;
        for (auto point : row) {
            Quaternion z = Quaternion(0, point.x, point.y, 0);
            auto evaluation = evaluate(z);
            evaluatedRow.push_back(evaluation);
        }
        evaluatedPoints.push_back(evaluatedRow);
    }

    std::string filename = "quad2.png";
    Plotter plotter = Plotter(evaluatedPoints);
    plotter.histogramSmoothing = false;
    plotter.makeImage(filename, "viridis");*/
    return;

    while (testPoints.size() < 20) {
        double height = Y0 + (1-Y0) * rand01();
        double x0 = rand01()*10.0 - 5.0;
        double x1 = rand01()*10.0 - 5.0;
        Quaternion point = Quaternion(x0, x1, height, 0);
        Quaternion pullback = Quaternion(point);
        pullback.reduce(d);
        if ((point - pullback).abs() > 0.000001) {
            testPoints.push_back(point);
            testPointPullbacks.push_back(pullback);
        }
    }

    std::vector<double> differences;
    for (int i = 0; i < testPoints.size(); i++) {
        std::complex<double> testPointEvaluation = evaluate(testPoints[i]);
        std::complex<double> pullbackEvaluation = evaluate(testPointPullbacks[i]);
        differences.push_back(std::abs(testPointEvaluation - pullbackEvaluation));
    }
    std::cout << "abs of differences" << std::endl;
    for (auto diff : differences) {
        std::cout << std::setprecision(16) << diff << std::endl;
    }
}

double CoefficientComputer::rand01() {
    return std::rand()/ (RAND_MAX + 1.0);
}

void CoefficientComputer::checkRamanujanPetersson(double r) {
    //check that coefficients for r have been computed
    //verify that |a(p)| \leq 2*N(p)^(-1/2) <--ideal norm

    //make a list of ideals which are within the bounds of computed coefficients
    //for indexM0
    //if prime
    //check |a(p)| \leq 2 * N(p)^(-1/2)

    //we will compute prime elements up to the bound M0
    std::vector<int> rationalPrimes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};
    std::vector<Index> primes;

    for (auto index : indicesM0) {
        double indexNorm = std::real(index.getComplex()*conj(index.getComplex()));
        if (indexNorm > M0) {
            continue;
        }
    }


}

void CoefficientComputer::computeFourierCoefficients(double r) {
    //set r
    //initialize indices and testpoints
    //populate matrix
    //check conditioning, adjusting Y if necessary
    //solve system
    //continue computing coefficients... up to M(Y)?
    //save them in the pointwise coefficient map
    //indicate that this process has been done for this value of r in a private class variable
}

void CoefficientComputer::checkSatoTate(double r) {
    //check that coefficients for r have been computed
    //check sato tate
    //output a text file, probably?
}

void CoefficientComputer::produceLFunction(double r) {
    //check that coefficients for r have been computed
    //idk
}

void CoefficientComputer::printCoefficients() {
    for (auto index : indexTransversal) {
        std::cout << index << ", " << coefficientMaps[Y1R1][index] << "\n";
    }
}

void CoefficientComputer::computeMoreCoefficients() {
    Y1 *= 0.99;
    computeMY();
    computeIndexData();
    computeTestPointData();
    std::cout << "started here" << std::endl;

#pragma omp parallel for
    for (int i = 0; i < indexTransversalMY.size(); i++) {
        auto m = indexTransversalMY[i];
        if (m.getAbs() <= M0) {
            continue;
        }

        double errorDenominator = Y1 * kappaExact(2*pi/A * m.getAbs() * Y1);
        double error = 2.0*truncation/errorDenominator;
        if (error > 0.1) {
            std::cout << m << " skipped" << std::endl;
            continue;
        }
        double sum = 0;
        for (auto n : indexTransversal) {
            sum += coefficientMaps[Y1R1][n] * computeEntry(Y1R1, m, n);
        }
        sum /= errorDenominator;
        coefficientMaps[Y1R1][m] = sum;
        if (i % 100 == 0) {
            std::cout << i << "/" << indexTransversalMY.size() << std::endl;
        }

    }

    for (int i = 0; i < indexTransversalMY.size(); i++) {
        Index index = indexTransversalMY[i];
        if (index.getAbs() <= M0) {
            continue;
        }
        auto orbit = indexOrbitData[index];
        auto coeff = coefficientMaps[Y1R1][index];
        for (auto tup : orbit) {
            double mult = std::get<1>(tup);
            coefficientMaps[Y1R1][std::get<0>(tup)] = coeff * mult;
        }
    }

    for (auto itr : coefficientMaps[Y1R1]) {
        outputFile << std::setprecision(16) << itr.first << ", " << itr.second << "\n";
    }
}

void CoefficientComputer::secantSearch(double startR1, double startR2) {

    double r1 = startR1;
    std::cout << std::setprecision(16) << "r: " << r1;
    this->leftR = r1;
    this->rightR = r1;
    K->updateR(r1);
    computeM0();


    double Y = 0.9*r1/(2*pi/A*M0);
    watch(Y);
    Y1 = Y;
    Y2 = Y;

    computeMY();
    watch(MY);
    computeIndexData();
    watch(indexTransversal.size());

    Eigen::MatrixXd matrix;
    matrix.resize(indexTransversal.size(), indexTransversal.size());

    for (int j = 0; j < indexTransversal.size(); j++) {
        Index m = indexTransversal[j];
        auto points = pointwisePointsAndPullbacks(m, Y);

        if (j == 0) {
            std::cout << "teehee" << std::endl;
            K->updateRAndPrecompute(r1, 2*pi/A*1*Y, 2*pi/A*M0*maxYStar);
            std::cout << "teehee" << std::endl;
        } else {
            K->extendPrecomputedRange(2*pi/A*M0*maxYStar);
        }

        auto row = pointwiseGenerateMatrixRow(m, Y, points);
        for (int k = 0; k < indexTransversal.size(); k++) {
            matrix(j,k) = row[k];
        }
    }

    matrices[Y1R1] = matrix;
    solveMatrixAndComputeCoefficientMap(Y1R1);
    for (auto itr : coefficientMaps[Y1R1]) {
        outputFile << std::setprecision(16) << itr.first << ", " << itr.second << "\n";
    }
    double cost1 = heckeCheck(Y1R1);
    std::cout << std::setprecision(16) << " cost: " << cost1 << std::endl;

    //
    //
    //
    //

    double r2 = startR2;
    std::cout << std::setprecision(16) << "r: " << r2;
    this->leftR = r2;
    this->rightR = r2;
    K->updateR(r2);
    computeM0();

    Y = 0.9*r2/(2*pi/A*M0);
    watch(Y);
    Y1 = Y;
    Y2 = Y;

    computeMY();
    computeIndexData();

    matrix.resize(indexTransversal.size(), indexTransversal.size());


    for (int j = 0; j < indexTransversal.size(); j++) {
        Index m = indexTransversal[j];
        auto points = pointwisePointsAndPullbacks(m, Y);

        if (j == 0) {
            K->updateRAndPrecompute(r2, 2*pi/A*1*Y, 2*pi/A*M0*maxYStar);
        } else {
            K->extendPrecomputedRange(2*pi/A*M0*maxYStar);
        }

        auto row = pointwiseGenerateMatrixRow(m, Y, points);
        for (int k = 0; k < indexTransversal.size(); k++) {
            matrix(j,k) = row[k];
        }
    }
    matrices[Y1R1] = matrix;
    solveMatrixAndComputeCoefficientMap(Y1R1);
    double cost2 = heckeCheck(Y1R1);
    std::cout << std::setprecision(16) << " cost: " << cost2 << std::endl;

    double stop = std::pow(10,-16);
    int count = 2;
    while (abs(r1 - r2) > stop) {
        double nextR = findZeroOfLinearInterpolation(r1, cost1, r2, cost2);
        std::cout << std::setprecision(16) << "r: " << nextR;
        this->leftR = nextR;
        this->rightR = nextR;
        K->updateR(nextR);
        computeM0();

        Y = 0.9*r2/(2*pi/A*M0);
        watch(Y);
        Y1 = Y;
        Y2 = Y;

        computeMY();
        computeIndexData();

        matrix.resize(indexTransversal.size(), indexTransversal.size());

        //  for m in indexTransversal:
        for (int j = 0; j < indexTransversal.size(); j++) {
            Index m = indexTransversal[j];
            auto points = pointwisePointsAndPullbacks(m, Y);
            //      if (lowest Y)
            if (j == 0) {
                //          precompute K-bessel function
                K->updateRAndPrecompute(nextR, 2*pi/A*1*Y, 2*pi/A*M0*maxYStar);
            } else {
                K->extendPrecomputedRange(2*pi/A * M0 * maxYStar);
            }
            //      matrix.row = generateRow(points)
            auto row = pointwiseGenerateMatrixRow(m, Y, points);
            for (int k = 0; k < indexTransversal.size(); k++) {
                matrix(j,k) = row[k];
            }
        }
        matrices[Y1R1] = matrix;
        solveMatrixAndComputeCoefficientMap(Y1R1);
        double nextCost = heckeCheck(Y1R1);
        std::cout << std::setprecision(16) << " cost: " << nextCost << std::endl;
        /////
        r1 = r2;
        cost1 = cost2;

        r2 = nextR;
        cost2 = nextCost;
        count++;
        if (count > 20) {
            break;
        }
    }


}

double CoefficientComputer::heckeCheck(short matrixID) {
    if (d == 1) {
        return 0;
        //return coefficientMap[Index(1,1,d)]*coefficientMap[Index(3,0,d)] - coefficientMap[Index(3,3,d)];
    } else if (d == 2) {
        return 0;
        //return coefficientMap[Index(0,1,d)]*coefficientMap[Index(1,1,d)] - coefficientMap[Index(-2,1,d)];
    }
    else if (d == 19) {
        auto temp = coefficientMaps[matrixID];
        auto hecke =  temp[Index(2,0,d)]*temp[Index(3,0,d)] - temp[Index(6,0,d)];
        return hecke;
    } else {
        throw(std::invalid_argument("Hecke check not implemented."));
    }
}

std::vector<std::vector<Quaternion>> CoefficientComputer::pointwisePointsAndPullbacks(Index m, double Y) {
    maxYStar = 0;


    int Q0 = 0;
    int Q1 = (MY + abs(m.getComplex().real()))/2.0;
    if (EricNT::mod(-d,4) == 1) {
        Q0 = std::max(MY + abs(m.getComplex().real()), (MY + abs(m.getComplex().imag()))/(2.0*A));
    } else {
        Q0 = (MY + abs(m.getComplex().imag()))/(2*A);
    }
    Q0 += 2;
    Q1 += 2;

    int totalNumberOfPoints = (2*Q0)*(2*Q1)/2;

    std::vector<Quaternion> points;
    points.reserve(totalNumberOfPoints);

    for (int l0 = 1-Q0; l0 <= Q0; l0++) {
        for (int l1 = 1; l1 <= Q1; l1++) {
        //for (int l1 = 1; l1 <= Q1; l1++) {
            std::complex<double> x = ((double)(l0-0.5))/(2.0*Q0) + theta*((double)(l1-0.5))/((double)2.0*Q1);
            double x0 = x.real();
            double x1 = x.imag();

            Quaternion point = {x0, x1, Y, 0};
            points.push_back(point);
        }
    }

    std::vector<std::vector<Quaternion>> pointAndPullbackPairs;
    pointAndPullbackPairs.resize(totalNumberOfPoints);

#pragma omp parallel for
    for (int i = 0; i < points.size(); i++) {
        Quaternion point = points[i];
        Quaternion pullback = Quaternion(point);
        pullback.reduce(d);
        assert(pullback.getJ() > Y0);
        std::vector<Quaternion> pair = {point,pullback};
        pointAndPullbackPairs[i] = pair;
    }

    for (int i = 0; i < points.size(); i++) {
        Quaternion pullback = pointAndPullbackPairs[i][1];
        if (pullback.getJ() > maxYStar) {
            maxYStar = pullback.getJ();
        }
    }

    return pointAndPullbackPairs;
}

std::vector<double> CoefficientComputer::pointwiseGenerateMatrixRow(const Index &m,
                                                                    double Y,
                                                                    const std::vector<std::vector<Quaternion>> &pointsAndPullbacks) {
    std::vector<double> row;
    row.resize(indexTransversal.size());
#pragma omp parallel for
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index n = indexTransversal[i];

        double answer = 0;

        for (auto pair : pointsAndPullbacks) {
            auto point = pair[0];
            auto pullback = pair[1];

            std::complex<double> x = point.getComplex();
            std::complex<double> xStar = pullback.getComplex();

            double yStar = pullback.getJ();
            assert(yStar > Y0);

            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double term = yStar * kappaApprox(2*pi/A*n.getAbs()*yStar);

            if (symClass== 'D' || symClass == 'G' || d == 1) {
                term *= std::cos(pi/A * traceProduct(-I*m.getComplex(),x));
            } else {
                term *= -std::sin(pi/A * traceProduct(-I*m.getComplex(),x));
            }


            //e(symClass, n, x^*)
            //sum_{l \in S_n} +/- 1 * exp(pi I /A * <I*l, xstar>)
            double cs = 0;
            if (symClass == 'D' || symClass == 'G' || d == 1) {
                for (auto tup : indexOrbitDataModMinusOne[n]) {
                    Index l = std::get<0>(tup);
                    double ellTerm = std::get<1>(tup);
                    ellTerm *= std::cos(pi/A * traceProduct(I*l.getComplex(),xStar));
                    cs += ellTerm;
                }
            } else {
                for (auto tup : indexOrbitDataModMinusOne[n]) {
                    Index l = std::get<0>(tup);
                    double ellTerm = std::get<1>(tup);
                    ellTerm *= std::sin(pi/A * traceProduct(I*l.getComplex(),xStar));
                    cs += ellTerm;
                }
            }
            term *= cs;

            answer += term;
        }

        //Right now answer = sum_(points) Y^* kappa(2*pi/A*|n|Y^*)e(symClass, n, x^*) exp(-pi *I/A *<I*m, x>)
        answer *= -4.0/(2.0*pointsAndPullbacks.size());

        //Add on the kronecker delta term
        if (m == n) {
            //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
            double deltaTerm = Y * kappaApprox(2*pi/A * m.getAbs() * Y);

            answer += deltaTerm;
        }
        row[i] = answer;
    }
    return row;
}

unsigned long long int CoefficientComputer::pointwisePointCount(Index m, double Y) {


    int Q0 = 0;
    int Q1 = (MY + abs(m.getComplex().real()))/2.0;
    if (EricNT::mod(-d,4) == 1) {
        Q0 = std::max(MY + abs(m.getComplex().real()), (MY + abs(m.getComplex().imag()))/(2.0*A));
    } else {
        Q0 = (MY + abs(m.getComplex().imag()))/(2*A);
    }
    Q0 += 2;
    Q1 += 2;

    return 2*Q0*2*Q1;
}

std::vector<TestPointOrbitData> CoefficientComputer::getPointPullbackOrbit(const Index &m, const double &Y) {
    std::vector<TestPointOrbitData> answer;
    maxYStar = 0;

    //These formulas provide bounds to guarantee that the DFT result is valid
    int Q0 = 0;
    int Q1 = std::ceil((MY + abs(m.getComplex().real()))/2.0);
    if (EricNT::mod(-d,4) == 1) {
        int Q0option1 = std::ceil(MY + abs(m.getComplex().real()));
        int Q0option2 = std::ceil((MY + abs(m.getComplex().imag()))/(2.0*A));
        Q0 = std::max(Q0option1, Q0option2);
    } else {
        Q0 = std::ceil((MY + abs(m.getComplex().imag()))/(2*A));
    }


    //Tweak Q0 and Q1 to be exactly what we need for exploiting symmetry
    if (EricNT::mod(-d, 4) == 1 && d != 3) {
        //If -d = 1 mod 4 then we need Q0/Q1 to be an even integer for the map x -> -bar(x) to be defined on test points
        if (EricNT::mod(Q0, Q1) == 0 && EricNT::mod(Q0/Q1, 2) == 0) {
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
                int firstOptionQ0 = std::ceil(Q1 * 2 * (k+1));
                while (!EricNT::mod(firstOptionQ0, 2 * (k+1) == 0)) {
                    firstOptionQ0++;
                }
                int firstOptionQ1 = firstOptionQ0/(2*(k+1));

                int secondOptionQ0 = std::ceil(Q1 * 2 * k);
                while (!EricNT::mod(secondOptionQ0, 2 * k == 0)) {
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
        Q0 = std::max(Q0, Q1) + 1;
        Q1 = Q0;
    } else {
        //We don't need to fiddle with the divisibility of Q0 and Q1
        //Add some numerical wiggle room
        Q0 += 1;
        Q1 += 1;
    }
    numberOfPoints = 2*Q0*2*Q1;

    char nonsquareUnitSign = (symClass == 'D' || symClass == 'G') ? 1 : -1;
    char reflectionSign = (symClass == 'D' || symClass == 'C') ? 1 : -1;

    //Now generate the representatives
    //then generate the orbits
    if (EricNT::mod(-d, 4) == 1 && d != 3) {
        //iterate to make orbit representatives
        for (int l1 = 1; l1 <= Q1; l1++) {
            int lowerBound = std::ceil(0.5 - Q0*(l1-0.5)/(2.0*Q1));
            int upperBound = std::floor(0.5 + Q0 - Q0*(l1 - 0.5)/(2.0*Q1));
            for (int l0 = lowerBound; l0 <= upperBound; l0++) {
                //build the orbit
                TestPointOrbitData orbit;
                std::complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
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

                    answer.push_back(orbit);
                } else {
                    //The orbit in this case is {x, -x, -bar(x), bar(x)}
                    //So orbit/{+/-1} = {[x], [-bar(x)]}
                    //So the only proper translate is -bar(x)
                    std::tuple<std::complex<double>, char> tup (-conj(x), reflectionSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    answer.push_back(orbit);
                }
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
                std::complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                std::tuple<std::complex<double>, char> tup (-conj(x), reflectionSign);
                orbit.properTranslatesModSign.push_back(tup);

                answer.push_back(orbit);
            }
        }
    } else if (d == 1) {
        //Iterate through orbit representatives
        for (int l0 = 1; l0 <= Q0; l0++) {
            for (int l1 = 1; l1 <= Q1; l1++) {
                TestPointOrbitData orbit;

                std::complex<double> x = (l0 - 0.5)/(2.0*Q0) + theta*(l1 - 0.5)/(2.0*Q1);
                Quaternion pullback = Quaternion(x.real(), x.imag(), Y, 0);
                pullback.reduce(d);
                maxYStar = maxYStar < pullback.getJ() ? pullback.getJ() : maxYStar;

                orbit.representativeComplex = x;
                orbit.representativePullback = pullback;

                if (abs(l0) == abs(l1)) {
                    //The orbit is {x, ix, -x, -ix}
                    //So orbit/{+/-1} is {[x], [ix]}
                    //So the only proper translate is ix

                    std::tuple<std::complex<double>, char> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);
                } else {
                    //The orbit is {x, ix, -x, -ix, bar(x), ibar(x), -bar(x), -ibar(x)}
                    //So orbit/{+/-1} is {[x], [ix], [-bar(x)], [-ibar(x)]}
                    //So the proper translates are ix, -bar(x), and -ibar(x)

                    std::tuple<std::complex<double>, char> tup (I*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    std::get<0>(tup) = -conj(x);
                    std::get<1>(tup) = reflectionSign;
                    orbit.properTranslatesModSign.push_back(tup);

                    std::get<0>(tup) = -I*conj(x);
                    std::get<1>(tup) = reflectionSign * nonsquareUnitSign;
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
        int l1UpperBound = std::floor(2.0*Q1/3.0);
        for (int l1 = l1LowerBound; l1 <= l1UpperBound; l1++) {
            int l0LowerBound = l1; //should really be std::ceil(Q0*l1/(double)Q1), but Q0=Q1
            int l0UpperBound = std::floor(Q0 - l1/2.0); //should really be std::floor(Q0 - l1*Q0/(2*Q1)), but Q0=Q1
            for (int l0 = l0LowerBound; l0 <= l0UpperBound; l0++) {
                TestPointOrbitData orbit;

                std::complex<double> x = l0 / (2.0 * Q0) + theta * (double) l1 / (2.0 * Q1);
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

                    std::tuple<std::complex<double>, char> tup (theta*x, nonsquareUnitSign);
                    orbit.properTranslatesModSign.push_back(tup);

                    std::get<0>(tup) = theta*theta*x;
                    std::get<1>(tup) = 1; //theta*theta is a square unit
                    orbit.properTranslatesModSign.push_back(tup);

                    //magic quantity that tells me about the real part of x
                    int discriminant = 2*Q1*l0 + Q1*l1 - 2*Q0*Q1;
                    bool orbitIsHalfSize = (l0 == 0 || l1 == 0 || l0 + l1 == 0 || discriminant == 0);

                    if (!orbitIsHalfSize) {
                        //For this comment let zeta=e^(2pi i /6) be represented by w
                        //The orbit is {x, wx, w^2x, -x, -wx, -w^2x, -bar(x), -wbar(x), -w^2bar(x), bar(x), wbar(x), w^2bar(x)}
                        //So orbit/{+/-1} is {[x], [wx], [w^2x], [-bar(x)], [-wbar(x)], [-w^2bar(x)]}
                        //So the proper translates are wx, w^2x, -bar(x), -wbar(x), -w^2bar(x)
                        std::get<0>(tup) = -conj(x);
                        std::get<1>(tup) = reflectionSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        std::get<0>(tup) = -theta*conj(x);
                        std::get<1>(tup) = reflectionSign * nonsquareUnitSign;
                        orbit.properTranslatesModSign.push_back(tup);

                        std::get<0>(tup) = -theta*theta*conj(x);
                        std::get<1>(tup) = reflectionSign; //theta*theta is a square unit
                        orbit.properTranslatesModSign.push_back(tup);
                    }
                }
            }
        }
    }

    return answer;
}

void CoefficientComputer::populateMatrix2(double Y) {


    matrix2.resize(indexTransversal.size(), indexTransversal.size());

    for (int i = 0; i < indexTransversal.size(); i++) {
        Index m = indexTransversal[i];
        testPointOrbits = getPointPullbackOrbit(m, Y);
        #pragma omp parallel for default(none) shared(matrix2, indexTransversal, i, m)
        for (int j = 0; j < indexTransversal.size(); j++) {

            Index n = indexTransversal[j];
            double entry = computeEntry2(m, n);
            //std::cout << entry << std::endl;

            matrix2(i,j) = entry;
        }
    }

}

double CoefficientComputer::computeEntry2(const Index &m, const Index &n) {
    double answer = 0;

    if (symClass == 'D' || symClass == 'G' || d == 1) {
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;


            //toAdd =  ystar * kappa(2*pi/A*|n|*ystar) * e(symClass, n, x^*) * exp(-pi*I/A * <I*m,x>)
            double yStar = pullback.getJ();
            double term = yStar * kappaApprox(2*pi/A*n.getAbs()*yStar);

            double indexTerm = 0;
            for (auto tup : indexOrbitDataModMinusOne[n]) {
                Index l = std::get<0>(tup);
                std::complex<double> xStar = pullback.getComplex();
                double ellTerm = std::get<1>(tup); //This is a_l/a_n
                ellTerm *= std::cos(pi/A * traceProduct(I*l.getComplex(),xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = std::cos(pi/A * traceProduct(-I*m.getComplex(), x));
            for (auto tup : testPointOrbit.properTranslatesModSign) {
                std::complex<double> etaX = std::get<0>(tup);
                char sign = std::get<1>(tup);
                testPointTerm += sign * std::cos(pi/A * traceProduct(-I*m.getComplex(), etaX));
            }
            term *= testPointTerm;

            answer += term;
        }
        answer *= -4.0/numberOfPoints;
    } else {
        for (const auto& testPointOrbit : testPointOrbits) {
            auto x = testPointOrbit.representativeComplex;
            auto pullback = testPointOrbit.representativePullback;

            double yStar = pullback.getJ();
            double term = yStar * kappaApprox(2*pi/A*n.getAbs()*yStar);

            double indexTerm = 0;
            for (auto tup : indexOrbitDataModMinusOne[n]) {
                Index l = std::get<0>(tup);
                std::complex<double> xStar = pullback.getComplex();
                double ellTerm = std::get<1>(tup); //This is a_l/a_n
                ellTerm *= std::sin(pi/A * traceProduct(I*l.getComplex(),xStar));
                indexTerm += ellTerm;
            }
            term *= indexTerm;

            double testPointTerm = std::sin(pi/A * traceProduct(-I*m.getComplex(), x));
            for (auto tup : testPointOrbit.properTranslatesModSign) {
                std::complex<double> etaX = std::get<0>(tup);
                char sign = std::get<1>(tup);
                testPointTerm += sign * std::sin(pi/A * traceProduct(-I*m.getComplex(), etaX));
            }
            term *= testPointTerm;

            answer += term;
        }
        answer *= +4.0/numberOfPoints;
    }

    //Add on the kronecker delta term
    if (m == n) {
        //deltaTerm = Y * kappa(2*pi/A*|m|*Y)
        double deltaTerm = Y * kappaExact(2*pi/A *m.getAbs() * Y);

        answer += deltaTerm;
    }
    return answer;
}

void CoefficientComputer::checkSingleEigenvalue2(double r) {
    K->updateR(r);
    leftR = r;
    rightR = r;

    computeM0();

    /*std::vector<double> rConditionNumbers;

    Y1 = YGrid[0];
    Y2 = Y0;
    watch(Y1);
    computeMY();
    computeIndexData();
    computeTestPointData();
    refreshKBesselPrecompute(Y1R1);
    populateMatrix(Y1R1);
    watch(conditionNumbers[Y1R1]);
    rConditionNumbers.push_back(conditionNumbers[Y1R1]);
    for (int i = 1; i < YGrid.size(); i ++) {
        Y1 = YGrid[i];
        watch(Y1);
        computeMY();
        computeIndexData();
        computeTestPointData();
        populateMatrix(Y1R1);
        rConditionNumbers.push_back(conditionNumbers[Y1R1]);
        watch(conditionNumbers[Y1R1]);
    }
    rToConditionMap[r] = rConditionNumbers;

    double minCondition = pow(10,10);
    double minConditionY = Y0;
    for (int i = 0; i < YGrid.size(); i++) {
        if (rConditionNumbers[i] < minCondition) {
            minCondition = rConditionNumbers[i];
            minConditionY = YGrid[i];
        }
    }

    std::cout << std::setprecision(16) << "minCondition is " << minCondition << " for Y " << minConditionY << std::endl;*/

    Y1 = 0.065;
    Y2 = Y1;
    Index largestM = Index(std::ceil(MY), std::ceil(MY/abs(theta.imag())), d);
    computeMY();
    computeIndexData();
    getPointPullbackOrbit(largestM, Y1);
    refreshKBesselPrecompute(Y1R1);
    //populateMatrix(Y1R1);
    populateMatrix2(Y1);
    solveMatrixAndComputeCoefficientMap2();

    //generate a bunch of random points above Y
    //pull back
    //evaluate at paired points
    //compute differences
    std::vector<Quaternion> testPoints;
    std::vector<Quaternion> testPointPullbacks;

    /*PlotWindow window = PlotWindow(16.0/9, 1.5, {0,Y1},1080);
    auto domainPoints = window.getPoints();
    FunctionToEvaluate func = FunctionToEvaluate();

    std::vector<std::vector<double>> evaluatedPoints;
    for (auto row : domainPoints) {
        std::vector<double> evaluatedRow;
        for (auto point : row) {
            Quaternion z = Quaternion(0, point.x, point.y, 0);
            auto evaluation = evaluate(z);
            evaluatedRow.push_back(evaluation);
        }
        evaluatedPoints.push_back(evaluatedRow);
    }

    std::string filename = "quad2.png";
    Plotter plotter = Plotter(evaluatedPoints);
    plotter.histogramSmoothing = false;
    plotter.makeImage(filename, "viridis");
    return;*/

    while (testPoints.size() < 20) {
        double height = Y0 + (1-Y0) * rand01();
        double x0 = rand01()*10.0 - 5.0;
        double x1 = rand01()*10.0 - 5.0;
        Quaternion point = Quaternion(x0, x1, height, 0);
        Quaternion pullback = Quaternion(point);
        pullback.reduce(d);
        if ((point - pullback).abs() > 0.000001) {
            testPoints.push_back(point);
            testPointPullbacks.push_back(pullback);
        }
    }

    std::vector<double> differences;
    for (int i = 0; i < testPoints.size(); i++) {
        std::complex<double> testPointEvaluation = evaluate2(testPoints[i]);
        std::complex<double> pullbackEvaluation = evaluate2(testPointPullbacks[i]);
        differences.push_back(std::abs(testPointEvaluation - pullbackEvaluation));
    }
    std::cout << "abs of differences" << std::endl;
    for (auto diff : differences) {
        std::cout << std::setprecision(16) << diff << std::endl;
    }
}

void CoefficientComputer::solveMatrixAndComputeCoefficientMap2() {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> x;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> b;
    Eigen::MatrixXcd subMatrix;

    auto tempCoefficientMap = coefficientMap2;
    tempCoefficientMap.clear();
    long size = indexTransversal.size() - 1;
    x.resize(size);
    b.resize(size);
    subMatrix.resize(size, size);

    //set b to negative of the indexOfNormalization column
    //set tempA to be matrixA deleting row and column indexOfNormalization
    auto tempMatrix = matrix2;
    for (int row  = 0; row < indexTransversal.size(); row++) {
        for (int column = 0; column < indexTransversal.size(); column++) {
            if (column == indexOfNormalization && row != indexOfNormalization) {
                //add it to b
                if (row > indexOfNormalization) {
                    b(row - 1) = tempMatrix(row, column);
                } else {
                    b(row) = tempMatrix(row, column);
                }
            } else if (column != indexOfNormalization && row != indexOfNormalization) {
                //add it to tempA somehow
                int modifiedRow = (row > indexOfNormalization) ? row - 1 : row;
                int modifiedColumn = (column > indexOfNormalization) ? column - 1 : column;
                subMatrix(modifiedRow, modifiedColumn) = tempMatrix(row, column);
            }
        }
    }


    //make b negative
    b = -b;

    //solve the equation tempA * x = b for x
    x = subMatrix.colPivHouseholderQr().solve(b);

    /*std::cout << "computed coefficients: " << std::endl;
    std::cout << x << std::endl;
    std::cout << "--------------" << std::endl;*/

    //put the coefficients in the coefficient map, remembering to add back in the coefficients for indexOfNormalization
    for (int i = 0; i < indexTransversal.size(); i++) {
        Index index = indexTransversal[i];
        if (i == indexOfNormalization) {
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = std::get<1>(tup);
                tempCoefficientMap[std::get<0>(tup)] = coeff;
            }
        } else {
            //look through the orbit and assign the same value to all of them using the +/- factor
            int modifiedi = (i > indexOfNormalization) ? i - 1 : i;
            auto orbit = indexOrbitData[index];
            for (auto tup : orbit) {
                double coeff = x(modifiedi).real();
                coeff *= std::get<1>(tup);
                tempCoefficientMap[std::get<0>(tup)] = coeff;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    double coeffAverage = 0;
    for (auto itr : tempCoefficientMap) {
        coeffAverage += std::abs(itr.second);
    }
    //std::cout << "coeff average " << matrixID << ": " << coeffAverage/tempCoefficientMap.size() << std::endl;
    coefficientMap2 = tempCoefficientMap;
}

double CoefficientComputer::evaluate2(const Quaternion &z) {
    //check that complex is above Y
    assert(!coefficientMap2.empty());
    auto coeffMap = coefficientMap2;

    double answer = 0;

    for (auto n : indexTransversal) {
        auto a_n = coeffMap[n];

        double term = a_n*z.getJ();

        //kappa(2*pi/A*|n|*y)
        term *= kappaApprox(2*pi/A * n.getAbs() * z.getJ());

        //exp(pi*I/A <I*n, x>)
        double cs = 0;
        for (auto tup : indexOrbitData[n]) {
            auto l = std::get<0>(tup);
            int sign = std::get<1>(tup);

            if (symClass == 'D' || symClass == 'G' || d == 1) {
                cs += sign * 2.0 * std::cos(pi/A * traceProduct(I*l.getComplex(), z.getComplex()));
            } else {
                cs += sign * 2.0 * std::sin(pi/A * traceProduct(I*l.getComplex(), z.getComplex()));
            }
        }
        term *= cs;

        answer += term;
    }
    return answer;
}


