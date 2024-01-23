//
// Created by Eric Moss on 6/27/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_COEFFICIENTCOMPUTER_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_COEFFICIENTCOMPUTER_H

#include <iostream>
#include <complex>
#include <map>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <fstream>
#include <omp.h>

#include "Index.h"
#include "Quaternion.h"
#include "EricNT.h"
#include "KBesselApproximator.h"
#include "TrigApproximator.h"
#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Core>

struct TestPointOrbitData {
    std::complex<double> representativeComplex;
    Quaternion representativePullback;
    std::vector<std::tuple<std::complex<double>, char>> properTranslatesModSign;
};

class CoefficientComputer {
public:
    CoefficientComputer(int d, int D, char symClass);
    ~CoefficientComputer();

    double evaluate(const Quaternion &z);
    void searchForEigenvalues(double leftR, double rightR);
    void checkSingleEigenvalue(double r);
    void computeFourierCoefficients(double r);
    void checkSatoTate(double r);
    void checkRamanujanPetersson(double r);
    void produceLFunction(double r);
    void recursiveSearchForEigenvalues(double localLeftR, double localRightR);
    void printTime();
    void printCoefficients();
    void computeMoreCoefficients();
    void secantSearch(double r1, double r2);

    std::map<double, long long> histogram;
    std::vector<std::vector<double>> probableEigenvalues;

    //test zone
    void checkSingleEigenvalue2(double r);
    //
private:
    int d;
    int D;
    char symClass;
    double Y1Parameter;

    std::ofstream outputFile;

    double A;
    double M0;
    double Y0;
    double MY;
    double Y1;
    double Y2;
    double pi = 3.141592653589793238462643383279502884197;
    std::complex<double> theta;
    std::complex<double> I = {0,1};
    int indexOfNormalization;
    double tolerance;
    double truncation;
    double volumeOfFD;
    double leftR;
    double rightR;

    KBesselApproximator* K;
    TrigApproximator* T;

    static double SUBINTERVAL_SEARCH_WIDTH;
    static double LOWER_BOUND_FOR_EIGENVALUE;
    static double MINIMUM_SPACING;
    static double SCALE_WEYL_LAW;
    static double BINARY_SEARCH_GRID_SPACING;
    static double RATIO_OF_Y1_TO_Y0;
    static double RATIO_OF_Y2_TO_Y0;
    static double EIGENVALUE_INTERVAL_WIDTH_CUTOFF;
    static double CONDITION_CUTOFF;

    std::vector<Eigen::MatrixXd> matrices;
    std::vector<double> conditionNumbers;
    std::vector<std::unordered_map<Index, double>> coefficientMaps;

    Eigen::MatrixXd matrixR1Y1;
    Eigen::MatrixXd matrixR2Y1;
    Eigen::MatrixXd matrixR1Y2;
    Eigen::MatrixXd matrixR2Y2;

    std::unordered_map<Index, double> coefficientMap;
    std::vector<double> YGrid;
    std::unordered_map<double, std::vector<double>> rToConditionMap;

    struct heightConditionPair {
        double Y;
        double condition;
    };


    /*Store index and test point data.*/
    std::vector<Index> indicesM0;
    std::vector<Index> indexTransversal;
    std::vector<Index> indicesMY;
    std::vector<Index> indexTransversalMY;
    std::unordered_map<Index, std::vector<std::tuple<Index, int>>> indexOrbitData;
    std::unordered_map<Index, std::vector<std::tuple<Index, int>>> indexOrbitDataModMinusOne;
    std::unordered_map<Index, std::vector<std::vector<Quaternion>>> mToTestPointPairsY1;
    std::unordered_map<Index, std::vector<std::vector<Quaternion>>> mToTestPointPairsY2;

    //test zone
    std::vector<TestPointOrbitData> getPointPullbackOrbit(const Index& m, const double& Y);
    std::vector<TestPointOrbitData> testPointOrbits;
    double Y;
    int numberOfPoints;
    Eigen::MatrixXd matrix2;
    double computeEntry2(const Index &m, const Index &n);
    void solveMatrixAndComputeCoefficientMap2();
    std::unordered_map<Index, double> coefficientMap2;
    double evaluate2(const Quaternion &z);
    //test zone

    std::vector<std::vector<Quaternion>> getPointAndPullbackPairsY1(const Index& m);
    std::vector<std::vector<Quaternion>> getPointAndPullbackPairsY2(const Index& m);

    /*Helper functions for initialization.*/
    void computeA();
    void computeTheta();
    void computeY0();
    void computeY1();
    void computeY2();
    void computeM0();
    void computeMY();
    void computeY1Y2();
    double computeMTime;
    //void setKnownEigenvalue();
    void computeIndexData();
    void computeIndexTransversal();
    double computeIndexDataTime;
    void computeTestPointData();
    double maxYStar;
    double computeTestPointDataTime;
    std::string heckeCheck();
    double heckeCheck(short matrixID);
    void adjustY1();
    void adjustY2();
    void refreshKBesselPrecompute(short matrixID);
    short Y1R1 = 0;
    short Y1R2 = 1;
    short Y2R1 = 2;
    short Y2R2 = 3;


    /*Functions for use after initialization.*/
    void populateMatrix(short matrixID);
    void populateMatrix2(double Y);
    double populateMatrixTime;
    void solveMatrixAndComputeCoefficientMap(short matrixID);
    int solveComputeZerosVectorsCountSignChanges();
    double solveMatrixAndComputeCoefficientMapTime;
    std::vector<double> setNewRAndComputeZeroesVector(double newR);
    double computeEntry(short matrixID, const Index &m, const Index &n);

    double kappaExact(const double &x);
    double kappaApprox(const double &x);

    int countSignChanges(const std::vector<double> &v1, const std::vector<double> &v2);

    bool signChangeVectorIsIncreasing(std::vector<int> &v);
    void nudgeYAndComputeMatrix();

    double traceProduct(const std::complex<double> &z, const std::complex<double> &w);
    double findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1);
    bool areDifferentSign(const double &a, const double &b);

    double rand01();

    //These are parameters for pointwise computation, meaning we're not searching for eigenvalues
    //You give me a spectral parameter, I give you the Fourier coefficients


    void extendedPopulateMatrix(Index& m);
    std::vector<std::vector<Quaternion>> pointwisePointsAndPullbacks(Index m, double Y);
    unsigned long long int pointwisePointCount(Index m, double Y);
    std::vector<double> pointwiseGenerateMatrixRow(const Index& m, double Y, const std::vector<std::vector<Quaternion>>& pointsAndPullbacks);
};




#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_COEFFICIENTCOMPUTER_H
