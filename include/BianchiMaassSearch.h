#pragma once

#include <iostream>
#include <complex>
#include <map>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <utility>
#include <omp.h>

#include "Index.h"
#include "Quaternion.h"
#include "Auxiliary.h"
#include "KBessel.h"
#include "ImaginaryQuadraticIntegers.h"

#include <eigen3/Eigen/Dense>

//Containers
using std::vector, std::unordered_map, std::tuple, std::map, std::pair;

//Math
using std::complex, Eigen::MatrixXd;

//IO
using std::ofstream, std::string;



class BianchiMaassSearch {
public:
    BianchiMaassSearch(string mode, int d, double D, char symClass);

    void coarseSearchForEigenvalues(const double leftR, const double rightR);
    void mediumSearchForEigenvalues();
    void fineSearchForEigenvalues();
    void testForModularity();
    void testConjectures();
    void sandbox(const double leftR, const double rightR);
    void sandbox2(const double r);
    void sandbox3(const double r);
    void sandbox4(const double r);

private:

    /***************************************************
     * Private members used in ALL calculations.
     ***************************************************/
    int d;
    double D;
    char symClass;
    string mode;
    int finalPrecision;

    ImaginaryQuadraticIntegers Od;

    complex<double> theta;
    double A;
    double Y0;

    ofstream coarseOutputFile;
    ofstream mediumOutputFile;
    ofstream fineOutputFile;

    double coarseComplete = -0.4;
    double mediumComplete = -0.4;
    double fineComplete = -0.4;

    int secondsToComputeMatrix = 400;

    double pi = 3.141592653589793238462643383279502884197;
    double twoPiOverA;
    double piOverA;
    complex<double> I = {0,1};

    double tolerance;
    double truncation;
    double volumeOfFD;
    map<int, vector<double>> dToDMap;


    double EIGENVALUE_INTERVAL_CUTOFF = pow(2,-16);
    double nanosecondsPerTerm = 5;
    int maxSecondsPerMatrix = 180;
    bool autoPrecision = true;

    /***************************************************
     * Private methods used in ALL calculations.
     ***************************************************/

    void computeMaximumD(double r, int timeLimitSeconds);

    double computeM0General(double r);
    double computeMYGeneral(double M0, double Y);

    int getIndexOfNormalization(const vector<Index> & indexTransversal);

    vector<TestPointOrbitData> getPointPullbackOrbits(const Index &m, double Y, double MY);
    unsigned long long int countPointPullbackOrbits(const Index &m, double Y, double MY);

    double traceProduct(complex<double> z, complex<double> w);
    bool areDifferentSign(double a, double b);

    double evaluate(const Quaternion& z, KBessel& K, const vector<Index>& indexTransversal, const vector<double>& coeffs);

    Auxiliary Aux;

    map<int, vector<Index>> dToPrimes;

    void setUpOutputLogFiles();
    static void createOutputDirectory(const std::string& directory);
    static bool isFileEmpty(const string& filename);

    vector<pair<double,double>> getIntervalsForCoarseSearch(double startR, double endR);
    vector<pair<double,double>> getIntervalsForMediumSearch();
    vector<pair<double,double>> getIntervalsForFineSearch();
    void saveSearchResult(vector<double> spectralParameters, vector<Index> &indexTransversal,
                          vector<double> &coeffs, double Y);

    bool possiblyContainsEigenvalue(double leftR, double rightR, KBessel *leftRK, KBessel *rightRK);
    void medianIllinoisSearch(double leftR, double rightR);

    double minBess(KBessel *K, const vector<Index> &indexTransversal, double Y);
    double computeWellConditionedY(KBessel *K, double r, double M0, vector<Index> &indexTransversal);
    double computeWellConditionedY(KBessel *leftRK, KBessel *rightRK, double leftR, double rightR, double M0, const vector<Index>& indexTransversal);
    pair<double, double> computeTwoWellConditionedY(KBessel *leftRK, KBessel *rightRK, double leftR, double rightR, double M0, const vector<Index>& indexTransversal);

    MatrixXd produceMatrix(const vector<Index> &indexTransversal,
                           map<Index, vector<TestPointOrbitData>> &mToTestPointData,
                           map<Index, vector<pair<Index, int>>> &ntoIndexOrbitData,
                           double Y,
                           KBessel &K);
    pair<MatrixXd, double> solveMatrix(const MatrixXd &matrix, const vector<Index> &indexTransversal,
                         int indexOfNormalization);

    double computeEntry(const Index &m, const Index &n, KBessel &K,
                        const vector<TestPointOrbitData> &mTestPointOrbits, double Y,
                        const vector<pair<Index, int>> &nIndexOrbitDataModSign);

    int countSignChanges(const vector<double> &v1, const vector<double> &v2);
};
