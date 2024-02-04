#ifndef BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H
#define BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H

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
#include "Auxilliary.h"
#include "KBesselApproximator.h"
#include "ImaginaryQuadraticIntegers.h"

#include <eigen3/Eigen/Dense>

//Containers
using std::vector, std::unordered_map, std::tuple, std::map;

//Math
using std::complex, Eigen::MatrixXd;

//IO
using std::ofstream, std::string;

struct TestPointOrbitData {
    complex<double> representativeComplex;
    Quaternion representativePullback;
    vector<tuple<complex<double>, char>> properTranslatesModSign;
};

class BianchiMaassSearch {
public:
    BianchiMaassSearch(int d, int D, char symClass);
    ~BianchiMaassSearch();

    /***************************************************
     * These functions are for SEARCHING for eigenvalues. They only
     * spit out intervals where they think an eigenvalue lives.
     ***************************************************/

    void searchForEigenvalues(const double leftR, const double rightR);
    void clearSearchData();

    /***************************************************
     * These functions NARROW an already small interval containing a likely eigenvalue.
     ***************************************************/

    void narrowLikelyInterval(const double leftR, const double rightR);
private:

    /***************************************************
     * Private members used in ALL calculations.
     ***************************************************/
    int d;
    int D;
    char symClass;

    ImaginaryQuadraticIntegers Od;

    complex<double> theta;
    double A;
    double Y0;

    ofstream outputFile;

    double pi = 3.141592653589793238462643383279502884197;
    complex<double> I = {0,1};

    double tolerance;
    double truncation;
    double volumeOfFD;

    /***************************************************
     * Private methods used in ALL calculations.
     ***************************************************/

    double computeM0General(const double r);
    double computeMYGeneral(const double M0, const double Y);

    vector<TestPointOrbitData> getPointPullbackOrbits(const Index &m, const double Y, const double MY);

    double traceProduct(const complex<double> &z, const complex<double> &w);
    bool areDifferentSign(const double &a, const double &b);

    /***************************************************
     * Private members used in SEARCH calculations.
     ***************************************************/

    enum MatrixID { Y1, Y2 };

    map<Index, vector<TestPointOrbitData>> searchMToY1TestPointOrbits;
    map<Index, vector<TestPointOrbitData>> searchMToY2TestPointOrbits;

    int searchTestPointCount;

    vector<MatrixXd> searchMatrices;
    vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> searchMatrixSolutions;
    vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> searchMatrixG;
    vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> searchMatrixB;

    int searchPossibleSignChanges;

    KBesselApproximator searchK;

    int searchIndexOfNormalization;

    double searchM0;

    vector<Index> searchIndicesM0;
    vector<Index> searchIndexTransversal;
    unordered_map<Index, vector<tuple<Index, int>>> searchIndexOrbitData;
    unordered_map<Index, vector<tuple<Index, int>>> searchIndexOrbitDataModMinusOne;

    double searchY1;
    double searchY2;

    double searchMaxYStar;

    /***************************************************
     * Private methods used in SEARCH calculations.
     ***************************************************/

    void recursiveSearchForEigenvalues(const double leftR, const double rightR,
                                       vector<double>& leftG,
                                       vector<double>& rightG);
    void computeIndexData();
    void computeTestPointData();
    void populateMatrix(MatrixID matrixId);
    void solveMatrix(MatrixID matrixId);
    vector<double> getGVector();
    double computeEntry(const Index& m, const Index& n, MatrixID matrixId);

    void secantSearch(double r1, double r2);

    int countSignChanges(const vector<double> &v1, const vector<double> &v2);

    bool signChangeVectorIsIncreasing(vector<int> &v);

    double findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1);
};


#endif //BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H
