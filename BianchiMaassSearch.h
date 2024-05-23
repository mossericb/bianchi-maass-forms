#ifndef BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H
#define BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H

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

    double EIGENVALUE_INTERVAL_CUTOFF = pow(2,-16);

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

    bool computeAllConditionNumbers;
    double conditionGoal;
    double argminY1;
    double argminY2;

    Auxiliary Aux;

    map<double, KBessel> rToKBess;
    /***************************************************
     * Private methods used in SEARCH calculations.
     ***************************************************/

    int findMaxFileNumber(const string &directory, const string &prefix);
    void createOutputDirectory(const std::string& directory);
    /*void recursiveSearchForEigenvalues(const double leftR, const double rightR,
                                       vector<double>& leftG,
                                       vector<double>& rightG);
    */
    vector<pair<double, double>> conditionedSearchForEigenvalues(const double leftR, const double rightR);
    void computeIndexData();
    void computeTestPointData();
    vector<double> computeAndGetGVector();

    double minBess(KBessel &K, const vector<Index> &indexTransversal, double Y);

    MatrixXd produceMatrix(const vector<Index> &indexTransversal,
                           map<Index, vector<TestPointOrbitData>> &mToTestPointData,
                           map<Index, vector<pair<Index, int>>> &ntoIndexOrbitData,
                           double Y,
                           KBessel &K);
    pair<MatrixXd, double> solveMatrix(const MatrixXd &matrix, const vector<Index> &indexTransversal,
                         const int indexOfNormalization);
    vector<double> getGVector();
    double computeEntry(const Index &m, const Index &n, KBessel &K,
                        const vector<TestPointOrbitData> &mTestPointOrbits, const double Y,
                        const vector<pair<Index, int>> &nIndexOrbitDataModSign);

    void secantMethod(double r1, double r2);

    int countSignChanges(const vector<double> &v1, const vector<double> &v2);

    bool signChangeVectorIsIncreasing(vector<int> &v);

    double findZeroOfLinearInterpolation(double x0, double y0, double x1, double y1);

    double heckeCheck(map<Index, double>& coeffMap);
    bool sufficientSignChanges(const vector<double>& v1, const vector<double>& v2);

    vector<double> mergeToVector(const MatrixXd& v1, const MatrixXd& v2);
};


#endif //BIANCHI_MAASS_FORMS_BIANCHIMAASSSEARCH_H
