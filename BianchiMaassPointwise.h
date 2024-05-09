//
// Created by Eric Moss on 2/3/24.
//

#ifndef BIANCHI_MAASS_FORMS_BIANCHIMAASSPOINTWISE_H
#define BIANCHI_MAASS_FORMS_BIANCHIMAASSPOINTWISE_H

#include <iostream>
#include <complex>
#include <map>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <omp.h>

#include "Index.h"
#include "Quaternion.h"
#include "Auxiliary.h"
#include "KBessel.h"
#include "ImaginaryQuadraticIntegers.h"

#include <eigen3/Eigen/Dense>

//Containers
using std::vector, std::unordered_map, std::map;

//Math
using std::complex, Eigen::MatrixXd;

//IO
using std::ofstream, std::string;

#ifndef BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
#define BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
struct TestPointOrbitData {
    complex<double> representativeComplex;
    Quaternion representativePullback;
    vector<tuple<complex<double>, char>> properTranslatesModSign;
};
#endif

class BianchiMaassPointwise {
public:
    BianchiMaassPointwise(int d, int D, char symClass);
    ~BianchiMaassPointwise();

    void checkSingleEigenvalue(const double r, const double Y = 0);
    unordered_map<Index, double> computeFourierCoefficients(const double r, const double Y = 0);
    double evaluate(const double r, const Quaternion& z, const double Y = 0);
    void checkSatoTate(const double r, const double Y = 0);
    bool checkRamanujanPetersson(const double r, const double Y = 0);
    void produceLFunction(const double r, const double Y = 0);
    void clearData();

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

    static constexpr double PI = 3.141592653589793238462643383279502884197;
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
     * Private members used in POINTWISE calculations.
     ***************************************************/

    double Y;
    double M0;
    double MY;

    MatrixXd matrix;
    map<Index, double> coefficientMap;

    KBessel K;
    Auxiliary Aux;

    double maxYStar;

    int testPointCount;

    int indexOfNormalization;

    vector<TestPointOrbitData> testPointOrbits;

    vector<Index> indicesM0;
    vector<Index> indexTransversal;
    vector<Index> indicesMY;
    vector<Index> indexTransversalMY;
    unordered_map<Index, vector<pair<Index, int>>> indexOrbitData;
    unordered_map<Index, vector<pair<Index, int>>> indexOrbitDataModMinusOne;

    double rand01();

    /***************************************************
     * Private methods used in POINTWISE calculations.
     ***************************************************/

    void populateMatrix();
    double computeEntry(const Index& m, const Index& n);
    void computeIndexData();
    void solveMatrix();
    double computeK(int d, int D, char symClass);
};


#endif //BIANCHI_MAASS_FORMS_BIANCHIMAASSPOINTWISE_H
