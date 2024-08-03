//
// Created by Eric Moss on 2/3/24.
//

#ifndef BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
#define BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
#include <complex>
#include "Index.h"
#include <vector>
#include "Quaternion.h"
#include <complex>
#include <map>

using std::complex;

using std::vector, std::pair, std::map;

#ifndef BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
#define BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
struct TestPointOrbitData {
    complex<double> representativeComplex;
    Quaternion representativePullback;
    vector<pair<complex<double>, short>> properTranslatesModSign;
};
#endif

class ImaginaryQuadraticIntegers {
public:
    explicit ImaginaryQuadraticIntegers(int d);
    ImaginaryQuadraticIntegers();

    [[nodiscard]] int getd() const { return d; };
    [[nodiscard]] std::complex<double> getTheta() const { return theta; };
    [[nodiscard]] double getA() const { return A; };
    [[nodiscard]] double getY0() const { return Y0; };
    [[nodiscard]] double getVolumeOfFD() const { return volumeOfFD; };

    double weylLaw(double r);
    double eigenvalueIntervalRightEndpoint(double leftEndpoint, double numEigenvalues);

    bool isPrime(const Index& index);
    static bool isRationalPrime(long n);

    vector<Index> indicesUpToM(const double M);
    pair<vector<Index>, map<Index, vector<pair<Index, int>>>>
    indexOrbitQuotientData(vector<Index> indices, const char symClass);

private:
    int d;
    int disc;
    std::complex<double> theta;
    int normTheta;
    int twoRealTheta;
    double A;
    double Y0;
    double volumeOfFD;
};


#endif //BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
