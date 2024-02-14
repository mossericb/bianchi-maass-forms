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

using std::vector, std::pair, std::tuple, std::map;

#ifndef BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
#define BIANCHI_MAASS_FORMS_TESTPOINTORBITDATA
struct TestPointOrbitData {
    complex<double> representativeComplex;
    Quaternion representativePullback;
    vector<pair<complex<double>, char>> properTranslatesModSign;
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

    vector<Index> indicesUpToM(const double M);
    tuple<vector<Index>, map<Index, vector<pair<Index, int>>>>
    indexOrbitQuotientData(vector<Index> indices, const char symClass);

private:
    int d;
    std::complex<double> theta;
    double A;
    double Y0;
};


#endif //BIANCHI_MAASS_FORMS_IMAGINARYQUADRATICINTEGERS_H
