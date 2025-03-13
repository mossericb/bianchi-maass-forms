#pragma once

#include <iostream>
#include <complex>
#include "Auxiliary.h"
#include "SL2C.h"
#include <vector>
#include <map>

using std::vector, std::map, std::complex;

class Quaternion {
public:
    Quaternion(double x, double y, double z, double w);
    Quaternion(complex<double> z);
    Quaternion(const Quaternion &q);
    Quaternion();
    double x, y, z, w;

    double abs() const;
    double squareAbs() const;
    [[nodiscard]] Quaternion inverse() const;
    [[nodiscard]] Quaternion conj() const;
    complex<double> getComplex() const;
    double getJ() const;

    void setEqual(Quaternion &q);
    void reduce(int d);

    friend Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
    friend Quaternion operator+(const Quaternion& lhs, const Quaternion& rhs);
    friend Quaternion operator-(const Quaternion& lhs, const Quaternion& rhs);
    friend Quaternion operator/(const Quaternion& q1, const Quaternion& q2);
    friend Quaternion operator-(const Quaternion& q);
    friend bool operator==(const Quaternion& q1, const Quaternion& q2);
    void operator+=(const Quaternion& rhs);
    void operator-=(const Quaternion& rhs);
    void operator*=(const Quaternion& rhs);
    void operator/=(const Quaternion& rhs);
    friend std::ostream& operator<<(std::ostream&, const Quaternion&);
    friend Quaternion operator*(const double &c, const Quaternion& q);
    friend Quaternion operator*(const SL2C &gamma, const Quaternion& q);

private:
    void reduceModT();
    bool reduceModInversion(int d);
    void reduceThetaGeneral(int d);
    void reduceUnits(int d);

    void vectorReduce(const complex<double> v);
    complex<double> vectorReduce(const complex<double> &u, const complex<double> &v);
    static std::complex<double> getTheta(int d);
    static int mod(int x, int n);

    const static map<int, vector<complex<double>>> alphas;
    const static map<int, vector<SL2C>> MAlphas;
    const static map<int, vector<complex<double>>> translators;

    const static SL2C S;
    const static vector<complex<double>> d19Alphas;
    const static vector<SL2C> d19MAlphas;
    const static vector<complex<double>> d19Translators;

    const static vector<complex<double>> d43Alphas;
    const static vector<SL2C> d43MAlphas;
    const static vector<complex<double>> d43Translators;

    const static vector<complex<double>> d67Alphas;
    const static vector<SL2C> d67MAlphas;
    const static vector<complex<double>> d67Translators;

    const static vector<complex<double>> d163Alphas;
    const static vector<SL2C> d163MAlphas;
    const static vector<complex<double>> d163Translators;

    [[nodiscard]] bool isInInversionFundamentalDomain() const;
};