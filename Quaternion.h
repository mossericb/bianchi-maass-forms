//
// Created by Eric Moss on 6/26/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_QUATERNION_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_QUATERNION_H

#include <iostream>
#include <complex>
#include "Auxiliary.h"
#include "SL2C.h"
#include <vector>

class Quaternion {
public:
    Quaternion(double x, double y, double z, double w);
    Quaternion(std::complex<double> z);
    Quaternion(const Quaternion &q);
    Quaternion();
    double x, y, z, w;

    double abs() const;
    double squareAbs() const;
    [[nodiscard]] Quaternion inverse() const;
    [[nodiscard]] Quaternion conj() const;
    std::complex<double> getComplex() const;
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

    const static SL2C S;
    const static std::vector<std::complex<double>> d19Alphas;
    const static std::vector<SL2C> d19MAlphas;
    const static std::vector<std::complex<double>> d19Translators;

    const static std::vector<std::complex<double>> d43Alphas;
    const static std::vector<SL2C> d43MAlphas;
    const static std::vector<std::complex<double>> d43Translators;

    const static std::vector<std::complex<double>> d67Alphas;
    const static std::vector<SL2C> d67MAlphas;
    const static std::vector<std::complex<double>> d67Translators;

    const static std::vector<std::complex<double>> d163Alphas;
    const static std::vector<SL2C> d163MAlphas;
    const static std::vector<std::complex<double>> d163Translators;

private:
    void reduceModT();
    bool reduceModInversion(int d);
    void reduceThetaGeneral(int d);
    void reduceUnits(int d);

    void vectorReduce(const std::complex<double> v);
    std::complex<double> vectorReduce(const std::complex<double> &u, const std::complex<double> &v);
    static std::complex<double> getTheta(int d);
    static int mod(int x, int n);

    /*[[nodiscard]] bool isInUnitsFundamentalDomain(int d) const;
    [[nodiscard]] bool isInParabolicFundamentalDomain(int d) const;*/
    [[nodiscard]] bool isInInversionFundamentalDomain() const;
    /*[[nodiscard]] bool eisensteinIsInNEDomain() const;
    [[nodiscard]] bool eisensteinIsInNWDomain() const;
    [[nodiscard]] bool eisensteinIsInSEDomain() const;
    [[nodiscard]] bool eisensteinIsInSWDomain() const;*/
};




#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_QUATERNION_H
