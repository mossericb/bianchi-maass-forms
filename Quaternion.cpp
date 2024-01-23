//
// Created by Eric Moss on 6/26/23.
//

#include "Quaternion.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

const SL2C Quaternion::S = {0,-1,1,0};
const std::vector<std::complex<double>> Quaternion::d19Alphas = {{0,0},
                                                                 {0.25, getTheta(19).imag()/2},
                                                                 {-0.25, getTheta(19).imag()/2}};
const std::vector<SL2C> Quaternion::d19MAlphas = {{0,-1,1,0},
                                                  {-getTheta(19)+ 1.0, -2, -2, getTheta(19)},
                                                  {-getTheta(19), -2, -2, getTheta(19) - 1.0}};
const std::vector<std::complex<double>> Quaternion::d19Translators = {{0,0},
                                                                      {1,0},
                                                                      {-1,0},
                                                                      getTheta(19),
                                                                      -getTheta(19),
                                                                      getTheta(19)-1.0,
                                                                      -getTheta(19)+1.0};



Quaternion::Quaternion(double x, double y, double z, double w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

Quaternion::Quaternion(const Quaternion &q) {
    this->x = q.x;
    this->y = q.y;
    this->z = q.z;
    this->w = q.w;
}

std::ostream &operator<<(std::ostream &strm, const Quaternion &q) {
    std::ostringstream out;
    out << std::setprecision(16);
    out << q.x << " + ";
    out << q.y << "*i + ";
    out << q.z << "*j + ";
    out << q.w << "*ij";
    return strm << std::move(out).str();
}

Quaternion operator*(const Quaternion &q1, const Quaternion &q2) {
    double c1 = q1.x*q2.x - q1.y*q2.y - q1.z*q2.z - q1.w*q2.w;
    double c2 = q1.x*q2.y + q1.y*q2.x + q1.z*q2.w - q1.w*q2.z;
    double c3 = q1.x*q2.z - q1.y*q2.w + q1.z*q2.x + q1.w*q2.y;
    double c4 = q1.x*q2.w + q1.y*q2.z - q1.z*q2.y + q1.w*q2.x;
    return {c1, c2, c3, c4};
}

Quaternion operator+(const Quaternion &lhs, const Quaternion& rhs) {
    Quaternion answer = Quaternion(lhs);
    answer += rhs;
    return answer;
}

Quaternion operator-(const Quaternion &lhs, const Quaternion &rhs) {
    Quaternion answer = Quaternion(lhs);
    answer -= rhs;
    return answer;
}

Quaternion operator/(const Quaternion &q1, const Quaternion &q2) {
    return q1*(q2.inverse());
}

Quaternion operator*(const double &c, const Quaternion& q) {
    return {c*q.x, c*q.y, c*q.z, c*q.w};
}

std::complex<double> Quaternion::getComplex() const {
    return {this->x, this->y};
}

double Quaternion::squareAbs() const {
    return pow(this->x,2) + pow(this->y,2) + pow(this->z,2) + pow(this->w,2);
}

double Quaternion::abs() const {
    return sqrt(squareAbs());
}

Quaternion Quaternion::inverse() const {
    return (1/squareAbs())*(this->conj());
}

Quaternion Quaternion::conj() const {
    return {this->x, -this->y, -this->z, -this->w};
}

Quaternion operator-(const Quaternion &q) {
    return {-q.x, -q.y, -q.z, -q.w};
}

void Quaternion::operator+=(const Quaternion& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
}

void Quaternion::operator-=(const Quaternion& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
}

void Quaternion::operator*=(const Quaternion& rhs) {
    *this = (*this) * rhs;
}

void Quaternion::reduceModT() {
    this->vectorReduce(std::complex<double> {1,0});
}

bool Quaternion::reduceModInversion(int d) {
    if (d < 19) {
        double temp = x*x + y*y + z*z;
        if (temp < 1) {
            x = -x;
            *this = (1/temp) * (*this);
            return true;
        } else {
            return false;
        }
    } else if (d == 19) {
        /*[ 0  1]  [-w + 1     -2]  [   -w    -2]
        [-1  0], [    -2      w], [   -2    w - 1]*/
        //compute the image wrt each of these matrices
        //compare the j part of each 3 images to the j part of doing nothing
        //the answer is the one with greatest j part
        Quaternion largestSoFar = Quaternion(*this);
        for (int i = 0; i < d19MAlphas.size(); i++) {
            auto alpha = d19Alphas[i];
            auto M_alpha = d19MAlphas[i];
            auto tempThis = Quaternion(*this);
            auto tempThisComplex = this->getComplex();

            //find nearest element of OK to *this and translate by that
            double minDistance = 100;
            std::complex<double> toTranslate;

            for (auto n : d19Translators) {
                if (std::abs((tempThisComplex - n) - alpha) < minDistance) {
                    minDistance = std::abs((tempThisComplex - n) - alpha);
                    toTranslate = n;
                }
            }

            auto translated = tempThisComplex - toTranslate;
            tempThis.x = translated.real();
            tempThis.y = translated.imag();

            tempThis = M_alpha * tempThis;

            if (tempThis.getJ() > largestSoFar.getJ()) {
                largestSoFar = tempThis;
            }
        }
        if (*this == largestSoFar) {
            return false;
        } else {
            *this = largestSoFar;
            return true;
        }

    } else {
        //throw(std::invalid_argument("not implemented for this d"));
    }
}

void Quaternion::reduceThetaGeneral(int d) {
    if (d == 3) {
        double minDistance = 100;

        auto theta = getTheta(d);
        std::complex<double> othertheta = {-theta.real(), theta.imag()};

        std::complex<double> temp1 = this->getComplex();

        std::complex<double> translated = {0,0};

        int horizbound = std::ceil(std::abs(temp1.real()));
        int verbound = std::ceil(std::abs(temp1.imag()/theta.imag()));
        for (int a = -horizbound; a <= horizbound; a++) {
            for (int b = -verbound; b <= verbound; b++) {
                std::complex<double> latticePoint = (double)a + (double)b*theta;
                double thisDistance = std::abs(temp1 - latticePoint);
                if (thisDistance < minDistance) {
                    minDistance = thisDistance;
                    translated = latticePoint;
                }
            }
        }
        x = temp1.real();
        y = temp1.imag();
    } else {
        auto theta = getTheta(d);
        auto complex = this->getComplex();
        //double n = -std::floor((y + theta.imag()/2)/(theta.imag()));
        double n = -std::floor(y/theta.imag() + 0.5);
        complex += n*theta;
        x = complex.real();
        y = complex.imag();
    }
}

void Quaternion::reduceUnits(int d) {
    if (d == 1) {
        if (x < 0) {
            Quaternion gaussianUnit = {0,1,0,0};
            *this = gaussianUnit * (*this) * gaussianUnit;
        }
    } else if (d == 3) {
        Quaternion quatTheta = Quaternion(getTheta(d).real(), getTheta(d).imag(), 0, 0);

        while (x < 0 || y < -1.0/sqrt(3)*x) {
            *this = quatTheta * (*this) * quatTheta;
        }
    }
}

double Quaternion::getJ() const {
    return this->z;
}

void Quaternion::reduce(int d) {
    reduceModInversion(d);
    do {
        reduceThetaGeneral(d);
        reduceModT();
        reduceUnits(d);
    } while (reduceModInversion(d));
}

bool Quaternion::isInInversionFundamentalDomain() const {
    assert(z > 0);
    return this->squareAbs() >= 1.0;
}

void Quaternion::setEqual(Quaternion &q) {
    this->x = q.x;
    this->y = q.y;
    this->z = q.z;
    this->w = q.w;
}

int Quaternion::mod(int a, int n) {
    int answer = a % n;
    if (answer < 0) {
        return answer + n;
    } else {
        return answer;
    }
}

std::complex<double> Quaternion::getTheta(int d) {
    std::complex<double> theta;
    if (mod(-d,4) == 1) {
        theta = std::complex<double> {1.0/2, sqrt(d)/2};
    } else if (mod(-d,4) == 0) {
        throw(std::invalid_argument("d should be squarefree"));
    } else {
        theta = std::complex<double> {0, sqrt(d)};
    }
    return theta;
}

/**
 * This reduces the Quaternion<T> via the action of addition/subtraction by v.
 * @tparam T
 * @param v
 */
void Quaternion::vectorReduce(const std::complex<double> v) {
    std::complex<double> complex = this->getComplex();
    double modifiedScalarProjection = complex.real()*v.real() + complex.imag()*v.imag();
    modifiedScalarProjection /= pow(v.real(),2) + pow(v.imag(),2);
    modifiedScalarProjection = floor(modifiedScalarProjection + 1.0/2);
    this->x -= modifiedScalarProjection*v.real();
    this->y -= modifiedScalarProjection*v.imag();
}

void Quaternion::operator/=(const Quaternion &rhs) {
    *this = (*this)*(rhs.inverse());
}

Quaternion operator*(const SL2C &gamma, const Quaternion &q) {
    Quaternion a = Quaternion(gamma.a);
    Quaternion b = Quaternion(gamma.b);
    Quaternion c = Quaternion(gamma.c);
    Quaternion d = Quaternion(gamma.d);

    Quaternion answer = a*q + b;
    answer /= (c*q + d);
    //answer.w = 0; //mitigate error propagation
    return answer;
}

Quaternion::Quaternion(std::complex<double> z) {
    this->x = z.real();
    this->y = z.imag();
    this->z = 0.0;
    this->w = 0.0;
}

bool operator==(const Quaternion& q1, const Quaternion& q2) {
    if ((q1.x == q2.x) && (q1.y == q2.y) && (q1.z == q2.z) && (q1.w == q2.w)) {
        return true;
    } else {
        return false;
    }
}

std::complex<double> Quaternion::vectorReduce(const std::complex<double> &u, const std::complex<double> &v) {
    double modifiedScalarProjection = u.real()*v.real() + u.imag()*v.imag();
    modifiedScalarProjection /= pow(v.real(),2) + pow(v.imag(),2);
    modifiedScalarProjection = floor(modifiedScalarProjection + 1.0/2);
    return u - modifiedScalarProjection*v;
}

Quaternion::Quaternion() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}


//
// Created by Eric Moss on 8/17/23.
//
