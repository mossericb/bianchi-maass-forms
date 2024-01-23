//
// Created by Eric Moss on 12/13/23.
//

#ifndef EUCLIDEAN_BIANCHI_MAASS_FORMS_TRIGAPPROXIMATOR_H
#define EUCLIDEAN_BIANCHI_MAASS_FORMS_TRIGAPPROXIMATOR_H
#include <arb.h>
#include <vector>

class TrigApproximator {
public:
    TrigApproximator();

    double sinApprox(const double& x);
    double cosApprox(const double& x);

    double sinExact(const double& x);
    double cosExact(const double& x);

private:
    int bits = 300;
    double spacing = pow(2,-3);
    int D = -3;

    double piOverFour = 0.785398163397448309615660845819;
    double piOverTwo = 1.57079632679489661923132169163;
    double pi = 3.1415926535897932384626433;
    double threePiOverTwo = 4.71238898038468985769396507491925432629;
    double twoPi = 6.283185307179586476925286766;

    double interpolate(const double& x);

    double rand01();
    double randab(double a, double b);

    double MAX_ERROR = pow(10,-11);

    std::vector<double> precomputedValues;
};


#endif //EUCLIDEAN_BIANCHI_MAASS_FORMS_TRIGAPPROXIMATOR_H
