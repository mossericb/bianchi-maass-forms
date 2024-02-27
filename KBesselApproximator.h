//
// Created by Eric Moss on 8/9/23.
//

#ifndef COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
#define COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
#include <map>
#include <string>
#include <chrono>
#include <vector>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

class KBesselApproximator {
public:
    KBesselApproximator(int bitsOfPrecision = 53);

    void setRAndPrecompute(double r, double precomputeLowerBound, double precomputeUpperBound);
    void extendPrecomputedRange(double newLowerBound, double newUpperBound);
    void setRAndClear(double newR);
    double approxKBessel(const double& x);
    double getR() { return this->r; }
    void printTiming();
    std::vector<double> maximize();

private:
    int prec;

    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;

    double r;

    double derivativeShift = pow(2,-25);
    double fineSplineSpacing = pow(2,-19);
    double coarseSplineSpacing = pow(2,-12);
    int fineSplineKnotCount = 100;
    int coarseSplineKnotCount = 100;

    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> fineSplines;
    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> coarseSplines;

    double fineIntervalWidth;
    double coarseIntervalWidth;

    void runTest();

};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
