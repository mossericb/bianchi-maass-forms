//
// Created by Eric Moss on 8/9/23.
//

#ifndef COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
#define COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
#include <map>
#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>
#include <string>
#include <chrono>
#include <vector>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

class KBesselApproximator {
public:
    KBesselApproximator(int bitsOfPrecision);
    ~KBesselApproximator();

    void setRAndPrecompute(double r, double precomputeLowerBound, double precomputeUpperBound);
    void extendPrecomputedRange(double newLowerBound, double newUpperBound);
    void setRAndClearPrecompute(double newR);
    double approxKBessel(const double& x);
    double exactKBessel(const double& x);
    void printTiming();
    std::vector<double> maximize();

private:
    int threads;
    int prec;
    double fineSpacing = pow(2,-3);
    double coarseSpacing = pow(2,-3);
    int fineD = -3;
    int coarseD = -3;
    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;

    double MAX_ERROR = std::pow(10,-10);

    double r;

    std::vector<arb_struct> realTemp;
    std::vector<arb_struct> arbR;
    std::vector<arb_struct> besselScale;
    std::vector<arb_struct> preComputedBesselScale;
    std::vector<acb_struct> complexTemp;
    std::vector<acb_struct> acbNu;


    std::chrono::duration<double> firstCheckDuration;
    std::chrono::duration<double> nearestThreeDuration;
    std::chrono::duration<double> immediatelyComputeDuration;
    std::chrono::duration<double> computeBecausePointsNotCloseEnoughDuration;
    std::chrono::duration<double> interpolateDuration;

    std::vector<double> finePrecomputedValues;
    std::vector<double> coarsePrecomputedValues;

    double fineInterpolate(const double& x);
    double coarseInterpolate(const double& x);
    std::string dtos(const double& x);

    struct Point {double x; double y; Point(double x_, double y_) {x = x_; y = y_;}};

    double neville(std::vector<Point>& points, short i, short j, const double &x);
    double rand01();
    double randab(double a, double b);
    double fineLagrange(const Point& p1, const Point& p2, const Point& p3, const double& x);
    double coarseLagrange(const Point& p1, const Point& p2, const Point& p3, const double& x);

    std::vector<boost::math::chebyshev_transform<double>> fineChebyshevApproximators;
    std::vector<boost::math::chebyshev_transform<double>> coarseChebyshevApproximators;
    double applyChebyshev(const double& x);

    int coarseChebyshevIndex(const double& x);
    int fineChebyshevIndex(const double& x);

    double fineSplineSpacing = pow(2,-17);
    double coarseSplineSpacing = pow(2,-12);
    int fineSplineKnotCount = 100;
    int coarseSplineKnotCount = 100;
    std::vector<double> fineSplinePrecompute;
    std::vector<double> coarseSplinePrecompute;

    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> fineSplines;
    std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> coarseSplines;

    double fineIntervalWidth;
    double coarseIntervalWidth;

    std::vector<double> fineSplineEndpoints;
    std::vector<double> coarseSplineEndpoints;

};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
