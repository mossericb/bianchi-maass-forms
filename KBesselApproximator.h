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


class KBesselApproximator {
public:
    KBesselApproximator(double r);

    void updateRAndPrecompute(double r, double precomputeLowerBound, double precomputeUpperBound);
    void extendPrecomputedRange(double precomputeUpperBound);
    void updateR(double newR);
    double approxKBessel(const double& x);
    double computeKBessel(const double& x);
    void printTiming();
    std::vector<double> maximize();

private:
    int bits = 300;
    double fineSpacing = pow(2,-3);
    double coarseSpacing = pow(2,-3);
    int fineD = -3;
    int coarseD = -3;
    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;

    double MAX_ERROR = std::pow(10,-11);

    double r;

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
};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
