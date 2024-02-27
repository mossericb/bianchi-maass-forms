//
// Created by Eric Moss on 8/9/23.
//

#include "KBesselApproximator.h"
#include <flint/flint.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include "omp.h"
#include "KBesselExact.h"

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl

/*
 * Idea: Instead of dividing the range into two parts: coarse and fine, separate it into many more spacings
 * Could be like, 10 different levels of fineness, the finest being a point every 2^-19, and the coarsest being 2^-9
 * Could be something even more than that, something more continuous.
 * Either way, it should be something dynamic, so that once I start using large spectral parameters, it still works.
 *
 *
 */

KBesselApproximator::KBesselApproximator(int bitsOfPrecision) {
    prec = bitsOfPrecision;
    r = 0;

    fineIntervalWidth = fineSplineSpacing * (fineSplineKnotCount - 1);
    coarseIntervalWidth = coarseSplineSpacing * (coarseSplineKnotCount - 1);

    precomputedRegionLeftBound = r;
    precomputedRegionRightBound = r;
}

double KBesselApproximator::approxKBessel(const double& x) {
    if (x < precomputedRegionLeftBound || x > precomputedRegionRightBound) {
        watch(x);
        watch(r);
        watch(precomputedRegionLeftBound);
        watch(precomputedRegionRightBound);
        throw std::invalid_argument("KBessel approximation out of bounds");
    }

    if (x < r) {
        int leftIndex = std::floor((r - x)/fineIntervalWidth);
        return fineSplines[leftIndex](x);
    } else {
        int leftIndex = std::floor((x - r)/coarseIntervalWidth);
        return coarseSplines[leftIndex](x);
    }
}

void KBesselApproximator::setRAndPrecompute(double newR, double precomputeLowerBound, double precomputeUpperBound) {
    setRAndClear(newR);

    extendPrecomputedRange(precomputeLowerBound, precomputeUpperBound);
}

void KBesselApproximator::setRAndClear(double newR) {
    this->r = newR;

    fineSplines.clear();
    coarseSplines.clear();

    precomputedRegionLeftBound = r;
    precomputedRegionRightBound = r;
}

std::vector<double> KBesselApproximator::maximize() {
    double start = r;
    double left = r - 0.01;

    KBesselExact K = KBesselExact(r);
    double startValue = K.exactKBessel(start);
    double leftValue = K.exactKBessel(left);

    double shift = 0;
    if (leftValue > startValue) {
        shift = -0.1;
    } else {
        shift = 0.1;
    }

    double here = start;
    double hereValue = startValue;
    for (int a = 0; a <= 5; a++) {
        double next = here + shift/std::pow(10,a);
        double nextValue = K.exactKBessel(next);
        while (nextValue > hereValue) {
            if (next < 0) {
                break;
            }
            here = next;
            next += shift/std::pow(10,a);


            hereValue = nextValue;
            nextValue = K.exactKBessel(next);
        }
    }

    return {here, hereValue};
}

void KBesselApproximator::extendPrecomputedRange(double newLowerBound, double newUpperBound) {
    if (newLowerBound < precomputedRegionLeftBound) {
        //Compute number of new intervals (to be inserted on the left)
        int existingFineIntervals = fineSplines.size();
        int newIntervalCount = std::ceil((r - newLowerBound)/fineIntervalWidth);
        newIntervalCount = newIntervalCount - existingFineIntervals;
        newIntervalCount = std::max(0, newIntervalCount);

        std::vector<std::vector<double>> precomputes;
        precomputes.resize(newIntervalCount);

        //Allocate memory for random access in parallelized loop
        for (int i = 0; i < newIntervalCount; i++) {
            precomputes[i].resize(fineSplineKnotCount);
        }

        KBesselExact K = KBesselExact(r);
        //Compute K-Bessel values at all the new points
#pragma omp parallel for collapse(2) default(none) shared(newIntervalCount, precomputes, existingFineIntervals, K)
        for (int i = 0; i < newIntervalCount; i++) {
            for (int j = 0; j < fineSplineKnotCount; j++) {
                double left = r - (existingFineIntervals + i + 1) * fineIntervalWidth;
                double x = left + j*fineSplineSpacing;
                double y = K.exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        //Build the new spline objects for each interval...
        //std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> newFineSplines;
        //newFineSplines.resize(newIntervalCount);
        fineSplines.resize(existingFineIntervals + newIntervalCount);
        //...in parallel
#pragma omp parallel for default(none) shared(newIntervalCount, precomputes, existingFineIntervals, K)
        for (int i = 0; i < newIntervalCount; i++) {
            double left = r - (existingFineIntervals + i + 1) * fineIntervalWidth;
            double derivativeLeft = K.estimateDerivativeKBessel(left);
            double right = left + fineSplineSpacing * (fineSplineKnotCount - 1);
            double derivativeRight = K.estimateDerivativeKBessel(right);
            boost::math::interpolators::cardinal_cubic_b_spline spline(precomputes[i].begin(),
                                                                       precomputes[i].end(),
                                                                       left,
                                                                       fineSplineSpacing,
                                                                       derivativeLeft,
                                                                       derivativeRight);
            fineSplines[existingFineIntervals + i] = spline;
        }

        //Update the new left bound for the precomputed region and insert all the new spline objects
        precomputedRegionLeftBound = r - (existingFineIntervals + newIntervalCount) * fineIntervalWidth;
    }

    if (precomputedRegionRightBound < newUpperBound) {
        //
        //Now to compute the new stuff to the right
        //

        //Compute number of new intervals (to be inserted on the right)
        int existingCoarseIntervals = coarseSplines.size();
        int newIntervalCount = std::ceil((newUpperBound - r)/coarseIntervalWidth);
        newIntervalCount = newIntervalCount - existingCoarseIntervals;
        newIntervalCount = std::max(0, newIntervalCount);

        std::vector<std::vector<double>> precomputes;
        precomputes.resize(newIntervalCount);

        //Allocate memory for random access in parallelized loop
        for (int i = 0; i < newIntervalCount; i++) {
            precomputes[i].resize(coarseSplineKnotCount);
        }
        KBesselExact K = KBesselExact(r);
        //Compute K-Bessel values at all the new points
#pragma omp parallel for collapse(2) default(none) shared(newIntervalCount, precomputes, existingCoarseIntervals, K)
        for (int i = 0; i < newIntervalCount; i++) {
            for (int j = 0; j < coarseSplineKnotCount; j++) {
                double left = r + (existingCoarseIntervals + i) * coarseIntervalWidth;
                double x = left + j * coarseSplineSpacing;
                double y = K.exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        //Build the new spline objects for each interval...
        //std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> newCoarseSplines;
        //newCoarseSplines.resize(newIntervalCount);
        coarseSplines.resize(existingCoarseIntervals + newIntervalCount);

        //...in parallel
#pragma omp parallel for default(none) shared(newIntervalCount, precomputes, existingCoarseIntervals, K)
        for (int i = 0; i < newIntervalCount; i++) {
            double left = r + (existingCoarseIntervals + i) * coarseIntervalWidth;
            double derivativeLeft = K.estimateDerivativeKBessel(left);
            double right = left + coarseSplineSpacing * (coarseSplineKnotCount - 1);
            double derivativeRight = K.estimateDerivativeKBessel(right);
            boost::math::interpolators::cardinal_cubic_b_spline spline(precomputes[i].begin(),
                                                                       precomputes[i].end(),
                                                                       left,
                                                                       coarseSplineSpacing,
                                                                       derivativeLeft,
                                                                       derivativeRight);
            coarseSplines[existingCoarseIntervals + i] = spline;
        }

        //Update the new right bound for the precomputed region and insert all the new spline objects
        precomputedRegionRightBound = r + (existingCoarseIntervals + newIntervalCount) * coarseIntervalWidth;
    }
    //runTest();
}

void KBesselApproximator::runTest() {
    double maxErrorSoFar = 0;
    KBesselExact K = KBesselExact(r);
    for (int i = 0; i < 100000; i++) {
        double x = r - i*(r-precomputedRegionLeftBound)/100000.0;
        double exact = K.exactKBessel(x);
        double approx = approxKBessel(x);
        double error = abs((approx - exact)/exact);
        if (error > maxErrorSoFar && exact != 0) {
            maxErrorSoFar = error;
            watch(maxErrorSoFar);
            std::cout << x << ", " << exact << std::endl;
        }
    }

    maxErrorSoFar = 0;
    for (int i = 0; i < 100000; i++) {
        double x = r + i*(precomputedRegionRightBound - r)/100000.0;
        double exact = K.exactKBessel(x);
        double approx = approxKBessel(x);
        double error = abs((approx - exact)/exact);
        if (error > maxErrorSoFar && exact != 0) {
            maxErrorSoFar = error;
            watch(maxErrorSoFar);
            std::cout << x << ", " << exact << std::endl;
        }
    }
}

/*
 * double chebyWidth = pow(2,-7);

double rand01() {
    return std::rand() / (RAND_MAX + 1.0);
}

double randab(double a, double b) {
    return a + rand01()*(b-a);
}

int transformerIndex(double x, double a) {
    int floor = std::floor((x-a)/chebyWidth);
    return floor;
}
 *
 * #include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions//chebyshev_transform.hpp>
#include <fftw3.h>
 *
 * KBesselApproximator K = KBesselApproximator(9);

    auto kBesselFunction = [&K](double x) {
        return K.computeKBessel(x);
    };



    double a = 8.5;
    double b = 50;
    K.setRAndPrecompute(9, a, b);

    double tol = std::numeric_limits<double>::epsilon();

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<boost::math::chebyshev_transform<double>> chebyshevTransformers;

    double left = a;

    while (left < b) {
        boost::math::chebyshev_transform<double> chebyshevTransformer([&](double x) { return kBesselFunction(x); }, left, left + chebyWidth, tol);
        chebyshevTransformers.push_back(chebyshevTransformer);
        left += chebyWidth;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "step 1 ";
    watch(duration);


    std::chrono::duration<double> exactDuration = std::chrono::duration<double> (0);
    std::chrono::duration<double> approxDuration = std::chrono::duration<double> (0);
    std::chrono::duration<double> lagrangeDuration = std::chrono::duration<double> (0);
    int numPoints = pow(10,6);
    double maxDiff = 0;
    for (int i = 0; i < numPoints; i++) {
        double x = a + i*(b-a)/(double)numPoints;
        start = std::chrono::high_resolution_clock::now();
        double exact = K.exactKBessel(x);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        exactDuration += duration;

        start = std::chrono::high_resolution_clock::now();
        double approx = chebyshevTransformers[transformerIndex(x,a)](x);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        approxDuration += duration;

        start = std::chrono::high_resolution_clock::now();
        double lagrangeApprox = K.approxKBessel(x);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        approxDuration += duration;

        double diff = exact - approx;
        diff = abs(diff);
        if (diff > maxDiff) {
            maxDiff = diff;
            watch(maxDiff);
            watch(x);
        }
    }

    watch(approxDuration);
    watch(exactDuration);
    watch(lagrangeDuration);
 */
