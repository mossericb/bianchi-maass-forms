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

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl

KBesselApproximator::KBesselApproximator(int bitsOfPrecision) {
    prec = bitsOfPrecision;
    r = 0;

    #pragma omp parallel default(none)
    {
        #pragma omp single
        threads = omp_get_num_threads();
    }

    realTemp.resize(threads);
    arbR.resize(threads);
    besselScale.resize(threads);
    preComputedBesselScale.resize(threads);
    complexTemp.resize(threads);
    acbNu.resize(threads);

    for (int i = 0; i < threads; i++) {
        arb_init(&realTemp[i]);
        arb_init(&arbR[i]);
        arb_init(&besselScale[i]);
        arb_init(&preComputedBesselScale[i]);
        acb_init(&complexTemp[i]);
        acb_init(&acbNu[i]);
    }

    fineIntervalWidth = fineSplineSpacing * (fineSplineKnotCount - 1);
    coarseIntervalWidth = coarseSplineSpacing * (coarseSplineKnotCount - 1);
}

double KBesselApproximator::approxKBessel(const double& x) {
    if (x < precomputedRegionLeftBound || x > precomputedRegionRightBound) {
        throw std::invalid_argument("KBessel approximation out of bounds");
    }

    if (x < r) {
        int leftIndex = std::floor((x - precomputedRegionLeftBound)/fineIntervalWidth);
        return fineSplines[leftIndex](x);
        double a = fineSplineEndpoints[leftIndex];
        double b = fineSplineEndpoints[leftIndex + 1];
        if (a <= x && x <= b) {
            return fineSplines[leftIndex](x);
        } else if (x < a) {
            return fineSplines[leftIndex - 1](x);
        } else {
            return fineSplines[leftIndex + 1](x);
        }
    } else {
        int leftIndex = std::floor((x - r)/coarseIntervalWidth);
        return coarseSplines[leftIndex](x);
        double a = coarseSplineEndpoints[leftIndex];
        double b = coarseSplineEndpoints[leftIndex + 1];
        if (a <= x && x <= b) {
            return coarseSplines[leftIndex](x);
        } else if (x < a) {
            return coarseSplines[leftIndex - 1](x);
        } else {
            return coarseSplines[leftIndex + 1](x);
        }
    }
}

void KBesselApproximator::setRAndPrecompute(double newR, double precomputeLowerBound, double precomputeUpperBound) {
    setRAndClear(r);

    finePrecomputedValues.clear();
    coarsePrecomputedValues.clear();
    fineChebyshevApproximators.clear();
    coarseChebyshevApproximators.clear();
    fineSplinePrecompute.clear();
    coarseSplinePrecompute.clear();
    fineSplines.clear();
    coarseSplines.clear();

    if (precomputeLowerBound < r) {
        //Guaranteed to initialize the interval [precomputeLowerBound, r]
        precomputedRegionLeftBound = precomputeLowerBound;
        std::vector<std::vector<double>> precomputes;

        int intervalCount = std::ceil((r - precomputeLowerBound) / fineIntervalWidth);

        precomputes.resize(intervalCount);
        for (int i = 0; i < intervalCount; i++) {
            precomputes[i].resize(fineSplineKnotCount);
        }

#pragma omp parallel for collapse(2) default(none) shared(precomputes, precomputeLowerBound)
        for (int i = 0; i < precomputes.size(); i++) {
            for (int j = 0; j < fineSplineKnotCount; j++) {
                double x = precomputeLowerBound + i * fineIntervalWidth + j * fineSplineSpacing;
                double y = exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        fineSplines.resize(intervalCount);
#pragma omp parallel for default(none) shared(precomputes, precomputeLowerBound)
        for (int i = 0; i < precomputes.size(); i++) {
            boost::math::interpolators::cardinal_cubic_b_spline<double> spline(precomputes[i].begin(),
                                                                               precomputes[i].end(),
                                                                               precomputeLowerBound +
                                                                               i * fineIntervalWidth,
                                                                               fineSplineSpacing);

            fineSplines[i] = spline;
        }

    } else {
        precomputedRegionLeftBound = r;
    }


    if (r <= precomputeUpperBound) {
        //Guaranteed to initialize the interval [r, precomputeUpperBound]
        //Will in fact go a little further up than that
        //Goes up to r + coarseIntervalWidth*intervalCount
        std::vector<std::vector<double>> precomputes;

        int intervalCount = std::ceil((precomputeUpperBound - r) / coarseIntervalWidth);
        precomputedRegionRightBound = r + coarseIntervalWidth * intervalCount;

        precomputes.resize(intervalCount);
        for (int i = 0; i < intervalCount; i++) {
            precomputes[i].resize(coarseSplineKnotCount);
        }

#pragma omp parallel for collapse(2) default(none) shared(precomputes)
        for (int i = 0; i < precomputes.size(); i++) {
            for (int j = 0; j < coarseSplineKnotCount; j++) {
                double x = r + i * coarseIntervalWidth + j * coarseSplineSpacing;
                double y = exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        coarseSplines.resize(intervalCount);
#pragma omp parallel for default(none) shared(precomputes)
        for (int i = 0; i < precomputes.size(); i++) {
            boost::math::interpolators::cardinal_cubic_b_spline<double> spline(precomputes[i].begin(),
                                                                               precomputes[i].end(),
                                                                               r + i * coarseIntervalWidth,
                                                                               coarseSplineSpacing);

            coarseSplines[i] = spline;
        }

    } else {
        precomputedRegionRightBound = precomputeLowerBound + fineIntervalWidth * fineSplines.size();
    }
}

void KBesselApproximator::setRAndClear(double newR) {
    this->r = newR;
    for (int i = 0; i < threads; i++) {
        arb_set_d(&arbR[i],r);
        acb_set_d_d(&acbNu[i], 0, r);

        arb_const_pi(&preComputedBesselScale[i], 2*prec);
        arb_mul(&preComputedBesselScale[i], &preComputedBesselScale[i], &arbR[i], 2*prec);
        arb_div_ui(&preComputedBesselScale[i], &preComputedBesselScale[i], 2, 2*prec);
        arb_exp(&preComputedBesselScale[i], &preComputedBesselScale[i], 2*prec);
    }
}

double KBesselApproximator::exactKBessel(const double& x) {
    if (x == 0) {
        return 0;
    }

    int threadNum = omp_get_thread_num();

    acb_set_d_d(&complexTemp[threadNum], x, 0);
    //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    //probably this precision will be enough, if not...
    int tryOnceBits = 2 * prec;

    acb_hypgeom_bessel_k(&complexTemp[threadNum],
                         &acbNu[threadNum],
                         &complexTemp[threadNum],
                         tryOnceBits);
    //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    acb_get_real(&realTemp[threadNum],&complexTemp[threadNum]);
    //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    arb_mul(&realTemp[threadNum],
            &realTemp[threadNum],
            &preComputedBesselScale[threadNum],
            tryOnceBits);
    //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
    //std::cout << std::endl;

    long accuracy = arb_rel_accuracy_bits(&realTemp[threadNum]);
    if (accuracy >= prec || accuracy >= 53) {
        double answer;
        try {
            answer = std::stod(arb_get_str(&realTemp[threadNum], 50, ARB_STR_NO_RADIUS));
        } catch (...) {
            answer = 0;
        }
        return answer;
    } else {
        //...then it enters this loop which will run until it reaches the desired precision
        for (int bits = 4 * prec; ; bits *= 2) {
            acb_set_d_d(&complexTemp[threadNum], x, 0);

            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_const_pi(&besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_mul(&besselScale[threadNum], &arbR[threadNum], &besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_div_ui(&besselScale[threadNum], &besselScale[threadNum], 2, bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_exp(&besselScale[threadNum], &besselScale[threadNum], bits);
            //arb_printn(&besselScale[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            acb_hypgeom_bessel_k(&complexTemp[threadNum],
                                 &acbNu[threadNum],
                                 &complexTemp[threadNum],
                                 bits);
            //acb_printn(&complexTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            acb_get_real(&realTemp[threadNum],&complexTemp[threadNum]);
            //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;
            arb_mul(&realTemp[threadNum],
                    &realTemp[threadNum],
                    &besselScale[threadNum],
                    bits);
            //arb_printn(&realTemp[threadNum], 20, ARB_STR_NO_RADIUS);
            //std::cout << std::endl;

            accuracy = arb_rel_accuracy_bits(&realTemp[threadNum]);
            if (accuracy >= prec || accuracy >= 53) {
                break;
            }
        }
    }

    double answer;
    try {
        answer = std::stod(arb_get_str(&realTemp[threadNum], 50, ARB_STR_NO_RADIUS));
    } catch (...) {
        answer = 0;
    }

    return answer;
}

std::string KBesselApproximator::dtos(const double& x) {
    std::stringstream out;
    out << std::setprecision(100) << x;
    return out.str();
}

void KBesselApproximator::printTiming() {
    std::cout << "First check: " << firstCheckDuration.count() << std::endl;
    std::cout << "Nearest three search: " << nearestThreeDuration.count() << std::endl;
    std::cout << "Compute immediately: " << immediatelyComputeDuration.count() << std::endl;
    std::cout << "Try to interpolate, but not close enough: " << computeBecausePointsNotCloseEnoughDuration.count() << std::endl;
    std::cout << "Interpolate: " << interpolateDuration.count() << std::endl << std::flush;
}

double KBesselApproximator::neville(std::vector<Point> &points, short i, short j, const double &x) {
    if (i == j) {
        return points[i].y;
    } else {
        double answer = (x-points[i].x)*neville(points, i+1, j,x);
        answer -= (x-points[j].x)*neville(points, i, j-1, x);
        answer /= (points[j].x - points[i].x);
        return answer;
    }
}

std::vector<double> KBesselApproximator::maximize() {
    double start = r;
    double left = r - 0.01;

    double startValue = exactKBessel(start);
    double leftValue = exactKBessel(left);

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
        double nextValue = exactKBessel(next);
        while (nextValue > hereValue) {
            if (next < 0) {
                break;
            }
            here = next;
            next += shift/std::pow(10,a);


            hereValue = nextValue;
            nextValue = exactKBessel(next);
        }
    }

    return {here, hereValue};
}

double KBesselApproximator::rand01() {
    return std::rand() / (RAND_MAX + 1.0);
}

double KBesselApproximator::randab(double a, double b) {
    return a + rand01()*(b-a);
}

double KBesselApproximator::fineLagrange(const KBesselApproximator::Point &p1, const KBesselApproximator::Point &p2,
                                     const KBesselApproximator::Point &p3, const double &x) {
    double x1 = p1.x;
    double x2 = p2.x;
    double x3 = p3.x;
    //watch(x1);
    //watch(x2);
    //watch(x3);

    double y1 = p1.y;
    double y2 = p2.y;
    double y3 = p3.y;
    //watch(y1);
    //watch(y2);
    //watch(y3);

    double diff1 = x - x1;
    double diff2 = x - x2;
    double diff3 = x - x3;

    double answer = y1*diff2*diff3/2;
    answer += -y2*diff1*diff3;
    answer += y3*diff1*diff2/2;
    answer /= pow(fineSpacing,2);
    return answer;
}

double KBesselApproximator::fineInterpolate(const double& x) {
    double measure = (x - precomputedRegionLeftBound)/fineSpacing;
    int left = std::floor(measure);
    int index1;
    int index2;
    int index3;

    if (left == 0) {
        index1 = 0;
        index2 = 1;
        index3 = 2;
    } else if (left == finePrecomputedValues.size() - 2) {
        index1 = finePrecomputedValues.size() - 3;
        index2 = finePrecomputedValues.size() - 2;
        index3 = finePrecomputedValues.size() - 1;
    } else if  (x - (precomputedRegionLeftBound + left*fineSpacing) < fineSpacing/2) {
        index1 = left - 1;
        index2 = left;
        index3 = left + 1;
    } else {
        index1 = left;
        index2 = left + 1;
        index3 = left + 2;
    }

    Point p1 = Point(precomputedRegionLeftBound + index1 * fineSpacing, finePrecomputedValues[index1]);
    Point p2 = Point(precomputedRegionLeftBound + index2 * fineSpacing, finePrecomputedValues[index2]);
    Point p3 = Point(precomputedRegionLeftBound + index3 * fineSpacing, finePrecomputedValues[index3]);


    return fineLagrange(p1, p2, p3, x);
}

double KBesselApproximator::coarseInterpolate(const double& x) {
    double measure = (x - r)/coarseSpacing;
    int left = std::floor(measure);
    int index1;
    int index2;
    int index3;

    if (left == 0) {
        index1 = 0;
        index2 = 1;
        index3 = 2;
    } else if (left == coarsePrecomputedValues.size() - 2) {
        index1 = coarsePrecomputedValues.size() - 3;
        index2 = coarsePrecomputedValues.size() - 2;
        index3 = coarsePrecomputedValues.size() - 1;
    } else if  (x - (r + left * coarseSpacing) < coarseSpacing/2) {
        index1 = left - 1;
        index2 = left;
        index3 = left + 1;
    } else {
        index1 = left;
        index2 = left + 1;
        index3 = left + 2;
    }

    Point p1 = Point(r + index1 * coarseSpacing, coarsePrecomputedValues[index1]);
    Point p2 = Point(r + index2 * coarseSpacing, coarsePrecomputedValues[index2]);
    Point p3 = Point(r + index3 * coarseSpacing, coarsePrecomputedValues[index3]);


    return coarseLagrange(p1, p2, p3, x);
}

double KBesselApproximator::coarseLagrange(const KBesselApproximator::Point &p1, const KBesselApproximator::Point &p2,
                                           const KBesselApproximator::Point &p3, const double &x) {
    double x1 = p1.x;
    double x2 = p2.x;
    double x3 = p3.x;
    //watch(x1);
    //watch(x2);
    //watch(x3);

    double y1 = p1.y;
    double y2 = p2.y;
    double y3 = p3.y;
    //watch(y1);
    //watch(y2);
    //watch(y3);

    double diff1 = x - x1;
    double diff2 = x - x2;
    double diff3 = x - x3;

    double answer = y1*diff2*diff3/2;
    answer += -y2*diff1*diff3;
    answer += y3*diff1*diff2/2;
    answer /= pow(coarseSpacing,2);
    return answer;
}

void KBesselApproximator::extendPrecomputedRange(double newLowerBound, double newUpperBound) {
    if (precomputedRegionLeftBound <= newLowerBound && newUpperBound <= precomputedRegionRightBound) {
        return; //do nothing
    }

    if (r <= newLowerBound) {
        //this would be unusual!
        //the whole interval lies in the "coarse" range

    } else if (newUpperBound < r) {
        //this would be unusual!
        //the whole interval lies in the "fine" range

    } else {
        //the typical case
        //interval straddles the x=r point

        //Compute number of new intervals (to be inserted on the left)
        int newIntervalCount = std::ceil((precomputedRegionLeftBound - newLowerBound)/fineIntervalWidth);

        std::vector<std::vector<double>> precomputes;
        precomputes.resize(newIntervalCount);

        //Allocate memory for random access in parallelized loop
        for (int i = 0; i < newIntervalCount; i++) {
            precomputes[i].resize(fineSplineKnotCount);
        }

        //Compute K-Bessel values at all the new points
        #pragma omp parallel for collapse(2) default(none) shared(newIntervalCount, precomputes)
        for (int i = 0; i < newIntervalCount; i++) {
            for (int j = 0; j < fineSplineKnotCount; j++) {
                double left = precomputedRegionLeftBound - (i+1)*fineIntervalWidth;
                double x = left + j*fineSplineSpacing;
                double y = exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        //Build the new spline objects for each interval...
        std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> newFineSplines;
        newFineSplines.resize(newIntervalCount);
        //...in parallel
        #pragma omp parallel for default(none) shared(newIntervalCount, precomputes, newFineSplines)
        for (int i = 0; i < newIntervalCount; i++) {
            double left = precomputedRegionLeftBound - (i+1)*fineIntervalWidth;
            boost::math::interpolators::cardinal_cubic_b_spline spline(precomputes[i].begin(),
                                                                       precomputes[i].end(),
                                                                       left,
                                                                       fineSplineSpacing);
            newFineSplines[i] = spline;
        }

        //Update the new left bound for the precomputed region and insert all the new spline objects
        precomputedRegionLeftBound = precomputedRegionLeftBound - newIntervalCount*fineIntervalWidth;
        for (int i = 0; i < newIntervalCount; i++) {
            fineSplines.insert(fineSplines.begin(), newFineSplines[i]);
        }

        //
        //Now to compute the new stuff to the right
        //

        //Compute number of new intervals (to be inserted on the right)
        newIntervalCount = std::ceil((newUpperBound - precomputedRegionRightBound)/coarseIntervalWidth);

        precomputes.resize(newIntervalCount);

        //Allocate memory for random access in parallelized loop
        for (int i = 0; i < newIntervalCount; i++) {
            precomputes[i].resize(coarseSplineKnotCount);
        }

        //Compute K-Bessel values at all the new points
#pragma omp parallel for collapse(2) default(none) shared(newIntervalCount, precomputes)
        for (int i = 0; i < newIntervalCount; i++) {
            for (int j = 0; j < coarseSplineKnotCount; j++) {
                double left = precomputedRegionRightBound + i*coarseIntervalWidth;
                double x = left + j*coarseSplineSpacing;
                double y = exactKBessel(x);
                precomputes[i][j] = y;
            }
        }

        //Build the new spline objects for each interval...
        std::vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> newCoarseSplines;
        newCoarseSplines.resize(newIntervalCount);
        //...in parallel
#pragma omp parallel for default(none) shared(newIntervalCount, precomputes, newCoarseSplines)
        for (int i = 0; i < newIntervalCount; i++) {
            double left = precomputedRegionLeftBound + i*coarseIntervalWidth;
            boost::math::interpolators::cardinal_cubic_b_spline spline(precomputes[i].begin(),
                                                                       precomputes[i].end(),
                                                                       left,
                                                                       coarseSplineSpacing);
            newCoarseSplines[i] = spline;
        }

        //Update the new right bound for the precomputed region and insert all the new spline objects
        precomputedRegionRightBound = precomputedRegionRightBound + newIntervalCount*coarseIntervalWidth;
        for (int i = 0; i < newIntervalCount; i++) {
            coarseSplines.push_back(newCoarseSplines[i]);
        }
    }
    double maxErrorSoFar = 0;
    for (int i = 0; i < 100000; i++) {
        double x = precomputedRegionLeftBound + i*(r-precomputedRegionLeftBound)/100000.0;
        double exact = exactKBessel(x);
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
        double exact = exactKBessel(x);
        double approx = approxKBessel(x);
        double error = abs((approx - exact)/exact);
        if (error > maxErrorSoFar && exact != 0) {
            maxErrorSoFar = error;
            watch(maxErrorSoFar);
            std::cout << x << ", " << exact << std::endl;
        }
    }
}

double KBesselApproximator::applyChebyshev(const double &x) {
    if (x < r) {
        //use fine cheby
        int index = fineChebyshevIndex(x);
        return fineChebyshevApproximators[index](x);
    } else {
        //use coarse cheby
        int index = coarseChebyshevIndex(x);
        return coarseChebyshevApproximators[index](x);
    }
}

int KBesselApproximator::coarseChebyshevIndex(const double &x) {
    return std::floor(std::sqrt(x) - std::sqrt(r));
}

int KBesselApproximator::fineChebyshevIndex(const double &x) {
    return std::floor(r/x - 0.5);
}

KBesselApproximator::~KBesselApproximator() {
    for (int i = 0; i < threads; i++) {
        arb_clear(&realTemp[i]);
        arb_clear(&arbR[i]);
        arb_clear(&besselScale[i]);
        arb_clear(&preComputedBesselScale[i]);
        acb_clear(&complexTemp[i]);
        acb_clear(&acbNu[i]);
    }

    flint_cleanup();
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
