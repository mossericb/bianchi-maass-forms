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

KBesselApproximator::KBesselApproximator(double r) {
    this->r = r;
}

double KBesselApproximator::approxKBessel(const double& x) {
    if (x < precomputedRegionLeftBound || x > precomputedRegionRightBound) {
        throw std::invalid_argument("KBessel approximation out of bounds");
    }

    if (x <= r) {
        return fineInterpolate(x);
    } else {
        return coarseInterpolate(x);
    }
}

void KBesselApproximator::updateRAndPrecompute(double newR, double precomputeLowerBound, double precomputeUpperBound) {
    updateR(newR);
    precomputedRegionLeftBound = precomputeLowerBound;
    precomputedRegionRightBound = precomputeUpperBound;

    finePrecomputedValues.clear();
    coarsePrecomputedValues.clear();

    for (int D = fineD; ; D--) {
        bool errorIsSatisfactory = true;
        fineSpacing = pow(2,D);
        fineD = D;
        finePrecomputedValues.clear();
        int numberOfFinePrecomputeValues = std::floor((r - precomputeLowerBound)/fineSpacing) + 1;
        finePrecomputedValues.resize(numberOfFinePrecomputeValues);

#pragma omp parallel for default(none) shared(fineSpacing, numberOfFinePrecomputeValues, precomputeLowerBound)
        for (int i = 0; i < numberOfFinePrecomputeValues; i++) {
            double x = precomputeLowerBound + i * fineSpacing;
            finePrecomputedValues[i] = computeKBessel(x);
        }

        //Check at 10 random points to see if the error is satisfactory
        for (int i = 0; i < 10; i++) {
            double x = randab(precomputedRegionLeftBound, r);
            double exact = computeKBessel(x);
            double approx = approxKBessel(x);
            double error  = abs((exact - approx)/exact);
            if (error > MAX_ERROR) {
                errorIsSatisfactory = false;
                break;
            }
        }
        if (!errorIsSatisfactory) {
            continue;
        }
        //Now we check 100 more times
        for (int i = 0; i < 100; i++) {
            double x = randab(precomputedRegionLeftBound, r);
            double exact = computeKBessel(x);
            double approx = approxKBessel(x);
            double error  = abs((exact - approx)/exact);
            if (error > MAX_ERROR) {
                errorIsSatisfactory = false;
                break;
            }
        }
        if (!errorIsSatisfactory) {
            continue;
        } else {
            break;
        }
    }

    for (int D = coarseD; ; D--) {
        bool errorIsSatisfactory = true;
        coarseSpacing = pow(2,D);
        coarseD = D;
        coarsePrecomputedValues.clear();
        int numberOfCoarsePrecomputeValues = std::floor((precomputedRegionRightBound - r)/coarseSpacing) + 1;
        coarsePrecomputedValues.resize(numberOfCoarsePrecomputeValues);

#pragma omp parallel for default(none) shared(coarseSpacing, numberOfCoarsePrecomputeValues, r)
        for (int i = 0; i < numberOfCoarsePrecomputeValues; i++) {
            double x = r + i * coarseSpacing;
            coarsePrecomputedValues[i] = computeKBessel(x);
        }

        //Check at 10 random points to see if the error is satisfactory
        for (int i = 0; i < 10; i++) {
            double x = randab(r, precomputedRegionRightBound);
            double exact = computeKBessel(x);
            double approx = approxKBessel(x);
            double error  = abs((exact - approx)/exact);
            if (error > MAX_ERROR) {
                errorIsSatisfactory = false;
                break;
            }
        }
        if (!errorIsSatisfactory) {
            continue;
        }
        //Now we check 100 more times
        for (int i = 0; i < 100; i++) {
            double x = randab(r, precomputedRegionRightBound);
            double exact = computeKBessel(x);
            double approx = approxKBessel(x);
            double error  = abs((exact - approx)/exact);
            if (error > MAX_ERROR) {
                errorIsSatisfactory = false;
                break;
            }
        }
        if (!errorIsSatisfactory) {
            continue;
        } else {
            break;
        }
    }
}

void KBesselApproximator::updateR(double newR) {
    this->r = newR;
}

double KBesselApproximator::computeKBessel(const double& x) {
    if (x == 0) {
        return 0;
    }
    arb_t localRealTemp;
    arb_t localR;
    acb_t localComplexTemp;
    arb_t localBesselScale;
    acb_t localNu;
    arb_init(localRealTemp);
    acb_init(localComplexTemp);
    arb_init(localBesselScale);
    arb_init(localR);
    acb_init(localNu);

    arb_set_d(localR, r);
    acb_set_d_d(localNu, 0, r);

    for (int prec = 64; ; prec *= 2) {

        arb_const_pi(localBesselScale, prec);
        arb_mul(localBesselScale, localBesselScale, localR, prec);
        arb_div_ui(localBesselScale, localBesselScale, 2, prec);
        arb_exp(localBesselScale, localBesselScale, prec);

         acb_set_d_d(localComplexTemp, x, 0);

        acb_hypgeom_bessel_k(localComplexTemp, localNu, localComplexTemp, prec);

        acb_get_real(localRealTemp,localComplexTemp);

        arb_mul(localRealTemp, localRealTemp, localBesselScale, prec);
        if (arb_rel_accuracy_bits(localRealTemp) >= 53) {
            break;
        }
    }


    double answer;
    try {
        answer = std::stod(arb_get_str(localRealTemp, 50, ARB_STR_NO_RADIUS));
    } catch (...) {
        answer = 0;
    }

    arb_clear(localRealTemp);
    acb_clear(localComplexTemp);
    arb_clear(localR);
    arb_clear(localBesselScale);
    acb_clear(localNu);

    flint_cleanup();

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

    double startValue = computeKBessel(start);
    double leftValue = computeKBessel(left);

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
        double nextValue = computeKBessel(next);
        while (nextValue > hereValue) {
            if (next < 0) {
                break;
            }
            here = next;
            next += shift/std::pow(10,a);


            hereValue = nextValue;
            nextValue = computeKBessel(next);
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

void KBesselApproximator::extendPrecomputedRange(double precomputeUpperBound) {
    if (precomputeUpperBound < precomputedRegionRightBound) {
        return; //do nothing
    }

    if (precomputeUpperBound < r) {
        updateRAndPrecompute(r, precomputedRegionLeftBound, precomputeUpperBound);
    } else {
        while (precomputeUpperBound > precomputedRegionRightBound) {
            precomputedRegionRightBound += fineSpacing;
            double exact = computeKBessel(precomputedRegionRightBound);
            finePrecomputedValues.push_back(exact);
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
    K.updateRAndPrecompute(9, a, b);

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
        double exact = K.computeKBessel(x);
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
