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

    precomputedRegionLeftBound = r;
    precomputedRegionRightBound = r;
}

double KBesselApproximator::approxKBessel(const double x) {
    //TODO: When everything is working properly, delete this.
    if (x < precomputedRegionLeftBound || x > precomputedRegionRightBound) {
        watch(x);
        watch(r);
        watch(precomputedRegionLeftBound);
        watch(precomputedRegionRightBound);
        throw std::invalid_argument("KBessel approximation out of bounds");
    }

    if (x < CHUNK_WIDTH) {
        /*
         * The scheme here is that there are "shrinking chunks" starting at CHUNK_WIDTH on downward.
         * Let C = CHUNK_WIDTH
         * The first chunk is [C/2, C)
         * The second chunk is [C/4, C/2)
         * and so on.
         * So the nth (zero-indexed) shrinking chunk is the interval [C/2^(n+1), C/2^n)
         *
         * If log_2(C/x) <= n+1 then C/2^(n+1) <= x.
         */
        int shrinkingChunkIndex = ceil(log2(CHUNK_WIDTH/x)) - 1;
        int subIndex = floor(
                (x - CHUNK_WIDTH/pow(2,shrinkingChunkIndex)) //distance from x to the left endpoint of its shrinking chunk
                /(SPLINE_KNOT_COUNT * shrinkingChunkStepSize[shrinkingChunkIndex]) //this is the width of the spline interval
                );
        return shrinkingChunks[shrinkingChunkIndex][subIndex](x);
    } else {
        int chunkIndex = floor((x - CHUNK_WIDTH)/CHUNK_WIDTH);
        int subIndex = floor(
                (x - (CHUNK_WIDTH + chunkIndex * CHUNK_WIDTH)) // distance from x to left endpoint of its chunk
                /(SPLINE_KNOT_COUNT * chunkStepSize[chunkIndex]) //this is the width of the spline interval
                );

        return chunks[chunkIndex][subIndex](x);
    }
}

void KBesselApproximator::setRAndPrecompute(double newR, double precomputeLowerBound, double precomputeUpperBound) {
    setRAndClear(newR);

    extendPrecomputedRange(precomputeLowerBound, precomputeUpperBound);
}

void KBesselApproximator::setRAndClear(double newR) {
    this->r = newR;

    chunkStepSize.clear();
    shrinkingChunkStepSize.clear();

    chunks.clear();
    shrinkingChunks.clear();

    precomputedRegionLeftBound = r;
    precomputedRegionRightBound = r;
}

vector<double> KBesselApproximator::maximize() {
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

/**
 * Precomputes values of the K-Bessel function on the interval [newLowerBound, newUpperBound]. This interval must be
 * contained in the interval (0, infty).
 *
 * The generic usage of this method should have the lower bound be quite small, say on the order of 10^-1.
 * The upper bound can be as large as desired.
 *
 * This entire class is optimized for this usage, so in fact this routine will always precompute the interval
 * [CHUNK_WIDTH, newUpperBound] since we expect to always need at least this much. Precomputing less than CHUNK_WIDTH
 * is more expensive and is done as requested.
 *
 * The user may call this function multiple times for a fixed parameter r. It will only precomputes more if the range
 * expands, i.e. it never recomputes things.
 *
 * This method checks for accuracy automatically up to the relative error which is a static private member of this class.
 * The accuracy check is not rigorously verified, only empirically at some problem points in each chunk.
 *
 * @param newLowerBound Positive number, usually small, on the order 10^-1
 * @param newUpperBound Positive number, usually large, on the order 10^3 to 10^4
 */
void KBesselApproximator::extendPrecomputedRange(double newLowerBound, double newUpperBound) {
    KBesselExact kbe = KBesselExact(r);

    /* ********************************************************************
     *
     * Precompute chunks, meaning values in the range [CHUNK_WIDTH, newUpperBound]
     *
     **********************************************************************/

    //Figure what new chunks need to be computed
    int indexOfLastNewChunk = ceil((newUpperBound - CHUNK_WIDTH)/CHUNK_WIDTH) - 1;

    //Add new slots in the chunkStepSize and chunks vectors
    if (indexOfLastNewChunk >= chunkStepSize.size()) {
        chunkStepSize.resize(indexOfLastNewChunk + 1);
        chunks.resize(indexOfLastNewChunk + 1);
    }

    //Compute spacing requirement for each new chunk
    //Do it in reverse since earlier chunks require spacing at least as fine as later chunks
    double spacing = pow(2.0, -1);
    for (int i = indexOfLastNewChunk; i >= computedChunkCount; i--) {
        bool errorIsSatisfactory = false;
        while (!errorIsSatisfactory) {
            //generate one interpolant
            int knots = std::min(SPLINE_KNOT_COUNT, (int)ceil(CHUNK_WIDTH/spacing) + 1);
            knots = std::max(3, knots); //3 knots required for cubic interpolation
            vector<double> precompute(knots, 0.0);

            double intervalLeft = CHUNK_WIDTH + i * CHUNK_WIDTH;
            double intervalRight = intervalLeft + (knots - 1) * spacing;
            for (int j = 0; j < knots; j++) {
                double x = intervalLeft + j * spacing;
                double y = kbe.exactKBessel(x);
                precompute[j] = y;
            }

            CubicSpline testSpline(precompute.begin(),
                                   precompute.end(),
                                   intervalLeft,
                                   spacing,
                                   kbe.estimateDerivativeKBessel(intervalLeft),
                                   kbe.estimateDerivativeKBessel(intervalRight));

            //check it for accuracy
            errorIsSatisfactory = true;
            for (int j = 1; j < 10; j++) {
                double x = intervalLeft + j * (spacing/10);
                double exact = kbe.exactKBessel(x);
                double approx = testSpline(x);
                if (exact == 0) {
                    continue;
                } else {
                    double relativeError = (approx - exact)/exact;
                    relativeError = abs(relativeError);

                    if (relativeError > ABS_ERROR_CUTOFF) {
                        errorIsSatisfactory = false;
                        break;
                    }
                }
            }

            if (errorIsSatisfactory) {
                break;
            }

            spacing /= 2.0;
        }
        chunkStepSize[computedChunkCount + i] = spacing;
    }
    
    //Compute all the splines, make it parallel
#pragma omp parallel for default(none) shared(indexOfLastNewChunk, kbe)
    for (int i = computedChunkCount; i <= indexOfLastNewChunk; i++) {
        //get left endpoint for this chunk
        //get spacing for this chunk
        //figure out how many intervals happen in this chunk
        //for each interval
        //  precompute K-Bessel values
        //  make spline
        //  add splines to the chunks vector

        //get left and right endpoints for this chunk
        double left = CHUNK_WIDTH + i * CHUNK_WIDTH;

        //get spacing for this chunk
        double spacing = chunkStepSize[i];

        //figure out how many intervals happen in this chunk
        int knots = std::min(SPLINE_KNOT_COUNT, (int)ceil(CHUNK_WIDTH/spacing) + 1);
        knots = std::max(3, knots);
        double intervalLength = (knots - 1) * spacing;
        int numIntervals = ceil(CHUNK_WIDTH/intervalLength);
        chunks[i].resize(numIntervals);

        //for each interval
        for (int j = 0; j < numIntervals; j++) {
            //  precompute K-Bessel values
            vector<double> precompute(knots, 0.0);
            double intervalLeft = left + j * intervalLength;
            double intervalRight = intervalLeft + (knots - 1) * spacing;
            for (int k = 0; k < knots; k++) {
                double x = intervalLeft + k * spacing;
                double y = kbe.exactKBessel(x);
                precompute[k] = y;
            }

            //  make spline
            CubicSpline spline(precompute.begin(),
                               precompute.end(),
                               intervalLeft,
                               spacing,
                               kbe.estimateDerivativeKBessel(intervalLeft),
                               kbe.estimateDerivativeKBessel(intervalRight));

            //  add splines to the chunks vector
            chunks[i][j] = spline;
        }
    }
    computedChunkCount = chunks.size();
    precomputedRegionRightBound = CHUNK_WIDTH + CHUNK_WIDTH * (computedChunkCount);


    /* ********************************************************************
     *
     * Precompute shrinking chunks, meaning values in the range [newLowerBound, CHUNK_WIDTH)
     *
     **********************************************************************/

    int indexOfLastNewShrinkingChunk = ceil(log2(CHUNK_WIDTH / newLowerBound)) - 1;
    indexOfLastNewShrinkingChunk = std::max(indexOfLastNewShrinkingChunk, 0);

    if (indexOfLastNewShrinkingChunk >= shrinkingChunks.size()) {
        shrinkingChunks.resize(indexOfLastNewShrinkingChunk + 1);
        shrinkingChunkStepSize.resize(indexOfLastNewShrinkingChunk + 1);
    }

    for (int newShrinkingChunk = computedShrinkingChunkCount; newShrinkingChunk <= indexOfLastNewShrinkingChunk; newShrinkingChunk++) {

        if (newShrinkingChunk == 0) {
            spacing = chunkStepSize[0];
        } else {
            spacing = shrinkingChunkStepSize[newShrinkingChunk - 1];
        }
        double left = CHUNK_WIDTH/pow(2,newShrinkingChunk + 1);

        bool errorIsSatisfactory = false;
        while (!errorIsSatisfactory) {
            //generate one interpolant
            vector<double> precompute(SPLINE_KNOT_COUNT, 0.0);

            for (int i = 0; i < SPLINE_KNOT_COUNT; i++) {
                double x = left + i * spacing;
                double y = kbe.exactKBessel(x);
                precompute[i] = y;
            }
            double right = left + SPLINE_KNOT_COUNT * spacing;

            CubicSpline testSpline(precompute.begin(),
                                   precompute.end(),
                                   left,
                                   spacing,
                                   kbe.estimateDerivativeKBessel(left),
                                   kbe.estimateDerivativeKBessel(right));

            //check it for accuracy
            errorIsSatisfactory = true;
            for (int j = 1; j < 10; j++) {
                double x = left + spacing * j / 10;
                double exact = kbe.exactKBessel(x);
                double approx = testSpline(x);
                if (exact == 0) {
                    continue;
                } else {
                    double relativeError = (approx - exact)/exact;
                    relativeError = abs(relativeError);

                    if (relativeError > ABS_ERROR_CUTOFF) {
                        errorIsSatisfactory = false;
                        break;
                    }
                }
            }

            if (errorIsSatisfactory) {
                break;
            }

            spacing /= 2.0;
        }
        shrinkingChunkStepSize[newShrinkingChunk] = spacing;
    }

#pragma omp parallel for default(none) shared(indexOfLastNewShrinkingChunk, spacing, kbe)
    for (int newShrinkingChunk = computedShrinkingChunkCount; newShrinkingChunk <= indexOfLastNewShrinkingChunk; newShrinkingChunk++) {

        double left = CHUNK_WIDTH/pow(2,newShrinkingChunk + 1);
        double right = CHUNK_WIDTH/pow(2, newShrinkingChunk);

        double intervalWidth = SPLINE_KNOT_COUNT * spacing;
        int numberOfIntervals = ceil((right - left)/intervalWidth);

        shrinkingChunks[newShrinkingChunk].resize(numberOfIntervals);
        for (int i = 0; i < numberOfIntervals; i++) {
            vector<double> precompute(SPLINE_KNOT_COUNT, 0.0);
            double intervalLeft = left + i * intervalWidth;
            double intervalRight = left + (SPLINE_KNOT_COUNT - 1) * intervalWidth;
            for (int j = 0; j < SPLINE_KNOT_COUNT; j++) {
                double x = intervalLeft + j * spacing;
                double y = kbe.exactKBessel(x);
                precompute[j] = y;
            }

            CubicSpline spline(precompute.begin(),
                               precompute.end(),
                               intervalLeft,
                               spacing,
                               kbe.estimateDerivativeKBessel(intervalLeft),
                               kbe.estimateDerivativeKBessel(intervalRight));

            shrinkingChunks[newShrinkingChunk][i] = spline;
        }
    }
    computedShrinkingChunkCount = shrinkingChunks.size();
    precomputedRegionLeftBound = CHUNK_WIDTH/pow(2.0, computedShrinkingChunkCount + 1);

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

double KBesselApproximator::getSpacing(double x) {
    int chunkIndex = floor((x - CHUNK_WIDTH)/CHUNK_WIDTH);
    return chunkStepSize[chunkIndex];
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
    vector<boost::math::chebyshev_transform<double>> chebyshevTransformers;

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
