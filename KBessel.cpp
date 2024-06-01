//
// Created by Eric Moss on 8/9/23.
//

#include "KBessel.h"
#include <flint/flint.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <boost/math/differentiation/finite_difference.hpp>
#include "omp.h"
#include "Auxiliary.h"

#define watch(x) std::cout << (#x) << " is " << (x) << std::endl
/*
 * Idea: Instead of dividing the range into two parts: coarse and fine, separate it into many more spacings
 * Could be like, 10 different levels of fineness, the finest being a point every 2^-19, and the coarsest being 2^-9
 * Could be something even more than that, something more continuous.
 * Either way, it should be something dynamic, so that once I start using large spectral parameters, it still works.
 *
 *
 */

KBessel::KBessel(double precomputeLowerBound, double r) {
    precomputedRegionLeftBound = precomputeLowerBound;
    precomputedRegionRightBound = precomputeLowerBound;

#pragma omp parallel default(none) shared(threads)
    {
#pragma omp single
        threads = omp_get_max_threads();
    }

    K.reserve(threads);
    for (int i = 0; i < threads; i++) {
        K.push_back(new archtKBessel(r));
    }

    setRAndClear(r);
}

KBessel::~KBessel() {
    for (auto bess : K) {
        delete bess;
    }
    K.clear();
}

double KBessel::exactKBessel(double x) {
    int threadNum = omp_get_thread_num();
    double ans = K[threadNum]->evaluate(x);
    return ans;
}

double KBessel::approxKBessel(const double x) {
    if (x >= zeroCutoff) {
        return 0;
    } else if (x < precomputedRegionLeftBound) {
        throw std::invalid_argument("x cannot be below 2*pi/A*1*Y0");
    } else if (x < firstChunkLeftEndpoint) { /* 2*pi/A*1*Y0 <= x < firstChunkLeftEndpoint */
        /*
         * The scheme here is that the interval [2*pi/A*1*Y0, firstChunkLeftEndpoint) = [a,b) is divided
         * into k shrinking pieces in the following way
         * a = a + (b-a)(0/k)^2, a + (b-a)(1/k)^2, a + (b-a)(2/k)^2, ... , a + (b-a)(k/k)^2 = b
         *
         * we have
         * a + (b-a)(l/k)^2 <= x < a + (b-a)((l+1)/k)^2
         * iff
         * l <= k sqrt((x-a)/(b-a)) < l + 1
         */

        int shrinkingChunkIndex = floor(
            numberOfShrinkingChunks
            * sqrt(
                (x - precomputedRegionLeftBound)
                /(firstChunkLeftEndpoint - precomputedRegionLeftBound)
            )
        );

        double distToLeftEndpoint = x
            - (precomputedRegionLeftBound
                + (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                * pow(((1. * shrinkingChunkIndex)/numberOfShrinkingChunks), 2)
              );

        double widthOfSpline = (SPLINE_KNOT_COUNT - 1) * shrinkingChunkStepSize[shrinkingChunkIndex];
        int subIndex = floor(distToLeftEndpoint/widthOfSpline);

        return shrinkingChunks[shrinkingChunkIndex][subIndex](x);
    } else { /* firstChunkLeftEndpoint <= x < zeroCutoff */
        int chunkIndex = floor((x - firstChunkLeftEndpoint) / chunkWidth);
        int subIndex = floor(
                (x - (firstChunkLeftEndpoint + chunkIndex * chunkWidth)) // distance from x to left endpoint of its chunk
                / ((chunkKnotCount[chunkIndex] - 1) * chunkStepSize[chunkIndex]) //this is the width of the spline interval
                );

        return chunks[chunkIndex][subIndex](x);
    }
}

void KBessel::setRAndPrecompute(double newR, double precomputeUpperBound) {
    setRAndClear(newR);

    extendPrecomputedRange(precomputeUpperBound);
}

void KBessel::setRAndClear(double newR) {
    r = newR;

    for (auto A : K) {
        A->setR(r);
    }

    chunkWidth = 1.0;

    //In general, I expect oscillation to start happening around r - 1
    //However, I always want to make sure that this value is always past the left bound
    firstChunkLeftEndpoint = std::max(r - 1, precomputedRegionLeftBound + 1);

    //set what I want the last width (this is the widest shrinking chunk)
    double shrinkingChunkMaxLastWidth = std::min(1.0, (firstChunkLeftEndpoint - precomputedRegionLeftBound)/2);

    int n = std::max(threads * 2, 2);

    while (true) {
        double lastWidth = firstChunkLeftEndpoint
                - (precomputedRegionLeftBound
                    + (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                    * pow((n - 1)/(1.*n), 2)
                 );
        if (lastWidth > shrinkingChunkMaxLastWidth) {
            ++n;
        } else {
            break;
        }
    }

    numberOfShrinkingChunks = n;

    zeroCutoff = (1136 + PI*r/2.0*log2(E) - 0.5*log2(E) + 0.5*log2(PI/2))/log2(E);

    chunkStepSize.clear();
    shrinkingChunkStepSize.clear();

    chunks.clear();
    shrinkingChunks.clear();

    chunkKnotCount.clear();

    precomputedRegionRightBound = precomputedRegionLeftBound;
}


/**
 * Precomputes values of the K-Bessel function on the interval [newLowerBound, newUpperBound]. This interval must be
 * contained in the interval (0, infty).
 *
 * The generic usage of this method should have the lower bound be quite small, say on the order of 10^-1.
 * The upper bound can be as large as desired.
 *
 * This entire class is optimized for this usage, so in fact this routine will always precompute the interval
 * [CHUNK_WIDTH, newUpperBound] since we expect to always need at least this much. Precomputing less than chunkWidth
 * is more expensive and is done as requested.
 *
 * The user may call this function multiple times for a fixed parameter r. It will only precomputes more if the range
 * expands, i.e. it never recomputes things.
 *
 * This method checks for accuracy automatically up to the relative relativeError which is a static private member of this class.
 * The accuracy check is not rigorously verified, only empirically at some problem points in each chunk.
 *
 * @param newLowerBound Positive number, usually small, on the order 10^-1
 * @param newUpperBound Positive number, usually large, on the order 10^3 to 10^4
 */
void KBessel::extendPrecomputedRange(double newUpperBound) {
    /* ********************************************************************
     *
     * Precompute chunks, meaning values in the range [chunkWidth, newUpperBound]
     *
     **********************************************************************/

    if (newUpperBound < precomputedRegionLeftBound && !shrinkingChunks.empty()) {
        return;
    } else if (newUpperBound > zeroCutoff) {
        newUpperBound = zeroCutoff;
    }

    int previouslyComputedChunks = chunks.size();

    //Figure what new chunks need to be computed
    int indexOfLastNewChunk = floor((newUpperBound - firstChunkLeftEndpoint) / chunkWidth);

    //Medium lift
    setUpChunkSplineComputation(indexOfLastNewChunk, previouslyComputedChunks);

    //Heavy lift
    computeChunkSplines(indexOfLastNewChunk, previouslyComputedChunks);


    /* ********************************************************************
     *
     * Precompute shrinking chunks, meaning values in the range [newLowerBound, chunkWidth)
     * This is only done once per update of spectral parameter r.
     *
     **********************************************************************/

    if (!shrinkingChunks.empty()) {
        return;
    }

    //Light lift
    setUpShrinkingChunkSplineComputation();

    //Heavy lift
    computeShrinkingChunkSplines();
}

void KBessel::runTest() {
    double maxErrorSoFar = 0;

    for (int i = 0; i < 1000; i++) {
        double x = r - i*(r-precomputedRegionLeftBound)/1000.0;
        double exact = exactKBessel(x);
        double approx = approxKBessel(x);
        double error = relativeError(exact, approx);
        if (error > maxErrorSoFar && exact != 0 && error > ABS_ERROR_CUTOFF) {
            maxErrorSoFar = error;
            watch(maxErrorSoFar);
            std::cout << x << ", " << exact << std::endl;
        }
    }

    maxErrorSoFar = 0;
    for (int i = 0; i < 1000; i++) {
        double x = r + i*(precomputedRegionRightBound - r)/1000.0;
        double exact = exactKBessel(x);
        double approx = approxKBessel(x);
        double error = relativeError(exact, approx);
        if (error > maxErrorSoFar && exact != 0 && error > ABS_ERROR_CUTOFF) {
            maxErrorSoFar = error;
            watch(maxErrorSoFar);
            std::cout << x << ", " << exact << std::endl;
        }
    }
}

double KBessel::relativeError(double exact, double approx, double x) {
    return abs(exact - approx)/(pow(10,-30) + abs(exact));
}

double KBessel::approxDerivativeKBessel(double x) {
    auto f = [this](double x) { return this->exactKBessel(x); };
    double dfdx = boost::math::differentiation::finite_difference_derivative(f, x);

    return dfdx;
}

void KBessel::setUpChunkSplineComputation(int indexOfLastNewChunk, int previouslyComputedChunks) {
    //Add new slots in the chunkStepSize and chunks vectors
    if (indexOfLastNewChunk >= chunks.size()) {
        chunkStepSize.resize(indexOfLastNewChunk + 1);
        chunkKnotCount.resize(indexOfLastNewChunk + 1);
        chunks.resize(indexOfLastNewChunk + 1);
    }

    //Compute spacing requirement for each new chunk
    //Do it in reverse since earlier chunks require spacing at least as fine as later chunks

#pragma omp parallel for schedule(dynamic) default(none) shared(indexOfLastNewChunk, previouslyComputedChunks)
    for (int i = indexOfLastNewChunk; i >= previouslyComputedChunks; i--) {

        double spacing = pow(2,-1);
        for (int j = i + 1; j < chunkStepSize.size(); j++) {
            if (chunkStepSize[j] > 0) {
                spacing = chunkStepSize[j];
                break;
            }
        }

        int knots = 0;
        bool errorIsSatisfactory = false;

        CubicSpline testSpline;

        while (!errorIsSatisfactory) {
            //generate one interpolant
            knots = std::min(SPLINE_KNOT_COUNT, (int)ceil(chunkWidth / spacing) + 1);
            knots = std::max(3, knots); //3 knots required for cubic interpolation
            vector<double> precompute(knots, 0.0);

            double intervalLeft = firstChunkLeftEndpoint + i * chunkWidth;
            double intervalRight = intervalLeft + (knots - 1) * spacing;
            for (int j = 0; j < knots; j++) {
                double x = intervalLeft + j * spacing;
                double y = exactKBessel(x);
                precompute[j] = y;
            }

            testSpline = CubicSpline (precompute.begin(),
                                      precompute.end(),
                                      intervalLeft,
                                      spacing,
                                      approxDerivativeKBessel(intervalLeft),
                                      approxDerivativeKBessel(intervalRight));

            //check it for accuracy
            errorIsSatisfactory = true;
            for (int j = 1; j < 10; j++) {
                double x = intervalLeft + j * (spacing/10);
                double exact = exactKBessel(x);
                double approx = testSpline(x);

                double error = relativeError(exact, approx);
                if (error > ABS_ERROR_CUTOFF) {
                    errorIsSatisfactory = false;
                    break;
                }
            }

            if (errorIsSatisfactory) {
                break;
            }

            spacing /= 2.0;
        }
        chunkStepSize[i] = spacing;
        chunkKnotCount[i] = knots;

        //Set up the vector of splines for this chunk
        double intervalLength = (knots - 1) * spacing;
        int numIntervals = ceil(chunkWidth / intervalLength);
        chunks[i].resize(numIntervals);

        //The testSpline is the first spline in the chunk
        chunks[i][0] = testSpline;
    }
}

void KBessel::computeChunkSplines(int indexOfLastNewChunk, int previouslyComputedChunks) {
    //Compute all the splines, make it parallel
#pragma omp parallel for schedule(dynamic) default(none) shared(indexOfLastNewChunk, previouslyComputedChunks)
    for (int i = previouslyComputedChunks; i <= indexOfLastNewChunk; i++) {
        //get left endpoint for this chunk
        //get spacing for this chunk
        //figure out how many intervals happen in this chunk
        //for each interval
        //  precompute K-Bessel values
        //  make spline
        //  add splines to the chunks vector

        //get left and right endpoints for this chunk
        double left = firstChunkLeftEndpoint + i * chunkWidth;

        //get spacing for this chunk
        double spacing = chunkStepSize[i];

        //figure out how many intervals happen in this chunk
        int knots = chunkKnotCount[i];
        double intervalLength = (knots - 1) * spacing;
        int numIntervals = ceil(chunkWidth / intervalLength);

        //start at j = 1 because j = 0 was already computed in the previous step
        for (int j = 1; j < numIntervals; j++) {
            //  precompute K-Bessel values
            vector<double> precompute(knots, 0.0);
            double intervalLeft = left + j * intervalLength;
            double intervalRight = intervalLeft + (knots - 1) * spacing;
            for (int k = 0; k < knots; k++) {
                double x = intervalLeft + k * spacing;
                double y = exactKBessel(x);
                precompute[k] = y;
            }

            //  make spline
            CubicSpline spline(precompute.begin(),
                               precompute.end(),
                               intervalLeft,
                               spacing,
                               approxDerivativeKBessel(intervalLeft),
                               approxDerivativeKBessel(intervalRight));

            //  add splines to the chunks vector
            chunks[i][j] = spline;
        }
    }
    precomputedRegionRightBound = firstChunkLeftEndpoint + chunkWidth * chunks.size();
}

void KBessel::setUpShrinkingChunkSplineComputation() {
    shrinkingChunks.resize(numberOfShrinkingChunks);
    shrinkingChunkStepSize.resize(numberOfShrinkingChunks);


#pragma omp parallel for schedule(dynamic) default(none)
    for (int newShrinkingChunk = numberOfShrinkingChunks - 1; newShrinkingChunk >= 0; newShrinkingChunk--) {

        double spacing = chunkStepSize[0];
        for (int j = newShrinkingChunk; j < shrinkingChunkStepSize.size(); j++) {
            if (shrinkingChunkStepSize[j] > 0) {
                spacing = shrinkingChunkStepSize[j];
                break;
            }
        }

        double left = precomputedRegionLeftBound
                + (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                * pow((1. * newShrinkingChunk)/numberOfShrinkingChunks, 2);

        CubicSpline testSpline;

        bool errorIsSatisfactory = false;
        while (!errorIsSatisfactory) {
            //generate one interpolant
            vector<double> precompute(SPLINE_KNOT_COUNT, 0.0);

            for (int i = 0; i < SPLINE_KNOT_COUNT; i++) {
                double x = left + i * spacing;
                double y = exactKBessel(x);
                precompute[i] = y;
            }
            double right = left + (SPLINE_KNOT_COUNT - 1) * spacing;

            testSpline = CubicSpline(precompute.begin(),
                                   precompute.end(),
                                   left,
                                   spacing,
                                   approxDerivativeKBessel(left),
                                   approxDerivativeKBessel(right));

            //check it for accuracy
            errorIsSatisfactory = true;
            for (int j = 1; j < 10; j++) {
                double x = left + spacing * j / 10;
                double exact = exactKBessel(x);
                double approx = testSpline(x);
                double error = relativeError(exact, approx);

                if (error > ABS_ERROR_CUTOFF) {
                    errorIsSatisfactory = false;
                    break;
                }
            }

            if (errorIsSatisfactory) {
                break;
            }

            spacing /= 2.0;
        }
        shrinkingChunkStepSize[newShrinkingChunk] = spacing;

        double subIntervalWidth = (SPLINE_KNOT_COUNT - 1) * spacing;
        int numberOfSubIntervals = ceil(
                (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                * (pow((1. * newShrinkingChunk + 1)/numberOfShrinkingChunks, 2)
                - pow((1. * newShrinkingChunk)/numberOfShrinkingChunks, 2)
                )
                /subIntervalWidth);


        shrinkingChunks[newShrinkingChunk].resize(numberOfSubIntervals);
        shrinkingChunks[newShrinkingChunk][0] = testSpline;
    }
}

void KBessel::computeShrinkingChunkSplines() {
//#pragma omp parallel for schedule(dynamic) default(none)
    for (int newShrinkingChunk = 0; newShrinkingChunk < numberOfShrinkingChunks; newShrinkingChunk++) {

        double left = precomputedRegionLeftBound + (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                * pow((1. * newShrinkingChunk)/numberOfShrinkingChunks, 2);

        double spacing = shrinkingChunkStepSize[newShrinkingChunk];

        double subIntervalWidth = (SPLINE_KNOT_COUNT - 1) * spacing;
        int numberOfSubIntervals = ceil(
                (firstChunkLeftEndpoint - precomputedRegionLeftBound)
                * (pow((1. * newShrinkingChunk + 1)/numberOfShrinkingChunks, 2)
                 - pow((1. * newShrinkingChunk)/numberOfShrinkingChunks, 2)
                )
                /subIntervalWidth);


        //i = 0 spline was already computed
#pragma omp parallel for default(none) shared(numberOfSubIntervals, left, subIntervalWidth, spacing, newShrinkingChunk)
        for (int i = 1; i < numberOfSubIntervals; i++) {
            vector<double> precompute(SPLINE_KNOT_COUNT, 0.0);
            double subIntervalLeft = left + i * subIntervalWidth;
            double subIntervalRight = left + (i + 1) * subIntervalWidth;
            for (int j = 0; j < SPLINE_KNOT_COUNT; j++) {
                double x = subIntervalLeft + j * spacing;
                double y = exactKBessel(x);
                precompute[j] = y;
            }

            CubicSpline spline(precompute.begin(),
                               precompute.end(),
                               subIntervalLeft,
                               spacing,
                               approxDerivativeKBessel(subIntervalLeft),
                               approxDerivativeKBessel(subIntervalRight));

            shrinkingChunks[newShrinkingChunk][i] = spline;
        }
    }

}

