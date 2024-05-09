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
    this->r = r;

    /*This is put here so that it can be initialized dynamically depending on r,
     * in case that becomes useful later.
     */
    firstChunkLeftEndpoint = 1.0;
    chunkWidth = 5.0;

    numberOfShrinkingChunks = 12;
    shrinkingChunkFirstWidth = (firstChunkLeftEndpoint - precomputeLowerBound)/pow(2, numberOfShrinkingChunks - 1);

    precomputedRegionLeftBound = precomputeLowerBound;
    precomputedRegionRightBound = precomputeLowerBound;

    if (precomputeLowerBound >= chunkWidth) {
        throw std::invalid_argument("precomputeLowerBound must be less than "
        + std::to_string(chunkWidth));
    }


    /**
     * Initialize the exact functionality
     */
#pragma omp parallel default(none) shared(threads)
    {
#pragma omp single
        threads = omp_get_max_threads();
    }

    K.reserve(threads);
    for (int i = 0; i < threads; i++) {
        K.push_back(new archtKBessel(r));
    }
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
    //TODO: Bake in logic for when K evaluates to 0, and also throw errors when it is evaluated below
    //the specified bound. Make this logic reflect in the precompute routine.0

    if (x >= zeroCutoff) {
        return 0;
    } else if (x < precomputedRegionLeftBound) {
        throw std::invalid_argument("x cannot be below 2*pi/A*1*Y0");
    } else if (x < firstChunkLeftEndpoint) { /* 2*pi/A*1*Y0 <= x < firstChunkLeftEndpoint */
        /*
         * The scheme here is that the interval [2*pi/A*1*Y0, firstChunkLeftEndpoint) = [a,b) is divided
         * into n exponentially shrinking pieces. The value w is chosen so that endpoints
         * a, a + w, a + 2w, a + 2^2w, ..., a + 2^(n-1)w = b
         * define the "shrinking chunks".
         * The 0th shrinking chunk is [a,a+w)
         * The (n-1)st shrinking chunk is [a + 2^(n-2)w, a + 2^(n-1)w = b)
         *
         * y is in the kth shrinking chunk if only if
         * y < a + 2^k w,  where k>= 0 is the smallest such k
         * This is equivalent to
         * log_2((y-a)/w) < k
         *
         * So we calculate k = next(log_2((y-a)/w)), and if k is less than 0 we make it 0
         */

        /*Do this branch first to avoid weirdness with log being negative*/
        if (x < precomputedRegionLeftBound + shrinkingChunkFirstWidth) {
            int subIndex = floor(
                    (x - precomputedRegionLeftBound) //distance from x to the left endpoint of it shrinking chunk
                    /((SPLINE_KNOT_COUNT - 1) * shrinkingChunkStepSize[0])); //width of spline interval
            return shrinkingChunks[0][subIndex](x);
        }

        int shrinkingChunkIndex = Auxiliary::next(log2((x - precomputedRegionLeftBound)/shrinkingChunkFirstWidth));

        double distToLeftEndpoint = x - (precomputedRegionLeftBound + pow(2, shrinkingChunkIndex - 1) * shrinkingChunkFirstWidth);
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


    //firstChunkLeftEndpoint = std::max(r/2.0, precomputedRegionLeftBound + 1);
    firstChunkLeftEndpoint = 1;
    shrinkingChunkFirstWidth = (firstChunkLeftEndpoint - precomputedRegionLeftBound)/pow(2, numberOfShrinkingChunks - 1);

    zeroCutoff = (1136 + PI*r/2.0*log2(E) - 0.5*log2(E) + 0.5*log2(PI/2))/log2(E);

    chunkStepSize.clear();
    shrinkingChunkStepSize.clear();

    chunks.clear();
    shrinkingChunks.clear();

    computedChunkCount = 0;

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

    if (newUpperBound < precomputedRegionLeftBound) {
        return;
    } else if (newUpperBound > zeroCutoff) {
        newUpperBound = zeroCutoff;
    }

    //Figure what new chunks need to be computed
    int indexOfLastNewChunk = ceil((newUpperBound - firstChunkLeftEndpoint) / chunkWidth) - 1;

    //Add new slots in the chunkStepSize and chunks vectors
    if (indexOfLastNewChunk >= chunkStepSize.size()) {
        chunkStepSize.resize(indexOfLastNewChunk + 1);
        chunkKnotCount.resize(indexOfLastNewChunk + 1);
        chunks.resize(indexOfLastNewChunk + 1);
    }

    //Compute spacing requirement for each new chunk
    //Do it in reverse since earlier chunks require spacing at least as fine as later chunks

    /*TODO: I need to parallelize this process since it takes quite a long time in a single thread.
     * Keep the logical structure, but give each thread its own continguous sequence of chunks.
     * Knowing that spacing decreases monotonically is the biggest speedup that we can't afford to lose.
     * Within each contiguous sequence of chunks we would expect the spacing to not change much, except for the
     * leftmost one.
     */
    int threads;

#pragma omp parallel default(none) shared(threads)
    {
#pragma omp single
        threads = omp_get_num_threads();
    }

#pragma omp parallel default(none) shared(indexOfLastNewChunk, threads)
    {
        int thread = omp_get_thread_num();
        int skip = (indexOfLastNewChunk - computedChunkCount)/threads;
        int start = indexOfLastNewChunk - thread * skip;
        int end = indexOfLastNewChunk - (thread + 1) * skip + 1;
        if (thread == threads - 1) {
            end = computedChunkCount;
        }

        double spacing = pow(2.0, -1);
        int knots = 0;
        for (int i = start; i >= end; i--) {
            bool errorIsSatisfactory = false;
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

                CubicSpline testSpline(precompute.begin(),
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
            chunkStepSize[computedChunkCount + i] = spacing;
            chunkKnotCount[computedChunkCount + i] = knots;
        }
    }
    
    //Compute all the splines, make it parallel
#pragma omp parallel for default(none) shared(indexOfLastNewChunk)
    for (int i = computedChunkCount; i <= indexOfLastNewChunk; i++) {
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
        chunks[i].resize(numIntervals);

        //for each interval
        for (int j = 0; j < numIntervals; j++) {
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
    computedChunkCount = chunks.size();
    precomputedRegionRightBound = firstChunkLeftEndpoint + chunkWidth * (computedChunkCount);


    /* ********************************************************************
     *
     * Precompute shrinking chunks, meaning values in the range [newLowerBound, chunkWidth)
     * This is only done once per update of spectral parameter r.
     *
     **********************************************************************/

    if (!shrinkingChunks.empty()) {
        return;
    }

    shrinkingChunks.resize(numberOfShrinkingChunks);
    shrinkingChunkStepSize.resize(numberOfShrinkingChunks);


    double spacing;
    for (int newShrinkingChunk = numberOfShrinkingChunks - 1; newShrinkingChunk >= 0; newShrinkingChunk--) {

        if (newShrinkingChunk == numberOfShrinkingChunks - 1) {
            spacing = chunkStepSize[0];
        } else {
            spacing = shrinkingChunkStepSize[newShrinkingChunk + 1];
        }

        double left;
        if (newShrinkingChunk == 0) {
            left = precomputedRegionLeftBound;
        } else {
            left = precomputedRegionLeftBound + pow(2, newShrinkingChunk - 1) * shrinkingChunkFirstWidth;
        }

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

            CubicSpline testSpline(precompute.begin(),
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
    }

#pragma omp parallel for default(none)
    for (int newShrinkingChunk = 0; newShrinkingChunk < numberOfShrinkingChunks; newShrinkingChunk++) {

        double left;
        double right;
        if (newShrinkingChunk == 0) {
            left = precomputedRegionLeftBound;
            right = precomputedRegionLeftBound + shrinkingChunkFirstWidth;
        } else {
            left = precomputedRegionLeftBound + pow(2, newShrinkingChunk - 1) * shrinkingChunkFirstWidth;
            right = precomputedRegionLeftBound + pow(2, newShrinkingChunk) * shrinkingChunkFirstWidth;
        }

        double spacing = shrinkingChunkStepSize[newShrinkingChunk];

        double subIntervalWidth = (SPLINE_KNOT_COUNT - 1) * spacing;
        int numberOfSubIntervals = ceil((right - left)/subIntervalWidth);

        shrinkingChunks[newShrinkingChunk].resize(numberOfSubIntervals);
        for (int i = 0; i < numberOfSubIntervals; i++) {
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

double KBessel::relativeError(double exact, double approx) {
    return abs(approx-exact)/(1 + abs(exact));
}

double KBessel::approxDerivativeKBessel(double x) {
    auto f = [this](double x) { return this->exactKBessel(x); };
    double dfdx = boost::math::differentiation::finite_difference_derivative(f, x);

    if (isnan(dfdx)) {
        std::cout << x << std::endl;
    }
    return dfdx;
}


