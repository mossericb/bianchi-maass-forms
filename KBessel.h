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
#include <map>
#include "archtKBessel.h"
#include "KBesselReal.h"

using std::map, std::vector;
using CubicSpline = boost::math::interpolators::cardinal_cubic_b_spline<double>;

class KBessel {
public:
    KBessel(double precomputeLowerBound, double lambda);
    ~KBessel();

    void setLambdaAndPrecompute(double lambda, double precomputeUpperBound);
    void extendPrecomputedRange(double newUpperBound);
    void setLambdaAndClear(double lambda);
    double approxKBessel(double x);
    double exactKBessel(double x);
    double approxDerivativeKBessel(double x);

    double getLambda() { return this->lambda; }

private:
    bool withinError(double exact, double approx);

    void setUpChunkSplineComputation(int indexOfLastNewChunk, int previouslyComputedChunks);
    void computeChunkSplines(int indexOfLastNewChunk, int previouslyComputedChunks);

    void setUpShrinkingChunkSplineComputation();
    void computeShrinkingChunkSplines();

    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;
    double zeroCutoff;

    bool small;
    double lambda;

    int numberOfShrinkingChunks;
    double chunkWidth;
    double firstChunkLeftEndpoint;
    static constexpr int SPLINE_KNOT_COUNT = 128;
    static constexpr double ABS_ERROR_CUTOFF = 5.0e-15;

    vector<double> chunkStepSize;
    vector<double> shrinkingChunkStepSize;

    vector<vector<CubicSpline>> chunks;
    vector<int> chunkKnotCount;
    vector<vector<CubicSpline>> shrinkingChunks;

    int threads;
    std::vector<archtKBessel*> K;
    std::vector<KBesselReal*> KReal;

    static constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
    static constexpr double E = 2.7182818284590452353602874713526624977572470936999595749669676277;
};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
