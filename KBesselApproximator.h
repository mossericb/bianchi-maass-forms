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

using std::map, std::vector;
using CubicSpline = boost::math::interpolators::cardinal_cubic_b_spline<double>;

class KBesselApproximator {
public:
    KBesselApproximator(double precomputeLowerBound);

    void setRAndPrecompute(double r, double precomputeUpperBound);
    void extendPrecomputedRange(double newUpperBound);
    void setRAndClear(double r);
    double approxKBessel(const double x);
    double getR() { return this->r; }
    void printTiming();
    std::vector<double> maximize();
    void runTest();

private:
    /*
     * TODO New computation paradigm.
     * Let lb = 2*pi/A * 1 * Y_0
     * def approx(x):
     *     if (x < lb):
     *         return exact(x)
     *     else if (x >= lb and x < zeroCutoff):
     *         return splines(x)
     *     else:
     *         return 0;
     *
     * def splines(x):
     *     if (x < r):
     *         return shrinkingChunks(x)
     *     else (x >= r):
     *         return chunks(x)
    */

    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;
    double zeroCutoff;

    int computedChunkCount;

    double r;

    vector<double> chunkStepSize;
    vector<double> shrinkingChunkStepSize;

    int numberOfShrinkingChunks;
    double shrinkingChunkFirstWidth;
    double chunkWidth;
    double firstChunkLeftEndpoint;
    static constexpr int SPLINE_KNOT_COUNT = 128;
    static constexpr double ABS_ERROR_CUTOFF = 5e-16;
    static constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
    static constexpr double E = 2.7182818284590452353602874713526624977572470936999595749669676277;

    vector<vector<CubicSpline>> chunks;
    vector<int> chunkKnotCount;
    vector<vector<CubicSpline>> shrinkingChunks;

    double relativeError(double exact, double approx);


};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
