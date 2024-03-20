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
    KBesselApproximator(int bitsOfPrecision = 53);

    void setRAndPrecompute(double r, double precomputeLowerBound, double precomputeUpperBound);
    void extendPrecomputedRange(double newLowerBound, double newUpperBound);
    void setRAndClear(double newR);
    double approxKBessel(const double x);
    double getR() { return this->r; }
    void printTiming();
    std::vector<double> maximize();

private:
    /*
     *
     * TODO Make it so that I can initialize the class so that the minimum value approximated is
     * 2*pi/A * 1 * Y_0
     * All other K-Bessel function calls can be computed as a 1-off exactly!
    */
    int prec;

    double precomputedRegionLeftBound;
    double precomputedRegionRightBound;
    int computedChunkCount;
    int computedShrinkingChunkCount;

    double r;

    vector<double> chunkStepSize;
    vector<double> shrinkingChunkStepSize;
    double getSpacing(double x);

    static constexpr double CHUNK_WIDTH = 1.0;
    static constexpr int SPLINE_KNOT_COUNT = 128;
    static constexpr double ABS_ERROR_CUTOFF = 5.0e-15;

    vector<vector<CubicSpline>> chunks;
    vector<vector<CubicSpline>> shrinkingChunks;

    double fineIntervalWidth;
    double coarseIntervalWidth;

    void runTest();
};


#endif //COMPUTINGKBESSEL_KBESSELAPPROXIMATOR_H
