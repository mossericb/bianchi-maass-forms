//
// Created by Eric Moss on 12/13/23.
//

#include "TrigApproximator.h"

TrigApproximator::TrigApproximator() {

    for ( ; ; D--) {
        bool errorIsSatisfactory = true;
        spacing = pow(2,D);
        precomputedValues.clear();

        int numberOfPrecomputedValues = piOverTwo/spacing + 1;
        precomputedValues.resize(numberOfPrecomputedValues);

        for (int i = 0; ; i++) {
            double x = spacing * i;
            if (x > piOverTwo) {
                break;
            }
            precomputedValues[i] = sinExact(x);
        }
        precomputedValues[precomputedValues.size() - 1] = 1.00000000000000000000;

        //Check at 10 random points to see if the error is satisfactory
        for (int i = 0; i < 10; i++) {
            double x = randab(-100, 100);
            double sinExactTest = sinExact(x);
            double sinApproxText = sinApprox(x);
            double sinError  = abs((sinExactTest - sinApproxText)/sinExactTest);

            double cosExactTest = cosExact(x);
            double cosApproxText = cosApprox(x);
            double cosError  = abs((cosExactTest - cosApproxText)/cosExactTest);

            if (sinError > MAX_ERROR || cosError > MAX_ERROR) {
                errorIsSatisfactory = false;
                break;
            }
        }
        if (!errorIsSatisfactory) {
            continue;
        }
        //Now we check 100 more times
        for (int i = 0; i < 100; i++) {
            double x = randab(-100, 100);
            double sinExactTest = sinExact(x);
            double sinApproxText = sinApprox(x);
            double sinError  = abs((sinExactTest - sinApproxText)/sinExactTest);

            double cosExactTest = cosExact(x);
            double cosApproxText = cosApprox(x);
            double cosError  = abs((cosExactTest - cosApproxText)/cosExactTest);

            if (sinError > MAX_ERROR || cosError > MAX_ERROR) {
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

double TrigApproximator::sinExact(const double &x) {
    arb_t theta;
    arb_t ans;
    arb_init(theta);
    arb_init(ans);
    arb_set_d(theta,x);

    for (int i = 1; ; i++) {
        arb_sin(ans, theta, bits*i);
        if (arb_rel_accuracy_bits(ans) > 53) {
            break;
        }
    }

    double answer;
    try {
        answer = std::stod(arb_get_str(ans, 50, ARB_STR_NO_RADIUS));
    } catch (...) {
        answer = 0;
    }

    arb_clear(theta);
    arb_clear(ans);

    flint_cleanup();

    return answer;
}

double TrigApproximator::cosExact(const double &x) {
    arb_t theta;
    arb_t ans;
    arb_init(theta);
    arb_init(ans);
    arb_set_d(theta,x);

    for (int i = 1; ; i++) {
        arb_cos(ans, theta, bits*i);
        if (arb_rel_accuracy_bits(ans) > 53) {
            break;
        }
    }

    double answer;
    try {
        answer = std::stod(arb_get_str(ans, 50, ARB_STR_NO_RADIUS));
    } catch (...) {
        answer = 0;
    }

    arb_clear(theta);
    arb_clear(ans);

    flint_cleanup();

    return answer;
}

double TrigApproximator::sinApprox(const double &x) {
    double reduce = x - twoPi * floor(x/twoPi); //Do I need to do this using arb?

    if (reduce < pi) {
        if (reduce < piOverTwo) {
            // 0 <= reduce < pi/2
            return interpolate(reduce);
        } else {
            // pi/2 <= reduce < pi
            return interpolate(pi - reduce);
        }
    } else {
        if (reduce < threePiOverTwo) {
            //pi <= reduce < 3pi/2
            return -interpolate(reduce - pi);
        } else {
            //3pi/2 <= reduce < 2pi
            return -interpolate(twoPi - reduce);
        }
    }
}

double TrigApproximator::cosApprox(const double &x) {
    double reduce = x - twoPi * floor(x/twoPi); //Do I need to do this using arb?

    if (reduce < pi) {
        if (reduce < piOverTwo) {
            // 0 <= reduce < pi/2
            return interpolate(piOverTwo - reduce);
        } else {
            // pi/2 <= reduce < pi
            return -interpolate(reduce - piOverTwo);
        }
    } else {
        if (reduce < threePiOverTwo) {
            //pi <= reduce < 3pi/2
            return -interpolate(threePiOverTwo - reduce);
        } else {
            //3pi/2 <= reduce < 2pi
            return interpolate(reduce - threePiOverTwo);
        }
    }
}

/**
 *
 * @param x Is in the interval [0, pi/2)
 * @return An approximation to sin(x) using Lagrange interpolation on precomputed values.
 */
double TrigApproximator::interpolate(const double &x) {
    double measure = x/spacing;
    int left = std::floor(measure);

    int index1;
    int index2;
    int index3;

    if (left == 0) {
        index1 = 0;
        index2 = 1;
        index3 = 2;
    } else if (left == precomputedValues.size() - 2) {
        index1 = left - 1;
        index2 = left;
        index3 = left + 1;
    } else if (x - left*spacing < spacing/2) { //closer to left half of the interval it's in
        index1 = left - 1;
        index2 = left;
        index3 = left + 1;
    } else {
        index1 = left;
        index2 = left + 1;
        index3 = left + 2;
    }

    double x1 = index1 * spacing;
    double x2 = index2 * spacing;
    double x3;
    if (index3 == precomputedValues.size() - 1) {
        x3 = piOverTwo;
    } else {
        x3 = index3 * spacing;
    }

    double y1 = precomputedValues[index1];
    double y2 = precomputedValues[index2];
    double y3 = precomputedValues[index3];

    double diff1 = x - x1;
    double diff2 = x - x2;
    double diff3 = x - x3;

    double answer = y1*diff2*diff3/2;
    answer += -y2*diff1*diff3;
    answer += y3*diff1*diff2/2;
    answer /= pow(spacing,2);
    return answer;
}

double TrigApproximator::rand01() {
    return std::rand() / (RAND_MAX + 1.0);
}

double TrigApproximator::randab(double a, double b) {
    return a + rand01()*(b-a);
}


