//
// Created by Eric Moss on 6/27/23.
//

#include <stdexcept>

#include "Auxilliary.h"
//#include <math.h>
#include "arb.h"
//#include "pari/pari.h"
//#include "flint/flint.h"
//#include "BianchiMaassComputer.h"


int Auxilliary::mod(int x, int modulus) {
    int answer = x % modulus;
    if (answer < 0) {
        return answer + abs(modulus);
    } else {
        return answer;
    }
}

int Auxilliary::next(long double x) {
    int answer = floor(x) + 1;
    return answer;
}

/*void Auxilliary::computeTheta(acb_t theta, int d, int bits) {
    acb_init(theta);
    if (Auxilliary::mod(-d,4) == 1) {
        //theta = 1/2 + I*sqrt(d)/2
        arb_t real, imag;
        arb_init(real);
        arb_init(imag);

        arb_set_si(real,1);
        arb_div_si(real,real,2, bits);

        arb_sqrt_ui(imag, d, bits);
        arb_div_si(imag, imag, 2, bits);

        acb_set_arb_arb(theta, real, imag);

        arb_clear(real);
        arb_clear(imag);
    } else if (Auxilliary::mod(-d,4) == 2 || Auxilliary::mod(-d,4) == 3) {
        //theta = I*sqrt(d)
        arb_t real, imag;
        arb_init(real);
        arb_init(imag);

        arb_zero(real);

        arb_sqrt_ui(imag, d, bits);

        acb_set_arb_arb(theta, real, imag);

        arb_clear(real);
        arb_clear(imag);
    } else {
        throw std::invalid_argument("d is incorrect");
    }
}

*//*
 * If the interval [-scale, scale] centered at 0 is translated to the rest of the real line,
 * where in the interval does x fall?
 * The answer is set to the variable answer.
 *
 * If x is positive, the right endpoint is favored.
 * If x is negative, the left endpoint is favored.
 * If x overlaps with an endpoint, we try to deal with it in a natural way, but there is no guarantee
 * that answer is a ball totally within [-scale, scale].
 *//*
void Auxilliary::intervalReduce(arb_struct *answer, const arb_struct *x, const arb_struct *scale, const int bits) {
    arb_t localCopyOfX;
    arb_init(localCopyOfX);
    arb_set(localCopyOfX,x);

    arb_t copyOfScale;
    arb_init(copyOfScale);
    arb_set(copyOfScale, scale);

    arb_t negScale;
    arb_init(negScale);
    arb_neg(negScale, scale);

    arb_t temp;
    arb_init(temp);

    if (arb_gt(localCopyOfX, copyOfScale)) {
        //use g = x - 2*scale*ceil(x/(2*scale) - 1/2)
        arb_sub(temp, localCopyOfX, copyOfScale, bits);
        arb_div_si(temp, temp, 2, bits);
        arb_div(temp, temp, copyOfScale, bits);
        arb_ceil(temp, temp, bits);
        arb_mul(temp, temp, copyOfScale, bits);
        arb_mul_si(temp,temp, 2, bits);
        arb_sub(answer, localCopyOfX, temp, bits);
    } else if (arb_lt(localCopyOfX, negScale)) {
        //use f = =x - 2*scale*floor(x/(2*scale) + 1/2)
        arb_add(temp, localCopyOfX, copyOfScale, bits);
        arb_div_si(temp, temp, 2, bits);
        arb_div(temp, temp, copyOfScale, bits);
        arb_floor(temp, temp, bits);
        arb_mul(temp, temp, copyOfScale, bits);
        arb_mul_si(temp,temp, 2, bits);
        arb_sub(answer, localCopyOfX, temp, bits);
    } else {
        //do something??
    }
    arb_clear(negScale);
    arb_clear(temp);
    arb_clear(localCopyOfX);
    arb_clear(copyOfScale);
}

*//*
 * Like interval reduce, but for the region defined by [-1,1]*v
 *
 * Project toReduce onto v to get the scalar projection arb_t s
 * final answer is toReduce - 2ceil((s-1)/2)*v
 *//*
void Auxilliary::vectorReduce(acb_struct *answer, const acb_struct *v, const acb_t toReduce, int bits) {
    acb_t copyOfToReduce;
    acb_init(copyOfToReduce);
    acb_set(copyOfToReduce, toReduce);

    acb_t copyOfV;
    acb_init(copyOfV);
    acb_set(copyOfV, v);

    arb_t vx, vy, rx, ry, s, denominator, temp;
    arb_init(vx);
    arb_init(vy);
    arb_init(rx);
    arb_init(ry);
    arb_init(s);
    arb_init(denominator);
    arb_init(temp);

    acb_get_real(vx,v);
    acb_get_imag(vy, v);
    acb_get_real(rx, copyOfToReduce);
    acb_get_imag(ry, copyOfToReduce);

    //Scalar projection of toReduce on to v
    arb_mul(s, vx, rx, bits);
    arb_addmul(s, vy, ry, bits);
    arb_sqr(denominator, vx, bits);
    arb_addmul(denominator, vy, vy, bits);
    arb_div(s, s, denominator, bits);

    //Do the rounding thing
    //temp = -2ceil((s-1)/2)
    arb_sub_si(temp, s, 1, bits);
    arb_div_si(temp, temp, 2, bits);
    arb_ceil(temp, temp, bits);
    arb_mul_si(temp, temp, 2, bits);
    arb_neg(temp, temp);

    //std::cout << "vector scaling is " << arb_get_str(temp, bits, 0) << std::endl;

    //Subtract the proper multiple of v from toReduce
    acb_mul_arb(answer, copyOfV, temp, bits);
    acb_add(answer, answer, copyOfToReduce, bits);

    acb_clear(copyOfToReduce);
    arb_clear(vx);
    arb_clear(vy);
    arb_clear(rx);
    arb_clear(ry);
    arb_clear(s);
    arb_clear(denominator);
    arb_clear(temp);
    acb_clear(copyOfV);
}

void Auxilliary::acbPlusArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_add(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxilliary::acbTimesArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_mul(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxilliary::acbDivArb(acb_struct *answer, const acb_struct *z, const arb_struct *x, int bits) {
    arb_t zero;
    arb_init(zero);

    arb_t copyOfX;
    arb_init(copyOfX);
    arb_set(copyOfX, x);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_set_arb_arb(complexEmbedding, copyOfX, zero);
    acb_div(answer, copyOfZ, complexEmbedding, bits);

    arb_clear(zero);
    acb_clear(complexEmbedding);
    arb_clear(copyOfX);
    acb_clear(copyOfZ);
}

void Auxilliary::acbDivArb(acb_struct *answer, const acb_struct *z, int x, int bits) {
    acb_t complexEmbedding;
    acb_init(complexEmbedding);

    acb_t copyOfZ;
    acb_init(copyOfZ);
    acb_set(copyOfZ, z);

    acb_set_si_si(complexEmbedding, x, 0);
    acb_div(answer, copyOfZ, complexEmbedding, bits);

    acb_clear(complexEmbedding);
    acb_clear(copyOfZ);
}

void Auxilliary::traceProduct(acb_struct *answer, const acb_struct *z, const acb_struct *w, int bits) {
    acb_t copyOfZ, copyOfW;
    acb_init(copyOfZ);
    acb_init(copyOfW);

    acb_set(copyOfZ, z);
    acb_set(copyOfW, w);

    acb_mul(answer, copyOfZ, copyOfW, bits);
    acb_conj(answer, answer);
    acb_addmul(answer, copyOfZ, copyOfW, bits);

    acb_clear(copyOfW);
    acb_clear(copyOfZ);
}*/

double Auxilliary::imagTheta(int d) {
    if (mod(-d, 4) == 1) {
        return sqrt(1.0*d)/2;
    } else {
        return sqrt(1.0*d);
    }
}

double Auxilliary::multiPrecisionSummation(const std::vector<double> &numbers) {
    arb_t sum;
    arb_t temp;
    arb_init(sum);
    arb_init(temp);

    for (double number : numbers) {
        arb_set_d(temp, number);
        arb_add(sum, sum, temp, 500);
    }
    double answer = std::stod(arb_get_str(sum, 60, ARB_STR_NO_RADIUS));
    arb_clear(sum);
    arb_clear(temp);
    flint_cleanup();
    return answer;
}



/*double Auxilliary::computeVolumeOfFD(int d) {
    pari_init(1000000,2);
    std::string arg = "L = lfuncreate(x^2+";
    arg.append(std::to_string(d));
    arg.append(");");
    gp_read_str(arg.c_str());
    arg = "lfun(L,2)";
    GEN res = gp_read_str(arg.c_str());
    char* out = GENtostr(res);
    pari_close();
    double zeta = std::stod(std::string(out));
    return pow(abs(d), 3.0/2)*zeta/(4*pow(Auxilliary::pi, 2));
}*/

//double Auxilliary::pi = std::numbers::pi_v<double>;


