//
// Created by Eric Moss on 6/26/23.
//

#include "Quaternion.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

const SL2C Quaternion::S = {0,-1,1,0};
const vector<complex<double>> Quaternion::d19Alphas = {(0.0)/(1.0),
                                                                 (getTheta(19))/(2.0),
                                                                 (getTheta(19) - 1.0)/(2.0)};
const vector<SL2C> Quaternion::d19MAlphas = {{0.0, 1.0, -1.0, 0.0},
                                                  {-getTheta(19) + 1.0, -2.0, -2.0, getTheta(19)},
                                                  {-getTheta(19), -2.0, -2.0, getTheta(19) - 1.0}};
const vector<complex<double>> Quaternion::d19Translators = {0.0,
                                                                      1.0,
                                                                      -1.0,
                                                                      getTheta(19),
                                                                      -getTheta(19),
                                                                      getTheta(19) - 1.0,
                                                                      -getTheta(19) + 1.0};

const vector<complex<double>> Quaternion::d43Alphas = {(0.0)/(1.0),
                                                                 (getTheta(43))/(2.0),
                                                                 (getTheta(43) - 1.0)/(2.0),
                                                                 (getTheta(43))/(3.0),
                                                                 (-getTheta(43))/(3.0),
                                                                 (-getTheta(43) + 1.0)/(3.0),
                                                                 (getTheta(43) - 1.0)/(3.0),
                                                                 (getTheta(43) + 1.0)/(3.0),
                                                                 (-getTheta(43) - 1.0)/(3.0)};
const vector<SL2C> Quaternion::d43MAlphas = {{0.0, 1.0, -1.0, 0.0},
                                                  {-getTheta(43) + 1.0, -5.0, -2.0, getTheta(43)},
                                                  {-getTheta(43), -5.0, -2.0, getTheta(43) - 1.0},
                                                  {getTheta(43) - 1.0, 4.0, -3.0, getTheta(43)},
                                                  {-getTheta(43) + 1.0, 4.0, -3.0, -getTheta(43)},
                                                  {-getTheta(43), 4.0, -3.0, -getTheta(43) + 1.0},
                                                  {getTheta(43), 4.0, -3.0, getTheta(43) - 1.0},
                                                  {-getTheta(43) - 1.0, getTheta(43) - 3.0, -3.0, getTheta(43) + 1.0},
                                                  {getTheta(43) + 1.0, getTheta(43) - 3.0, -3.0, -getTheta(43) - 1.0}};
const vector<complex<double>> Quaternion::d43Translators = {0.0,
                                                                      1.0,
                                                                      -1.0,
                                                                      getTheta(43),
                                                                      -getTheta(43),
                                                                      getTheta(43) - 1.0,
                                                                      -getTheta(43) + 1.0};

const vector<complex<double>> Quaternion::d67Alphas = {(0.0)/(1.0),
                                                                 (getTheta(67))/(2.0),
                                                                 (getTheta(67) - 1.0)/(2.0),
                                                                 (getTheta(67))/(3.0),
                                                                 (-getTheta(67))/(3.0),
                                                                 (-getTheta(67) + 1.0)/(3.0),
                                                                 (getTheta(67) - 1.0)/(3.0),
                                                                 (getTheta(67) + 1.0)/(3.0),
                                                                 (-getTheta(67) - 1.0)/(3.0),
                                                                 (getTheta(67) + 1.0)/(4.0),
                                                                 (-getTheta(67) - 1.0)/(4.0),
                                                                 (-getTheta(67) + 2.0)/(4.0),
                                                                 (getTheta(67) - 2.0)/(4.0),
                                                                 (getTheta(67) - 1.0)/(4.0),
                                                                 (-getTheta(67) + 1.0)/(4.0),
                                                                 (getTheta(67))/(4.0),
                                                                 (-getTheta(67))/(4.0),
                                                                 (getTheta(67) - 3.0)/(getTheta(67) + 2.0),
                                                                 (-getTheta(67) + 3.0)/(getTheta(67) + 2.0),
                                                                 (getTheta(67) - 7.0)/(getTheta(67) + 2.0),
                                                                 (-getTheta(67) + 7.0)/(getTheta(67) + 2.0),
                                                                 (-getTheta(67) - 2.0)/(-getTheta(67) + 3.0),
                                                                 (getTheta(67) + 2.0)/(-getTheta(67) + 3.0),
                                                                 (-getTheta(67) - 6.0)/(-getTheta(67) + 3.0),
                                                                 (getTheta(67) + 6.0)/(-getTheta(67) + 3.0)};
const vector<SL2C> Quaternion::d67MAlphas = {{0.0, 1.0, -1.0, 0.0},
                                                  {-getTheta(67) + 1.0, -8.0, -2.0, getTheta(67)},
                                                  {-getTheta(67), -8.0, -2.0, getTheta(67) - 1.0},
                                                  {getTheta(67) - 1.0, 6.0, -3.0, getTheta(67)},
                                                  {-getTheta(67) + 1.0, 6.0, -3.0, -getTheta(67)},
                                                  {-getTheta(67), 6.0, -3.0, -getTheta(67) + 1.0},
                                                  {getTheta(67), 6.0, -3.0, getTheta(67) - 1.0},
                                                  {-getTheta(67) - 1.0, getTheta(67) - 5.0, -3.0, getTheta(67) + 1.0},
                                                  {getTheta(67) + 1.0, getTheta(67) - 5.0, -3.0, -getTheta(67) - 1.0},
                                                  {getTheta(67) - 2.0, 5.0, -4.0, getTheta(67) + 1.0},
                                                  {-getTheta(67) + 2.0, 5.0, -4.0, -getTheta(67) - 1.0},
                                                  {-getTheta(67) - 1.0, 5.0, -4.0, -getTheta(67) + 2.0},
                                                  {getTheta(67) + 1.0, 5.0, -4.0, getTheta(67) - 2.0},
                                                  {-getTheta(67), -4.0, -4.0, getTheta(67) - 1.0},
                                                  {getTheta(67), -4.0, -4.0, -getTheta(67) + 1.0},
                                                  {-getTheta(67) + 1.0, -4.0, -4.0, getTheta(67)},
                                                  {getTheta(67) - 1.0, -4.0, -4.0, -getTheta(67)},
                                                  {-getTheta(67) + 7.0, -getTheta(67) - 6.0, -getTheta(67) - 2.0, getTheta(67) - 3.0},
                                                  {getTheta(67) - 7.0, -getTheta(67) - 6.0, -getTheta(67) - 2.0, -getTheta(67) + 3.0},
                                                  {-getTheta(67) + 3.0, -getTheta(67) - 6.0, -getTheta(67) - 2.0, getTheta(67) - 7.0},
                                                  {getTheta(67) - 3.0, -getTheta(67) - 6.0, -getTheta(67) - 2.0, -getTheta(67) + 7.0},
                                                  {-getTheta(67) - 6.0, -getTheta(67) + 7.0, -getTheta(67) + 3.0, getTheta(67) + 2.0},
                                                  {getTheta(67) + 6.0, -getTheta(67) + 7.0, -getTheta(67) + 3.0, -getTheta(67) - 2.0},
                                                  {-getTheta(67) - 2.0, -getTheta(67) + 7.0, -getTheta(67) + 3.0, getTheta(67) + 6.0},
                                                  {getTheta(67) + 2.0, -getTheta(67) + 7.0, -getTheta(67) + 3.0, -getTheta(67) - 6.0}};
const vector<complex<double>> Quaternion::d67Translators = {0.0,
                                                                      1.0,
                                                                      -1.0,
                                                                      getTheta(67),
                                                                      -getTheta(67),
                                                                      getTheta(67) - 1.0,
                                                                      -getTheta(67) + 1.0};

const vector<complex<double>> Quaternion::d163Alphas = {(0.0)/(1.0),
                                                                  (getTheta(163))/(2.0),
                                                                  (getTheta(163) - 1.0)/(2.0),
                                                                  (getTheta(163))/(3.0),
                                                                  (-getTheta(163))/(3.0),
                                                                  (-getTheta(163) + 1.0)/(3.0),
                                                                  (getTheta(163) - 1.0)/(3.0),
                                                                  (getTheta(163) + 1.0)/(3.0),
                                                                  (-getTheta(163) - 1.0)/(3.0),
                                                                  (3.0*getTheta(163) + 2.0)/(7.0),
                                                                  (-3.0*getTheta(163) - 2.0)/(7.0),
                                                                  (getTheta(163) + 1.0)/(4.0),
                                                                  (-getTheta(163) - 1.0)/(4.0),
                                                                  (-getTheta(163) + 2.0)/(4.0),
                                                                  (getTheta(163) - 2.0)/(4.0),
                                                                  (getTheta(163) - 1.0)/(4.0),
                                                                  (-getTheta(163) + 1.0)/(4.0),
                                                                  (getTheta(163))/(4.0),
                                                                  (-getTheta(163))/(4.0),
                                                                  (2.0*getTheta(163) - 3.0)/(5.0),
                                                                  (-2.0*getTheta(163) + 3.0)/(5.0),
                                                                  (getTheta(163) - 2.0)/(5.0),
                                                                  (-getTheta(163) + 2.0)/(5.0),
                                                                  (getTheta(163) + 2.0)/(5.0),
                                                                  (-getTheta(163) - 2.0)/(5.0),
                                                                  (-2.0*getTheta(163) + 1.0)/(5.0),
                                                                  (2.0*getTheta(163) - 1.0)/(5.0),
                                                                  (2.0*getTheta(163) - 2.0)/(5.0),
                                                                  (-2.0*getTheta(163) + 2.0)/(5.0),
                                                                  (-2.0*getTheta(163))/(5.0),
                                                                  (2.0*getTheta(163))/(5.0),
                                                                  (2.0*getTheta(163) + 1.0)/(5.0),
                                                                  (-2.0*getTheta(163) - 1.0)/(5.0),
                                                                  (getTheta(163) + 1.0)/(5.0),
                                                                  (-getTheta(163) - 1.0)/(5.0),
                                                                  (getTheta(163) - 1.0)/(5.0),
                                                                  (-getTheta(163) + 1.0)/(5.0),
                                                                  (getTheta(163))/(5.0),
                                                                  (-getTheta(163))/(5.0),
                                                                  (getTheta(163) + 2.0)/(6.0),
                                                                  (-getTheta(163) - 2.0)/(6.0),
                                                                  (-getTheta(163) + 3.0)/(6.0),
                                                                  (getTheta(163) - 3.0)/(6.0),
                                                                  (getTheta(163) - 2.0)/(6.0),
                                                                  (-getTheta(163) + 2.0)/(6.0),
                                                                  (getTheta(163) + 1.0)/(6.0),
                                                                  (-getTheta(163) - 1.0)/(6.0),
                                                                  (getTheta(163) - 1.0)/(6.0),
                                                                  (-getTheta(163) + 1.0)/(6.0),
                                                                  (-getTheta(163))/(6.0),
                                                                  (getTheta(163))/(6.0),
                                                                  (-17.0)/(-getTheta(163) + 1.0),
                                                                  (17.0)/(-getTheta(163) + 1.0),
                                                                  (-12.0)/(-getTheta(163) + 1.0),
                                                                  (12.0)/(-getTheta(163) + 1.0),
                                                                  (17.0)/(getTheta(163)),
                                                                  (-17.0)/(getTheta(163)),
                                                                  (12.0)/(getTheta(163)),
                                                                  (-12.0)/(getTheta(163)),
                                                                  (-getTheta(163) - 16.0)/(-getTheta(163) + 2.0),
                                                                  (getTheta(163) + 16.0)/(-getTheta(163) + 2.0),
                                                                  (12.0)/(-getTheta(163) + 2.0),
                                                                  (-12.0)/(-getTheta(163) + 2.0),
                                                                  (12.0)/(getTheta(163) + 1.0),
                                                                  (-12.0)/(getTheta(163) + 1.0),
                                                                  (getTheta(163) - 17.0)/(getTheta(163) + 1.0),
                                                                  (-getTheta(163) + 17.0)/(getTheta(163) + 1.0),
                                                                  (-getTheta(163) - 10.0)/(-getTheta(163) + 3.0),
                                                                  (getTheta(163) + 10.0)/(-getTheta(163) + 3.0),
                                                                  (-getTheta(163) - 15.0)/(-getTheta(163) + 3.0),
                                                                  (getTheta(163) + 15.0)/(-getTheta(163) + 3.0),
                                                                  (-7.0)/(-getTheta(163) + 3.0),
                                                                  (7.0)/(-getTheta(163) + 3.0),
                                                                  (-getTheta(163) - 17.0)/(-getTheta(163) + 3.0),
                                                                  (getTheta(163) + 17.0)/(-getTheta(163) + 3.0),
                                                                  (7.0)/(getTheta(163) + 2.0),
                                                                  (-7.0)/(getTheta(163) + 2.0),
                                                                  (-getTheta(163) + 18.0)/(getTheta(163) + 2.0),
                                                                  (getTheta(163) - 18.0)/(getTheta(163) + 2.0),
                                                                  (getTheta(163) - 11.0)/(getTheta(163) + 2.0),
                                                                  (-getTheta(163) + 11.0)/(getTheta(163) + 2.0),
                                                                  (getTheta(163) - 16.0)/(getTheta(163) + 2.0),
                                                                  (-getTheta(163) + 16.0)/(getTheta(163) + 2.0),
                                                                  (getTheta(163) + 3.0)/(7.0),
                                                                  (-getTheta(163) - 3.0)/(7.0),
                                                                  (2.0*getTheta(163) - 1.0)/(7.0),
                                                                  (-2.0*getTheta(163) + 1.0)/(7.0),
                                                                  (2.0*getTheta(163) - 3.0)/(7.0),
                                                                  (-2.0*getTheta(163) + 3.0)/(7.0),
                                                                  (-2.0*getTheta(163) - 1.0)/(7.0),
                                                                  (2.0*getTheta(163) + 1.0)/(7.0),
                                                                  (-getTheta(163) - 4.0)/(-getTheta(163) + 4.0),
                                                                  (getTheta(163) + 4.0)/(-getTheta(163) + 4.0),
                                                                  (getTheta(163) + 16.0)/(-getTheta(163) + 4.0),
                                                                  (-getTheta(163) - 16.0)/(-getTheta(163) + 4.0),
                                                                  (getTheta(163) - 5.0)/(getTheta(163) + 3.0),
                                                                  (-getTheta(163) + 5.0)/(getTheta(163) + 3.0),
                                                                  (-getTheta(163) + 17.0)/(getTheta(163) + 3.0),
                                                                  (getTheta(163) - 17.0)/(getTheta(163) + 3.0)};
const vector<SL2C> Quaternion::d163MAlphas = {{0.0, 1.0, -1.0, 0.0},
                                                   {-getTheta(163) + 1.0, -20.0, -2.0, getTheta(163)},
                                                   {-getTheta(163), -20.0, -2.0, getTheta(163) - 1.0},
                                                   {getTheta(163) - 1.0, 14.0, -3.0, getTheta(163)},
                                                   {-getTheta(163) + 1.0, 14.0, -3.0, -getTheta(163)},
                                                   {-getTheta(163), 14.0, -3.0, -getTheta(163) + 1.0},
                                                   {getTheta(163), 14.0, -3.0, getTheta(163) - 1.0},
                                                   {-getTheta(163) - 1.0, getTheta(163) - 13.0, -3.0, getTheta(163) + 1.0},
                                                   {getTheta(163) + 1.0, getTheta(163) - 13.0, -3.0, -getTheta(163) - 1.0},
                                                   {-3.0*getTheta(163) - 2.0, 3.0*getTheta(163) - 52.0, -7.0, 3.0*getTheta(163) + 2.0},
                                                   {3.0*getTheta(163) + 2.0, 3.0*getTheta(163) - 52.0, -7.0, -3.0*getTheta(163) - 2.0},
                                                   {getTheta(163) - 2.0, 11.0, -4.0, getTheta(163) + 1.0},
                                                   {-getTheta(163) + 2.0, 11.0, -4.0, -getTheta(163) - 1.0},
                                                   {-getTheta(163) - 1.0, 11.0, -4.0, -getTheta(163) + 2.0},
                                                   {getTheta(163) + 1.0, 11.0, -4.0, getTheta(163) - 2.0},
                                                   {-getTheta(163), -10.0, -4.0, getTheta(163) - 1.0},
                                                   {getTheta(163), -10.0, -4.0, -getTheta(163) + 1.0},
                                                   {-getTheta(163) + 1.0, -10.0, -4.0, getTheta(163)},
                                                   {getTheta(163) - 1.0, -10.0, -4.0, -getTheta(163)},
                                                   {-getTheta(163) + 2.0, -getTheta(163) - 15.0, -5.0, 2.0*getTheta(163) - 3.0},
                                                   {getTheta(163) - 2.0, -getTheta(163) - 15.0, -5.0, -2.0*getTheta(163) + 3.0},
                                                   {-2.0*getTheta(163) + 3.0, -getTheta(163) - 15.0, -5.0, getTheta(163) - 2.0},
                                                   {2.0*getTheta(163) - 3.0, -getTheta(163) - 15.0, -5.0, -getTheta(163) + 2.0},
                                                   {2.0*getTheta(163) - 1.0, -getTheta(163) + 17.0, -5.0, getTheta(163) + 2.0},
                                                   {-2.0*getTheta(163) + 1.0, -getTheta(163) + 17.0, -5.0, -getTheta(163) - 2.0},
                                                   {-getTheta(163) - 2.0, -getTheta(163) + 17.0, -5.0, -2.0*getTheta(163) + 1.0},
                                                   {getTheta(163) + 2.0, -getTheta(163) + 17.0, -5.0, 2.0*getTheta(163) - 1.0},
                                                   {2.0*getTheta(163), 33.0, -5.0, 2.0*getTheta(163) - 2.0},
                                                   {-2.0*getTheta(163), 33.0, -5.0, -2.0*getTheta(163) + 2.0},
                                                   {-2.0*getTheta(163) + 2.0, 33.0, -5.0, -2.0*getTheta(163)},
                                                   {2.0*getTheta(163) - 2.0, 33.0, -5.0, 2.0*getTheta(163)},
                                                   {-getTheta(163) - 1.0, getTheta(163) - 16.0, -5.0, 2.0*getTheta(163) + 1.0},
                                                   {getTheta(163) + 1.0, getTheta(163) - 16.0, -5.0, -2.0*getTheta(163) - 1.0},
                                                   {-2.0*getTheta(163) - 1.0, getTheta(163) - 16.0, -5.0, getTheta(163) + 1.0},
                                                   {2.0*getTheta(163) + 1.0, getTheta(163) - 16.0, -5.0, -getTheta(163) - 1.0},
                                                   {-getTheta(163), -8.0, -5.0, getTheta(163) - 1.0},
                                                   {getTheta(163), -8.0, -5.0, -getTheta(163) + 1.0},
                                                   {-getTheta(163) + 1.0, -8.0, -5.0, getTheta(163)},
                                                   {getTheta(163) - 1.0, -8.0, -5.0, -getTheta(163)},
                                                   {getTheta(163) - 3.0, 8.0, -6.0, getTheta(163) + 2.0},
                                                   {-getTheta(163) + 3.0, 8.0, -6.0, -getTheta(163) - 2.0},
                                                   {-getTheta(163) - 2.0, 8.0, -6.0, -getTheta(163) + 3.0},
                                                   {getTheta(163) + 2.0, 8.0, -6.0, getTheta(163) - 3.0},
                                                   {-getTheta(163) - 1.0, -7.0, -6.0, getTheta(163) - 2.0},
                                                   {getTheta(163) + 1.0, -7.0, -6.0, -getTheta(163) + 2.0},
                                                   {-getTheta(163) + 2.0, -7.0, -6.0, getTheta(163) + 1.0},
                                                   {getTheta(163) - 2.0, -7.0, -6.0, -getTheta(163) - 1.0},
                                                   {getTheta(163), 7.0, -6.0, getTheta(163) - 1.0},
                                                   {-getTheta(163), 7.0, -6.0, -getTheta(163) + 1.0},
                                                   {-getTheta(163) + 1.0, 7.0, -6.0, -getTheta(163)},
                                                   {getTheta(163) - 1.0, 7.0, -6.0, getTheta(163)},
                                                   {-12.0, -5.0*getTheta(163), -getTheta(163) + 1.0, 17.0},
                                                   {12.0, -5.0*getTheta(163), -getTheta(163) + 1.0, -17.0},
                                                   {-17.0, -5.0*getTheta(163), -getTheta(163) + 1.0, 12.0},
                                                   {17.0, -5.0*getTheta(163), -getTheta(163) + 1.0, -12.0},
                                                   {-12.0, -5.0*getTheta(163) + 5.0, -getTheta(163), 17.0},
                                                   {12.0, -5.0*getTheta(163) + 5.0, -getTheta(163), -17.0},
                                                   {-17.0, -5.0*getTheta(163) + 5.0, -getTheta(163), 12.0},
                                                   {17.0, -5.0*getTheta(163) + 5.0, -getTheta(163), -12.0},
                                                   {12.0, 5.0*getTheta(163) - 7.0, -getTheta(163) + 2.0, getTheta(163) + 16.0},
                                                   {-12.0, 5.0*getTheta(163) - 7.0, -getTheta(163) + 2.0, -getTheta(163) - 16.0},
                                                   {-getTheta(163) - 16.0, 5.0*getTheta(163) - 7.0, -getTheta(163) + 2.0, -12.0},
                                                   {getTheta(163) + 16.0, 5.0*getTheta(163) - 7.0, -getTheta(163) + 2.0, 12.0},
                                                   {-getTheta(163) + 17.0, 5.0*getTheta(163) + 2.0, -getTheta(163) - 1.0, 12.0},
                                                   {getTheta(163) - 17.0, 5.0*getTheta(163) + 2.0, -getTheta(163) - 1.0, -12.0},
                                                   {-12.0, 5.0*getTheta(163) + 2.0, -getTheta(163) - 1.0, getTheta(163) - 17.0},
                                                   {12.0, 5.0*getTheta(163) + 2.0, -getTheta(163) - 1.0, -getTheta(163) + 17.0},
                                                   {-getTheta(163) - 15.0, -4.0*getTheta(163) + 18.0, -getTheta(163) + 3.0, getTheta(163) + 10.0},
                                                   {getTheta(163) + 15.0, -4.0*getTheta(163) + 18.0, -getTheta(163) + 3.0, -getTheta(163) - 10.0},
                                                   {-getTheta(163) - 10.0, -4.0*getTheta(163) + 18.0, -getTheta(163) + 3.0, getTheta(163) + 15.0},
                                                   {getTheta(163) + 10.0, -4.0*getTheta(163) + 18.0, -getTheta(163) + 3.0, -getTheta(163) - 15.0},
                                                   {-getTheta(163) - 17.0, -3.0*getTheta(163) + 1.0, -getTheta(163) + 3.0, 7.0},
                                                   {getTheta(163) + 17.0, -3.0*getTheta(163) + 1.0, -getTheta(163) + 3.0, -7.0},
                                                   {-7.0, -3.0*getTheta(163) + 1.0, -getTheta(163) + 3.0, getTheta(163) + 17.0},
                                                   {7.0, -3.0*getTheta(163) + 1.0, -getTheta(163) + 3.0, -getTheta(163) - 17.0},
                                                   {getTheta(163) - 18.0, -3.0*getTheta(163) + 2.0, -getTheta(163) - 2.0, 7.0},
                                                   {-getTheta(163) + 18.0, -3.0*getTheta(163) + 2.0, -getTheta(163) - 2.0, -7.0},
                                                   {-7.0, -3.0*getTheta(163) + 2.0, -getTheta(163) - 2.0, -getTheta(163) + 18.0},
                                                   {7.0, -3.0*getTheta(163) + 2.0, -getTheta(163) - 2.0, getTheta(163) - 18.0},
                                                   {-getTheta(163) + 16.0, -4.0*getTheta(163) - 14.0, -getTheta(163) - 2.0, getTheta(163) - 11.0},
                                                   {getTheta(163) - 16.0, -4.0*getTheta(163) - 14.0, -getTheta(163) - 2.0, -getTheta(163) + 11.0},
                                                   {-getTheta(163) + 11.0, -4.0*getTheta(163) - 14.0, -getTheta(163) - 2.0, getTheta(163) - 16.0},
                                                   {getTheta(163) - 11.0, -4.0*getTheta(163) - 14.0, -getTheta(163) - 2.0, -getTheta(163) + 16.0},
                                                   {-2.0*getTheta(163) + 1.0, getTheta(163) - 12.0, -7.0, getTheta(163) + 3.0},
                                                   {2.0*getTheta(163) - 1.0, getTheta(163) - 12.0, -7.0, -getTheta(163) - 3.0},
                                                   {-getTheta(163) - 3.0, getTheta(163) - 12.0, -7.0, 2.0*getTheta(163) - 1.0},
                                                   {getTheta(163) + 3.0, getTheta(163) - 12.0, -7.0, -2.0*getTheta(163) + 1.0},
                                                   {2.0*getTheta(163) + 1.0, 24.0, -7.0, 2.0*getTheta(163) - 3.0},
                                                   {-2.0*getTheta(163) - 1.0, 24.0, -7.0, -2.0*getTheta(163) + 3.0},
                                                   {-2.0*getTheta(163) + 3.0, 24.0, -7.0, -2.0*getTheta(163) - 1.0},
                                                   {2.0*getTheta(163) - 3.0, 24.0, -7.0, 2.0*getTheta(163) + 1.0},
                                                   {getTheta(163) + 16.0, 2.0*getTheta(163) - 15.0, -getTheta(163) + 4.0, getTheta(163) + 4.0},
                                                   {-getTheta(163) - 16.0, 2.0*getTheta(163) - 15.0, -getTheta(163) + 4.0, -getTheta(163) - 4.0},
                                                   {-getTheta(163) - 4.0, 2.0*getTheta(163) - 15.0, -getTheta(163) + 4.0, -getTheta(163) - 16.0},
                                                   {getTheta(163) + 4.0, 2.0*getTheta(163) - 15.0, -getTheta(163) + 4.0, getTheta(163) + 16.0},
                                                   {getTheta(163) - 17.0, 2.0*getTheta(163) + 13.0, -getTheta(163) - 3.0, getTheta(163) - 5.0},
                                                   {-getTheta(163) + 17.0, 2.0*getTheta(163) + 13.0, -getTheta(163) - 3.0, -getTheta(163) + 5.0},
                                                   {-getTheta(163) + 5.0, 2.0*getTheta(163) + 13.0, -getTheta(163) - 3.0, -getTheta(163) + 17.0},
                                                   {getTheta(163) - 5.0, 2.0*getTheta(163) + 13.0, -getTheta(163) - 3.0, getTheta(163) - 17.0}};
const vector<complex<double>> Quaternion::d163Translators = {0.0,
                                                             1.0,
                                                             -1.0,
                                                             getTheta(163),
                                                             -getTheta(163),
                                                             getTheta(163) - 1.0,
                                                             -getTheta(163) + 1.0};

const map<int, vector<complex<double>>> Quaternion::alphas = {{19, d19Alphas},
                                                              {43, d43Alphas},
                                                              {67, d67Alphas},
                                                              {163, d163Alphas}};

const map<int, vector<SL2C>> Quaternion::MAlphas = {{19, d19MAlphas},
                                                    {43, d43MAlphas},
                                                    {67, d67MAlphas},
                                                    {163, d163MAlphas}};

const map<int, vector<complex<double>>> Quaternion::translators = {{19, d19Translators},
                                                                   {43, d43Translators},
                                                                   {67, d67Translators},
                                                                   {163, d163Translators}};

Quaternion::Quaternion(double x, double y, double z, double w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

Quaternion::Quaternion(const Quaternion &q) {
    this->x = q.x;
    this->y = q.y;
    this->z = q.z;
    this->w = q.w;
}

std::ostream &operator<<(std::ostream &strm, const Quaternion &q) {
    std::ostringstream out;
    out << std::setprecision(16);
    out << q.x << " + ";
    out << q.y << "*i + ";
    out << q.z << "*j + ";
    out << q.w << "*ij";
    return strm << std::move(out).str();
}

Quaternion operator*(const Quaternion &q1, const Quaternion &q2) {
    double c1 = q1.x*q2.x - q1.y*q2.y - q1.z*q2.z - q1.w*q2.w;
    double c2 = q1.x*q2.y + q1.y*q2.x + q1.z*q2.w - q1.w*q2.z;
    double c3 = q1.x*q2.z - q1.y*q2.w + q1.z*q2.x + q1.w*q2.y;
    double c4 = q1.x*q2.w + q1.y*q2.z - q1.z*q2.y + q1.w*q2.x;
    return {c1, c2, c3, c4};
}

Quaternion operator+(const Quaternion &lhs, const Quaternion& rhs) {
    Quaternion answer = Quaternion(lhs);
    answer += rhs;
    return answer;
}

Quaternion operator-(const Quaternion &lhs, const Quaternion &rhs) {
    Quaternion answer = Quaternion(lhs);
    answer -= rhs;
    return answer;
}

Quaternion operator/(const Quaternion &q1, const Quaternion &q2) {
    return q1*(q2.inverse());
}

Quaternion operator*(const double &c, const Quaternion& q) {
    return {c*q.x, c*q.y, c*q.z, c*q.w};
}

complex<double> Quaternion::getComplex() const {
    return {this->x, this->y};
}

double Quaternion::squareAbs() const {
    return pow(this->x,2) + pow(this->y,2) + pow(this->z,2) + pow(this->w,2);
}

double Quaternion::abs() const {
    return sqrt(squareAbs());
}

Quaternion Quaternion::inverse() const {
    return (1/squareAbs())*(this->conj());
}

Quaternion Quaternion::conj() const {
    return {this->x, -this->y, -this->z, -this->w};
}

Quaternion operator-(const Quaternion &q) {
    return {-q.x, -q.y, -q.z, -q.w};
}

void Quaternion::operator+=(const Quaternion& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
}

void Quaternion::operator-=(const Quaternion& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
}

void Quaternion::operator*=(const Quaternion& rhs) {
    *this = (*this) * rhs;
}

void Quaternion::reduceModT() {
    this->vectorReduce(complex<double> {1,0});
}

bool Quaternion::reduceModInversion(int d) {
    if (d < 19) {
        double temp = x*x + y*y + z*z;
        if (temp < 1) {
            x = -x;
            *this = (1/temp) * (*this);
            return true;
        } else {
            return false;
        }
    } else {
        //compute the image wrt each MAlpha
        //compare the j part of each image to the j part of doing nothing
        //the answer is the one with greatest j part

        auto dAlphas = alphas.at(d);
        auto dMAlphas = MAlphas.at(d);
        auto dTranslators = translators.at(d);

        Quaternion largestSoFar = Quaternion(*this);
        SL2C inverter = {1,0,0,1};
        complex<double> translator = 0.0;

        auto tempComplex = this->getComplex();
        double leastDenominator = 1.0;
        for (int i = 0; i < dAlphas.size(); i++) {
            auto alpha = dAlphas[i];
            auto mAlpha = dMAlphas[i];

            //find nearest element of OK to *this and translate by that
            auto minDistance = (double)+INFINITY;
            complex<double> toTranslate;

            //TODO: Should be possible to calculate this directly. Is it faster than the loop?
            for (auto n : dTranslators) {
                double dist = std::abs((tempComplex - n) - alpha);
                if (dist < minDistance) {
                    minDistance = dist;
                    toTranslate = n;
                }
            }

            //height of point under the map MAlpha is height/(|cz+d|^2 + |c|^2*height^2)
            //The point is only raised if the denominator is less than 1
            double denominator = norm(mAlpha.c * (tempComplex - toTranslate) + mAlpha.d) + norm(mAlpha.c * z);

            if (denominator < leastDenominator) {
                translator = toTranslate;
                inverter = mAlpha;
                leastDenominator = denominator;
            }
        }
        if (leastDenominator == 1) {
            return false;
        } else {
            *this  = inverter * (*this - Quaternion(translator));
            return true;
        }
    }
}

void Quaternion::reduceThetaGeneral(int d) {
    if (d == 3) {
        double minDistance = 100;

        auto theta = getTheta(d);
        complex<double> othertheta = {-theta.real(), theta.imag()};

        complex<double> temp1 = this->getComplex();

        complex<double> translated = {0,0};

        int horizbound = std::ceil(std::abs(temp1.real()));
        int verbound = std::ceil(std::abs(temp1.imag()/theta.imag()));
        for (int a = -horizbound; a <= horizbound; a++) {
            for (int b = -verbound; b <= verbound; b++) {
                complex<double> latticePoint = (double)a + (double)b*theta;
                double thisDistance = std::abs(temp1 - latticePoint);
                if (thisDistance < minDistance) {
                    minDistance = thisDistance;
                    translated = latticePoint;
                }
            }
        }
        x = temp1.real();
        y = temp1.imag();
    } else {
        auto theta = getTheta(d);
        auto complex = this->getComplex();
        //double n = -std::floor((y + theta.imag()/2)/(theta.imag()));
        double n = -std::floor(y/theta.imag() + 0.5);
        complex += n*theta;
        x = complex.real();
        y = complex.imag();
    }
}

void Quaternion::reduceUnits(int d) {
    if (d == 1) {
        if (x < 0) {
            Quaternion gaussianUnit = {0,1,0,0};
            *this = gaussianUnit * (*this) * gaussianUnit;
        }
    } else if (d == 3) {
        Quaternion quatTheta = Quaternion(getTheta(d).real(), getTheta(d).imag(), 0, 0);

        while (x < 0 || y < -1.0/sqrt(3)*x) {
            *this = quatTheta * (*this) * quatTheta;
        }
    }
}

double Quaternion::getJ() const {
    return this->z;
}

void Quaternion::reduce(int d) {
    reduceModInversion(d);
    do {
        reduceThetaGeneral(d);
        reduceModT();
        reduceUnits(d);
    } while (reduceModInversion(d));
}

bool Quaternion::isInInversionFundamentalDomain() const {
    assert(z > 0);
    return this->squareAbs() >= 1.0;
}

void Quaternion::setEqual(Quaternion &q) {
    this->x = q.x;
    this->y = q.y;
    this->z = q.z;
    this->w = q.w;
}

int Quaternion::mod(int a, int n) {
    int answer = a % n;
    if (answer < 0) {
        return answer + n;
    } else {
        return answer;
    }
}

complex<double> Quaternion::getTheta(int d) {
    complex<double> theta;
    if (mod(-d,4) == 1) {
        theta = complex<double> {1.0/2, sqrt(d)/2};
    } else if (mod(-d,4) == 0) {
        throw(std::invalid_argument("d should be squarefree"));
    } else {
        theta = complex<double> {0, sqrt(d)};
    }
    return theta;
}

/**
 * This reduces the Quaternion<T> via the action of addition/subtraction by v.
 * @tparam T
 * @param v
 */
void Quaternion::vectorReduce(const complex<double> v) {
    complex<double> complex = this->getComplex();
    double modifiedScalarProjection = complex.real()*v.real() + complex.imag()*v.imag();
    modifiedScalarProjection /= pow(v.real(),2) + pow(v.imag(),2);
    modifiedScalarProjection = floor(modifiedScalarProjection + 1.0/2);
    this->x -= modifiedScalarProjection*v.real();
    this->y -= modifiedScalarProjection*v.imag();
}

void Quaternion::operator/=(const Quaternion &rhs) {
    *this = (*this)*(rhs.inverse());
}

Quaternion operator*(const SL2C &gamma, const Quaternion &q) {
    complex<double> z = q.getComplex();

    double denominator = norm(gamma.c * z + gamma.d) + norm(gamma.c * q.z);

    complex<double> newComplex = (gamma.a * z + gamma.b) * conj((gamma.c * z + gamma.d)) + gamma.a * conj((gamma.c)) * (q.z) * (q.z);

    Quaternion answer = {newComplex.real()/denominator, newComplex.imag()/denominator, q.z/denominator, 0};

    return answer;
}

Quaternion::Quaternion(complex<double> z) {
    this->x = z.real();
    this->y = z.imag();
    this->z = 0.0;
    this->w = 0.0;
}

bool operator==(const Quaternion& q1, const Quaternion& q2) {
    if ((q1.x == q2.x) && (q1.y == q2.y) && (q1.z == q2.z) && (q1.w == q2.w)) {
        return true;
    } else {
        return false;
    }
}

complex<double> Quaternion::vectorReduce(const complex<double> &u, const complex<double> &v) {
    double modifiedScalarProjection = u.real()*v.real() + u.imag()*v.imag();
    modifiedScalarProjection /= pow(v.real(),2) + pow(v.imag(),2);
    modifiedScalarProjection = floor(modifiedScalarProjection + 1.0/2);
    return u - modifiedScalarProjection*v;
}

Quaternion::Quaternion() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}




//
// Created by Eric Moss on 8/17/23.
//
