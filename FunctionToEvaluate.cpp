//
// Created by Eric Moss on 8/26/23.
//

#include "FunctionToEvaluate.h"
#include <complex>

FunctionToEvaluate::FunctionToEvaluate() {

}

std::vector<std::vector<double>> FunctionToEvaluate::cubicRoots(const std::vector<std::vector<Point>> &domainPoints) {
    std::vector<std::complex<double>> roots;
    int bound = 10;
    for (int a = -bound; a <= bound; a++) {
        for (int b = -bound; b <= bound; b++) {
            for (int c = -bound; c <= bound; c++) {
                double disc = (double)(b*b - 4*a*c);
                if (disc < 0) {
                    //add the root -b/2 + I*sqrt(-disc)/2
                    double real = -b / (2.0 * a);
                    double imag = sqrt(-disc) / (2.0 * a);

                    roots.emplace_back(real, imag);
                }
            }
        }
    }


    std::vector<std::vector<double>> output;
    for (int i = 0; i < domainPoints.size(); i++){
        std::vector<double> row;
        for (int j = 0; j < domainPoints[0].size(); j++) {
            row.push_back(0.0);
        }
        output.push_back(row);
    }

    double xSpacing = domainPoints[0][1].x - domainPoints[0][0].x;
    double ySpacing = domainPoints[0][0].y - domainPoints[1][0].y;
    double xMin = domainPoints[0][0].x - xSpacing/2;
    double yMax = domainPoints[0][0].y + ySpacing/2;
    double xMax = domainPoints[0][domainPoints[0].size() - 1].x + xSpacing/2;
    double yMin = domainPoints[domainPoints.size() - 1][0].y - ySpacing/2;

    for (auto root : roots) {
        if (root.real() > xMax || root.real() < xMin) {
            continue;
        }
        if (root.imag() > yMax || root.imag() < yMin) {
            continue;
        }
        int i  = floor((yMax + ySpacing/2 - root.imag())/ySpacing);
        int j = floor((root.real() - xMin - xSpacing/2)/xSpacing);
        output[i][j]++;
    }

    return output;
}

std::vector<std::vector<double>> FunctionToEvaluate::quadraticRoots(const std::vector<std::vector<Point>> &domainPoints) {
    double xSpacing = domainPoints[0][1].x - domainPoints[0][0].x;
    double ySpacing = domainPoints[0][0].y - domainPoints[1][0].y;
    double xMin = domainPoints[0][0].x - xSpacing/2;
    double yMax = domainPoints[0][0].y + ySpacing/2;
    double xMax = domainPoints[0][domainPoints[0].size() - 1].x + xSpacing/2;
    double yMin = domainPoints[domainPoints.size() - 1][0].y - ySpacing/2;

    std::vector<std::vector<double>> output;
    for (int i = 0; i < domainPoints.size(); i++){
        std::vector<double> row;
        for (int j = 0; j < domainPoints[0].size(); j++) {
            row.push_back(0.0);
        }
        output.push_back(row);
    }

    int bound = 10;
    for (int a = -bound; a <= bound; a++) {
        for (int b = -bound; b <= bound; b++) {
            for (int c = -bound; c <= bound; c++) {
                double disc = (double)(b*b - 4*a*c);
                if (disc < 0) {
                    //add the root -b/2 + I*sqrt(-disc)/2
                    double real = -b / (2.0 * a);
                    double imag = sqrt(-disc) / (2.0 * a);
                    std::complex<double> root = {real, imag};
                    if (root.real() > xMax || root.real() < xMin) {
                        continue;
                    }
                    if (root.imag() > yMax || root.imag() < yMin) {
                        continue;
                    }
                    int i  = floor((yMax - root.imag())/ySpacing);
                    int j = floor((root.real() - xMin)/xSpacing);
                    output[i][j]++;
                }
            }
        }
    }

    return output;
}


double FunctionToEvaluate::f(Point point) {
    double x = point.x;
    double y = point.y;
    std::complex<double> z = {x,y};

    std::complex<double> ans = {0,0};
    for (int i = 1; i < 10; i++) {
        ans += exp(2*3.14159*std::complex<double>{0,1}*(double)i*z);
    }

    return arg(ans);
}

std::vector<std::vector<double>> FunctionToEvaluate::gridEvaluate(const std::vector<std::vector<Point>> &domainPoints) {
    std::vector<std::vector<double>> answer;

    for (const auto& itr1 : domainPoints) {
        std::vector<double> row;
        row.reserve(itr1.size());
        for (auto itr2 : itr1) {
            row.push_back(f(itr2));
        }
        answer.push_back(row);
    }
    return answer;
}

std::vector<std::vector<double>>
FunctionToEvaluate::indexPoints(const std::vector<std::vector<Point>> &domainPoints, const std::vector<Index> &indices) {
    std::vector<std::vector<double>> answer;

    std::vector<std::complex<double>> latticePoints;
    for (auto itr : indices) {
        latticePoints.push_back(itr.getComplex());
    }

    for (auto itr1 : domainPoints) {
        std::vector<double> row;
        row.resize(itr1.size());
        answer.push_back(row);
    }

    double DETECTION_RADIUS = 0.25;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < domainPoints.size(); i++) {
        for (int j = 0; j < domainPoints[0].size(); j++) {
            std::complex<double> pixel = {domainPoints[i][j].x, domainPoints[i][j].y};
            bool closeToLatticePoint = false;
            for (auto latticePoint : latticePoints) {
                if (abs(latticePoint - pixel) < DETECTION_RADIUS) {
                    closeToLatticePoint = true;
                    break;
                }
            }
            if (closeToLatticePoint) {
                //black
                answer[i][j] = 0.0;
            } else {
                //white
                answer[i][j] = 1.0;
            }
        }
    }

    return answer;
}

std::vector<std::vector<double>> FunctionToEvaluate::latticePoints(const std::vector<std::vector<Point>> &domainPoints,
                                                                   const std::vector<std::complex<double>> &latticePoints) {
    return std::vector<std::vector<double>>();
}


