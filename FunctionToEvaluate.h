//
// Created by Eric Moss on 8/26/23.
//

#ifndef IMAGEGENERATOR_FUNCTIONTOEVALUATE_H
#define IMAGEGENERATOR_FUNCTIONTOEVALUATE_H
#include "PlotWindow.h"
#include "Index.h"

class FunctionToEvaluate {
public:
    FunctionToEvaluate();
    std::vector<std::vector<double>> indexPoints(const std::vector<std::vector<Point>> &domainPoints, const std::vector<Index> &indices);
    std::vector<std::vector<double>> latticePoints(const std::vector<std::vector<Point>> &domainPoints, const std::vector<std::complex<double>> &latticePoints);
    std::vector<std::vector<double>> cubicRoots(const std::vector<std::vector<Point>> &domainPoints);
    std::vector<std::vector<double>> quadraticRoots(const std::vector<std::vector<Point>> &domainPoints);
    std::vector<std::vector<double>> gridEvaluate(const std::vector<std::vector<Point>> &domainPoints);
    double f(Point point);
};


#endif //IMAGEGENERATOR_FUNCTIONTOEVALUATE_H
