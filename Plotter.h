//
// Created by Eric Moss on 8/26/23.
//

#ifndef IMAGEGENERATOR_PLOTTER_H
#define IMAGEGENERATOR_PLOTTER_H
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include <map>

enum color {greyscale, viridis};

class Plotter {
public:
    Plotter(std::vector<std::vector<double>> &values);
    void makeImage(const std::string &filename, const std::string& color);
    bool histogramSmoothing = false;
    void drawAxes(int originRow, int originColumn);
private:
    std::vector<std::vector<double>> values;
    cv::Vec3b doubleToPixel(const double &x);

    color colorType;
    int verticalPixels;
    int horizontalPixels;
    int numPixels;
    double dynamicRange;
    double histogramBinSize;
    double minValue;
    bool willDrawAxes = false;
    int originRow;
    int originColumn;
    int doubleToBinNumber(const double &x);
    std::map<int, int> histogram;
    std::map<int, double> integratedHistogram;
};


#endif //IMAGEGENERATOR_PLOTTER_H
