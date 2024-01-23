//
// Created by Eric Moss on 8/26/23.
//

#include "Plotter.h"
#include "tinycolormap.hpp"
#include <opencv2/opencv.hpp>
#include <iostream>

Plotter::Plotter(std::vector<std::vector<double>> &values) {
    this->values = values;

    verticalPixels = values.size();
    horizontalPixels = values[0].size();
    numPixels = verticalPixels * horizontalPixels;

    double max, min;

    for (const auto& row : values) {
        for (auto value : row) {
            if (value > max) {
                max = value;
            }
            if (value < min) {
                min = value;
            }
        }
    }
    dynamicRange = max - min;
    histogramBinSize = dynamicRange/1024;
    minValue = min;

    for (const auto& row : values) {
        for (auto value : row) {
            int bin = doubleToBinNumber(value);
            histogram[bin]++;
        }
    }

    int soFar = 0;
    for (auto itr : histogram) {
        integratedHistogram[itr.first] = (double)soFar/numPixels;
        soFar += itr.second;
    }
}


void
Plotter::makeImage(const std::string &filename, const std::string& color) {
    if (color != "viridis" && color != "greyscale") {
        std::cout << "color string is invalid\n";
        return;
    }

    if (color == "greyscale") {
        colorType = greyscale;
    } else if (color == "viridis") {
        colorType = viridis;
    }

    cv::Mat image = cv::Mat(verticalPixels, horizontalPixels, CV_8UC3, cv::Scalar(0,0,0));


    for (int i = 0; i < verticalPixels; ++i) {
        for (int j = 0; j < horizontalPixels; ++j) {
            cv::Vec3b pixel = doubleToPixel(values[i][j]);
            image.at<cv::Vec3b>(i,j) = pixel;
        }
    }

    if (willDrawAxes) {
        cv::Vec3b black = {0,0,0};
        for (int i = 0; i < verticalPixels; i++) {
            image.at<cv::Vec3b>(i,originColumn) = black;
        }
        for (int i = 0; i < horizontalPixels; i++) {
            image.at<cv::Vec3b>(originRow,i) = black;
        }
    }

    cv::imwrite(filename, image);
    std::cout << "Image saved." << std::endl;
}

cv::Vec3b Plotter::doubleToPixel(const double &x) {
    double scale = integratedHistogram[doubleToBinNumber(x)];

    if (histogramSmoothing == false) {
        scale = (x-minValue)/dynamicRange;
    }
    if (colorType == viridis) {
        tinycolormap::Color color = tinycolormap::GetColor(scale, tinycolormap::ColormapType::Viridis);
        return cv::Vec3b(255*color.b(), 255*color.g(), 255*color.r());
    } else {
        return cv::Vec3b(255*scale, 255*scale, 255*scale);
    }
}

int Plotter::doubleToBinNumber(const double &x) {
    double bin = (x - minValue)/histogramBinSize;
    return floor(bin);
}

void Plotter::drawAxes(int originRow, int originColumn) {
    this->originRow = originRow;
    this->originColumn = originColumn;
    willDrawAxes = true;
}
