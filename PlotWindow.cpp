//
// Created by Eric Moss on 8/26/23.
//

#include "PlotWindow.h"

PlotWindow::PlotWindow(double aspectRatio, double windowHeight, Point bottomCenter, int verticalPixels) {
    yMin = bottomCenter.y;
    yMax = yMin + windowHeight;

    double verticalLength = yMax - yMin;
    double horizontalLength = aspectRatio * verticalLength;
    xMin = bottomCenter.x - horizontalLength / 2;
    xMax = bottomCenter.x + horizontalLength / 2;

    this->verticalPixels = verticalPixels;
    horizontalPixels = round(aspectRatio * verticalPixels);

    xSpacing = horizontalLength / horizontalPixels;
    ySpacing = verticalLength / verticalPixels;

    points.clear();
    for (int i = 0; i < verticalPixels; i++) {
        std::vector<Point> row;
        for (int j = 0; j < horizontalPixels; j++) {
            Point point;
            point.x = xMin + xSpacing/2 + j*xSpacing;
            point.y = yMax - ySpacing/2 - i*ySpacing;
            row.push_back(point);
        }
        points.push_back(row);
    }
}

PlotWindow::PlotWindow(double xMin, double xMax, double yMin, double yMax, int horizontalPixels) {
    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->horizontalPixels = horizontalPixels;
    double aspectRatio = (xMax - xMin)/(yMax - yMin);
    verticalPixels = round(horizontalPixels/aspectRatio);
    xSpacing = (xMax - xMin)/horizontalPixels;
    ySpacing = (yMax - yMin)/verticalPixels;

    points.clear();
    for (int i = 0; i < verticalPixels; i++) {
        std::vector<Point> row;
        for (int j = 0; j < horizontalPixels; j++) {
            Point point;
            point.x = xMin + xSpacing/2 + j*xSpacing;
            point.y = yMax - ySpacing/2 - i*ySpacing;
            row.push_back(point);
        }
        points.push_back(row);
    }
}

std::vector<std::vector<Point>> &PlotWindow::getPoints() {
    return points;
}

std::vector<int> PlotWindow::getOrigin() {
    int originRow = (0 - yMin)/ySpacing;
    int originColumn = (0 - xMin)/xSpacing;
    return {originRow, originColumn};
}



