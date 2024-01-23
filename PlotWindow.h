//
// Created by Eric Moss on 8/26/23.
//

#ifndef IMAGEGENERATOR_PLOTWINDOW_H
#define IMAGEGENERATOR_PLOTWINDOW_H
#include <vector>

struct Point {double x; double y;};

/***
 * PlotWindow is a rectangular region in the x-y plane, along with a number of pixels in each the horizontal
 * and vertical directions. This can be specified in a number of ways. Make another constructor to suit
 * your preference.
 *
 * A PlotWindow provides points where we will evaluate whatever function we want. It provides the center point
 * of each pixel.
 *
 * Final usage for plotting the function f might look like
 *
 * PlotWindow window = ...
 * for (auto point : window.points) {
 *      auto evaluation = f(point);
 *      evaluations.push_back(evaluation);
 * }
 * auto evaluation = f(window
 */
class PlotWindow {
public:
    PlotWindow(double aspectRatio, double windowHeight, Point bottomCenter, int verticalPixels);
    PlotWindow(double xMin, double xMax, double yMin, double yMax, int horizontalPixels);

    std::vector<std::vector<Point>> &getPoints();
    std::vector<int> getOrigin();
private:
    double xMin;
    double xMax;
    double yMin;
    double yMax;

    int verticalPixels;
    int horizontalPixels;

    double xSpacing;
    double ySpacing;
    std::vector<std::vector<Point>> points;
};


#endif //IMAGEGENERATOR_PLOTWINDOW_H
