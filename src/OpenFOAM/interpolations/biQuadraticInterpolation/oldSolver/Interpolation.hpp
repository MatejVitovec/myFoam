#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <math.h>
#include <vector>

class Interpolation
{
public:
    enum Transformation{NONE, LOG, LOG10, LOGINV};
    Interpolation() : gridSizeX(std::vector<int>()), gridSizeY(std::vector<int>()),
    boundaryX(std::vector<double>()), boundaryY(std::vector<double>()),
    transformationX(NONE), transformationY(NONE) {}

    Interpolation(std::vector<int> gridSizeX_, std::vector<int> gridSizeY_,
                  std::vector<double> boundaryX_, std::vector<double> boundaryY_,
                  Transformation transformationX_, Transformation transformationY_) :
                  gridSizeX(gridSizeX_), gridSizeY(gridSizeY_),
                  boundaryX(boundaryX_), boundaryY(boundaryY_),
                  transformationX(transformationX_), transformationY(transformationY_) {}

                  virtual ~Interpolation() {}

                  virtual double calc(double xx, double yy) const = 0;
                  virtual double calcFastFind(double xx, double yy) const = 0;

                  virtual double calcFastFindDiffX(double xx, double yy) const = 0;
                  virtual double calcFastFindDiffY(double xx, double yy) const = 0;

                  virtual double calcInverseX(double zz, double yy, double guessXX) const = 0;
                  virtual double calcInverseY(double xx, double zz, double guessYY) const = 0;

protected:
    std::vector<int> gridSizeX;
    std::vector<int> gridSizeY;
    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    Transformation transformationX;
    Transformation transformationY;

    std::vector<double> createInterpolationAxis(std::vector<int> gridSize, std::vector<double> boundary, Transformation transformation) const;

    double transform(double x, Transformation transformation) const;
    double backTransform(double x, Transformation transformation) const;


};

#endif // INTERPOLATION_HPP
