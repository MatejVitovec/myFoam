#include <algorithm>
#include "Interpolation.hpp"


std::vector<double> Interpolation::createInterpolationAxis(std::vector<int> gridSize, std::vector<double> boundary, Transformation transformation) const
{
    for (int i = 0; i < boundary.size(); i++)
    {
        boundary[i] = transform(boundary[i], transformation);
    }

    if(boundary[0] > boundary[1])
    {
        std::reverse(boundary.begin(), boundary.end());
    }

    double size = 0.0;
    for (int i = 0; i < gridSize.size(); i++)
    {
        size += gridSize[i];
    }

    std::vector<double> out(size+1);

    std::vector<double> dx = std::vector<double>(gridSize.size());

    for (int i = 0; i < dx.size(); i++)
    {
        dx[i] = (boundary[i+1] - boundary[i])/gridSize[i];
    }

    int idx = 0;
    out[idx] = backTransform(boundary[0], transformation);
    idx++;

    for (int ii = 0; ii < gridSize.size(); ii++)
    {
        for (int i = 0; i < gridSize[ii]; i++)
        {
            out[idx] = backTransform(boundary[ii] + dx[ii]*(i+1), transformation);
            idx++;
        }
    }

    if(out[0] > out[1])
    {
        std::reverse(out.begin(), out.end());
    }

    return out;
}

double Interpolation::transform(double x, Transformation transformation) const
{
    switch(transformation)
    {
        case LOG:
            return log(x);
        case LOG10:
            return log10(x);
        case LOGINV:
            return log(1/x);
    }

    return x;
}

double Interpolation::backTransform(double x, Transformation transformation) const
{
    switch(transformation)
    {
        case LOG:
            return exp(x);
        case LOG10:
            return pow(x, 10.0);
        case LOGINV:
            return 1.0/exp(x);
    }

    return x;
}
