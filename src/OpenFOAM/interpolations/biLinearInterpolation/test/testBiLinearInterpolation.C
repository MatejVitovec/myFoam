#include "biLinearInterpolation.H"

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>


namespace
{

Foam::dictionary makeSettings
(
    const Foam::word& transformX,
    const Foam::word& transformY
)
{
    Foam::labelList gridSizeX(2);
    gridSizeX[0] = 8;
    gridSizeX[1] = 7;

    Foam::labelList gridSizeY(2);
    gridSizeY[0] = 6;
    gridSizeY[1] = 9;

    Foam::scalarField boundaryX(3);
    boundaryX[0] = 0.1;
    boundaryX[1] = 1.0;
    boundaryX[2] = 1.5;

    Foam::scalarField boundaryY(3);
    boundaryY[0] = 0.1;
    boundaryY[1] = 1.0;
    boundaryY[2] = 1.5;

    Foam::dictionary settings("interpolationSettings");
    settings.add("gridSizeX", gridSizeX);
    settings.add("gridSizeY", gridSizeY);
    settings.add("boundaryX", boundaryX);
    settings.add("boundaryY", boundaryY);
    settings.add("transformationX", transformX);
    settings.add("transformationY", transformY);
    return settings;
}

}


int main()
{
    const struct
    {
        const char* name;
        const char* transformX;
        const char* transformY;
    } cases[] =
    {
        {"none", "none", "none"},
        {"log", "log", "log"},
        {"log10", "log10", "log10"},
        {"loginv", "loginv", "loginv"},
        {"mixed", "log", "loginv"}
    };

    const std::function<Foam::scalar(Foam::scalar, Foam::scalar)> function
    (
        [](const Foam::scalar x, const Foam::scalar y)
        {
            return 1.0 + 2.0*x - 3.0*y + 0.5*x*y;
        }
    );
    const std::function<Foam::scalar(Foam::scalar, Foam::scalar)>
        plotFunction
    (
        [](const Foam::scalar x, const Foam::scalar y)
        {
            return Foam::sin(x*y);
        }
    );

    const Foam::scalar tolerance = 1e-10;
    const Foam::scalar roundTripTolerance = 1e-4;
    bool passed = true;

    for (const auto& testCase : cases)
    {
        const Foam::dictionary settings
        (
            makeSettings(testCase.transformX, testCase.transformY)
        );
        const Foam::biLinearInterpolation calculated(settings, function);
        const Foam::biLinearInterpolation plotted(settings, plotFunction);
        const Foam::fileName coefficientFile
        (
            std::string("test_coefficients_") + testCase.name + ".dict"
        );
        calculated.write(coefficientFile);
        const Foam::biLinearInterpolation loaded(coefficientFile);

        const std::string dataFileName
        (
            std::string("test_") + testCase.name + ".dat"
        );
        std::ofstream data(dataFileName.c_str());
        if (!data.good())
        {
            std::cerr << "Cannot open output file " << dataFileName << '\n';
            return 1;
        }
        data << "# x y exact interpolated absError\n";
        data << std::setprecision(16);

        Foam::scalar maxError = 0;
        Foam::scalar maxRoundTripError = 0;
        for (Foam::label j = 0; j <= 100; ++j)
        {
            const Foam::scalar y = 0.1 + 1.4*j/100;
            for (Foam::label i = 0; i <= 100; ++i)
            {
                const Foam::scalar x = 0.1 + 1.4*i/100;
                const Foam::scalar exactPlotValue = plotFunction(x, y);
                const Foam::scalar interpolatedPlotValue = plotted(x, y);
                data << x << ' ' << y << ' ' << exactPlotValue << ' '
                     << interpolatedPlotValue << ' '
                     << Foam::mag(interpolatedPlotValue - exactPlotValue)
                     << '\n';

                maxError = Foam::max
                (
                    maxError,
                    Foam::mag(calculated(x, y) - function(x, y))
                );
                maxRoundTripError = Foam::max
                (
                    maxRoundTripError,
                    Foam::mag(loaded(x, y) - calculated(x, y))
                );
            }
            data << '\n';
        }
        if (!data.good())
        {
            std::cerr << "Error writing output file " << dataFileName << '\n';
            return 1;
        }

        std::cout << testCase.name
                  << ": max interpolation error = " << maxError
                  << ", max save/load difference = " << maxRoundTripError
                  << '\n';
        if (maxError > tolerance || maxRoundTripError > roundTripTolerance)
        {
            passed = false;
        }
    }

    if (!passed)
    {
        std::cerr << "Bilinear interpolation checks failed\n";
        return 1;
    }

    std::cout << "All bilinear interpolation checks passed\n";
    return 0;
}
