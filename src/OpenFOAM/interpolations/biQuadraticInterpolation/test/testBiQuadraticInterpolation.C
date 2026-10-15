#include "biQuadraticInterpolation.H"

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>


namespace
{

struct TransformCase
{
    const char* name;
    const char* transformX;
    const char* transformY;
};


Foam::dictionary makeSettings
(
    const Foam::word& transformX,
    const Foam::word& transformY
)
{
    Foam::labelList gridSizeX(2);
    gridSizeX[0] = 14;
    gridSizeX[1] = 14;

    Foam::labelList gridSizeY(2);
    gridSizeY[0] = 13;
    gridSizeY[1] = 13;

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


bool runCase
(
    const TransformCase& testCase,
    const std::function<Foam::scalar(Foam::scalar, Foam::scalar)>& function,
    Foam::scalar& maxInterpolationError,
    Foam::scalar& rmsInterpolationError,
    Foam::scalar& maxRoundTripError
)
{
    const Foam::word transformX(testCase.transformX);
    const Foam::word transformY(testCase.transformY);
    const Foam::dictionary settings(makeSettings(transformX, transformY));

    const Foam::biQuadraticInterpolation calculated(settings, function);

    const std::string coefficientFile
    (
        std::string("test_coefficients_") + testCase.name + ".dict"
    );
    calculated.write(Foam::fileName(coefficientFile));

    const Foam::biQuadraticInterpolation loaded
    {
        Foam::fileName(coefficientFile)
    };

    const std::string dataFile
    (
        std::string("test_") + testCase.name + ".dat"
    );
    std::ofstream data(dataFile.c_str());
    if (!data.good())
    {
        std::cerr << "Cannot open output file " << dataFile << '\n';
        return false;
    }

    data << "# x y exact computed loaded absError roundTripError\n";
    data << std::setprecision(16);

    const Foam::label sampleCountX = 60;
    const Foam::label sampleCountY = 60;
    const Foam::scalar minX = 0.1;
    const Foam::scalar maxX = 1.5;
    const Foam::scalar minY = 0.1;
    const Foam::scalar maxY = 1.5;

    Foam::scalar errorSquaredSum = 0;
    Foam::scalar maxError = 0;
    Foam::scalar maxReloadDifference = 0;
    Foam::label sampleCount = 0;

    for (Foam::label j = 0; j <= sampleCountY; j++)
    {
        const Foam::scalar y
        (
            minY + (maxY - minY)*j/sampleCountY
        );

        for (Foam::label i = 0; i <= sampleCountX; i++)
        {
            const Foam::scalar x
            (
                minX + (maxX - minX)*i/sampleCountX
            );
            const Foam::scalar exact = function(x, y);
            const Foam::scalar interpolated = calculated(x, y);
            const Foam::scalar reloaded = loaded(x, y);
            const Foam::scalar error = Foam::mag(interpolated - exact);
            const Foam::scalar reloadDifference
            (
                Foam::mag(reloaded - interpolated)
            );

            data << x << ' ' << y << ' ' << exact << ' '
                 << interpolated << ' ' << reloaded << ' '
                 << error << ' ' << reloadDifference << '\n';

            errorSquaredSum += error*error;
            maxError = Foam::max(maxError, error);
            maxReloadDifference = Foam::max
            (
                maxReloadDifference,
                reloadDifference
            );
            sampleCount++;
        }
        data << '\n';
    }

    maxInterpolationError = maxError;
    rmsInterpolationError = Foam::sqrt(errorSquaredSum/sampleCount);
    maxRoundTripError = maxReloadDifference;

    std::cout << testCase.name
              << ": max interpolation error = " << maxInterpolationError
              << ", RMS interpolation error = " << rmsInterpolationError
              << ", max save/load difference = " << maxRoundTripError
              << '\n';

    return data.good();
}

}


int main()
{
    const TransformCase cases[] =
    {
        {"none",   "none",   "none"},
        {"log",    "log",    "log"},
        {"log10",  "log10",  "log10"},
        {"loginv", "loginv", "loginv"},
        {"mixed",  "log",    "loginv"}
    };

    const std::function<Foam::scalar(Foam::scalar, Foam::scalar)> function
    (
        [](const Foam::scalar x, const Foam::scalar y)
        {
            return 1/(x*x + y*y);
        }
    );

    const Foam::scalar roundTripTolerance = 1e-4;
    bool passed = true;

    for (Foam::label i = 0; i < 5; i++)
    {
        Foam::scalar maxInterpolationError = 0;
        Foam::scalar rmsInterpolationError = 0;
        Foam::scalar maxRoundTripError = 0;

        if
        (
            !runCase
            (
                cases[i],
                function,
                maxInterpolationError,
                rmsInterpolationError,
                maxRoundTripError
            )
        )
        {
            passed = false;
            continue;
        }

        if (maxRoundTripError > roundTripTolerance)
        {
            std::cerr << cases[i].name
                      << ": save/load difference exceeds tolerance "
                      << roundTripTolerance << '\n';
            passed = false;
        }
    }

    if (!passed)
    {
        return 1;
    }

    std::cout << "All save/load checks passed. Gnuplot data files were written.\n";
    return 0;
}
