/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "biLinearInterpolation.H"
#include "IFstream.H"
#include "OFstream.H"
#include <cmath>


Foam::biLinearInterpolation::biLinearInterpolation(const fileName& fName)
:
    fileName_(fName),
    n_(0),
    m_(0),
    transformationX_(NONE),
    transformationY_(NONE)
{
    readTable();
}


Foam::biLinearInterpolation::biLinearInterpolation
(
    const dictionary& dict,
    const std::function<scalar(scalar, scalar)>& f
)
:
    fileName_(),
    n_(0),
    m_(0),
    transformationX_(NONE),
    transformationY_(NONE)
{
    readSettings(dict);
    calcCoeffs(f);
}


void Foam::biLinearInterpolation::readTable()
{
    IFstream is(fileName_);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open coefficient file " << fileName_ << nl
            << exit(FatalError);
    }

    const dictionary dict(is);
    readSettings(dict);
    coefficients_ = dict.get<List<FixedList<scalar, 4>>>("coefficients");

    if (coefficients_.size() != n_*m_)
    {
        FatalErrorInFunction
            << "Coefficient file " << fileName_ << " contains "
            << coefficients_.size() << " cell entries, expected "
            << n_*m_ << nl
            << exit(FatalError);
    }
}


Foam::biLinearInterpolation::Transformation
Foam::biLinearInterpolation::readTransformation(const word& name)
{
    if (name == "none") return NONE;
    if (name == "log") return LOG;
    if (name == "log10") return LOG10;
    if (name == "loginv") return LOGINV;

    FatalErrorInFunction
        << "Unknown coordinate transformation '" << name << "'. Valid values: "
        << "none, log, log10, loginv" << nl
        << exit(FatalError);
    return NONE;
}


Foam::word Foam::biLinearInterpolation::transformationName
(
    const Transformation transformation
)
{
    switch (transformation)
    {
        case NONE: return "none";
        case LOG: return "log";
        case LOG10: return "log10";
        case LOGINV: return "loginv";
    }
    return "none";
}


void Foam::biLinearInterpolation::check() const
{
    if
    (
        gridSizeX_.empty()
     || gridSizeY_.empty()
     || boundaryX_.size() != gridSizeX_.size() + 1
     || boundaryY_.size() != gridSizeY_.size() + 1
    )
    {
        FatalErrorInFunction
            << "Grid sizes must be non-empty and each boundary list must "
            << "contain one more entry than its grid-size list" << nl
            << exit(FatalError);
    }

    forAll(gridSizeX_, i)
    {
        if (gridSizeX_[i] <= 0 || boundaryX_[i + 1] <= boundaryX_[i])
        {
            FatalErrorInFunction
                << "The x grid sizes must be positive and x boundaries "
                << "strictly increasing" << nl
                << exit(FatalError);
        }
    }

    forAll(gridSizeY_, j)
    {
        if (gridSizeY_[j] <= 0 || boundaryY_[j + 1] <= boundaryY_[j])
        {
            FatalErrorInFunction
                << "The y grid sizes must be positive and y boundaries "
                << "strictly increasing" << nl
                << exit(FatalError);
        }
    }

    if (transformationX_ != NONE && boundaryX_[0] <= 0)
    {
        FatalErrorInFunction
            << "The x boundaries must be positive for logarithmic "
            << "transformations" << nl
            << exit(FatalError);
    }

    if (transformationY_ != NONE && boundaryY_[0] <= 0)
    {
        FatalErrorInFunction
            << "The y boundaries must be positive for logarithmic "
            << "transformations" << nl
            << exit(FatalError);
    }
}


void Foam::biLinearInterpolation::write(Ostream& os) const
{
    os.writeEntry("gridSizeX", gridSizeX_);
    os.writeEntry("gridSizeY", gridSizeY_);
    os.writeEntry("boundaryX", boundaryX_);
    os.writeEntry("boundaryY", boundaryY_);
    os.writeEntry("transformationX", transformationName(transformationX_));
    os.writeEntry("transformationY", transformationName(transformationY_));
    os.writeEntry("coefficients", coefficients_);
}


void Foam::biLinearInterpolation::write(const fileName& fName) const
{
    OFstream os(fName);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open coefficient file " << fName << " for writing" << nl
            << exit(FatalError);
    }

    os << token::BEGIN_BLOCK << nl;
    write(os);
    os << token::END_BLOCK << nl;
}


void Foam::biLinearInterpolation::createGrid()
{
    n_ = sum(gridSizeX_);
    m_ = sum(gridSizeY_);

    x_.setSize(n_ + 1);
    dx_.setSize(n_);
    sectionStepX_.setSize(gridSizeX_.size());

    label xIndex = 0;
    forAll(gridSizeX_, section)
    {
        const scalar transformedStart
        (
            transform(boundaryX_[section], transformationX_)
        );
        const scalar transformedEnd
        (
            transform(boundaryX_[section + 1], transformationX_)
        );
        const scalar transformedStep
        (
            (transformedEnd - transformedStart)/gridSizeX_[section]
        );
        sectionStepX_[section] = mag(transformedStep);

        for
        (
            label localNode = 0;
            localNode <= gridSizeX_[section];
            ++localNode
        )
        {
            x_[xIndex + localNode] = backTransform
            (
                transformedStart + localNode*transformedStep,
                transformationX_
            );
        }
        xIndex += gridSizeX_[section];
    }

    y_.setSize(m_ + 1);
    dy_.setSize(m_);
    sectionStepY_.setSize(gridSizeY_.size());

    label yIndex = 0;
    forAll(gridSizeY_, section)
    {
        const scalar transformedStart
        (
            transform(boundaryY_[section], transformationY_)
        );
        const scalar transformedEnd
        (
            transform(boundaryY_[section + 1], transformationY_)
        );
        const scalar transformedStep
        (
            (transformedEnd - transformedStart)/gridSizeY_[section]
        );
        sectionStepY_[section] = mag(transformedStep);

        for
        (
            label localNode = 0;
            localNode <= gridSizeY_[section];
            ++localNode
        )
        {
            y_[yIndex + localNode] = backTransform
            (
                transformedStart + localNode*transformedStep,
                transformationY_
            );
        }
        yIndex += gridSizeY_[section];
    }

    for (label i = 0; i < n_; ++i)
    {
        dx_[i] = x_[i + 1] - x_[i];
    }
    for (label j = 0; j < m_; ++j)
    {
        dy_[j] = y_[j + 1] - y_[j];
    }

    coefficients_.setSize(n_*m_);
}


void Foam::biLinearInterpolation::readSettings(const dictionary& dict)
{
    gridSizeX_ = dict.get<labelList>("gridSizeX");
    gridSizeY_ = dict.get<labelList>("gridSizeY");
    boundaryX_ = dict.get<scalarField>("boundaryX");
    boundaryY_ = dict.get<scalarField>("boundaryY");
    transformationX_ = readTransformation
    (
        dict.get<word>("transformationX")
    );
    transformationY_ = readTransformation
    (
        dict.get<word>("transformationY")
    );

    check();
    createGrid();
}


void Foam::biLinearInterpolation::calcCoeffs
(
    const std::function<scalar(scalar, scalar)>& f
)
{
    for (label j = 0; j < m_; ++j)
    {
        for (label i = 0; i < n_; ++i)
        {
            const scalar f00 = f(x_[i], y_[j]);
            const scalar f10 = f(x_[i + 1], y_[j]);
            const scalar f01 = f(x_[i], y_[j + 1]);
            const scalar f11 = f(x_[i + 1], y_[j + 1]);

            coefficients_[j*n_ + i] =
            {
                f00,
                f10 - f00,
                f01 - f00,
                f11 - f10 - f01 + f00
            };
        }
    }
}


Foam::scalar Foam::biLinearInterpolation::transform
(
    const scalar x,
    const Transformation transformation
) const
{
    switch (transformation)
    {
        case LOG: return log(x);
        case LOG10: return log10(x);
        case LOGINV: return -log(x);
        case NONE: return x;
    }
    return x;
}


Foam::scalar Foam::biLinearInterpolation::backTransform
(
    const scalar x,
    const Transformation transformation
) const
{
    switch (transformation)
    {
        case LOG: return exp(x);
        case LOG10: return pow(10.0, x);
        case LOGINV: return exp(-x);
        case NONE: return x;
    }
    return x;
}


Foam::scalar Foam::biLinearInterpolation::calc
(
    const scalar xx,
    const scalar yy
) const
{
    const std::pair<int, int> position = fastFindPosition(xx, yy);
    const label i = position.first;
    const label j = position.second;
    const scalar tx = (xx - x_[i])/dx_[i];
    const scalar ty = (yy - y_[j])/dy_[j];
    const FixedList<scalar, 4>& coeff = coefficients_[j*n_ + i];

    return coeff[0] + tx*(coeff[1] + ty*coeff[3]) + ty*coeff[2];
}


std::pair<int, int> Foam::biLinearInterpolation::fastFindPosition
(
    const scalar xx,
    const scalar yy
) const
{
    label sectionX = 0;
    label sectionY = 0;
    label shiftIdxX = 0;
    label shiftIdxY = 0;

    for (; sectionX < gridSizeX_.size() - 1; ++sectionX)
    {
        if (xx < boundaryX_[sectionX + 1])
        {
            break;
        }
        shiftIdxX += gridSizeX_[sectionX];
    }

    for (; sectionY < gridSizeY_.size() - 1; ++sectionY)
    {
        if (yy < boundaryY_[sectionY + 1])
        {
            break;
        }
        shiftIdxY += gridSizeY_[sectionY];
    }

    label localCellX = 0;
    if (xx >= boundaryX_.last())
    {
        localCellX = gridSizeX_.last() - 1;
    }
    else if (xx > boundaryX_[0])
    {
        const scalar distance
        (
            mag
            (
                transform(xx, transformationX_)
              - transform(boundaryX_[sectionX], transformationX_)
            )
        );
        localCellX = min
        (
            label(floor(distance/sectionStepX_[sectionX])),
            gridSizeX_[sectionX] - 1
        );
    }

    label localCellY = 0;
    if (yy >= boundaryY_.last())
    {
        localCellY = gridSizeY_.last() - 1;
    }
    else if (yy > boundaryY_[0])
    {
        const scalar distance
        (
            mag
            (
                transform(yy, transformationY_)
              - transform(boundaryY_[sectionY], transformationY_)
            )
        );
        localCellY = min
        (
            label(floor(distance/sectionStepY_[sectionY])),
            gridSizeY_[sectionY] - 1
        );
    }

    return std::pair<int, int>
    (
        min(shiftIdxX + localCellX, n_ - 1),
        min(shiftIdxY + localCellY, m_ - 1)
    );
}


std::pair<int, int> Foam::biLinearInterpolation::findPosition
(
    const scalar xx,
    const scalar yy
) const
{
    label lowerX = 0;
    label upperX = n_;
    while (upperX - lowerX > 1)
    {
        const label middle = (lowerX + upperX)/2;
        if (xx < x_[middle])
        {
            upperX = middle;
        }
        else
        {
            lowerX = middle;
        }
    }

    label lowerY = 0;
    label upperY = m_;
    while (upperY - lowerY > 1)
    {
        const label middle = (lowerY + upperY)/2;
        if (yy < y_[middle])
        {
            upperY = middle;
        }
        else
        {
            lowerY = middle;
        }
    }

    return std::pair<int, int>(lowerX, lowerY);
}


Foam::scalar Foam::biLinearInterpolation::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    return calc(valueX, valueY);
}


// ************************************************************************* //
