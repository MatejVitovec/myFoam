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

#include "biCubicInterpolation.H"
#include "IFstream.H"
#include "OFstream.H"
#include <cmath>


Foam::biCubicInterpolation::biCubicInterpolation(const fileName& fName)
:
    fileName_(fName),
    n_(0),
    m_(0),
    transformationX_(NONE),
    transformationY_(NONE)
{
    readTable();
}


Foam::biCubicInterpolation::biCubicInterpolation
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


void Foam::biCubicInterpolation::readTable()
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
    coefficients_ = dict.get<List<FixedList<scalar, 16>>>("coefficients");

    if (coefficients_.size() != n_*m_)
    {
        FatalErrorInFunction
            << "Coefficient file " << fileName_ << " contains "
            << coefficients_.size() << " cell entries, expected "
            << n_*m_ << nl
            << exit(FatalError);
    }
}


Foam::biCubicInterpolation::Transformation
Foam::biCubicInterpolation::readTransformation(const word& name)
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


Foam::word Foam::biCubicInterpolation::transformationName
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


void Foam::biCubicInterpolation::check() const
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

    label totalX = 0;
    forAll(gridSizeX_, i)
    {
        if (gridSizeX_[i] <= 0 || boundaryX_[i + 1] <= boundaryX_[i])
        {
            FatalErrorInFunction
                << "The x grid sizes must be positive and x boundaries "
                << "strictly increasing" << nl
                << exit(FatalError);
        }
        totalX += gridSizeX_[i];
    }

    label totalY = 0;
    forAll(gridSizeY_, j)
    {
        if (gridSizeY_[j] <= 0 || boundaryY_[j + 1] <= boundaryY_[j])
        {
            FatalErrorInFunction
                << "The y grid sizes must be positive and y boundaries "
                << "strictly increasing" << nl
                << exit(FatalError);
        }
        totalY += gridSizeY_[j];
    }

    if (totalX < 2 || totalY < 2)
    {
        FatalErrorInFunction
            << "At least two cells are required in each direction" << nl
            << exit(FatalError);
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


void Foam::biCubicInterpolation::write(Ostream& os) const
{
    os.writeEntry("gridSizeX", gridSizeX_);
    os.writeEntry("gridSizeY", gridSizeY_);
    os.writeEntry("boundaryX", boundaryX_);
    os.writeEntry("boundaryY", boundaryY_);
    os.writeEntry("transformationX", transformationName(transformationX_));
    os.writeEntry("transformationY", transformationName(transformationY_));
    os.writeEntry("coefficients", coefficients_);
}


void Foam::biCubicInterpolation::write(const fileName& fName) const
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


void Foam::biCubicInterpolation::createGrid()
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


void Foam::biCubicInterpolation::readSettings(const dictionary& dict)
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


void Foam::biCubicInterpolation::calcCoeffs
(
    const std::function<scalar(scalar, scalar)>& f
)
{
    scalarMatrix u(n_ + 1, m_ + 1);
    scalarMatrix p(n_ + 1, m_ + 1);
    scalarMatrix q(n_ + 1, m_ + 1);
    scalarMatrix r(n_ + 1, m_ + 1);

    for (label i = 0; i <= n_; ++i)
    {
        for (label j = 0; j <= m_; ++j)
        {
            u(i, j) = f(x_[i], y_[j]);
        }
    }

    for (label j = 0; j <= m_; ++j)
    {
        p(0, j) = (u(1, j) - u(0, j))/dx_[0];
        p(n_, j) = (u(n_, j) - u(n_ - 1, j))/dx_[n_ - 1];
    }
    for (label i = 0; i <= n_; ++i)
    {
        q(i, 0) = (u(i, 1) - u(i, 0))/dy_[0];
        q(i, m_) = (u(i, m_) - u(i, m_ - 1))/dy_[m_ - 1];
    }

    for (label iSide = 0; iSide < 2; ++iSide)
    {
        const label i = iSide ? n_ - 1 : 0;
        for (label jSide = 0; jSide < 2; ++jSide)
        {
            const label j = jSide ? m_ - 1 : 0;
            r(iSide ? n_ : 0, jSide ? m_ : 0) =
                (u(i + 1, j + 1) - u(i, j + 1) - u(i + 1, j) + u(i, j))
               /(dx_[i]*dy_[j]);
        }
    }

    for (label j = 0; j <= m_; ++j)
    {
        scalarField values(n_ + 1);
        for (label i = 0; i <= n_; ++i)
        {
            values[i] = u(i, j);
        }
        const scalarField derivatives
        (
            solveSplineDerivatives(values, dx_, p(0, j), p(n_, j))
        );
        for (label i = 0; i <= n_; ++i)
        {
            p(i, j) = derivatives[i];
        }
    }

    for (label i = 0; i <= n_; ++i)
    {
        scalarField values(m_ + 1);
        for (label j = 0; j <= m_; ++j)
        {
            values[j] = u(i, j);
        }
        const scalarField derivatives
        (
            solveSplineDerivatives(values, dy_, q(i, 0), q(i, m_))
        );
        for (label j = 0; j <= m_; ++j)
        {
            q(i, j) = derivatives[j];
        }
    }

    for (label side = 0; side < 2; ++side)
    {
        const label j = side ? m_ : 0;
        scalarField values(n_ + 1);
        for (label i = 0; i <= n_; ++i)
        {
            values[i] = q(i, j);
        }
        const scalarField derivatives
        (
            solveSplineDerivatives
            (
                values,
                dx_,
                r(0, j),
                r(n_, j)
            )
        );
        for (label i = 0; i <= n_; ++i)
        {
            r(i, j) = derivatives[i];
        }
    }

    for (label i = 0; i <= n_; ++i)
    {
        scalarField values(m_ + 1);
        for (label j = 0; j <= m_; ++j)
        {
            values[j] = p(i, j);
        }
        const scalarField derivatives
        (
            solveSplineDerivatives
            (
                values,
                dy_,
                r(i, 0),
                r(i, m_)
            )
        );
        for (label j = 0; j <= m_; ++j)
        {
            r(i, j) = derivatives[j];
        }
    }

    const Mat4x4 transformMatrix
    ({
        1, 0, 0, 0,
        0, 0, 1, 0,
       -3, 3,-2,-1,
        2,-2, 1, 1
    });

    for (label j = 0; j < m_; ++j)
    {
        for (label i = 0; i < n_; ++i)
        {
            const scalar hx = dx_[i];
            const scalar hy = dy_[j];

            const Mat4x4 data
            ({
                u(i, j),       u(i, j + 1),       q(i, j)*hy,
                q(i, j + 1)*hy,
                u(i + 1, j),   u(i + 1, j + 1),   q(i + 1, j)*hy,
                q(i + 1, j + 1)*hy,
                p(i, j)*hx,    p(i, j + 1)*hx,    r(i, j)*hx*hy,
                r(i, j + 1)*hx*hy,
                p(i + 1, j)*hx, p(i + 1, j + 1)*hx,
                r(i + 1, j)*hx*hy, r(i + 1, j + 1)*hx*hy
            });

            coefficients_[j*n_ + i] =
                (transformMatrix*data*transformMatrix.transpose()).data();
        }
    }
}


Foam::scalarField Foam::biCubicInterpolation::solveSplineDerivatives
(
    const scalarField& values,
    const scalarField& intervals,
    const scalar leftDerivative,
    const scalar rightDerivative
) const
{
    const label nValues = values.size();
    scalarField derivatives(nValues, Zero);
    derivatives[0] = leftDerivative;
    derivatives[nValues - 1] = rightDerivative;

    if (nValues <= 2)
    {
        return derivatives;
    }

    const label nInterior = nValues - 2;
    scalarField lower(nInterior - 1, Zero);
    scalarField diagonal(nInterior, Zero);
    scalarField upper(nInterior - 1, Zero);
    scalarField rhs(nInterior, Zero);

    for (label row = 0; row < nInterior; ++row)
    {
        const label i = row + 1;
        const scalar hLeft = intervals[i - 1];
        const scalar hRight = intervals[i];

        diagonal[row] = 2.0*(hLeft + hRight);
        rhs[row] =
            3.0
           *(
                hRight*(values[i] - values[i - 1])/hLeft
              + hLeft*(values[i + 1] - values[i])/hRight
            );

        if (row == 0)
        {
            rhs[row] -= hRight*leftDerivative;
        }
        else
        {
            lower[row - 1] = hRight;
        }

        if (row == nInterior - 1)
        {
            rhs[row] -= hLeft*rightDerivative;
        }
        else
        {
            upper[row] = hLeft;
        }
    }

    const scalarField interior
    (
        solveTridiagonal(lower, diagonal, upper, rhs)
    );
    for (label i = 0; i < nInterior; ++i)
    {
        derivatives[i + 1] = interior[i];
    }
    return derivatives;
}


Foam::scalarField Foam::biCubicInterpolation::solveTridiagonal
(
    const scalarField& lower,
    const scalarField& diagonal,
    const scalarField& upper,
    const scalarField& rhs
) const
{
    const label n = rhs.size();
    scalarField out(n, Zero);

    if (n == 1)
    {
        out[0] = rhs[0]/diagonal[0];
        return out;
    }

    scalarField upperStar(n - 1, Zero);
    scalarField rhsStar(n, Zero);

    upperStar[0] = upper[0]/diagonal[0];
    rhsStar[0] = rhs[0]/diagonal[0];

    for (label i = 1; i < n - 1; ++i)
    {
        const scalar inverse
        (
            1.0/(diagonal[i] - lower[i - 1]*upperStar[i - 1])
        );
        upperStar[i] = upper[i]*inverse;
        rhsStar[i] = (rhs[i] - lower[i - 1]*rhsStar[i - 1])*inverse;
    }

    rhsStar[n - 1] =
        (rhs[n - 1] - lower[n - 2]*rhsStar[n - 2])
       /(diagonal[n - 1] - lower[n - 2]*upperStar[n - 2]);
    out[n - 1] = rhsStar[n - 1];

    for (label i = n - 2; i >= 0; --i)
    {
        out[i] = rhsStar[i] - upperStar[i]*out[i + 1];
    }
    return out;
}


Foam::scalar Foam::biCubicInterpolation::transform
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


Foam::scalar Foam::biCubicInterpolation::backTransform
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


Foam::scalar Foam::biCubicInterpolation::calc
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
    const Mat4x4 coeff(coefficients_[j*n_ + i]);

    scalar value = Zero;
    for (label row = 3; row >= 0; --row)
    {
        scalar polynomial = coeff(row, 3);
        for (label col = 2; col >= 0; --col)
        {
            polynomial = polynomial*ty + coeff(row, col);
        }
        value = value*tx + polynomial;
    }
    return value;
}


std::pair<int, int> Foam::biCubicInterpolation::findPosition
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


std::pair<int, int> Foam::biCubicInterpolation::fastFindPosition
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


Foam::scalar Foam::biCubicInterpolation::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    return calc(valueX, valueY);
}


// ************************************************************************* //
