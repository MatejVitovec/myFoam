/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "biQuadraticInterpolation.H"
#include "IFstream.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biQuadraticInterpolation::biQuadraticInterpolation
(
    const fileName& fName
)
:
    fileName_(fName),
    n_(0),
    m_(0),
    transformationX_(NONE),
    transformationY_(NONE)
{
    readTable();
}


Foam::biQuadraticInterpolation::biQuadraticInterpolation
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


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::biQuadraticInterpolation::readTable()
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

    coefficients_ = dict.get<List<FixedList<scalar, 9>>>("coefficients");

    if (coefficients_.size() != n_*m_)
    {
        FatalErrorInFunction
            << "Coefficient file " << fileName_ << " contains "
            << coefficients_.size() << " cell entries, expected "
            << n_*m_ << nl
            << exit(FatalError);
    }
}


Foam::biQuadraticInterpolation::Transformation
Foam::biQuadraticInterpolation::readTransformation(const word& name)
{
    if (name == "none")
    {
        return NONE;
    }
    if (name == "log")
    {
        return LOG;
    }
    if (name == "log10")
    {
        return LOG10;
    }
    if (name == "loginv")
    {
        return LOGINV;
    }

    FatalErrorInFunction
        << "Unknown coordinate transformation '" << name << "'. Valid values: "
        << "none, log, log10, loginv" << nl
        << exit(FatalError);

    return NONE;
}


Foam::word Foam::biQuadraticInterpolation::transformationName
(
    const Transformation transformation
)
{
    switch (transformation)
    {
        case NONE:
            return "none";
        case LOG:
            return "log";
        case LOG10:
            return "log10";
        case LOGINV:
            return "loginv";
    }

    return "none";
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::biQuadraticInterpolation::check() const
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


void Foam::biQuadraticInterpolation::write(Ostream& os) const
{
    os.writeEntry("gridSizeX", gridSizeX_);
    os.writeEntry("gridSizeY", gridSizeY_);
    os.writeEntry("boundaryX", boundaryX_);
    os.writeEntry("boundaryY", boundaryY_);
    os.writeEntry("transformationX", transformationName(transformationX_));
    os.writeEntry("transformationY", transformationName(transformationY_));
    os.writeEntry("coefficients", coefficients_);
}


void Foam::biQuadraticInterpolation::write(const fileName& fName) const
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


void Foam::biQuadraticInterpolation::createGrid()
{
    n_ = sum(gridSizeX_);
    m_ = sum(gridSizeY_);

    x_.setSize(n_);
    z_.setSize(n_ + 1);
    dx_.setSize(n_ + 1);
    dz_.setSize(gridSizeX_.size());

    y_.setSize(m_);
    t_.setSize(m_ + 1);
    dy_.setSize(m_ + 1);
    dt_.setSize(gridSizeY_.size());

    label xIndex = 0;
    z_[0] = boundaryX_[0];

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

        dz_[section] = mag(transformedStep);

        for (label localCell = 0; localCell < gridSizeX_[section]; localCell++)
        {
            x_[xIndex] = backTransform
            (
                transformedStart + (localCell + 0.5)*transformedStep,
                transformationX_
            );
            z_[xIndex + 1] = backTransform
            (
                transformedStart + (localCell + 1)*transformedStep,
                transformationX_
            );
            xIndex++;
        }
    }

    label yIndex = 0;
    t_[0] = boundaryY_[0];

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

        dt_[section] = mag(transformedStep);

        for (label localCell = 0; localCell < gridSizeY_[section]; localCell++)
        {
            y_[yIndex] = backTransform
            (
                transformedStart + (localCell + 0.5)*transformedStep,
                transformationY_
            );
            t_[yIndex + 1] = backTransform
            (
                transformedStart + (localCell + 1)*transformedStep,
                transformationY_
            );
            yIndex++;
        }
    }

    for (label i = 1; i < n_; i++)
    {
        dx_[i] = x_[i] - x_[i - 1];
    }
    dx_[0] = dx_[1];
    dx_[n_] = dx_[n_ - 1];

    for (label j = 1; j < m_; j++)
    {
        dy_[j] = y_[j] - y_[j - 1];
    }
    dy_[0] = dy_[1];
    dy_[m_] = dy_[m_ - 1];

    coefficients_.setSize(n_*m_);
}


void Foam::biQuadraticInterpolation::readSettings(const dictionary& dict)
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


void Foam::biQuadraticInterpolation::calcCoeffs
(
    const std::function<scalar(scalar, scalar)>& f
)
{
    scalarMatrix u(n_,   m_);
    scalarMatrix p(n_+1, m_);
    scalarMatrix q(n_,   m_+1);
    scalarMatrix r(n_+1, m_+1);

    for (int i = 0; i < n_; i++)
    {
        for (int j = 0; j < m_; j++)
        {
            //u[i][j] = f(x_[i], y_[j]);
            u(i, j) = f(x_[i], y_[j]);
        }        
    }

    const scalar divider = 100.0;

    for (int j = 0; j < m_; j++)
    {
        p(0, j)  = calcDerivativeX(f, z_[0],  y_[j], dx_[0]/divider);
        p(n_, j) = calcDerivativeX(f, z_[n_], y_[j], dx_[n_]/divider);
    }

    for (int i = 0; i < n_; i++)
    {
        q(i, 0)  = calcDerivativeY(f, x_[i], t_[0], dy_[0]/divider);
        q(i, m_) = calcDerivativeY(f, x_[i], t_[m_], dy_[m_]/divider);
    }

    r(0, 0)   = calcDerivativeXY(f, z_[0],  t_[0],  dx_[0]/divider,  dy_[0]/divider);
    r(n_, 0)  = calcDerivativeXY(f, z_[n_], t_[0],  dx_[n_]/divider, dy_[0]/divider);
    r(0, m_)  = calcDerivativeXY(f, z_[0],  t_[m_], dx_[0]/divider,  dy_[m_]/divider);
    r(n_, m_) = calcDerivativeXY(f, z_[n_], t_[m_], dx_[n_]/divider, dy_[m_]/divider);
    

    // Solve the derivative spline independently along each grid line.
    for (int j = 0; j < m_; j++)
    {
        scalarField values(n_);
        for (int i = 0; i < n_; i++)
        {
            values[i] = u(i, j);
        }

        const scalarField derivatives
        (
            solveSplineDerivatives(values, dx_, p(0, j), p(n_, j))
        );

        for (int i = 0; i < n_ + 1; i++)
        {
            p(i, j) = derivatives[i];
        }
    }

    for (int i = 0; i < n_; i++)
    {
        scalarField values(m_);
        for (int j = 0; j < m_; j++)
        {
            values[j] = u(i, j);
        }

        const scalarField derivatives
        (
            solveSplineDerivatives(values, dy_, q(i, 0), q(i, m_))
        );

        for (int j = 0; j < m_ + 1; j++)
        {
            q(i, j) = derivatives[j];
        }
    }

    for (int side = 0; side < 2; side++)
    {
        const label j = side ? m_ : 0;
        scalarField values(n_);
        for (int i = 0; i < n_; i++)
        {
            values[i] = q(i, j);
        }

        const scalarField derivatives
        (
            solveSplineDerivatives(values, dx_, r(0, j), r(n_, j))
        );

        for (int i = 0; i < n_ + 1; i++)
        {
            r(i, j) = derivatives[i];
        }
    }

    for (int i = 0; i < n_ + 1; i++)
    {
        scalarField values(m_);
        for (int j = 0; j < m_; j++)
        {
            values[j] = p(i, j);
        }

        const scalarField derivatives
        (
            solveSplineDerivatives(values, dy_, r(i, 0), r(i, m_))
        );

        for (int j = 0; j < m_ + 1; j++)
        {
            r(i, j) = derivatives[j];
        }
    }

    for (label j = 0; j < m_; j++)
    {
        for (label i = 0; i < n_; i++)
        {
            Mat3x3 VxInv({0 , 1.0, 0,
                          dx_[i+1]/(dx_[i] + dx_[i+1]), 0, dx_[i]/(dx_[i] + dx_[i+1]),
                          -1.0/(dx_[i] + dx_[i+1]), 0, 1.0/(dx_[i] + dx_[i+1])});

            Mat3x3 VyInv({0 , 1.0, 0,
                          dy_[j+1]/(dy_[j] + dy_[j+1]), 0, dy_[j]/(dy_[j] + dy_[j+1]),
                          -1.0/(dy_[j] + dy_[j+1]), 0, 1.0/(dy_[j] + dy_[j+1])});

            Mat3x3 C({r(i, j), p(i, j), r(i, j+1),
                      q(i, j), u(i, j), q(i, j+1),
                      r(i+1, j), p(i+1, j), r(i+1, j+1)});

            coefficients_[j*n_ + i] = (VxInv*C*(VyInv.transpose())).data();
        }
    }
}


Foam::scalarField Foam::biQuadraticInterpolation::solveSplineDerivatives
(
    const scalarField& values,
    const scalarField& intervals,
    const scalar leftDerivative,
    const scalar rightDerivative
) const
{
    const label nValues = values.size();
    scalarField derivatives(nValues + 1, Zero);
    derivatives[0] = leftDerivative;
    derivatives[nValues] = rightDerivative;

    if (nValues <= 1)
    {
        return derivatives;
    }

    scalarField lower(nValues - 2, Zero);
    scalarField diagonal(nValues - 1, Zero);
    scalarField upper(nValues - 2, Zero);
    scalarField rhs(nValues - 1, Zero);

    for (label row = 0; row < nValues - 1; row++)
    {
        const scalar hLeft = intervals[row];
        const scalar h = intervals[row + 1];
        const scalar hRight = intervals[row + 2];

        diagonal[row] = (1.0/h)
           *(
                2.0
              + hLeft/(hLeft + h)
              + hRight/(h + hRight)
            );

        rhs[row] = (4.0/h)*(values[row + 1] - values[row])/h;

        if (row == 0)
        {
            rhs[row] -= leftDerivative/(hLeft + h);
        }
        else
        {
            lower[row - 1] = 1.0/(hLeft + h);
        }

        if (row == nValues - 2)
        {
            rhs[row] -= rightDerivative/(h + hRight);
        }
        else
        {
            upper[row] = 1.0/(h + hRight);
        }
    }

    const scalarField interior
    (
        solveTridiagonal(lower, diagonal, upper, rhs)
    );

    for (label i = 0; i < interior.size(); i++)
    {
        derivatives[i + 1] = interior[i];
    }

    return derivatives;
}


Foam::scalarField Foam::biQuadraticInterpolation::solveTridiagonal
(
    const Foam::scalarField& L,
    const Foam::scalarField& D,
    const Foam::scalarField& U,
    const Foam::scalarField& b
) const
{
    const label n = b.size();
    scalarField out(n);

    if (n == 1)
    {
        out[0] = b[0]/D[0];
        return out;
    }
                                                                                                                                                                       
    scalarField UStar(n-1, 0.0);
    scalarField bStar(n, 0.0);
                                                                                                                                                    
    UStar[0] = U[0] / D[0];
    bStar[0] = b[0] / D[0];
                                                                                                                                            
    for (int i = 1; i < n-1; i++)
    {
        scalar m = 1.0/(D[i] - L[i-1]*UStar[i-1]);
        UStar[i] = U[i] * m;
        bStar[i] = (b[i] - L[i-1]*bStar[i-1])*m;
    }
    bStar[n-1] = (b[n-1] - L[n-2]*bStar[n-2])*(1.0/(D[n-1] - L[n-2] * UStar[n-2]));

    out[n-1] = bStar[n-1];
                                                                                                                                            
    for (int i = n-2; i >= 0; i--)
    {
        out[i] = bStar[i] - UStar[i]*out[i+1];
    }

    return out;
}


Foam::scalar Foam::biQuadraticInterpolation::calcDerivativeX
(
    const std::function<scalar(scalar, scalar)>& f,
    const scalar x,
    const scalar y,
    const scalar step
) const
{
    return (f(x + step, y) - f(x - step, y))/(2.0*step);
}


Foam::scalar Foam::biQuadraticInterpolation::calcDerivativeY
(
    const std::function<scalar(scalar, scalar)>& f,
    const scalar x,
    const scalar y,
    const scalar step
) const
{
    return (f(x, y + step) - f(x, y - step))/(2.0*step);
}


Foam::scalar Foam::biQuadraticInterpolation::calcDerivativeXY
(
    const std::function<scalar(scalar, scalar)>& f,
    const scalar x,
    const scalar y,
    const scalar stepX,
    const scalar stepY
) const
{
    return
    (
        f(x + stepX, y + stepY)
      - f(x + stepX, y - stepY)
      - f(x - stepX, y + stepY)
      + f(x - stepX, y - stepY)
    )/(4.0*stepX*stepY);
}


Foam::scalar Foam::biQuadraticInterpolation::transform(Foam::scalar x, Transformation transformation) const
{
    switch(transformation) 
    {
        case LOG:
            return log(x);
        case LOG10:
            return log10(x);
        case LOGINV:
            return -log(x);
        case NONE:
            return x;
    }

    return x;
}

Foam::scalar Foam::biQuadraticInterpolation::backTransform(Foam::scalar x, Transformation transformation) const
{
    switch(transformation) 
    {
        case LOG:
            return exp(x);
        case LOG10:
            return pow(10.0, x);
        case LOGINV:
            return exp(-x);
        case NONE:
            return x;
    }

    return x;
}


Foam::scalar Foam::biQuadraticInterpolation::calc
(
    Foam::scalar xx,
    Foam::scalar yy
) const
{
    const std::pair<int, int> position = fastFindPosition(xx, yy);

    scalar v = xx - x_[position.first];
    scalar w = yy - y_[position.second];

    const Mat3x3& coeff = coefficients_[position.second*n_ + position.first];

    return coeff(0, 0) + w*(coeff(0, 1) + w*coeff(0, 2)) + v*(coeff(1, 0) + w*(coeff(1, 1) + w*coeff(1, 2)) + v*(coeff(2, 0) + w*(coeff(2, 1) + w*coeff(2, 2))));
}


std::pair<int, int> Foam::biQuadraticInterpolation::findPosition(Foam::scalar xx, Foam::scalar yy) const
{
    label cellX = n_ - 1;
    for (label i = 0; i < n_; i++)
    {
        if (xx < z_[i + 1] || i == n_ - 1)
        {
            cellX = i;
            break;
        }
    }

    label cellY = m_ - 1;
    for (label j = 0; j < m_; j++)
    {
        if (yy < t_[j + 1] || j == m_ - 1)
        {
            cellY = j;
            break;
        }
    }

    return std::pair<int, int>(cellX, cellY);
}


std::pair<int, int> Foam::biQuadraticInterpolation::fastFindPosition
(
    const scalar xx,
    const scalar yy
) const
{
    label sectionX = 0;
    label sectionY = 0;
    label shiftIdxX = 0;
    label shiftIdxY = 0;

    for (; sectionX < gridSizeX_.size() - 1; sectionX++)
    {
        if (xx < boundaryX_[sectionX + 1])
        {
            break;
        }
        shiftIdxX += gridSizeX_[sectionX];
    }

    for (; sectionY < gridSizeY_.size() - 1; sectionY++)
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
            label(floor(distance/dz_[sectionX])),
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
            label(floor(distance/dt_[sectionY])),
            gridSizeY_[sectionY] - 1
        );
    }

    return std::pair<int, int>
    (
        min(shiftIdxX + localCellX, n_ - 1),
        min(shiftIdxY + localCellY, m_ - 1)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::scalar Foam::biQuadraticInterpolation::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    return calc(valueX, valueY);
}


// ************************************************************************* //

