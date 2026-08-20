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
//#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

/*void Foam::biQuadraticInterpolation::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    check();
}*/


void Foam::biQuadraticInterpolation::calcCoeffs(std::function<scalar(scalar, scalar)> f)
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
    

    /////// solution for p
    for (int j = 0; j < m_; j++)
    {
        scalarField U(n_-2);
        scalarField D(n_-1);
        scalarField L(n_-2);
        scalarField b(n_-1);
        int i = 0;

        U[i] = 1.0/(dx_[i+1] + dx_[i+2]);
        D[i] = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));

        b[i] = (-1.0/(dx_[i] + dx_[i+1]))*p(0, j) + (4.0/dx_[i+1])*((u(i+1, j) - u(i, j))/dx_[i+1]);

        i++;
        for ( ; i < n_ - 2; i++)
        {
            U[i]   = 1.0/(dx_[i+1] + dx_[i+2]);
            D[i]   = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));
            L[i-1] = 1.0/(dx_[i] + dx_[i+1]);
            b[i]   = (4.0/dx_[i+1])*((u(i+1, j) - u(i, j))/dx_[i+1]);
        }

        D[i] = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));
        L[i-1] = 1.0/(dx_[i]+dx_[i+1]);
        b[i] = (-1.0/(dx_[i+1] + dx_[i+2]))*p(n_, j) + (4.0/dx_[i+1])*((u(i+1, j) - u(i, j))/dx_[i+1]);

        scalarField pj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < pj.size()+1; i++)
        {
            p(i, j) = pj[i-1];
        }        
    }

    ///////// solution for q
    for (int i = 0; i < n_; i++)
    {
        scalarField U(m_-2);
        scalarField D(m_-1);
        scalarField L(m_-2);
        scalarField b(m_-1);
        int j = 0;

        U[j] = 1.0/(dy_[j+1] + dy_[j+2]);
        D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));

        b[j] = (-1.0/(dy_[j] + dy_[j+1]))*q(i, 0) + (4.0/dy_[j+1])*((u(i, j+1) - u(i, j))/dy_[j+1]);

        j++;
        for ( ; j < m_-2; j++)
        {
            U[j] = 1.0/(dy_[j+1] + dy_[j+2]);
            D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));
            L[j-1] = 1.0/(dy_[j] + dy_[j+1]);
            b[j] = (4.0/dy_[j+1])*((u(i, j+1) - u(i, j))/dy_[j+1]);
        }

        D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));
        L[j-1] = 1.0/(dy_[j] + dy_[j+1]);
        b[j] = (-1.0/(dy_[j+1] + dy_[j+2]))*q(i, m_) + (4.0/dy_[j+1])*((u(i, j+1) - u(i, j))/dy_[j+1]);

        scalarField qi = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < qi.size()+1; j++)
        {
            q(i, j) = qi[j-1];
        }        
    }

    
    ///////// solution for ri0 and rim
    for (int j = 0; j < m_+1; j += m_)
    {
        scalarField U(n_-2);
        scalarField D(n_-1);
        scalarField L(n_-2);
        scalarField b(n_-1);
        int i = 0;

        U[i] = 1.0/(dx_[i+1] + dx_[i+2]);
        D[i] = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));

        b[i] = (-1.0/(dx_[i] + dx_[i+1]))*r(i, j) + (4.0/dx_[i+1])*((q(i+1, j) - q(i, j))/dx_[i+1]);

        i++;
        for ( ; i < n_-2; i++)
        {
            U[i] = 1.0/(dx_[i+1] + dx_[i+2]);
            D[i] = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));
            L[i-1] = 1.0/(dx_[i] + dx_[i+1]);
            b[i] = (4.0/dx_[i+1])*((q(i+1, j) - q(i, j))/dx_[i+1]);
        }

        D[i] = (1.0/dx_[i+1])*(2.0 + dx_[i]/(dx_[i] + dx_[i+1]) + dx_[i+2]/(dx_[i+1] + dx_[i+2]));
        L[i-1] = 1.0/(dx_[i]+dx_[i+1]);
        b[i] = (-1.0/(dx_[i+1] + dx_[i+2]))*r(i+2, j) + (4.0/dx_[i+1])*((q(i+1, j) - q(i, j))/dx_[i+1]);

        scalarField rj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < rj.size()+1; i++)
        {
            r(i, j) = rj[i-1];
        }
    }


    //////solution for r
    for (int i = 0; i < n_+1; i++)
    {
        scalarField U(m_-2);
        scalarField D(m_-1);
        scalarField L(m_-2);
        scalarField b(m_-1);
        int j = 0;

        U[j] = 1.0/(dy_[j+1] + dy_[j+2]);
        D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));

        b[j] = (-1.0/(dy_[j] + dy_[j+1]))*r(i, 0) + (4.0/dy_[j+1])*((p(i, j+1) - p(i, j))/dy_[j+1]);

        j++;
        for ( ; j < m_-2; j++)
        {
            U[j] = 1.0/(dy_[j+1] + dy_[j+2]);
            D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));
            L[j-1] = 1.0/(dy_[j] + dy_[j+1]);
            b[j] = (4.0/dy_[j+1])*((p(i, j+1) - p(i, j))/dy_[j+1]);
        }

        D[j] = (1.0/dy_[j+1])*(2.0 + dy_[j]/(dy_[j] + dy_[j+1]) + dy_[j+2]/(dy_[j+1] + dy_[j+2]));
        L[j-1] = 1.0/(dy_[j] + dy_[j+1]);
        b[j] = (-1.0/(dy_[j+1] + dy_[j+2]))*r(i, m_) + (4.0/dy_[j+1])*((p(i, j+1) - p(i, j))/dy_[j+1]);

        scalarField ri = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < ri.size()+1; j++)
        {
            r(i, j) = ri[j-1];
        }        
    }

    int idx_ = 0;
    for(int j = 0; j < m_; j++)
    {
        for(int i = 0; i < n_; i++)
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

            coefficients_[idx_] = (VyInv*C*(VyInv.transpose())).data(); //TODO pristup 
            idx_++;
        }
    }
}


Foam::scalarField Foam::biQuadraticInterpolation::solveTridiagonal
(
    const Foam::scalarField& L,
    const Foam::scalarField& D,
    const Foam::scalarField& U,
    const Foam::scalarField& b
) const
{
    int n = b.size();
    scalarField out(n);
                                                                                                                                                                       
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


Foam::scalar Foam::biQuadraticInterpolation::transform(Foam::scalar x, Transformation transformation) const
{
    switch(transformation) 
    {
        case LOG:
            return log(x);
        case LOG10:
            return log10(x);
        case LOGINV:
            return log(1/x);
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
            return pow(x, 10.0);
        case LOGINV:
            return 1.0/exp(x);
        case NONE:
            return x;
    }

    return x;
}

Foam::scalar Foam::biQuadraticInterpolation::calc(Foam::scalar xx, Foam::scalar yy)
{
    std::pair<int, int> position = findPosition(xx, yy);

    scalar v = xx - x_[position.first];
    scalar w = yy - y_[position.second];

    const Mat3x3& coeff = coefficients_[position.second*n_ + position.first];

    return coeff(0, 0) + w*(coeff(0, 1) + w*coeff(0, 2)) + v*(coeff(1, 0) + w*(coeff(1, 1) + w*coeff(1, 2)) + v*(coeff(2, 0) + w*(coeff(2, 1) + w*coeff(2, 2))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::biQuadraticInterpolation::biQuadraticInterpolation()
:
    List<value_type>(),
    bounding_(bounds::normalBounding::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}*/


/*Foam::biQuadraticInterpolation::biQuadraticInterpolation(const fileName& fName)
:
    fileName_(fName)
{
    readTable();
}*/


/*Foam::biQuadraticInterpolation::biQuadraticInterpolation(const dictionary& dict)
:
    fileName_(dict.get<fileName>("file"))
{
    readTable();
}*/


/*Foam::biQuadraticInterpolation::biQuadraticInterpolation
(
     const biQuadraticInterpolation& tbl
)
:
    fileName_(tbl.fileName_)
{}*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*scalar Foam::biQuadraticInterpolation::interpolateValue
(
    const List<Tuple2<scalar, Type>>& list,
    scalar lookupValue
) const
{
    return interpolationTable<Type>::interpolateValue
    (
        list,
        lookupValue,
        bounds::repeatableBounding(bounding_)
    );
}*/


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*void Foam::biQuadraticInterpolation::operator=
(
    const biQuadraticInterpolation& rhs
)
{
    if (this == &rhs)
    {
        return;
    }

    static_cast<List<value_type>&>(*this) = rhs;
    fileName_ = rhs.fileName_;
}*/


/*scalar Foam::biQuadraticInterpolation::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{

    //return (y0 + (y1 - y0)*(valueX - x0)/(x1 - x0));
}*/


/*void Foam::biQuadraticInterpolation::check() const
{

}*/


/*void Foam::biQuadraticInterpolation::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", bounds::normalBoundingNames[bounding_]);

    os  << *this;
}*/


std::pair<int, int> Foam::biQuadraticInterpolation::findPosition(Foam::scalar xx, Foam::scalar yy) const
{
    int shiftIdx_X = 0;
    int shiftIdx_Y = 0;

    int ii;
    for (ii = 0; ii < gridSizeX_.size(); ii++)
    {
        if (xx < boundaryX_[ii+1]) { break; }
        shiftIdx_X += gridSizeX_[ii];
    }
    int jj; 
    for (jj = 0; jj < gridSizeY_.size(); jj++)
    {
        if (yy < boundaryY_[jj+1]) { break; }
        shiftIdx_Y += gridSizeY_[jj];
    }
    
    return std::pair<int, int>(std::floor((abs(transform(xx, transformationX_) - transform(boundaryX_[ii], transformationX_)))/dz_[ii]) + shiftIdx_X,
                               std::floor((abs(transform(yy, transformationY_) - transform(boundaryY_[jj], transformationY_)))/dt_[jj]) + shiftIdx_Y);
}


// ************************************************************************* //

