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

#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::biQuadraticInterpolation::readTable()
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
}


void Foam::biQuadraticInterpolation::calcCoeffs();
{
    Eigen::MatrixXd u(n, m);
    Eigen::MatrixXd p(n+1, m);
    Eigen::MatrixXd q(n, m+1);
    Eigen::MatrixXd r(n+1, m+1);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            u[i][j] = f(x[i], y[j]);
        }        
    }

    const double divider = 100.0;

    for (int j = 0; j < m; j++)
    {
        p[0][j] = calcDerivativeX(f, z[0], y[j], dx[0]/divider);
        p[n][j] = calcDerivativeX(f, z[n], y[j], dx[n]/divider);
    }

    for (int i = 0; i < n; i++)
    {
        q[i][0] = calcDerivativeY(f, x[i], t[0], dy[0]/divider);
        q[i][m] = calcDerivativeY(f, x[i], t[m], dy[m]/divider);
    }

    r[0][0] = calcDerivativeXY(f, z[0], t[0], dx[0]/divider, dy[0]/divider);
    r[n][0] = calcDerivativeXY(f, z[n], t[0], dx[n]/divider, dy[0]/divider);
    r[0][m] = calcDerivativeXY(f, z[0], t[m], dx[0]/divider, dy[m]/divider);
    r[n][m] = calcDerivativeXY(f, z[n], t[m], dx[n]/divider, dy[m]/divider);
    

    /////// solution for p
    for (int j = 0; j < m; j++)
    {
        Eigen::VectorXd U(n-2);
        Eigen::VectorXd D(n-1);
        Eigen::VectorXd L(n-2);
        Eigen::VectorXd b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*p[0][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*p[n][j] + (4.0/dx[i+1])*((u[i+1][j] - u[i][j])/dx[i+1]);

        Eigen::VectorXd pj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < pj.size()+1; i++)
        {
            p[i][j] = pj[i-1];
        }        
    }

    ///////// solution for q
    for (int i = 0; i < n; i++)
    {
        Eigen::VectorXd U(m-2);
        Eigen::VectorXd D(m-1);
        Eigen::VectorXd L(m-2);
        Eigen::VectorXd b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*q[i][0] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*q[i][m] + (4.0/dy[j+1])*((u[i][j+1] - u[i][j])/dy[j+1]);

        Eigen::VectorXd qi = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < qi.size()+1; j++)
        {
            q[i][j] = qi[j-1];
        }        
    }

    
    ///////// solution for ri0 and rim
    for (int j = 0; j < m+1; j += m)
    {
        Eigen::VectorXd U(n-2);
        Eigen::VectorXd D(n-1);
        Eigen::VectorXd L(n-2);
        Eigen::VectorXd b(n-1);
        int i = 0;

        U[i] = 1.0/(dx[i+1] + dx[i+2]);
        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));

        b[i] = (-1.0/(dx[i] + dx[i+1]))*r[i][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        i++;
        for ( ; i < n-2; i++)
        {
            U[i] = 1.0/(dx[i+1] + dx[i+2]);
            D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
            L[i-1] = 1.0/(dx[i] + dx[i+1]);
            b[i] = (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);
        }

        D[i] = (1.0/dx[i+1])*(2.0 + dx[i]/(dx[i] + dx[i+1]) + dx[i+2]/(dx[i+1] + dx[i+2]));
        L[i-1] = 1.0/(dx[i]+dx[i+1]);
        b[i] = (-1.0/(dx[i+1] + dx[i+2]))*r[i+2][j] + (4.0/dx[i+1])*((q[i+1][j] - q[i][j])/dx[i+1]);

        Eigen::VectorXd rj = solveTridiagonal(L, D, U, b);
        for (int i = 1; i < rj.size()+1; i++)
        {
            r[i][j] = rj[i-1];
        }
    }


    //////solution for r
    for (int i = 0; i < n+1; i++)
    {
        Eigen::VectorXd U(m-2);
        Eigen::VectorXd D(m-1);
        Eigen::VectorXd L(m-2);
        Eigen::VectorXd b(m-1);
        int j = 0;

        U[j] = 1.0/(dy[j+1] + dy[j+2]);
        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));

        b[j] = (-1.0/(dy[j] + dy[j+1]))*r[i][0] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        j++;
        for ( ; j < m-2; j++)
        {
            U[j] = 1.0/(dy[j+1] + dy[j+2]);
            D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
            L[j-1] = 1.0/(dy[j] + dy[j+1]);
            b[j] = (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);
        }

        D[j] = (1.0/dy[j+1])*(2.0 + dy[j]/(dy[j] + dy[j+1]) + dy[j+2]/(dy[j+1] + dy[j+2]));
        L[j-1] = 1.0/(dy[j] + dy[j+1]);
        b[j] = (-1.0/(dy[j+1] + dy[j+2]))*r[i][m] + (4.0/dy[j+1])*((p[i][j+1] - p[i][j])/dy[j+1]);

        Eigen::VectorXd ri = solveTridiagonal(L, D, U, b);
        for (int j = 1; j < ri.size()+1; j++)
        {
            r[i][j] = ri[j-1];
        }        
    }

    int idx = 0;
    for(int j = 0; j < m; j++)
    {
        for(int i = 0; i < n; i++)
        {
            Mat3x3 VxInv({0 , 1.0, 0,
                          dx[i+1]/(dx[i] + dx[i+1]), 0, dx[i]/(dx[i] + dx[i+1]),
                          -1.0/(dx[i] + dx[i+1]), 0, 1.0/(dx[i] + dx[i+1])});

            Mat3x3 VyInv({0 , 1.0, 0,
                          dy[j+1]/(dy[j] + dy[j+1]), 0, dy[j]/(dy[j] + dy[j+1]),
                          -1.0/(dy[j] + dy[j+1]), 0, 1.0/(dy[j] + dy[j+1])});

            Mat3x3 C({r[i][j], p[i][j], r[i][j+1],
                      q[i][j], u[i][j], q[i][j+1],
                      r[i+1][j], p[i+1][j], r[i+1][j+1]});

            coeffs[idx] = VyInv*C*(VyInv.transpose());
            idx++;
        }
    }
}


Eigen::VectorXd Foam::biQuadraticInterpolation::solveTridiagonal
(
    const Eigen::VectorXd& L,
    const Eigen::VectorXd& D,
    const Eigen::VectorXd& U,
    const Eigen::VectorXd& b
) const
{
    int n = b.size();
    Eigen::VectorXd out(n);
                                                                                                                                                                       
    Eigen::VectorXd UStar(n-1, 0.0);
    Eigen::VectorXd bStar(n, 0.0);
                                                                                                                                                    
    UStar[0] = U[0] / D[0];
    bStar[0] = b[0] / D[0];
                                                                                                                                            
    for (int i = 1; i < n-1; i++)
    {
        double m = 1.0/(D[i] - L[i-1]*UStar[i-1]);
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


double Foam::biQuadraticInterpolation::transform(double x, Transformation transformation) const
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

double Foam::biQuadraticInterpolation::backTransform(double x, Transformation transformation) const
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

scalar Foam::biQuadraticInterpolation::calc(scalar xx, scalar yy)
{
    std::pair<int, int> position = findPosition(xx, yy);

    double v = xx - x[position.first];
    double w = yy - y[position.second];

    const Mat3x3& coeff = coeffs_[position.second*n + position.first];

    return coeff[0][0] + w*(coeff[0][1] + w*coeff[0][2]) + v*(coeff[1][0] + w*(coeff[1][1] + w*coeff[1][2]) + v*(coeff[2][0] + w*(coeff[2][1] + w*coeff[2][2])));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biQuadraticInterpolation::biQuadraticInterpolation()
:
    List<value_type>(),
    bounding_(bounds::normalBounding::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


Foam::biQuadraticInterpolation::biQuadraticInterpolation(const fileName& fName)
:
    fileName_(fName)
{
    readTable();
}


Foam::biQuadraticInterpolation::biQuadraticInterpolation(const dictionary& dict)
:
    fileName_(dict.get<fileName>("file"))
{
    readTable();
}


Foam::biQuadraticInterpolation::biQuadraticInterpolation
(
     const biQuadraticInterpolation& tbl
)
:
    fileName_(tbl.fileName_)
{}


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

void Foam::biQuadraticInterpolation::operator=
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
}


scalar Foam::biQuadraticInterpolation::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{

    //return (y0 + (y1 - y0)*(valueX - x0)/(x1 - x0));
}


void Foam::biQuadraticInterpolation::check() const
{

}


void Foam::biQuadraticInterpolation::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", bounds::normalBoundingNames[bounding_]);

    os  << *this;
}


std::pair<int, int> Foam::BiQuadraticInterpolation::findPosition(double xx, double yy) const
{
    int shiftIdxX = 0;
    int shiftIdxY = 0;

    int ii;
    for (ii = 0; ii < gridSizeX.size(); ii++)
    {
        if (xx < boundaryX[ii+1]) { break; }
        shiftIdxX += gridSizeX[ii];
    }
    int jj; 
    for (jj = 0; jj < gridSizeY.size(); jj++)
    {
        if (yy < boundaryY[jj+1]) { break; }
        shiftIdxY += gridSizeY[jj];
    }
    
    return std::pair<int, int>(std::floor((abs(transform(xx, transformationX) - transform(boundaryX[ii], transformationX)))/dz[ii]) + shiftIdxX,
                               std::floor((abs(transform(yy, transformationY) - transform(boundaryY[jj], transformationY)))/dt[jj]) + shiftIdxY);
}


// ************************************************************************* //

