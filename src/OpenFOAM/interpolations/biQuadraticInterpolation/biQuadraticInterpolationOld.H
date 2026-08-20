/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::biQuadraticInterpolation

Description
    biQuadraticInterpolation

SourceFiles
    biQuadraticInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef biQuadraticInterpolation_H
#define biQuadraticInterpolation_H

//#include "interpolationTable.H"
#include <eigen3/Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class biQuadraticInterpolation Declaration
\*---------------------------------------------------------------------------*/

//template<class Type>
class biQuadraticInterpolation
{
    // Private Data

        //- File name
        fileName fileName_;

        std::vector<int> gridSizeX;
        std::vector<int> gridSizeY;
        std::Eigen::VectorXd boundaryX;
        std::Eigen::VectorXd boundaryY;
        Transformation transformationX;
        Transformation transformationY;

        std::Eigen::VectorXd x;
        std::Eigen::VectorXd y;
        std::Eigen::VectorXd z;
        std::Eigen::VectorXd t;
        std::Eigen::VectorXd dx;
        std::Eigen::VectorXd dy;

        std::Eigen::VectorXd dz;
        std::Eigen::VectorXd dt;

        Eigen::Matrix<double, 9, Eigen::Dynamic> coeffs_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        void calcCoeffs();

        Eigen::VectorXd solveTridiagonal(const Eigen::VectorXd& L, const Eigen::VectorXd& D, const Eigen::VectorXd& U, const Eigen::VectorXd& b) const;
        scalar calcDerivativeX (std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar step) const;
        scalar calcDerivativeY (std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar step) const;
        scalar calcDerivativeXY(std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar stepX, scalar stepY) const;
        
        double transform(double x, Transformation transformation) const;
        
        double backTransform(double x, Transformation transformation) const;

        scalar calc(scalar xx, scalar yy);


public:
    enum Transformation{NONE, LOG, LOG10, LOGINV};

    // Constructors

        //- Default construct
        biQuadraticInterpolation();

        //- Construct from components
        biQuadraticInterpolation
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
            const bounds::normalBounding bounding,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        explicit biQuadraticInterpolation(const fileName& fName);

        //- Construct by reading file name and outOfBounds from dictionary
        explicit biQuadraticInterpolation(const dictionary& dict);

        //- Copy construct
        biQuadraticInterpolation(const biQuadraticInterpolation& tbl);


    // Member Functions

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

        void loadCoeffs();

        //- Write
        void write(Ostream& os) const;

        std::pair<int, int> findPosition(double xx, double yy) const;


    // Member Operators

        //- Copy assignment
        void operator=(const biQuadraticInterpolation<Type>& rhs);

        //- Return an interpolated value
        Type operator()(const scalar valueX, const scalar valueY) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

