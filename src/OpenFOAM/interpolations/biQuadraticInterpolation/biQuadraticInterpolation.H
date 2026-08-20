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
#include "IOobject.H"
#include "scalar.H"
#include "label.H"
#include "scalarField.H"
#include "FixedList.H"
#include <eigen3/Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class biQuadraticInterpolation Declaration
\*---------------------------------------------------------------------------*/

class biQuadraticInterpolation
{
    public:
        enum Transformation{NONE, LOG, LOG10, LOGINV};

    struct Mat3x3
    {
        FixedList<scalar, 9> v_;

        // Constructors
        Mat3x3() : v_(Zero) {}

        Mat3x3(const FixedList<scalar, 9>& values) : v_(values) {}

        // Element access
        inline scalar& operator()(const label i, const label j)
        {
            return v_[3*i + j];
        }

        inline const scalar& operator()(const label i, const label j) const
        {
            return v_[3*i + j];
        }

        // Transpose
        inline Mat3x3 transpose() const
        {
            Mat3x3 result;

            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    result(i, j) = (*this)(j, i);
                }
            }

            return result;
        }

        // Matrix multiplication
        inline Mat3x3 operator*(const Mat3x3& B) const
        {
            Mat3x3 result;

            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    result(i, j) = Zero;

                    for (label k = 0; k < 3; ++k)
                    {
                        result(i, j) += (*this)(i, k)*B(k, j);
                    }
                }
            }

            return result;
        }

        // Access to underlying data
        inline FixedList<scalar, 9>& data()
        {
            return v_;
        }

        inline const FixedList<scalar, 9>& data() const
        {
            return v_;
        }
    };


    class scalarMatrix
    {
        label nRows_;
        label nCols_;
        scalarField data_;

    public:

        // Constructors
        scalarMatrix() : nRows_(0), nCols_(0), data_() {}

        scalarMatrix ( const label nRows, const label nCols)
        :
            nRows_(nRows),
            nCols_(nCols),
            data_(nRows*nCols, Zero)
        {}


        inline label nRows() const
        {
            return nRows_;
        }

        inline label nCols() const
        {
            return nCols_;
        }

        inline label size() const
        {
            return data_.size();
        }

        // Element access
        inline scalar& operator()(const label i, const label j)
        {
            return data_[i*nCols_ + j];
        }

        inline const scalar& operator()(const label i, const label j) const
        {
            return data_[i*nCols_ + j];
        }

        // Access to underlying data
        inline scalarField& data()
        {
            return data_;
        }

        inline const scalarField& data() const
        {
            return data_;
        }
    };




    /////////////



    // Private Data

        //- File name
        fileName fileName_;

        int n_;
        int m_;

        labelList gridSizeX_;
        labelList gridSizeY_;
        scalarField boundaryX_;
        scalarField boundaryY_;
        Transformation transformationX_;
        Transformation transformationY_;

        scalarField x_;
        scalarField y_;
        scalarField z_;
        scalarField t_;
        scalarField dx_;
        scalarField dy_;

        scalarField dz_;
        scalarField dt_;

        List<FixedList<scalar, 9>> coefficients_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        void calcCoeffs(std::function<scalar(scalar, scalar)> f);

        scalarField solveTridiagonal(const scalarField& L, const scalarField& D, const scalarField& U, const scalarField& b) const;
        scalar calcDerivativeX (std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar step) const;
        scalar calcDerivativeY (std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar step) const;
        scalar calcDerivativeXY(std::function<scalar(scalar, scalar)> f, scalar x, scalar y, scalar stepX, scalar stepY) const;
        
        scalar transform(scalar x, Transformation transformation) const;
        
        scalar backTransform(scalar x, Transformation transformation) const;

        scalar calc(scalar xx, scalar yy);


public:

    // Constructors

        //- Default construct
        biQuadraticInterpolation();

        //- Construct from components
        /*biQuadraticInterpolation
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>>& values,
            const bounds::normalBounding bounding,
            const fileName& fName
        );*/

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

        std::pair<int, int> findPosition(scalar xx, scalar yy) const;


    // Member Operators

        //- Copy assignment
        void operator=(const biQuadraticInterpolation& rhs);

        //- Return an interpolated value
        scalar operator()(const scalar valueX, const scalar valueY) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

