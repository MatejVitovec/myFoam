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

Class
    Foam::biCubicInterpolation

Description
    Bicubic spline interpolation on a rectangular grid.

SourceFiles
    biCubicInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef biCubicInterpolation_H
#define biCubicInterpolation_H

#include "IOobject.H"
#include "dictionary.H"
#include "scalar.H"
#include "label.H"
#include "scalarField.H"
#include "FixedList.H"
#include "word.H"
#include <functional>
#include <utility>

namespace Foam
{

class biCubicInterpolation
{
public:
    enum Transformation { NONE, LOG, LOG10, LOGINV };

    struct Mat4x4
    {
        FixedList<scalar, 16> v_;

        Mat4x4() : v_(Zero) {}

        Mat4x4(const FixedList<scalar, 16>& values) : v_(values) {}

        inline scalar& operator()(const label i, const label j)
        {
            return v_[4*i + j];
        }

        inline const scalar& operator()(const label i, const label j) const
        {
            return v_[4*i + j];
        }

        inline Mat4x4 transpose() const
        {
            Mat4x4 result;
            for (label i = 0; i < 4; ++i)
            {
                for (label j = 0; j < 4; ++j)
                {
                    result(i, j) = (*this)(j, i);
                }
            }
            return result;
        }

        inline Mat4x4 operator*(const Mat4x4& B) const
        {
            Mat4x4 result;
            for (label i = 0; i < 4; ++i)
            {
                for (label j = 0; j < 4; ++j)
                {
                    for (label k = 0; k < 4; ++k)
                    {
                        result(i, j) += (*this)(i, k)*B(k, j);
                    }
                }
            }
            return result;
        }

        inline FixedList<scalar, 16>& data()
        {
            return v_;
        }

        inline const FixedList<scalar, 16>& data() const
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
        scalarMatrix() : nRows_(0), nCols_(0), data_() {}

        scalarMatrix(const label nRows, const label nCols)
        :
            nRows_(nRows),
            nCols_(nCols),
            data_(nRows*nCols, Zero)
        {}

        inline scalar& operator()(const label i, const label j)
        {
            return data_[i*nCols_ + j];
        }

        inline const scalar& operator()(const label i, const label j) const
        {
            return data_[i*nCols_ + j];
        }
    };

private:
    fileName fileName_;

    label n_;
    label m_;

    labelList gridSizeX_;
    labelList gridSizeY_;
    scalarField boundaryX_;
    scalarField boundaryY_;
    Transformation transformationX_;
    Transformation transformationY_;

    scalarField x_;
    scalarField y_;
    scalarField dx_;
    scalarField dy_;
    scalarField sectionStepX_;
    scalarField sectionStepY_;

    List<FixedList<scalar, 16>> coefficients_;

    void readTable();

    static Transformation readTransformation(const word& name);
    static word transformationName(Transformation transformation);

    void readSettings(const dictionary& dict);
    void createGrid();

    void calcCoeffs(const std::function<scalar(scalar, scalar)>& f);

    scalarField solveSplineDerivatives
    (
        const scalarField& values,
        const scalarField& intervals,
        scalar leftDerivative,
        scalar rightDerivative
    ) const;

    scalarField solveTridiagonal
    (
        const scalarField& lower,
        const scalarField& diagonal,
        const scalarField& upper,
        const scalarField& rhs
    ) const;

    scalar transform(scalar x, Transformation transformation) const;
    scalar backTransform(scalar x, Transformation transformation) const;

    std::pair<int, int> fastFindPosition(scalar xx, scalar yy) const;
    scalar calc(scalar xx, scalar yy) const;

public:
    explicit biCubicInterpolation(const fileName& fName);

    biCubicInterpolation
    (
        const dictionary& dict,
        const std::function<scalar(scalar, scalar)>& f
    );

    void check() const;

    void write(Ostream& os) const;
    void write(const fileName& fName) const;

    std::pair<int, int> findPosition(scalar xx, scalar yy) const;

    scalar operator()(scalar valueX, scalar valueY) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
