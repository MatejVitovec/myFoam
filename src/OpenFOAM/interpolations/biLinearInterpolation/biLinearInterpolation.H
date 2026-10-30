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
    Foam::biLinearInterpolation

Description
    Bilinear interpolation on a rectangular grid.

SourceFiles
    biLinearInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef biLinearInterpolation_H
#define biLinearInterpolation_H

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

class biLinearInterpolation
{
public:
    enum Transformation { NONE, LOG, LOG10, LOGINV };

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

    List<FixedList<scalar, 4>> coefficients_;

    void readTable();

    static Transformation readTransformation(const word& name);
    static word transformationName(Transformation transformation);

    void readSettings(const dictionary& dict);
    void createGrid();

    void calcCoeffs(const std::function<scalar(scalar, scalar)>& f);

    scalar transform(scalar x, Transformation transformation) const;
    scalar backTransform(scalar x, Transformation transformation) const;

    std::pair<int, int> fastFindPosition(scalar xx, scalar yy) const;
    scalar calc(scalar xx, scalar yy) const;

public:
    explicit biLinearInterpolation(const fileName& fName);

    biLinearInterpolation
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
