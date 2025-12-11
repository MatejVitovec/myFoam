/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 22 "/home/matej/myFoam/tutorials/waterCollum/0/alpha/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_efce25ed6e0bb9c813fe3bebc892f53ae1114b6e(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 39 "/home/matej/myFoam/tutorials/waterCollum/0/alpha/#codeStream"
const scalar epsilon = 1e-5;

        const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());
        scalarField alpha(mesh.nCells(), epsilon);

        const vector center(0, 0, 0);
        const scalar R = 0.0032;
        const scalar delta = 0.000025;

        forAll(mesh.C(), cellI)
        {
            vector pos = mesh.C()[cellI];
            scalar r = mag(pos - center);

            if (R - 2.0*delta <= r && r <= R + 2.0*delta)
            {
                scalar xi = (r - (R - 2.0*delta))/(4.0*delta);
                scalar G = -Foam::pow(xi, 2.0)*(2.0*xi - 3.0);
                alpha[cellI] = G*(epsilon) + (1.0 - G)*(1.0 - epsilon);
            }
            else if(r < R)
            {
                alpha[cellI] = 1.0 - epsilon;
            }
        }

        alpha.writeEntry("", os);
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

