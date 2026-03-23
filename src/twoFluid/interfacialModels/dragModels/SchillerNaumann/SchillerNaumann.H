/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::SchillerNaumann

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    SchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumann_H
#define SchillerNaumann_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{

class twoFluid;


/*---------------------------------------------------------------------------*\
                           Class SchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class SchillerNaumann
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        // -Reynolds number
        volScalarField Re_;

        //- Drag coefficient
        tmp<volScalarField> CdRe() const;

        scalar dCdRedx(const scalar Re, const scalar dRedx) const;
        vector dCdRedx(const scalar Re, const vector dRedx) const;

public:

    //- Runtime type information
    TypeName("SchillerNaumann");


    // Constructors

        //- Construct from a dictionary and a twoFluid
        SchillerNaumann
        (
            const dictionary& dict,
            const twoFluid& fluid,
            const volScalarField& dropletDiameter,
            const bool registerObject
        );


    //- Destructor
    virtual ~SchillerNaumann();


    // Member Functions

        virtual scalar dKidp(const label celli) const;
        virtual scalar dKidalpha(const label celli) const;
        virtual vector dKidU1(const label celli) const;
        virtual vector dKidU2(const label celli) const;
        virtual scalar dKidT1(const label celli) const;
        virtual scalar dKidT2(const label celli) const;

        virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
