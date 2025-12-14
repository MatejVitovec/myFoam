/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::transonicPressureOutletFvPatchScalarField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a generic outflow condition, with
    specified inflow for the case of return flow.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Flux field name         | no          | phi
        inletValue   | Inlet value for reverse flow | yes    |
        value        | initial field value          | optional |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            inletOutlet;
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

Note
    Sign conventions:
    - Positive flux (out of domain): apply zero-gradient condition
    - Negative flux (into domain): apply the "inletValue" fixed-value

See also
    Foam::mixedFvPatchField
    Foam::zeroGradientFvPatchField
    Foam::outletInletFvPatchField

SourceFiles
    transonicPressureOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_transonicPressureOutletFvPatchScalarField_H
#define Foam_transonicPressureOutletFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class transonicPressureOutletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class transonicPressureOutletFvPatchScalarField
:
    public mixedFvPatchField<Type>
{
protected:

    // Protected Data

        //- Fluid name
        word fluidName_;

        //- Pressure field name
        word pName_;

        //- Temperature field name
        word TName_;

        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("transonicPressureOutlet");


    // Constructors

        //- Construct from patch and internal field
        transonicPressureOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transonicPressureOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given transonicPressureOutletFvPatchScalarField onto a new patch
        transonicPressureOutletFvPatchScalarField
        (
            const transonicPressureOutletFvPatchScalarField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        transonicPressureOutletFvPatchScalarField
        (
            const transonicPressureOutletFvPatchScalarField<Type>&
        );

        //- Construct as copy setting internal field reference
        transonicPressureOutletFvPatchScalarField
        (
            const transonicPressureOutletFvPatchScalarField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        //- True: this patch field is altered by assignment.
        virtual bool assignable() const { return true; }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<Type>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transonicPressureOutletFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
