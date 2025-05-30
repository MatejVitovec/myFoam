/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 7b45920900006f2e0af41ee33a0c0d75de2ecece
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void mapFieldT1ToT2_7b45920900006f2e0af41ee33a0c0d75de2ecece(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    mapFieldT1ToT2FixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
mapFieldT1ToT2FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct mapFieldT1ToT2 : patch/DimensionedField");
    }
}


Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
mapFieldT1ToT2FixedValueFvPatchScalarField
(
    const mapFieldT1ToT2FixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct mapFieldT1ToT2 : patch/DimensionedField/mapper");
    }
}


Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
mapFieldT1ToT2FixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct mapFieldT1ToT2 : patch/dictionary");
    }
}


Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
mapFieldT1ToT2FixedValueFvPatchScalarField
(
    const mapFieldT1ToT2FixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct mapFieldT1ToT2");
    }
}


Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
mapFieldT1ToT2FixedValueFvPatchScalarField
(
    const mapFieldT1ToT2FixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct mapFieldT1ToT2 : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::
~mapFieldT1ToT2FixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy mapFieldT1ToT2");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
mapFieldT1ToT2FixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs mapFieldT1ToT2");
    }

//{{{ begin code
    #line 39 "/home/matej/myFoam/tutorials/SE1050condensation/0/T.2/boundaryField/inlet"
const volScalarField& sourceField = db().lookupObject<volScalarField>("T.1");
            operator== (sourceField.boundaryField()[patch().index()]);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

