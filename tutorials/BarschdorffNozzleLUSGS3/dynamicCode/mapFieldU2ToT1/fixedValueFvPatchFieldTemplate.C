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
// SHA1 = bf766773243118fdfeb54c6f9d87d057e8e10002
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void mapFieldU2ToT1_bf766773243118fdfeb54c6f9d87d057e8e10002(bool load)
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
    fvPatchVectorField,
    mapFieldU2ToT1FixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
mapFieldU2ToT1FixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct mapFieldU2ToT1 : patch/DimensionedField");
    }
}


Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
mapFieldU2ToT1FixedValueFvPatchVectorField
(
    const mapFieldU2ToT1FixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct mapFieldU2ToT1 : patch/DimensionedField/mapper");
    }
}


Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
mapFieldU2ToT1FixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct mapFieldU2ToT1 : patch/dictionary");
    }
}


Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
mapFieldU2ToT1FixedValueFvPatchVectorField
(
    const mapFieldU2ToT1FixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct mapFieldU2ToT1");
    }
}


Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
mapFieldU2ToT1FixedValueFvPatchVectorField
(
    const mapFieldU2ToT1FixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct mapFieldU2ToT1 : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::
~mapFieldU2ToT1FixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy mapFieldU2ToT1");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
mapFieldU2ToT1FixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs mapFieldU2ToT1");
    }

//{{{ begin code
    #line 32 "/home/matejv/myFoam/tutorials/BarschdorffNozzleLUSGS3/0/U.2/boundaryField/inlet"
const volVectorField& sourceField = db().lookupObject<volVectorField>("U.1");
            operator== (sourceField.boundaryField()[patch().index()]);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

