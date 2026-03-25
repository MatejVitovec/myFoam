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
// SHA1 = a7031561f0ee20f3ae7d277c310071fb5821f8ba
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void switchAlphaInlet_a7031561f0ee20f3ae7d277c310071fb5821f8ba(bool load)
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
    switchAlphaInletFixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
switchAlphaInletFixedValueFvPatchScalarField::
switchAlphaInletFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct switchAlphaInlet : patch/DimensionedField");
    }
}


Foam::
switchAlphaInletFixedValueFvPatchScalarField::
switchAlphaInletFixedValueFvPatchScalarField
(
    const switchAlphaInletFixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct switchAlphaInlet : patch/DimensionedField/mapper");
    }
}


Foam::
switchAlphaInletFixedValueFvPatchScalarField::
switchAlphaInletFixedValueFvPatchScalarField
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
        printMessage("Construct switchAlphaInlet : patch/dictionary");
    }
}


Foam::
switchAlphaInletFixedValueFvPatchScalarField::
switchAlphaInletFixedValueFvPatchScalarField
(
    const switchAlphaInletFixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct switchAlphaInlet");
    }
}


Foam::
switchAlphaInletFixedValueFvPatchScalarField::
switchAlphaInletFixedValueFvPatchScalarField
(
    const switchAlphaInletFixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct switchAlphaInlet : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
switchAlphaInletFixedValueFvPatchScalarField::
~switchAlphaInletFixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy switchAlphaInlet");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
switchAlphaInletFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs switchAlphaInlet");
    }

//{{{ begin code
    #line 36 "/home/matejv/myFoam/tutorials/SE1050BL2/0/alpha/boundaryField/inlet"
const scalar lowValue = 1e-10;
            const scalar highValue = 1e-5;
            const label switchIter = 4000;

            label iter = this->db().time().timeIndex();

            scalar val = (iter < switchIter) ? lowValue : highValue;

            operator==(val);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

