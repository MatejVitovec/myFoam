/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "transonicPressureOutletFvPatchScalarField.H"
#include "fluidThermo.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::transonicPressureOutletFvPatchScalarField<Type>::transonicPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    fluidName_(mesh.foundObject<volVectorField>("U.1") ? "1" : word::null),
    pName_("p"),
    TName_(IOobject::groupName("T", fluidName_)),
    UName_(IOobject::groupName("U", fluidName_))
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::transonicPressureOutletFvPatchScalarField<Type>::transonicPressureOutletFvPatchScalarField
(
    const transonicPressureOutletFvPatchScalarField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    fluidName_(ptf.fluidName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    UName_(ptf.UName_)
{}


template<class Type>
Foam::transonicPressureOutletFvPatchScalarField<Type>::transonicPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    fluidName_(mesh.foundObject<volVectorField>("U.1") ? "1" : word::null),
    pName_("p"),
    UName_(IOobject::groupName("U", fluidName_)),
    phiName_(IOobject::groupName("phi", fluidName_)),
{
    fvPatchFieldBase::readDict(dict);

    // Require inletValue (MUST_READ)
    this->refValue().assign("inletValue", dict, p.size());
    this->refGrad() = Zero;
    this->valueFraction() = 0;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<Type>::extrapolateInternal();
    }
}


template<class Type>
Foam::transonicPressureOutletFvPatchScalarField<Type>::transonicPressureOutletFvPatchScalarField
(
    const transonicPressureOutletFvPatchScalarField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    fluidName_(ptf.fluidName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    UName_(ptf.UName_)
{}


template<class Type>
Foam::transonicPressureOutletFvPatchScalarField<Type>::transonicPressureOutletFvPatchScalarField
(
    const transonicPressureOutletFvPatchScalarField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    fluidName_(ptf.fluidName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::transonicPressureOutletFvPatchScalarField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& p = 
        db().lookupObject<volScalarField>(pName_);

    const volScalarField& T = 
        db().lookupObject<volScalarField>(TName_);

    const volVectorField& U = 
        db().lookupObject<volScalarField>(UName_);

    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", fluidName_));

    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));

    scalarField& fildp = *this;

    scalarField Ma = scalarField(fildp.size(), 0.0);

    forAll(Ma, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];
        Ma[faceI] = U[faceCellI]/gasProps->c(p[faceCellI], T[faceCellI]);
    }
    
    //const Field<scalar>& phip =
    //    this->patch().template lookupPatchField<surfaceScalarField>(phiName_);

    this->valueFraction() = pos(1 - Ma);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::transonicPressureOutletFvPatchScalarField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("fluidName", "", fluidName_);
    os.writeEntryIfDifferent<word>(IOobject::groupName("p", fluidName_), "p", pName_);
    os.writeEntryIfDifferent<word>(IOobject::groupName("U", fluidName_), "U", UName_);
    os.writeEntryIfDifferent<word>(IOobject::groupName("phi", fluidName_), "phi", phiName_);
    this->refValue().writeEntry("outletValue", os);
    fvPatchField<Type>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::transonicPressureOutletFvPatchScalarField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //

