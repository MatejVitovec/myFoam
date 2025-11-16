/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "saturation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

/*namespace Foam
{
    defineTypeNameAndDebug(saturation, 0);
    defineRunTimeSelectionTable(saturation, dictionary);

}*/

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

saturation::saturation
(
    const dictionary& dict,
    const fluidThermo& thermo
)
:
    thermo_(thermo),
    pSaturation_(saturationCurve::New(dict)),
    saturationCurve_(pSaturation_()),
    ps_(
        IOobject
        (
            "ps",
            thermo.p().mesh().time().timeName(),
            thermo.p().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.p().mesh(),
        dimensionedScalar("initSaturPressure", dimPressure, 0.0)
    ),
    Ts_(
        IOobject
        (
            "Ts",
            thermo.p().mesh().time().timeName(),
            thermo.p().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.p().mesh(),
        dimensionedScalar("initSaturTemperature", dimTemperature, 0.0)
    ),
    L_(
        IOobject
        (
            "latentHeat",
            thermo.p().mesh().time().timeName(),
            thermo.p().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.p().mesh(),
        dimensionedScalar("initLatentHeat", dimEnergy/dimMass, 0.0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

/*autoPtr<saturation>
saturation::New
(
    const dictionary& dict
)
{
    word name
    (
        dict.lookup("saturation")
    );

    Info<< "Selecting saturation curve "
        << name << endl;

    auto cstrIter = dictConstructorTablePtr_->cfind(name);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown saturation "
            << name << endl << endl
            << "Valid saturation curves are : " << endl
            << dictConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<saturation>(cstrIter()(dict));
}*/

const volScalarField& saturation::ps() const
{
    return ps_;
}

const volScalarField& saturation::Ts() const
{
    return Ts_;
}

const volScalarField& saturation::L() const
{
    return L_;
}

tmp<volScalarField> saturation::dpsdT() const
{
    return saturationCurve_.dpsdT(Ts_);
}

tmp<volScalarField> saturation::rhosl() const
{
    return saturationCurve_.rhosl(Ts_);
}

tmp<volScalarField> saturation::rhosv() const
{
    return saturationCurve_.rhosv(Ts_);
}

tmp<volScalarField> saturation::hsl() const
{
    return saturationCurve_.hsl(Ts_);
}

tmp<volScalarField> saturation::hsv() const
{
    return saturationCurve_.hsv(Ts_);
}

tmp<volScalarField> saturation::esl() const
{
    return saturationCurve_.esl(Ts_);
}

tmp<volScalarField> saturation::esv() const
{
    return saturationCurve_.esv(Ts_);
}

void saturation::correct()
{
    ps_ = saturationCurve_.ps(thermo_.T());
    Ts_ = saturationCurve_.Ts(thermo_.p());
    //L_  = saturationCurve_.L(Ts_);

    //L_ = hsv() - hsl();
    L_ = (thermo_.he() + thermo_.p()/thermo_.rho()) - hsl();

    const auto& mesh = thermo_.p().mesh();
    scalar averageValue = gSum(mesh.V()*(hsv() - hsl() - L_))/gSum(mesh.V());

    Info << "latent heat mean diff: " << averageValue << endl;
}

}

