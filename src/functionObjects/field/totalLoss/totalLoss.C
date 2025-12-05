/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "totalLoss.H"
#include "fluidThermo.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "Pstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalLoss, 0);
    addToRunTimeSelectionTable(functionObject, totalLoss, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::totalLoss::totalLoss
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    xMin_(readScalar(dict.lookup("xMin"))),
    xMax_(readScalar(dict.lookup("xMax"))),
    Npoints_(readLabel(dict.lookup("Npoints"))),
    sigma_(readScalar(dict.lookup("sigma")))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::totalLoss::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::totalLoss::execute()
{
    return true;
}


bool Foam::functionObjects::totalLoss::write()
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volScalarField& T1 = mesh_.lookupObject<volScalarField>("T.1");
    const volScalarField& T2 = mesh_.lookupObject<volScalarField>("T.2");

    const fluidThermo& thermo1 =
        lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", "1"));
    
    const fluidThermo& thermo2 =
        lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", "2"));

    const volScalarField h1 = thermo1.he() + p/thermo1.rho();
    const volScalarField h2 = thermo2.he() + p/thermo2.rho();

    const volScalarField& mDot = mesh_.lookupObject<volScalarField>("mDotCondensation");

    const volScalarField zRelax
    (
        IOobject
        (
            scopedName("zRelax"),
            time_.timeName(),
            mesh_
        ),
        mDot*(h1 - h2)*(1.0/T1 - 1.0/T2)*(T1 + T2)/2.0
    );

    //Log << "    Relaxation losses field " << zRelax.name()
    //    << " to " << time_.timeName() << ", total losses: " << fvc::domainIntegrate(zRelax) << endl;


    List<scalar> cumulativeLoses = computeIntegral(zRelax);
    List<scalar> globalCumulativeLoses = cumulativeLoses;

    for (label i = 0; i < globalCumulativeLoses.size(); ++i)
    {
        scalar v = globalCumulativeLoses[i];
        reduce(v, sumOp<scalar>());
        globalCumulativeLoses[i] = v;
    }

    if (Pstream::master())
    {
        const scalar dx = (xMax_ - xMin_)/(Npoints_ - 1);

        OFstream out(fileName(mesh_.runTime().path()/"postProcessing"/"relaxationLoss"/mesh_.runTime().timeName()));
        out << "# X  Integral\n";

        for (label k = 0; k < Npoints_; k++)
        {
            out << xMin_ + dx*scalar(k) << " " << globalCumulativeLoses[k] << "\n";
        }
    }

    return true;
}


Foam::List<Foam::scalar> Foam::functionObjects::totalLoss::computeIntegral(const Foam::volScalarField& field) const
{
    Info << "GaussianVolumeIntegral: computing smooth cumulative integral…" << endl;

    const vectorField& C = mesh_.C(); // cell centres
    const scalarField& V = mesh_.V(); // cell volumes

    const scalar sqrt2sigma = Foam::sqrt(2.0)*sigma_;

    List<scalar> out(Npoints_, 0.0);

    for (label k = 0; k < Npoints_; k++)
    {
        const scalar X = xMin_ + (xMax_ - xMin_)*(scalar(k)/(Npoints_-1));
        scalar Ik = 0.0;

        forAll(C, i)
        {
            scalar dx = (X - C[i].x())/sqrt2sigma;

            // Smooth heaviside: 0.5*(1+erf())
            scalar w = 0.5*(1.0 + erf(dx));

            Ik += field[i]*V[i]*w;
        }
        reduce(Ik, sumOp<scalar>());

        out[k] = Ik;
    }

    return out;
}


// ************************************************************************* //

