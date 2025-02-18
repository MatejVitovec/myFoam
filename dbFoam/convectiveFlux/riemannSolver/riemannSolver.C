/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    riemannSolver

Description
    Basic class for of inviscid numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    riemannSolver.C

\*---------------------------------------------------------------------------*/

#include "riemannSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*namespace Foam
{
  defineTypeNameAndDebug(riemannSolver, 0);
  defineRunTimeSelectionTable(riemannSolver,dictionary);
}*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::riemannSolver> Foam::riemannSolver::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word key
)
{
    const word fluxName( dict.lookup(key) );

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(fluxName);

    if (!cstrIter.found())
    {
        FatalErrorIn
        (
            "riemannSolver::New(mesh, dict, key)"
        )   << "Unknown riemannSolver type "
            << fluxName << nl << nl
            << "Valid riemannSolver types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<riemannSolver>(cstrIter()(mesh, dict));
}*/


// ************************************************************************* //