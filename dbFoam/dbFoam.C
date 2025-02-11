/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    dbFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "twoFluidConvectiveFlux.H"
#include "twoFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " Riemann problem flux calculation"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        surfaceScalarField amaxSf("amaxSf", 
            mag(fvc::interpolate(U) & mesh.Sf()) +
            mesh.magSf() * fvc::interpolate(sqrt(thermo.Cp()/thermo.Cv()/thermo.psi())));

        #include "courantNo.H"
        #include "setDeltaT.H"

        //if (LTS)
        //{
        //    #include "setRDeltaT.H"
        //}

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
       
        dbFlux.computeFlux();

        // --- Solve density
        solve(fvm::ddt(alphaRho1) + twoFluid::fvc::div(dbFlux.alphaRhoFlux1_pos(), dbFlux.alphaRhoFlux1_neg()));
        solve(fvm::ddt(alphaRho2) + twoFluid::fvc::div(dbFlux.alphaRhoFlux2_pos(), dbFlux.alphaRhoFlux2_neg()));

        // --- Solve momentum
        solve(fvm::ddt(alphaRhoU1) + twoFluid::fvc::div(dbFlux.alphaRhoUFlux1_pos(), dbFlux.alphaRhoUFlux1_neg()));
        solve(fvm::ddt(alphaRhoU2) + twoFluid::fvc::div(dbFlux.alphaRhoUFlux2_pos(), dbFlux.alphaRhoUFlux2_neg()));

        // --- Solve energy
        solve(fvm::ddt(alphaRhoE1) + twoFluid::fvc::div(dbFlux.alphaRhoEFlux1_pos(), dbFlux.alphaRhoEFlux1_neg()));
        solve(fvm::ddt(alphaRhoE2) + twoFluid::fvc::div(dbFlux.alphaRhoEFlux2_pos(), dbFlux.alphaRhoEFlux2_neg()));


        #include "updateFields.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
