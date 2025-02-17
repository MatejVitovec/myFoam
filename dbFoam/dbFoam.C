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

#include "twoFluid.H"
#include "twoFluidConvectiveFlux.H"

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

        surfaceScalarField amaxSf
        (
            "amaxSf", 
            max(mag(fvc::interpolate(U1) & mesh.Sf()), mag(fvc::interpolate(U2) & mesh.Sf()))
            + max(mesh.magSf() * fvc::interpolate(sqrt(thermo1.Cp()/thermo1.Cv()/thermo1.psi())),
                  mesh.magSf() * fvc::interpolate(sqrt(thermo2.Cp()/thermo2.Cv()/thermo2.psi()))) //approx sound speed
        );

        #include "courantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
       
        twoFluidFlux.computeFlux();

        // --- Solve density
        solve(fvm::ddt(conservative.alphaRho1()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg()));
        solve(fvm::ddt(conservative.alphaRho2()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg()));

        // --- Solve momentum
        solve(fvm::ddt(conservative.alphaRhoU1()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
            + fluidSystem.pInt()*twoFluid::fvc::div(twoFluidFlux.alpha_pos()*mesh_.Sf(), twoFluidFlux.alpha_neg()*mesh_.Sf()));
        solve(fvm::ddt(conservative.alphaRhoU2()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
            + fluidSystem.pInt()*twoFluid::fvc::div((1.0 - twoFluidFlux.alpha_pos())*mesh_.Sf(), (1.0 - twoFluidFlux.alpha_neg())*mesh_.Sf()));

        // --- Solve energy
        solve(fvm::ddt(conservative.alphaRhoE1()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg()));
        solve(fvm::ddt(conservative.alphaRhoE2()) + twoFluid::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg()));

        fluidSystem.correct();
        fluidSystem.blendVanishingFluid();
        fluidSystem.correctBoundaryCondition();
        fluidSystem.correctThermo();
        fluidSystem.correctInterfacialPressure();
        fluidSystem.correctConservative();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
