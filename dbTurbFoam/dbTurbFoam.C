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
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "convectiveFlux.H"

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

        amaxSf = mag(fvc::interpolate(U) & mesh.Sf()) + mesh.magSf()*fvc::interpolate(a);
        #include "courantNo.H"
        #include "setDeltaT.H"
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        dimensionedScalar dt = runTime.deltaT();

        scalar coeff[3][3] = {{1.0, 0.0 , 1.0}, {3.0/4.0, 1.0/4.0, 1.0/4.0}, {1.0/3.0, 2.0/3.0, 2.0/3.0}};

        for (size_t i = 0; i < 3; i++)
        {
            flux.computeFlux();

            rho = coeff[i][0]*rho.oldTime()
                + coeff[i][1]*rho
                - coeff[i][2]*dt*(fvc::div(flux.rhoFlux()));

            U = coeff[i][0]*U.oldTime()
              + coeff[i][1]*U
              - coeff[i][2]*dt*fvc::div(flux.rhoUFlux());
            
            rhoE = coeff[i][0]*rhoE.oldTime()
                 + coeff[i][1]*rhoE
                 - coeff[i][2]*dt*fvc::div(flux.rhoEFlux());
            
            #include "updateFields.H"
        }

        runTime.write();

        //runTime.printExecutionTime(Info);
        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
