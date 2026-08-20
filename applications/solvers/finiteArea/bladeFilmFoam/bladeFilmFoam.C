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
    lusgsTwoFluidFoam

Group
    grpCompressibleSolvers

Description
    Density-based lusgs two-fluid compressible flow solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "localEulerDdtScheme.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "fvcSmooth.H"

#include "twoFluid.H"
#include "twoFluidConvectiveFlux.H"
#include "twoFluidFvc.H"
#include "dragModel.H"

#include <eigen3/Eigen/Dense>
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"

#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Blade surface film solver"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "../../../tools/film/createFilmFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;


    Info<< "\nStarting time loop\n" << endl;

    twoFluidFlux.computeFlux();

    while (runTime.run())
    {
        #include "readTimeControls.H"

        amaxSf = fluid.amaxSf();
        #include "courantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << nl << "Iter = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << endl;
        
        dimensionedScalar dt = runTime.deltaT();

        #include "../../../tools/film/solveFilm.H"

        h[270] = 0.0;
        h[271] = 0.0;

        runTime.write();

        //runTime.printExecutionTime(Info);
        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
