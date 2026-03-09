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
#include "fvcSmooth.H"

#include "twoFluid.H"
#include "twoFluidConvectiveFlux.H"
#include "twoFluidFvc.H"
#include "dragModel.H"

#include <eigen3/Eigen/Dense>

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

        amaxSf = fluid.amaxSf();
        #include "courantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Iter = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << endl;
        
        dimensionedScalar dt = runTime.deltaT();

        scalar coeff[3][3] = {{1.0, 0.0 , 1.0}, {3.0/4.0, 1.0/4.0, 1.0/4.0}, {1.0/3.0, 2.0/3.0, 2.0/3.0}};

        for (size_t i = 0; i < 3; i++)
        {
            twoFluidFlux.computeFlux();

            volScalarField dragK = drag.K(d);

            volScalarField rezAlphaRho1(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg())
            ));

            volScalarField rezAlphaRho2(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg())
            ));

            volVectorField rezAlphaRhoU1(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha1_pos()*mesh.Sf(), twoFluidFlux.alpha1_neg()*mesh.Sf())
                + dragK*(fluid.U1() - fluid.U2())
            ));

            volVectorField rezAlphaRhoU2(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                + dragK*(fluid.U2() - fluid.U1())
            ));

            volScalarField rezEpsilon1(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                + ((dragK*(fluid.U1() - fluid.U2())) & U1)
            ));

            volScalarField rezEpsilon2(-(
                TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                + ((dragK*(fluid.U2() - fluid.U1())) & U2)
            ));

            dp     = dimensionedScalar("pzero", dimPressure, 0.0);
            dalpha = dimensionedScalar("alphazero", dimless, 0.0);
            dU1    = dimensionedVector("U1zero", dimVelocity, vector(0, 0, 0));
            dU2    = dimensionedVector("U2zero", dimVelocity, vector(0, 0, 0));
            dT1    = dimensionedScalar("T1zero", dimTemperature, 0.0);
            dT2    = dimensionedScalar("T2zero", dimTemperature, 0.0);

            #include "calculateDeltaPrimitive.H"

            p     = coeff[i][0]*p.oldTime()     + coeff[i][1]*p     + coeff[i][2]*dp;
            alpha = coeff[i][0]*alpha.oldTime() + coeff[i][1]*alpha + coeff[i][2]*dalpha;
            U1    = coeff[i][0]*U1.oldTime()    + coeff[i][1]*U1    + coeff[i][2]*dU1;
            U2    = coeff[i][0]*U2.oldTime()    + coeff[i][1]*U2    + coeff[i][2]*dU2;
            T1    = coeff[i][0]*T1.oldTime()    + coeff[i][1]*T1    + coeff[i][2]*dT1;
            T2    = coeff[i][0]*T2.oldTime()    + coeff[i][1]*T2    + coeff[i][2]*dT2;

            fluid.blendVanishingFluid();
            fluid.correctBoundaryCondition();
            fluid.correctThermo();
            if(i == 2) fluid.correctInterfacialPressure();
            //fluid.correctConservative();
        }


        rho1.ref() = thermo1.rho();
        rho2.ref() = thermo2.rho();

        runTime.write();

        /*if (runTime.timeIndex() > 6930)
        {
            Info << ">>> Forcing write <<<" << endl;
            runTime.writeNow();
            runTime.write();
        }*/

        //runTime.printExecutionTime(Info);

        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
