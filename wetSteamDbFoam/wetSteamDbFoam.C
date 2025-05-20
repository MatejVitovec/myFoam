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
#include "twoFluidFvc.H"
#include "dragModel.H"
#include "saturationCurve.H"
#include "condensationModel.H"

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

        amaxSf = fluidSystem.amaxSf();
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

            volVectorField Uint = (alpha1*rho1*U1 + alpha2*rho2*U2)/(alpha1*rho1 + alpha2*rho2);

            volScalarField Ts = saturation.Ts(p);
            volScalarField Hvint = saturation.hsv(Ts) + (Uint & U1) + 0.5*magSqr(U1);
            volScalarField Hlint = saturation.hsl(Ts) + (Uint & U2) + 0.5*magSqr(U2);

            volVectorField dragSource = drag.K(condensation.dropletDiameter())*(fluidSystem.U1() - fluidSystem.U2());

            volScalarField condensationMassSource = condensation.condensationRateMassSource();
            volVectorField condensationMomentumSource = Uint*condensation.condensationRateMassSource();
            volScalarField condensationEnergyVaporSource = condensation.condensationRateMassSource()*(Hvint - condensation.L());
            volScalarField condensationEnergyLiquidSource = condensation.condensationRateMassSource()*Hlint;

            conservative.alphaRho1() = coeff[i][0]*conservative.alphaRho1().oldTime()
                                     + coeff[i][1]*conservative.alphaRho1()
                                     - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg())
                                        + condensationMassSource);
            
            conservative.alphaRho2() = coeff[i][0]*conservative.alphaRho2().oldTime()
                                     + coeff[i][1]*conservative.alphaRho2()
                                     - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg())
                                        - condensationMassSource);

            conservative.alphaRhoU1() = coeff[i][0]*conservative.alphaRhoU1().oldTime()
                                      + coeff[i][1]*conservative.alphaRhoU1()
                                      - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                                          - fluidSystem.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha1_pos()*mesh.Sf(), twoFluidFlux.alpha1_neg()*mesh.Sf())
                                          + dragSource + condensationMomentumSource);
          
            conservative.alphaRhoU2() = coeff[i][0]*conservative.alphaRhoU2().oldTime()
                                      + coeff[i][1]*conservative.alphaRhoU2()
                                      - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                                          - fluidSystem.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                                          - dragSource - condensationMomentumSource);

            conservative.epsilon1() = coeff[i][0]*conservative.epsilon1().oldTime()
                                    + coeff[i][1]*conservative.epsilon1()
                                    - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                                        + (dragSource & Uint) + condensationEnergyVaporSource);
            
            conservative.epsilon2() = coeff[i][0]*conservative.epsilon2().oldTime()
                                    + coeff[i][1]*conservative.epsilon2()
                                    - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                                        - (dragSource & Uint) - condensationEnergyLiquidSource);

            fluidSystem.correct();

            fluidSystem.blendVanishingFluid();
            fluidSystem.correctBoundaryCondition();
            fluidSystem.correctThermo();
            if(i == 2) fluidSystem.correctInterfacialPressure();
            fluidSystem.correctConservative();
        }

        //alphaRhoPhi2 = fvc::interpolate(thermo2.rho())*fvc::interpolate(alpha2)*fvc::flux(U2);
        twoFluidFlux.computeFlux();
        alphaRhoPhi2 = twoFluidFlux.alphaRhoFlux2_pos();
        condensation.correct();
        //Info << "condensation complete " << endl;

        /*if (runTime.timeIndex() > 4130)
        {
            Info << ">>> Forcing write <<<" << endl;
            runTime.writeNow();
            runTime.write();
        }*/

        Info <<endl;


        rho1.ref() = thermo1.rho();
        rho2.ref() = thermo2.rho();
        e2.ref() = thermo2.he();

        runTime.write();

        //runTime.printExecutionTime(Info);

        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
