{
    const volScalarField& alphaRho1  = fluid.conservative().alphaRho1();
    const volScalarField& alphaRho2  = fluid.conservative().alphaRho2();
    const volVectorField& alphaRhoU1 = fluid.conservative().alphaRhoU1();
    const volVectorField& alphaRhoU2 = fluid.conservative().alphaRhoU2();
    const volScalarField& epsilon1   = fluid.conservative().epsilon1();
    const volScalarField& epsilon2   = fluid.conservative().epsilon2();

    const volScalarField& a1 = fluid.a1();
    const volScalarField& a2 = fluid.a2();

    Field<scalar> D(mesh.V()/dt);
    /*if (LTS)
    {
        D *= dt * trDeltaTau();
    }*/
    
    /*const volScalarField nuMax( "nuMax", max(
        4./3.0 * turbulence->nuEff(),
        thermo.Cp()/thermo.Cv() * turbulence->alphaEff()/rho
    ));*/


    Field<scalar> dAlphaRho1Star (dAlphaRho1);
    Field<scalar> dAlphaRho2Star (dAlphaRho2);
    Field<vector> dAlphaRhoU1Star(dAlphaRhoU1);
    Field<vector> dAlphaRhoU2Star(dAlphaRhoU2);
    Field<scalar> dEpsilon1Star  (dEpsilon1);
    Field<scalar> dEpsilon2Star  (dEpsilon2);
    
    #include "forwardSweep.H"
    #include "backwardSweep.H"
}
