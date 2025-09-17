{
    U = rhoU/rho;

    e = rhoE/rho - 0.5*magSqr(U);
    
    forAll(p, cell)
    {
        p[cell] = gasProps.pEs(e[cell], rho[cell], p[cell]);
    }    

    //rho.correctBoundaryConditions();

    // Correct boundary conditions
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    e.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();

    // Update boundary field of rho
    //rho.boundaryFieldRef() = thermo.rho()().boundaryField();
    rho.boundaryFieldRef() = thermo.psi().boundaryField()*p.boundaryField();

    // Update boundary field of rhoU 
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() = rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

    forAll(a, cell)
    {
        a[cell] = gasProps.c(p[cell], T[cell]);
    }      
}
