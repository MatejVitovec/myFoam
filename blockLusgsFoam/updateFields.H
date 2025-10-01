{
    // Correct boundary conditions
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    T.correctBoundaryConditions();

    // Correct thermodynamics
    e = thermo.he(p, T);
    thermo.correct();

    rho = thermo.rho();

    forAll(a, cell)
    {
        a[cell] = gasProps.c(p[cell], T[cell]);
    }      
}
