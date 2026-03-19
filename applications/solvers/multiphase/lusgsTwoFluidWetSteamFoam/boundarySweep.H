{
    pStar.internalFieldRef()       = p.internalField() + dpStar;
    alphaStar.internalFieldRef()   = alpha.internalField() + dalphaStar;
    U1Star.internalFieldRef()      = U1.internalField() + dU1Star;
    U2Star.internalFieldRef()      = U2.internalField() + dU2Star;
    T1Star.internalFieldRef()      = T1.internalField() + dT1Star;
    T2Star.internalFieldRef()      = T2.internalField() + dT2Star;

    pStar.correctBoundaryConditions();
    alphaStar.correctBoundaryConditions();
    U1Star.correctBoundaryConditions();
    U2Star.correctBoundaryConditions();
    T1Star.correctBoundaryConditions();
    T2Star.correctBoundaryConditions();
}