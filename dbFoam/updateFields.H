

const volScalarField gamma = thermo.Cp()/thermo.Cv();

rho.boundaryFieldRef() = thermo.rho()().boundaryField();

U = rhoU/rho;
U.correctBoundaryConditions();

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();

p = (gamma - 1.0)*rho*e;
p.correctBoundaryConditions();


thermo.correct();

rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
rhoE.boundaryFieldRef() ==
    rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));
