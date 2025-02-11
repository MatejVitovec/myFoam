


twoFluid::primitiveFromConservative(p, alpha, U1, U2, T1, T2, alphaRho1, alphaRho2, alphaRhoU1, alphaRhoU2, alphaRhoE1, alphaRhoE2, thermo1, thermo2, mesh);

twoFluid::blend(p, alpha, U1, U2, T1, T2);

/*const volScalarField gamma = thermo.Cp()/thermo.Cv();

rho.boundaryFieldRef() = thermo.rho()().boundaryField();

U = rhoU/rho;
U.correctBoundaryConditions();

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();

p = (gamma - 1.0)*rho*e;
p.correctBoundaryConditions();*/

thermo1.he() = thermo1.he(p, T1);
thermo2.he() = thermo2.he(p, T2);

thermo1.psi() = thermo1.psi(p, T1);
thermo2.psi() = thermo2.psi(p, T2);

//TODO mu + alpha(thermal)

p.correctBoundaryConditions();
alpha.correctBoundaryConditions();
U1.correctBoundaryConditions();
U2.correctBoundaryConditions();
T1.correctBoundaryConditions();
T2.correctBoundaryConditions();

//thermo1.correct();
//thermo2.correct();

const auto& e1 = thermo1.he();
const auto& e2 = thermo2.he();


//TODO - pridat pInt celny, zkontrolavet poradi update thermo velicin
alphaRho1.boundaryFieldRef() == alpha.boundaryField()*thermo1.rho().boundaryField();
alphaRho2.boundaryFieldRef() == (1.0 - alpha).boundaryField()*thermo2.rho().boundaryField();
alphaRhoU1.boundaryFieldRef() == alphaRho1.boundaryField()*U1.boundaryField();
alphaRhoU2.boundaryFieldRef() == alphaRho2.boundaryField()*U2.boundaryField(); 
alphaRhoE1.boundaryFieldRef() == alphaRho1.boundaryField()*(e1.boundaryField() + 0.5*magSqr(U1.boundaryField()));
alphaRhoE2.boundaryFieldRef() == alphaRho2.boundaryField()*(e2.boundaryField() + 0.5*magSqr(U2.boundaryField()));

/*rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
rhoE.boundaryFieldRef() ==
    rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));*/
