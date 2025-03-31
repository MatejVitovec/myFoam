
const scalar dragCoeff = 2.0;

volScalarField drag = volScalarField
(
    "drag",
    dragCoeff*(3/4)*fluidSystem.thermo1.rho()*(1.0 - fluidSystem.alpha())*mag(fluidSystem.U1() - fluidSystem.U2())*(fluidSystem.U1() - fluidSystem.U2())
)