/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "makeGasProperties.H"

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "AungierRedlichKwongGas.H"
#include "pVirialGas.H"
#include "stiffenedGas.H"
#include "NASG.H"
#include "IAPWSIF97metaGas.H"
#include "IAPWSIF97reg1.H"

#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "IAPWSIF97metaThermo.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "IAPWSIF97Transport.H"

#include "hPolynomialThermo.H"
#include "hTabulatedThermo.H"
#include "polynomialTransport.H"

#ifdef COOLPROP
#include "CoolPropTransport.H"
#include "CoolPropThermo.H"
#include "CoolPropGas.H"
#endif

namespace Foam
{

makeGasProperties(
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeGasProperties(
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    stiffenedGas,
    specie
);

makeGasProperties(
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    NASG,
    specie
);

makeGasProperties(
    constTransport,
    sensibleInternalEnergy,
    IAPWSIF97metaThermo,
    IAPWSIF97metaGas,
    specie
);

makeGasProperties(
    IAPWSIF97Transport,
    sensibleInternalEnergy,
    IAPWSIF97metaThermo,
    IAPWSIF97metaGas,
    specie
);

makeGasProperties(
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    IAPWSIF97reg1,
    specie
);


makeGasProperties(
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    PengRobinsonGas,
    specie
);

makeGasProperties(
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeGasProperties(
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    AungierRedlichKwongGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    AungierRedlichKwongGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    AungierRedlichKwongGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    pVirialGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    pVirialGas,
    specie
);

makeGasProperties(
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    pVirialGas,
    specie
);

#ifdef COOLPROP
makeGasProperties(
    CoolPropTransport,
    sensibleEnthalpy,
    CoolPropThermo,
    CoolPropGas,
    specie
);
#endif

// - Optimized gas properties

// - perfectGas, hConst

#define optimizePerfectGas(Thermo, Transport,Type)                         \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<perfectGas<specie>>,Type>>        \
        >::c(scalar p, scalar T) const                                     \
    {                                                                      \
        scalar cp = Cp(p,T);                                               \
        scalar gamma = cp / (cp - R());                                    \
        return sqrt(gamma*R()*T);                                          \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<perfectGas<specie>>,Type>>        \
        >::beta_p(scalar p, scalar T) const                                \
    {                                                                      \
        return 1/T;                                                        \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<perfectGas<specie>>,Type>>        \
        >::beta_T(scalar p, scalar T) const                                \
    {                                                                      \
        return 1/p;                                                        \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<perfectGas<specie>>,Type>>        \
    >::THs(const scalar Hs, const scalar p, const scalar T0) const         \
    {                                                                      \
        scalar Href = this->Hs(0,0);                                       \
        return (Hs - Href)/Cp(p,T0);                                       \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<perfectGas<specie>>,Type>>        \
    >::pEs(const scalar Es, const scalar rho, const scalar p0) const       \
    {                                                                      \
	const scalar Tref = 288.15;                                            \
        scalar Eref = this->Es(p0,Tref);                                   \
        scalar T = Tref + (Es - Eref)/Cv(p0,Tref);                         \
        return rho*R()*T;                                                  \
    }                                                                      \
                                                                           \
    


optimizePerfectGas(hConstThermo, constTransport, sensibleEnthalpy);
optimizePerfectGas(hConstThermo, sutherlandTransport, sensibleEnthalpy);
optimizePerfectGas(eConstThermo, constTransport, sensibleInternalEnergy);


// =============== Optimization of Stiffened gas model

#define optimizeStiffenedGas(Thermo, Transport,Type)                       \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<stiffenedGas<specie>>,Type>>      \
        >::c(scalar p, scalar T) const                                     \
    {                                                                      \
        scalar cp = Cp(p,T);                                               \
        scalar gamma = cp / (cp - R());                                    \
        return sqrt(gamma*R()*T);                                          \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<stiffenedGas<specie>>,Type>>      \
        >::beta_p(scalar p, scalar T) const                                \
    {                                                                      \
        return 1/T;                                                        \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<stiffenedGas<specie>>,Type>>      \
        >::beta_T(scalar p, scalar T) const                                \
    {                                                                      \
        return Z(p, T)/p;                                                  \
    }                                                                      \
                                                                           \

optimizeStiffenedGas(eConstThermo, constTransport, sensibleInternalEnergy);


// =============== Optimization of modified NASG

/*#define optimizeModifiedNASG(Thermo, Transport,Type)                       \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<NASG<specie>>,Type>>      \
        >::c(scalar p, scalar T) const                                     \
    {                                                                      \
        scalar rho_ = transport_.rho(p, T);                                \
        scalar pPlusPInf = (R()*T)/(1 - (transport_.H(p, T)*rho_)/p);      \
        scalar cp = Cp(p,T);                                               \
        scalar gamma = cp/(cp - R());                                      \
        return (pPlusPInf/rho_)*sqrt(gamma/(R()*T));                       \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<NASG<specie>>,Type>>      \
        >::beta_p(scalar p, scalar T) const                                \
    {                                                                      \
        return (1.0 - (transport_.H(p, T)*rho(p, T))/p)/T;                 \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<NASG<specie>>,Type>>      \
        >::beta_T(scalar p, scalar T) const                                \
    {                                                                      \
        scalar pPlusPInf = (R()*T)/(1 - (transport_.H(p, T)*rho(p, T))/p); \
        return (rho(p, T)*T*R())/pow(pPlusPInf, 2);                        \
    }                                                                      \
                                                                           \
*/

//optimizeModifiedNASG(eConstThermo, constTransport, sensibleInternalEnergy);


// =============== Optimization of Aungier-Redlich-Kwong gas model

#define optimizeAungierRedlichKwongGas(Thermo,Transport,Type)             \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<AungierRedlichKwongGas<specie>>,Type>>  \
    >::pEs(const scalar Es, const scalar rho, const scalar p0) const       \
    {                                                                      \
        scalar p = p0;                                                     \
        scalar T = p0/(rho*R());                                           \
        const scalar tol = 1.e-8;                                          \
        const label maxIter = 100;                                         \
                                                                           \
        const scalar TTol = T*tol;                                         \
        const scalar pTol = p*tol;                                         \
        scalar dT, dp;                                                     \
        label iter = 0;                                                    \
        do                                                                 \
        {                                                                  \
            dT = (Es - this->Es(p, T))/this->Cv(p,T);                      \
            T += dT;                                                       \
            dp = transport_.p(rho,T) - p;                                  \
            p += dp;                                                       \
                                                                           \
            if (iter++ > maxIter)                                          \
            {                                                              \
                FatalErrorInFunction                                       \
                    << "Maximum number of iterations exceeded: " << maxIter \
                        << " when starting from p0:" << p0                 \
                        << " T  : " << T                                   \
                        << " p  : " << p                                   \
                        << " tol: " << tol                                 \
                        << abort(FatalError);                              \
            }                                                              \
                                                                           \
        } while( (mag(dT) > TTol) || (mag(dp) > pTol) );                   \
        return p;                                                          \
    }


optimizeAungierRedlichKwongGas(hConstThermo, constTransport, sensibleEnthalpy);
optimizeAungierRedlichKwongGas(hConstThermo, sutherlandTransport, sensibleEnthalpy);
optimizeAungierRedlichKwongGas(hPolynomialThermo, constTransport, sensibleEnthalpy);
optimizeAungierRedlichKwongGas(hPolynomialThermo, sutherlandTransport, sensibleEnthalpy);
optimizeAungierRedlichKwongGas(janafThermo, constTransport, sensibleEnthalpy);
optimizeAungierRedlichKwongGas(janafThermo, sutherlandTransport, sensibleEnthalpy);


// =============== Optimization of p-Virial gas model

#define optimizePVirialGas(Thermo,Transport,Type)                          \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<pVirialGas<specie>>,Type>>        \
        >::beta_p(scalar p, scalar T) const                                \
    {                                                                      \
        return 1/T;                                                        \
    }                                                                      \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<Thermo<pVirialGas<specie>>,Type>>        \
        >::beta_T(scalar p, scalar T) const                                \
    {                                                                      \
        scalar z  = this->Z(p, T);                                         \
        scalar z0 = this->Z(0.0, T);                                       \
        return z0/z / p;                                                   \
    }                                                                      \

#define optimizehConstPVirialGas(Transport,Type)                           \
                                                                           \
    template<>                                                             \
    scalar genericGasProperties<                                           \
        Transport<species::thermo<hConstThermo<pVirialGas<specie>>,Type>>  \
    >::THs(const scalar Hs, const scalar p, const scalar T0) const         \
    {                                                                      \
        scalar Href = this->Hs(0,0);                                       \
        return (Hs - Href)/Cp(p,T0);                                       \
    }                                                                      \
                                                                           \

optimizePVirialGas(hConstThermo, constTransport, sensibleEnthalpy);
optimizePVirialGas(hPolynomialThermo, sutherlandTransport, sensibleEnthalpy);

optimizehConstPVirialGas(constTransport, sensibleEnthalpy);
optimizehConstPVirialGas(sutherlandTransport, sensibleEnthalpy);

// =============== Optimization of CoolProp model =====================
#ifdef COOLPROP

template<> 
scalar genericGasProperties<CoolPropTransport<
    species::thermo<CoolPropThermo<CoolPropGas<specie>>,sensibleEnthalpy>>>::THs
    (const scalar Hs, const scalar p, const scalar T0) const         	
{                                                                      
    auto state = this->transportModel().state();
    state->update(CoolProp::HmassP_INPUTS, Hs, p);
    return state->T();   
} 


#endif

}

// ************************************************************************* //
