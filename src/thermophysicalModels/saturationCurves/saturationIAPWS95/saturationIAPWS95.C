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

#include "saturationIAPWS95.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(saturationIAPWS95, 0);
addToRunTimeSelectionTable(saturationCurve, saturationIAPWS95, dict);


saturationIAPWS95::saturationIAPWS95(const dictionary& dict)
:
    saturationCurve(dict)
{}

scalar saturationIAPWS95::alphaDalpha0(scalar T) const
{
    const scalar theta = 1 - T/Tc;
    
    return d[0] + d[1]*pow(theta, -19) + d[2]*theta + d[3]*pow(theta, 9/2)
        + d[4]*pow(theta, 5) + d[5]*pow(theta, 54.5);
}


scalar saturationIAPWS95::dpsdT(scalar T) const
{
    const scalar theta = 1 - T/Tc;
    const scalar theta05 = sqrt(theta);
    const scalar theta15 = theta*theta05;
    const scalar theta20 = theta*theta;
    const scalar theta25 = theta20*theta05;
    const scalar theta30 = theta20*theta;
    const scalar theta35 = theta30*theta05;
    const scalar theta40 = theta20*theta20;
    const scalar theta65 = theta40*theta25;
    const scalar theta75 = theta40*theta35;

    const scalar lnP = Tc/T*(a[0]*theta + a[1]*theta15 + a[2]*theta30 + a[3]*theta35 + a[4]*theta40 + a[5]*theta75);
    const scalar p = pc*exp(lnP);
    
    return -p/T*(lnP + a[0] + 1.5*a[1]*theta05 + 3*a[2]*theta20 +
        3.5*a[3]*theta25 + 4*a[4]*theta30 + 7.5*a[5]*theta65);
}


scalar saturationIAPWS95::ps(scalar T) const
{
    const scalar theta = 1 - T/Tc;
    const scalar theta05 = sqrt(theta);
    const scalar theta15 = theta*theta05;
    const scalar theta30 = theta*theta*theta;
    const scalar theta35 = theta30*theta05;
    const scalar theta40 = theta30*theta;
    const scalar theta75 = theta40*theta35;
    
    return pc*exp(Tc/T*(
        a[0]*theta + a[1]*theta15 + a[2]*theta30 + a[3]*theta35 + a[4]*theta40 + a[5]*theta75));
}


scalar saturationIAPWS95::Ts(scalar p) const
{
    scalar T = 647.096;
    
    // Newton method for f = ps(T) - p = 0    
    int iter = 0;
    while (iter < 10) 
    {
        scalar f = ps(T) - p;
        if (fabs(f) < 1e-4*p) break;
        T -= f / dpsdT(T);
        iter++;
    }

    return T;
}

scalar saturationIAPWS95::rhosl(scalar T) const
{
    const scalar theta = 1 - T/Tc;
    const scalar theta1d3 = pow(theta, 1/3);
    const scalar theta2d3 = theta1d3*theta1d3;
    const scalar theta5d3 = theta1d3*theta2d3*theta2d3;
    const scalar theta16d3 = theta5d3*theta5d3*theta5d3*theta1d3;
    const scalar theta43d3 = theta16d3*theta16d3*theta5d3*theta5d3*theta1d3;
    const scalar theta110d3 = theta43d3*theta43d3*theta16d3*theta5d3*theta2d3*theta1d3;

    return (1.0 + b[0]*theta1d3 + b[1]*theta2d3 + b[2]*theta5d3 + b[3]*theta16d3
        + b[4]*theta43d3 + b[5]*theta110d3)*rhoc;
}

scalar saturationIAPWS95::rhosv(scalar T) const
{
    const scalar theta = 1 - T/Tc;
    const scalar theta1d6 = pow(theta, 1/6);
    const scalar theta2d6 = theta1d6*theta1d6;
    const scalar theta4d6 = theta2d6*theta2d6;
    const scalar theta8d6 = theta4d6*theta4d6;
    const scalar theta18d6 = theta8d6*theta8d6*theta2d6;
    const scalar theta37d6 = theta18d6*theta18d6*theta1d6;
    const scalar theta71d6 = theta37d6*theta18d6*theta8d6*theta8d6;

    return exp(c[0]*theta2d6 + c[1]*theta4d6 + c[2]*theta8d6
        + c[3]*theta18d6 + c[4]*theta37d6 + c[5]*theta71d6)*rhoc;
}

scalar saturationIAPWS95::hsl(scalar T) const
{
    return alphaDalpha0(T) + 1e-3*(T/rhosl(T))*dpsdT(T);
}

scalar saturationIAPWS95::hsv(scalar T) const
{
    return alphaDalpha0(T) + 1e-3*(T/rhosv(T))*dpsdT(T);
}

scalar saturationIAPWS95::esl(scalar T) const
{
    return alphaDalpha0(T) + (1e-3*T*dpsdT(T) - ps(T))/rhosl(T);
}

scalar saturationIAPWS95::esv(scalar T) const
{
    return alphaDalpha0(T) + (1e-3*T*dpsdT(T) - ps(T))/rhosv(T);
}

}
