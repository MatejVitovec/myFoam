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

Class
    Foam::saturationIAPWS95

Description

SourceFiles
    saturationIAPWS95.C

\*---------------------------------------------------------------------------*/

#ifndef saturationIAPWS95_H
#define saturationIAPWS95_H

#include "saturationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class saturationCurve Declaration
\*---------------------------------------------------------------------------*/

class saturationIAPWS95
:
    public saturationCurve
{
private:

    // Critical values
    static constexpr scalar Tc = 647.096;
    static constexpr scalar pc = 22.064e6;
    static constexpr scalar rhoc = 322;

    // Saturation line pSat(T) (Wagner & Pruss, 2000)
    static constexpr scalar a[] = {-7.85951783, 1.84408259, -11.7866497,
        22.6807411, -15.9618719, 1.80122502 };

    static constexpr scalar b[] = { 1.99274064, 1.09965342, -0.510839303,
        -1.75493479, -45.5170352, -6.74694450e5 };

    static constexpr scalar c[] = { -2.03150240, -2.68302940, -5.38626492,
        -17.2991605, -44.7586581, -63.9201063 };

    static constexpr scalar d[] = { -1135.905627715, -5.65134998e-8, 2690.66631,
        127.287297, -135.003439, 0.981825814 };

    scalar alphaDalpha0(scalar T) const;

public:

    TypeName("IAPWS95");

    // Constructor
    saturationIAPWS95(const dictionary& dict);
    
    // Member Functions

    virtual scalar ps(scalar T) const;
    virtual scalar Ts(scalar p) const;
    virtual scalar dpsdT(scalar T) const;
    virtual scalar rhosl(scalar T) const;
    virtual scalar rhosv(scalar T) const;
    virtual scalar hsl(scalar T) const;
    virtual scalar hsv(scalar T) const;
    virtual scalar esl(scalar T) const;
    virtual scalar esv(scalar T) const;

    virtual ~saturationIAPWS95() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
