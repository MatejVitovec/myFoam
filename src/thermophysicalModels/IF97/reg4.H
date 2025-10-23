#ifndef IF97_reg4_H
#define IF97_reg4_H

namespace reg4
{

constexpr double pStar = 1.0e6;     // [Pa]
constexpr double TStar = 1.0;       // [K]

constexpr double n[] = {
     0.11670521452767e4,
    -0.72421316703206e6,
    -0.17073846940092e2,
     0.12020824702470e5,
    -0.32325550322333e7,
     0.14915108613530e2,
    -0.48232657361591e4,
     0.40511340542057e6,
    -0.23855557567849,
     0.65017534844798e3
    };

//======================================================================


// Saturation pressure [Pa]
inline double ps(double Ts)
{
    const double theta = Ts/TStar + n[8]/(Ts/TStar - n[9]);
    const double A =      theta*theta + n[0]*theta + n[1];
    const double B = n[2]*theta*theta + n[3]*theta + n[4];
    const double C = n[5]*theta*theta + n[6]*theta + n[7];

    return pStar*std::pow((2*C)/(-B + std::sqrt(std::pow(B,2) - 4*A*C)), 4);
}

// Saturation temperature [K]
inline double Ts(double ps)
{
    const double beta = std::pow(ps/pStar, 0.25);
    const double E =      beta*beta + n[2]*beta + n[5];
    const double F = n[0]*beta*beta + n[3]*beta + n[6];
    const double G = n[1]*beta*beta + n[4]*beta + n[7];
    const double D = 2*G/(-F - std::sqrt(F*F - 4*E*G));

    return 0.5*TStar*(n[9] + D - std::sqrt(std::pow(n[9] + D, 2) - 4*(n[8] + n[9]*D)));
}

inline double dpsdT(double T)
{
    double psat = ps(T);
    double beta = std::pow(psat*pStar, 0.25);
    double theta = T + n[8]/(T - n[9]);

    // Derivatives of implicit equation (29) with respect to beta and theta
    double dfdbeta = (2*beta + n[2])*std::pow(theta, 2) + (2*beta*n[0] + n[3])*theta + 2*n[1]*beta + n[4];
    double dfdtheta= (2*theta + n[0])*std::pow(beta, 2) + (2*theta*n[2] + n[3])*beta + 2*n[5]*theta + n[6]; 

    double dthetadT = (1 - n[8]/std::pow(T - n[9], 2));
    double dbetadtheta = -dfdtheta/dfdbeta;
    double dpdbeta = 4*std::pow(beta, 3)*pStar;

    return dpdbeta*dbetadtheta*dthetadT;
}


}

#endif
