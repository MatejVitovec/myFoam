#ifndef IF97_reg1_H
#define IF97_reg1_H

namespace reg1
{

constexpr double pStar = 16.53e6;     // [Pa]
constexpr double TStar = 1386.0;      // [K]

constexpr double n[] = {
     0.14632971213167,     -0.84548187169114,    -0.37563603672040e1,
     0.33855169168385e1,   -0.95791963387872,     0.15772038513228,
    -0.16616417199501e-1,   0.81214629983568e-3,  0.28319080123804e-3,
    -0.60706301565874e-3,  -0.18990068218419e-1, -0.32529748770505e-1,
    -0.21841717175414e-1,  -0.52838357969930e-4, -0.47184321073267e-3,
    -0.30001780793026e-3,   0.47661393906987e-4, -0.44141845330846e-5,
    -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5,
    -0.85205128120103e-9,  -0.22425281908000e-5, -0.65171222895601e-6,
    -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,
    -0.17424871230634e-9,  -0.68762131295531e-18, 0.14478307828521e-19,
     0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23,
    -0.93537087292458e-25};

constexpr int J[] = {
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4,
    0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};

constexpr int I[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4,
    4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};


inline double gamma(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma += n[i]*pow(7.1 - pi, I[i])*pow(tau - 1.222, J[i]);
    }
    return gamma;
}

inline double gamma_pi(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_pi = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_pi += -n[i]*I[i]*pow(7.1 - pi, I[i] - 1)*pow(tau - 1.222, J[i]);
    }
    return gamma_pi;
}

inline double gamma_pipi(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_pipi = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_pipi += n[i]*I[i]*(I[i] - 1)*pow(7.1 - pi, I[i] - 2)*pow(tau - 1.222, J[i]);
    }
    return gamma_pipi;
}

inline double gamma_tau(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_tau = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_tau += n[i]*pow(7.1 - pi, I[i])*J[i]*pow(tau - 1.222, J[i] - 1);
    }   
    return gamma_tau;
}

inline double gamma_tautau(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_tautau = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_tautau += n[i]*pow(7.1 - pi, I[i])*J[i]*(J[i] - 1)*pow(tau - 1.222, J[i] - 2);
    }
    return gamma_tautau;
}

inline double gamma_tautautau(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_tautautau = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_tautautau += n[i]*pow(7.1 - pi, I[i])*J[i]*(J[i] - 1)*(J[i] - 2)*pow(tau - 1.222, J[i] - 3);
    }
    return gamma_tautautau;
}

inline double gamma_pitau(double p, double T)
{
    const double pi = p/pStar;
    const double tau = TStar/T;
    double gamma_pitau = 0.0;
    for (int i = 0; i < 34; i++)
    {
        gamma_pitau += -n[i]*I[i]*pow(7.1 - pi, I[i] - 1)*J[i]*pow(tau - 1.222, J[i] - 1);
    }   
    return gamma_pitau;
}

//======================================================================


// Specific volume [m3/kg]
inline double v(double p, double T)
{
    return R*T*(1.0/pStar)*gamma_pi(p, T);
}

// Compression factor [-] (z=pv/RT)
inline double z(double p, double T)
{
    //return (p*v(p, T))/(R*T);
    return (p/pStar)*gamma_pi(p, T);
}

// Specific enthalpy [J/kg]
inline double h(double p, double T)
{
    return R*TStar*gamma_tau(p, T);
}

// Specific internal energy [J/kg]
inline double u(double p, double T)
{
    return R*(TStar*gamma_tau(p, T) - T*(p/pStar)*gamma_pi(p, T));
}

// Specific entropy [J/kg/K]
inline double s(double p, double T)
{
    return R*((TStar/T)*gamma_tau(p, T) - gamma(p, T));
}

// Specific heat capacity at const. pressure [J/kg/K]
inline double cp(double p, double T)
{
    return R*(-pow(TStar/T, 2)*gamma_tautau(p, T));
}

inline double dcpdT(double p, double T)
{
    const double tau = TStar/T;
    return (R/T)*pow(tau, 2)*(2*gamma_tautau(p, T) + tau*gamma_tautautau(p, T));
}


// Specific heat capacity at const. volume [J/kg/K]
inline double cv(double p, double T)
{
    const double tau = TStar/T;
    return R*(-pow(tau, 2)*gamma_tautau(p, T) + pow(gamma_pi(p, T) - tau*gamma_pitau(p, T), 2)/gamma_pipi(p, T));
}


inline double w2(double p, double T)
{
    const double tau = TStar/T;
    const double gp  = gamma_pi(p,T);
    const double gpp = gamma_pipi(p,T);
    const double gtt = gamma_tautau(p,T);
    const double gpt = gamma_pitau(p,T);

    return R*T*(pow(gp, 2)/(pow(gp - tau*gpt, 2)/(pow(tau, 2)*gtt) - gpp));
}


}

#endif
