#include <iostream>
#include <fstream>

#include "IF97.H"

using namespace std;
using namespace IF97;

double sutherland(double p, double T)
{
    return 1.823e-6*std::sqrt(T)/(1.0  + 673.0/T);
}

double sutherlandLambda(double p, double T)
{
    double cv = reg2meta::cv(p, T);
    return sutherland(p, T)*cv*(1.32 + 1.77*R/cv);
}

double lambdaPoly(double p, double T)
{
    return 7.341e-3 - 1.013e-5*T + 1.801e-7*T*T - 9.1e-11*T*T*T;
}

int main() {



    double p = 40000.0;

    ofstream f;
    f.open ("viscosity.txt");

    for (int i = 0; i < 101; i++)
    {
        double T = 300 + i;
        f << T << " " << mu(1/IF97::reg2meta::v(p, T), T) << "\n";
    }
    f << endl;
    f.close();

    /*f.open ("conductivity.txt");

    for (int i = 0; i < 101; i++)
    {
        double T = 300 + i;
        f << T << " " << lambda(1/IF97::reg2meta::v(p, T), T) << " " << sutherlandLambda(p, T) << " " << lambdaPoly(p, T) << "\n";
    }
    f << endl;
    f.close();*/
    return 0;
};
