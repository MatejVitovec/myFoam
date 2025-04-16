#include <iostream>

#include "IF97.H"

using namespace std;

int main() {

    double p = 10e4;
    double T = 350;

    cout << "v     = "   << IF97::reg1::v(p,T) << endl;
    cout << "rho   = " << 1/IF97::reg1::v(p,T) << endl;
    cout << endl;
    
    cout << "h     = "   << IF97::reg1::h(p,T) << endl;
    cout << endl;
    
    cout << "u     = "   << IF97::reg1::u(p,T) << endl;
    cout << endl;
    
    cout << "s     = "   << IF97::reg1::s(p,T) << endl;
    cout << endl;
    
    cout << "cp     = "   << IF97::reg1::cp(p,T) << endl;
    cout << endl;
    
    cout << "cv     = "   << IF97::reg1::cv(p,T) << endl;
    cout << endl;

    cout << "w      = "   << sqrt(IF97::reg1::w2(p,T)) << endl;
    cout << endl;

    return 0;
};
