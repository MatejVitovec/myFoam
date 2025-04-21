#include <iostream>
#include <fstream>
#include <vector>

#include "IF97.H"

using namespace std;

int main() {

    int n = 200;

    double p = 100e4;
    double TStart = 300;
    double TEnd = 550;
    double dT = (TEnd - TStart)/n;

    ofstream f;
    f.open ("datacp.txt");

    for (int i = 0; i < n; i++)
    {
        double T = TStart + i*dT;
        f << T << " " << IF97::reg1::cp(p, T) << "\n";
    }
    f << endl;
    f.close();
    

    /*double pStart = 10000;
    double pEnd = 200000;

    double TStart = 280;
    double TEnd = 600;

    int n = 100;

    double pDiff = (log(pEnd) - log(pStart))/(n-1);
    double TDiff = (TEnd - TStart)/(n-1);

    vector<double> pVec = vector<double>(n);
    vector<double> TVec = vector<double>(n);

    for (int i = 0; i < n; i++)
    {
        pVec[i] = exp(log(pStart) + i*pDiff);
        TVec[i] = TStart + i*TDiff;
    }
    
    vector<vector<double>> v = vector<vector<double>>(n);

    for (int i = 0; i < n; i++)
    {
        v[i] = vector<double>(n);
        for (int j = 0; j < n; j++)
        {
            v[i][j] = IF97::reg1::v(pVec[i], TVec[j]);
        }
    }

    ofstream f;
    f.open ("data.txt");

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            f << pVec[i] << " " << TVec[j] << " " << v[i][j] << "\n";
        }
    }
    f << endl;

    f.close();*/
    

    /*cout << "v     = "   << IF97::reg1::v(p,T) << endl;
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

    return 0;*/
};
