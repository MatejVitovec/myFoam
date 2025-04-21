#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {

    /*double p = 2e4;
    double T = 300;

    double cv = 3202.0;
    double gamma = 1.39;
    double pInf = 8899.0e5;
    double b = 4.78e-4;
    double q = -1244191.0;

    double u = cv*T + (((gamma - 1.0)*pInf)/(p + pInf))*cv*T + q;
    double h = gamma*cv*T + b*p + q;
    double v = ((gamma - 1.0)*cv*T)/(p + pInf) + b;

    cout << "v     = "   << v << endl;
    cout << "rho   = " << 1/v << endl;
    cout << endl;
    
    cout << "h     = "   << h << endl;
    cout << endl;
    
    cout << "u     = "   << u << endl;
    cout << endl;*/

    int n = 200;

    double p = 20e4;
    double TStart = 300;
    double TEnd = 600;
    double dT = (TEnd - TStart)/n;

    double cv = 3202.0;
    double gamma = 1.39;
    double pInf = 8899.0e5;
    double b = 4.78e-4;
    double q = -1244191.0;

    ofstream f;
    f.open ("nasg.txt");

    for (int i = 0; i < n; i++)
    {
        double T = TStart + i*dT;
        double v = ((gamma - 1.0)*cv*T)/(p + pInf) + b;
        //double u = cv*T + (((gamma - 1.0)*pInf)/(p + pInf))*cv*T + q;
        f << T << " " << v << "\n";
    }
    f << endl;
    f.close();

    return 0;
};
