#include <iostream>
#include <fstream>
#include <vector>

#include "IF97.H"

double mean(const std::vector<double>& u)
{
    double out = 0.0;
    for (int i = 0; i < u.size(); i++)
    {
        out += u[i];
    }
    
    return out/u.size();
}

double mean(const std::vector<std::vector<double>>& u)
{
    double out = 0.0;
    for (int j = 0; j < u.size(); j++)
    {
        double aux = 0.0;
        for (int i = 0; i < u[j].size(); i++)
        {
            aux += u[j][i];
        }
        out += aux/u[j].size();
    }
    
    return out/u.size();
}

double cpg(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& hg)
{
    const double meanT = mean(T);
    const double meanhg = mean(hg);

    double num = 0.0;
    double denom = 0.0;
    for (int j = 0; j < p.size(); j++)
    {
        for (int i = 0; i < T.size(); i++)
        {
            num += T[i]*(hg[j][i] - meanhg);
            denom += T[i]*(T[i] - meanT);
        } 
    }
    return num/denom;
}

double cpgMcvg(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& vg)
{
    double num = 0.0;
    double denom = 0.0;
    for (int j = 0; j < p.size(); j++)
    {
        for (int i = 0; i < T.size(); i++)
        {
            num += vg[j][i]*(T[i]/p[j]);
            denom += std::pow(T[i]/p[j], 2);
        } 
    }

    return num/denom;
}

double qg(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& hg)
{
    const double meanT = mean(T);
    const double meanhg = mean(hg);

    return meanhg - cpg(p, T, hg)*meanT;
}


double cpl(const std::vector<double>& p, const std::vector<double>& T,const std::vector<std::vector<double>>& hl, const double pInfl, const double bl)
{
    const double meanp = mean(p);
    const double meanT = mean(T);
    const double meanhl = mean(hl);

    double num1 = 0.0;
    double num2 = 0.0;
    double denom = 0.0;

    for (int j = 0; j < p.size(); j++)
    {
        for (int i = 0; i < T.size(); i++)
        {
            num1 += T[i]*(hl[j][i] - meanhl);
            num2 += T[i]*(p[j] - meanp);

            denom += T[i]*(T[i] - meanT);
        }        
    }

    return (num1 + bl*num2)/denom;
}

double meanTdP(const std::vector<double>& p, const std::vector<double>& T, double pInfl)
{
    double out = 0.0;
    for (int j = 0; j < p.size(); j++)
    {
        for (int i = 0; i < T.size(); i++)
        {
            out += T[i]/(p[j] + pInfl);
        }
    }
    
    return out/(p.size()*T.size());
}

double cplMcvl(const std::vector<double>& p, const std::vector<double>& T,const std::vector<std::vector<double>>& vl, double pInfl)
{
    const double meanp = mean(p);
    const double meanT = mean(T);
    const double meanvl = mean(vl);

    double num = 0.0;
    double denom = 0.0;

    double meanTdP_ = meanTdP(p, T, pInfl);
    

    for (int j = 0; j < p.size(); j++)
    {
        for (int i = 0; i < T.size(); i++)
        {
            num += (T[i]/(p[j] + pInfl))*(vl[j][i] - meanvl);
            denom += (T[i]/(p[j] + pInfl))*((T[i]/(p[j] + pInfl)) - meanTdP_);
        }        
    }
    
    return num/denom;
}


double bl(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& vl, double pInfl)
{
    return mean(vl) - cplMcvl(p, T, vl, pInfl)*meanTdP(p, T, pInfl);
}

double pInflFunc(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& vl, const std::vector<std::vector<double>>& hl, const double pInfl, const double p0, const double rhol0, const double cl0sqr)
{
    const double bl_ = bl(p, T, vl, pInfl);
    return p0 + pInfl - (1.0 - cplMcvl(p, T, vl, pInfl)/cpl(p, T, hl, pInfl, bl_))*rhol0*cl0sqr*(1.0 - bl_*rhol0);
}

double ql(const std::vector<double>& p, const std::vector<double>& T, const std::vector<std::vector<double>>& vl, const std::vector<std::vector<double>>& hl, double pInfl)
{
    const double meanp = mean(p);
    const double meanT = mean(T);
    const double meanhl = mean(hl);

    const double bl_ = bl(p, T, vl, pInfl);
    
    return meanhl - cpl(p, T, hl, pInfl, bl_)*meanT - bl_*meanp;
}





int main() {

    const int np = 800;
    const int nT = 800;

    const double pStart = 3e4;
    const double pEnd = 8e4;
    const double dp = (pEnd - pStart)/np;

    std::vector<double> p(np);
    for (int i = 0; i < nT; i++)
    {
        p[i] = pStart + i*dp;
    }

    const double TStart = 280;
    const double TEnd = 350;
    const double dT = (TEnd - TStart)/nT;

    std::vector<double> T(nT);
    for (int i = 0; i < nT; i++)
    {
        T[i] = TStart + i*dT;
    }

    std::vector<std::vector<double>> vg = std::vector<std::vector<double>>(np);
    std::vector<std::vector<double>> vl = std::vector<std::vector<double>>(np);
    std::vector<std::vector<double>> hg = std::vector<std::vector<double>>(np);
    std::vector<std::vector<double>> hl = std::vector<std::vector<double>>(np);

    for (int j = 0; j < p.size(); j++)
    {
        vg[j] = std::vector<double>(nT);
        vl[j] = std::vector<double>(nT);
        hg[j] = std::vector<double>(nT);
        hl[j] = std::vector<double>(nT);
        for (int i = 0; i < T.size(); i++)
        {
            vg[j][i] = IF97::reg2meta::v(p[j], T[i]);
            vl[j][i] = IF97::reg1::v(p[j], T[i]);
            hg[j][i] = IF97::reg2meta::h(p[j], T[i]);
            hl[j][i] = IF97::reg1::h(p[j], T[i]);
        }
    }
    

    const double p0 = 104530;
    const double rhol0 = 957.74;
    const double cl0sqr = std::pow(1542, 2);

    double pInfl = 8e8;
    double dpInfl = 1000.0;
    int iter = 0;
    while(dpInfl > 1e-10)
    {
        const double pInflFunc0 = pInflFunc(p, T, vl, hl, pInfl, p0, rhol0, cl0sqr);
        const double pInfl1 = (1.0 + 1e-3)*pInfl;
        const double pInflFunc1 = pInflFunc(p, T, vl, hl, pInfl1, p0, rhol0, cl0sqr);
        
        dpInfl = (pInflFunc0*(pInfl1 - pInfl))/(pInflFunc1 - pInflFunc0);

        pInfl -= dpInfl;
        iter++;
    }

    double bl_ = bl(p, T, vl, pInfl);
    double cpl_ = cpl(p, T, hl, pInfl, bl_);
    double cvl_ = -(cplMcvl(p, T, vl, pInfl) - cpl_);
    double gammal = cpl_/cvl_;
    double ql_ = ql(p, T, vl, hl, pInfl);

    double cpg_ = cpg(p, T, hg);
    double cvg_ = -(cpgMcvg(p, T, vg) - cpg_);
    double gammag = cpg_/cvg_;
    double qg_ = qg(p, T, hg);


    std::cout << std::endl;
    std::cout << "ITER: " << iter << std::endl;
    std::cout << std::endl;
    std::cout << "Cp_g: " << cpg_ << " Cv_g: " << cvg_ << " gamma_g: " << gammag << " q_g: " << qg_ << " q'_g: " << 0.0 <<  std::endl;
    std::cout << "Cp_l: " << cpl_ << " Cv_l: " << cvl_ << " gamma_l: " << gammal << " pInf_l: " << pInfl << " b_l: " << bl_ <<  " q_l: " << ql_ << " q'_l: " << 0.0 <<  std::endl;
    std::cout << std::endl;

    return 0;
};


