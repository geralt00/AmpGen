#ifndef _DALITZ_
#define _DALITZ_

#include <complex>
#include <vector>
#include <string>

typedef std::complex<double> dcomplex;

#define CH_AB 1
#define CH_AC 2
#define CH_BC 3

const double gPi = 3.1415926535897932385;

dcomplex pol(double a, double phi);

double log_poisson(double nu, double x);

class TDalitz {

  public:

    TDalitz(double md, double ma, double mb, double mc);

    dcomplex scalar_ampl(double x, double y, double m, double gamma, int chan);
    dcomplex f0_ampl(double x, double y, double m, int chan);
    dcomplex vector_ampl(double x, double y, double m, double gamma, int chan);
    dcomplex gs_ampl(double x, double y, double m, double gamma, int chan);
    dcomplex tensor_ampl(double x, double y, double m, double gamma, int chan);
 
    int kine_limits(double x, double y);

  private:
    
    double fMa;
    double fMb;
    double fMc;
    double fMd;
    double fM2a;
    double fM2b;
    double fM2c;
    double fM2d;
    
    double fM2sum;
    
    int fMabForVectors;
    int fMabForTensors;
    int fWidthDep;
    
    double fR2r;
    double fR2d;

};

#endif
