#include "extern/DtoKpipiAmplitude.h"
#include <iostream>
#include <random>

DtoKpipiAmplitude::DtoKpipiAmplitude(std::string kaon_type) :
_kaon_type(kaon_type),
_mD(1.8648399),
_mKs(0.49761401),
_mPi(0.13957017)
{
    if (_kaon_type == "KL" || _kaon_type == "KS" ) {
        std::cout << "Initialising amplitude for : " << _kaon_type << std::endl;
    } else {
        std::cout << "UNKNOWN KAON TYPE: " << _kaon_type << "! EXITING!\n";
        exit(-1);
    }


};

DtoKpipiAmplitude::~DtoKpipiAmplitude(){}

dcomplex DtoKpipiAmplitude::get_amplitude_factor(double r, double delta) const {
    if (r<0) r=0;
    dcomplex tmp = dcomplex(1);
    tmp -= dcomplex(2*r*cos(delta), 2*r*sin(delta));
    tmp *= dcomplex(-1.0, 0.);
    return tmp;
}

bool DtoKpipiAmplitude::is_valid_point(double x, double y) const {

    static const double K0_mass  = 497.61401-3;
    static const double PI_mass  = 139.57017e-3;
    static const double mD0  = 1864.8399e-3;

    double msquared01 = x; // Kpi(RS)
    double msquared02 = y; // Kpi(WS)
    double msquared12 = mD0*mD0 + K0_mass*K0_mass + 2*PI_mass*PI_mass - msquared01 - msquared02;

    double local_ma = K0_mass;
    double local_mb = PI_mass;
    double local_mc = PI_mass;

    double local_xmin = pow(local_ma + local_mb,2);
    double local_xmax = pow(mD0 - local_mc,2);

    // Find energy of b(c) in ab frame
    double ebab = (x - local_ma*local_ma + local_mb*local_mb)/(2.0*sqrt(x));
    double ecab = (mD0*mD0 - x - local_mc*local_mc)/(2.0*sqrt(x));

    double yhi = pow(ebab+ecab,2) - pow( sqrt(ebab*ebab-local_mb*local_mb)-sqrt(ecab*ecab-local_mc*local_mc) ,2);
    double ylo = pow(ebab+ecab,2) - pow( sqrt(ebab*ebab-local_mb*local_mb)+sqrt(ecab*ecab-local_mc*local_mc) ,2);

    // Initialize boolean variable as false.
    bool inDal = false;

    // Return true, if within the Dalitz-plot phase space.
    if ((local_xmin <= x) && (x <= local_xmax) && (ylo <= msquared12) && (msquared12 <= yhi)) { inDal = true; }

    return inDal;
}

