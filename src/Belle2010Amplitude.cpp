#include "extern/Belle2010Amplitude.h"
#include <iostream>
#include <random>
#include "extern/TDalitz.h"
#include "TComplex.h"
#include "TMath.h"
#include <fstream>
#include <sstream>
#include <algorithm>

Belle2010Amplitude::Belle2010Amplitude(std::string kaon_type, std::string cp_mult_dir) :
DtoKpipiAmplitude(kaon_type)
{

    std::ifstream fin;
    fin.open(cp_mult_dir);
    double r, delta;
    for (int i = 0; i < 8; i++){
        fin >> r >> delta;
        _r_value[i] = r;
        _delta[i] = delta;
    }
    set_amplitude_factors();
};


Belle2010Amplitude::~Belle2010Amplitude(){}


void Belle2010Amplitude::set_amplitude_factors(){
    if (_kaon_type == "KS"){
        std::cout << "Using default KS amplitude factors\n";
        _CF_factor = dcomplex(1.);
        _DCS_factor = dcomplex(1.);
        _nonRes_factor = dcomplex(1.);
        for (int i = 0; i < 8; i++){
            _pipi_factors.push_back(dcomplex(1.));
        }

    } else {
        std::cout << "Using default KL amplitude factors\n";
        
        _CF_factor = dcomplex(-1.);
        _DCS_factor = dcomplex(1.);
        _nonRes_factor = dcomplex(-1.);
        // _nonRes_factor = get_amplitude_factor(r, delta);
        for (int i = 0; i < 8; i++){
            double r = _r_value[i];
            double delta = _delta[i];
            _pipi_factors.push_back(get_amplitude_factor(r, delta));
        } 
    }



}



dcomplex Belle2010Amplitude::get_amplitude(double sKp, double sKm) const {

  

    // const double m_ks  = 497.61401-3;
    // const double m_pi  = 139.57017e-3;
    // const double m_d0  = 1864.8399e-3;
    //need to be consistent
    const double m_ks  = 497.672e-3;
    const double m_pi  = 139.56995e-3;
    const double m_d0  = 1864.6e-3;
    TDalitz* dkpp_dalitz;


    dkpp_dalitz = new TDalitz(m_d0, m_ks, m_pi, m_pi);

        if(!dkpp_dalitz->kine_limits(sKp, sKm)){
        return 0;
    }
    dcomplex dk2pires[19]; 

    dcomplex amp=0;

    auto mkp = sKm;
    auto mkm = sKp;

    dk2pires[0]=    pol(1.560805, 214.065284)*dkpp_dalitz->scalar_ampl(mkp, mkm, 0.522477, 0.453106, CH_BC) ; //pipi
    dk2pires[1]=    pol(1.,       0.        )*dkpp_dalitz->gs_ampl(    mkp, mkm, 0.771700, 0.13600, CH_BC) ;
    dk2pires[2]=    pol(0.491463,  64.390750)*dkpp_dalitz->gs_ampl(    mkp, mkm, 1.4590,   0.4550,   CH_BC) ;
    dk2pires[3]=    pol(0.034337, 111.974402)*dkpp_dalitz->vector_ampl(mkp, mkm, 0.782650, 0.008490, CH_BC) ;
    dk2pires[4]=    pol(0.385497, 207.278721)*dkpp_dalitz->f0_ampl(    mkp, mkm, 0.97700,            CH_BC) ;
    dk2pires[5]=    pol(0.203222, 212.128769)*dkpp_dalitz->scalar_ampl(mkp, mkm, 1.033172, 0.087984, CH_BC) ;
    dk2pires[6]=    pol(1.436933, 342.852060)*dkpp_dalitz->tensor_ampl(mkp, mkm, 1.275400, 0.185200, CH_BC) ;
    dk2pires[7]=    pol(1.561670, 109.586718)*dkpp_dalitz->scalar_ampl(mkp, mkm, 1.434000, 0.173000, CH_BC) ;

    dk2pires[8]=    pol(1.638345, 133.218073)*dkpp_dalitz->vector_ampl(mkp, mkm, 0.89370, 0.048400, CH_AB) ;//CF
    dk2pires[9]=    pol(0.651326, 119.929357)*dkpp_dalitz->vector_ampl(mkp, mkm, 1.414000, 0.232000, CH_AB) ;
    dk2pires[10]=    pol(2.209476, 358.855281)*dkpp_dalitz->scalar_ampl(mkp, mkm, 1.414000, 0.290000, CH_AB) ; 
    dk2pires[11]=    pol(0.890176, 314.789317)*dkpp_dalitz->tensor_ampl(mkp, mkm, 1.425600, 0.098500, CH_AB) ;
    dk2pires[12]=    pol(0.877231,  82.271754)*dkpp_dalitz->vector_ampl(mkp, mkm, 1.717000, 0.322000, CH_AB) ;

    dk2pires[13]=     pol(0.149579, 325.356816)*dkpp_dalitz->vector_ampl(mkp, mkm, 0.891660, 0.050800, CH_AC) ;//DCS
    dk2pires[14]=     pol(0.423211, 252.523919)*dkpp_dalitz->vector_ampl(mkp, mkm, 1.414000, 0.232000, CH_AC) ; 
    dk2pires[15]=     pol(0.364030,  87.118694)*dkpp_dalitz->scalar_ampl(mkp, mkm, 1.414000, 0.290000, CH_AC) ; 
    dk2pires[16]=     pol(0.228236, 275.203595)*dkpp_dalitz->tensor_ampl(mkp, mkm, 1.425600, 0.098500, CH_AC) ;
    dk2pires[17]=    pol(2.081620, 130.047574)*dkpp_dalitz->vector_ampl(mkp, mkm, 1.717000, 0.322000, CH_AC) ;
        
    dk2pires[18]=    pol(2.666864, 160.480162);
    //Actually this part is building
    //\mathcal{T}(\mathbf{x}) = \sum_{k} g_k A_k (\mathbf{x})


    for(int i=0;i<19;++i){
        amp+=dk2pires[i];
        if (isnan(dk2pires[i].real()) || isnan(dk2pires[i].imag())){
            std::cout << i << " " << mkp << " " << mkm << " " << amp << std::endl;
            amp = dcomplex(0, 0);
        }
    }
    if (isnan(amp.real()) || isnan(amp.imag())){
        std::cout << mkp << " " << mkm << " " << amp << std::endl;
        amp = dcomplex(0, 0);
    }

  return amp;
}

