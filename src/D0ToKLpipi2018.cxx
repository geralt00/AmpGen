#include "extern/D0ToKLpipi2018.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <math.h>
#include "TMath.h"
#include "TLorentzVector.h"



using namespace std;  //::endl;

D0ToKLpipi2018::~D0ToKLpipi2018() {}

void D0ToKLpipi2018::init(){
  std::cout << "D0ToKLpipi2018 ==> Initialization !" << std::endl;

  _nd = 3;
  tan2thetaC = (0.22650*0.22650)/(1.-(0.22650*0.22650)) ;  //sin(theta_C) = 0.22650 +/- 0.00048
  pi180inv = 1.0*3.1415926/180;

	double Pi = 3.1415926;

  mass_R[0]=   0.77155; width_R[0]=  0.13469; spin_R[0]=  1; ar[0]=  1;        phir[0]=  0;
  mass_R[1]=   0.78265; width_R[1]=  0.00849; spin_R[1]=  1; ar[1]=  0.038791; phir[1]=  (180./Pi)*2.1073;
  mass_R[2]=   1.27510; width_R[2]=  0.18420; spin_R[2]=  2; ar[2]=  1.42887;  phir[2]=  (180./Pi)*-0.633296;
  mass_R[3]=   1.46500; width_R[3]=  0.40000; spin_R[3]=  1; ar[3]=  2.85131;  phir[3]=  (180./Pi)*1.7820801;
  mass_R[4]=   0.89371; width_R[4]=  0.04719; spin_R[4]=  1; ar[4]=  1.72044;  phir[4]=  (180./Pi)*2.38835877;
  mass_R[5]=   1.42560; width_R[5]=  0.09850; spin_R[5]=  2; ar[5]=  1.27268;  phir[5]=  (180./Pi)*-0.769095;
  mass_R[6]=   1.71700; width_R[6]=  0.3220;  spin_R[6]=  1; ar[6]=  3.307642; phir[6]=  (180./Pi)*-2.062227;
  mass_R[7]=   1.41400; width_R[7]=  0.2320;  spin_R[7]=  1; ar[7]=  0.286927; phir[7]=  (180./Pi)*1.7346186;
  mass_R[8]=   0.89371; width_R[8]=  0.04719; spin_R[8]=  1; ar[8]=  0.1641792;phir[8]=  (180./Pi)*-0.735903;
  mass_R[9]=   1.42560; width_R[9]=  0.0985;  spin_R[9]=  2; ar[9]=  0.1025736;phir[9]=  (180./Pi)*-1.56397;
  mass_R[10]=  1.41400; width_R[10]= 0.2320;  spin_R[10]= 1; ar[10]= 0.2090326;phir[10]= (180./Pi)*2.6208986;
	////////////////////////////////////////
  mass_R[11]=  1.42500; width_R[11]= 0.2700;  spin_R[11]= 1; ar[11]= 2.36;  phir[11]= 99.4;//not found
	//		beta[kb] = EvtComplex(mag*cos(phase*pi180inv), mag*sin(phase*pi180inv)) ;
  //    fprod[kf] = EvtComplex(mag*cos(phase*pi180inv), mag*sin(phase*pi180inv)) ;
  beta[0] = complex<double>( 8.521486*cos( 1.195641 ), 8.521486*sin( 1.195641));//
  beta[1] = complex<double>( 12.1895 *cos( 0.41802), 12.1895  *sin( 0.41802));
  beta[2] = complex<double>(29.14616  *cos(-0.0018386   ), 29.14616 *sin(-0.0018386   ));
  beta[3] = complex<double>(10.745735  *cos(-0.9057014 ), 10.745735  *sin(-0.9057014 ));
  beta[4] = complex<double>(0., 0.);

  fprod[0] = complex<double>(8.04427*cos(-2.19847), 8.04427*sin(-2.19847));
  fprod[1] = complex<double>(26.2986*cos(-2.65853), 26.2986*sin(-2.65853));
  fprod[2] = complex<double>(33.0349*cos(-1.62714), 33.0349*sin(-1.62714));
  fprod[3] = complex<double>(26.1741*cos(-2.11891), 26.1741*sin(-2.11891));
  fprod[4] = complex<double>(0., 0.);
                                                                                                
	//CP_mult[w] = (EvtComplex(1., 0.) - 2.*tan2thetaC*EvtComplex(r*cos(delta), r*sin(delta))) ;
	CP_mult[0] = complex<double>(1.,0.)-2.*tan2thetaC*complex<double>(1.851 *cos(-94.07  *pi180inv),1.851 *sin(-94.07  *pi180inv));
	CP_mult[1] = complex<double>(1.,0.)-2.*tan2thetaC*complex<double>(6.332 *cos(2.103   *pi180inv),6.332 *sin(2.103   *pi180inv));
	CP_mult[2] = complex<double>(1.,0.)-2.*tan2thetaC*complex<double>(3.229 *cos(-60.05  *pi180inv),3.229 *sin(-60.05  *pi180inv));
	CP_mult[3] = complex<double>(1.,0.)-2.*tan2thetaC*complex<double>(12.75 *cos(73.62   *pi180inv),12.75 *sin(73.62   *pi180inv));
	CP_mult[4] = complex<double>(1.,0.)-2.*tan2thetaC*complex<double>(0.7116*cos(-177.149*pi180inv),0.7116*sin(-177.149*pi180inv));

  ma[0]= 0.651;   g[0][0]= 0.22889; g[0][1]= -0.55377; g[0][2]=  0;       g[0][3]= -0.39899; g[0][4]= -0.34639;
  ma[1]= 1.20360; g[1][0]= 0.94128; g[1][1]=  0.55095; g[1][2]=  0;       g[1][3]=  0.39065; g[1][4]=  0.31503;
  ma[2]= 1.55817; g[2][0]= 0.36856; g[2][1]=  0.23888; g[2][2]=  0.55639; g[2][3]=  0.18340; g[2][4]=  0.18681;
  ma[3]= 1.21000; g[3][0]= 0.33650; g[3][1]=  0.40907; g[3][2]=  0.85679; g[3][3]=  0.19906; g[3][4]= -0.00984;
  ma[4]= 1.82206; g[4][0]= 0.18171; g[4][1]= -0.17558; g[4][2]= -0.79658; g[4][3]= -0.00355; g[4][4]=  0.22358;
  
	// Hadronic parameters for tag modes: 0=no-specified, 1=Kpi, 2=Kpipi0, 3=K3pi
  rd[0] = 0.0;
  rd[1] = 0.0586;
  rd[2] = 0.0440;
  rd[3] = 0.0546;
  deltad[0] = 0.0;
  deltad[1] = 194.7*pi180inv;
  deltad[2] = 196.0*pi180inv;
  deltad[3] = 167.0*pi180inv;
  Rf[0] = 0.0;
  Rf[1] = 1.0;
  Rf[2] = 0.78;
  Rf[3] = 0.52;

	return;
}                                                                                                

complex<double> D0ToKLpipi2018::Amp_PFT(vector<double> k0l, vector<double> pip, vector<double> pim) {

	vector<double> pD;pD.clear();
	if(k0l.size()!=4||pip.size()!=4||pim.size()!=4)cout<<"ERROR in KSPIPI daughter 4 momentum"<<endl;
	for(int i=0;i<k0l.size();i++){
		pD.push_back(k0l[i] + pip[i] + pim[i]);
	}

	//EvtResonance2 DK2piRes0(pD, pip, pim, ar[0], phir[0], width_R[0], mass_R[0], spin[0]);  //ar, phir, width, mass, spin   //Rho770      
	//EvtResonance2 DK2piRes1(pD, pip, pim, ar[1], phir[1], width_R[1], mass_R[1], spin[1]);  //ar, phir, width, mass, spin   //Omega782       
	//EvtResonance2 DK2piRes2(pD, pip, pim, ar[2], phir[2], width_R[2], mass_R[2], spin[2]);  //ar, phir, width, mass, spin   //ftwo1270       
	//EvtResonance2 DK2piRes3(pD, pip, pim, ar[3], phir[3], width_R[3], mass_R[3], spin[3]);  //ar, phir, width, mass, spin   //Rho1450        
	//EvtResonance2 DK2piRes4(pD, k0l, pim, ar[4], phir[4], width_R[4], mass_R[4], spin[4]);  //ar, phir, width, mass, spin   //Kstar892minus  
	//EvtResonance2 DK2piRes5(pD, k0l, pim, ar[5], phir[5], width_R[5], mass_R[5], spin[5]);  //ar, phir, width, mass, spin   //K2star1430minus
	//EvtResonance2 DK2piRes6(pD, k0l, pim, ar[6], phir[6], width_R[6], mass_R[6], spin[6]);  //ar, phir, width, mass, spin   //Kstar1680minus 
	//EvtResonance2 DK2piRes7(pD, k0l, pim, ar[7], phir[7], width_R[7], mass_R[7], spin[7]);  //ar, phir, width, mass, spin   //Kstar1410minus 
	//EvtResonance2 DK2piRes8(pD, k0l, pip, ar[8], phir[8], width_R[8], mass_R[8], spin[8]);  //ar, phir, width, mass, spin   //Kstar892plus   
	//EvtResonance2 DK2piRes9(pD, k0l, pip, ar[9], phir[9], width_R[9], mass_R[9], spin[9]);  //ar, phir, width, mass, spin   //K2star1430plus 
	//EvtResonance2 DK2piRes10(pD, k0l, pip, ar[10], phir[10], width_R[10], mass_R[10], spin[10]);  //ar, phir, width, mass, spin //Kstar1410plus
  complex<double> DK2piRes0  = Resonance2(pD, pip, pim, ar[0],  phir[0],  width_R[0],  mass_R[0],  spin_R[0]);  //ar, phir, width, mass, spin Rho770
  complex<double> DK2piRes1  = Resonance2(pD, pip, pim, ar[1],  phir[1],  width_R[1],  mass_R[1],  spin_R[1]);  //ar, phir, width, mass, spin Omega782
  complex<double> DK2piRes2  = Resonance2(pD, pip, pim, ar[2],  phir[2],  width_R[2],  mass_R[2],  spin_R[2]);  //ar, phir, width, mass, spin ftwo1270
  complex<double> DK2piRes3  = Resonance2(pD, pip, pim, ar[3],  phir[3],  width_R[3],  mass_R[3],  spin_R[3]);  //ar, phir, width, mass, spin Rho1450
  complex<double> DK2piRes4  = Resonance2(pD, k0l, pim, ar[4],  phir[4],  width_R[4],  mass_R[4],  spin_R[4]);  //ar, phir, width, mass, spin Kstar892-
  complex<double> DK2piRes5  = Resonance2(pD, k0l, pim, ar[5],  phir[5],  width_R[5],  mass_R[5],  spin_R[5]);  //ar, phir, width, mass, spin K2star1430-
  complex<double> DK2piRes6  = Resonance2(pD, k0l, pim, ar[6],  phir[6],  width_R[6],  mass_R[6],  spin_R[6]);  //ar, phir, width, mass, spin Kstar1680-
  complex<double> DK2piRes7  = Resonance2(pD, k0l, pim, ar[7],  phir[7],  width_R[7],  mass_R[7],  spin_R[7]);  //ar, phir, width, mass, spin Kstar1410-
  complex<double> DK2piRes8  = Resonance2(pD, k0l, pip, ar[8],  phir[8],  width_R[8],  mass_R[8],  spin_R[8]);  //ar, phir, width, mass, spin Kstar892+
  complex<double> DK2piRes9  = Resonance2(pD, k0l, pip, ar[9],  phir[9],  width_R[9],  mass_R[9],  spin_R[9]);  //ar, phir, width, mass, spin K2star1430+
  complex<double> DK2piRes10 = Resonance2(pD, k0l, pip, ar[10], phir[10], width_R[10], mass_R[10], spin_R[10]); //ar, phir, width, mass, spin Kstar1410+

	complex<double> pipi_s_wave = K_matrix(pip, pim) ;
	if(pipi_s_wave == complex<double>(9999., 9999.)) return 1e-20 ;

	complex<double> kpi_s_wave = amplitude_LASS(k0l, pip, pim, "k0lpim", ar[11], phir[11]*pi180inv) ;
	//complex<double> kpi_s_wave_dcs = amplitude_LASS(k0l, pip, pim, "k0spip", ar[12], phir[12]*pi180inv) ; should be there but not observed yet GUESS

	complex<double> TOT_PFT_AMP   = DK2piRes0 * CP_mult[0]  
		+ DK2piRes1 * CP_mult[1]
		+ DK2piRes2 * CP_mult[2]
		+ DK2piRes3 * CP_mult[3]
		+ DK2piRes4
		+ DK2piRes5
		+ DK2piRes6
		+ DK2piRes7
		+ DK2piRes8 * (-1.)
		+ DK2piRes9 * (-1.)
		+ DK2piRes10* (-1.) 
		+ pipi_s_wave * CP_mult[4]
		+ kpi_s_wave ;


	return TOT_PFT_AMP;

}

complex<double> D0ToKLpipi2018::Resonance2(vector<double> p4_p, vector<double> p4_d1, vector<double> p4_d2, double mag, double theta, double gamma, double bwm, int spin) {

  complex<double> ampl;

  //EvtVector4R  p4_d3 = p4_p - p4_d1 - p4_d2;
	TLorentzVector _p4_p;_p4_p.SetX(p4_p[0]);_p4_p.SetY(p4_p[1]);_p4_p.SetZ(p4_p[2]);_p4_p.SetT(p4_p[3]); 
	TLorentzVector _p4_d1;_p4_d1.SetX(p4_d1[0]);_p4_d1.SetY(p4_d1[1]);_p4_d1.SetZ(p4_d1[2]);_p4_d1.SetT(p4_d1[3]); 
	TLorentzVector _p4_d2;_p4_d2.SetX(p4_d2[0]);_p4_d2.SetY(p4_d2[1]);_p4_d2.SetZ(p4_d2[2]);_p4_d2.SetT(p4_d2[3]); 
	TLorentzVector _p4_d3=_p4_p-_p4_d1-_p4_d2;


  double mAB= (_p4_d1 + _p4_d2).Mag();
  double mBC= (_p4_d2 + _p4_d3).Mag();
  double mAC= (_p4_d1 + _p4_d3).Mag();
  double mA = _p4_d1.Mag();
  double mB = _p4_d2.Mag();
  double mD = _p4_p.Mag();
  double mC = _p4_d3.Mag();


  double mR = bwm;
  double gammaR = gamma;
  double pAB = sqrt( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mAB*mAB));
  double pR  = sqrt( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0)     - mA*mA*mB*mB)/(mR*mR));

  double pD= (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) - mR*mR*mC*mC)/(mD*mD);
  if ( pD>0 ) { pD = sqrt(pD); }
  else        { pD = 0;}
  double pDAB=sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) - mAB*mAB*mC*mC)/(mD*mD));
  double fR = 1;
  double fD = 1;
  int power = 0;
  switch (spin) {
    case 0:
      fR = 1.0;
      fD = 1.0;
      power = 1;
      break;
    case 1:
      fR = sqrt(1.0+1.5*1.5*pR*pR)/sqrt(1.0+1.5*1.5*pAB*pAB);
      fD = sqrt(1.0+5.0*5.0*pD*pD)/sqrt(1.0+5.0*5.0*pDAB*pDAB);
      power = 3;
      break;
    case 2:
      fR = sqrt( (9+3*pow((1.5*pR),2)+pow((1.5*pR),4))/(9+3*pow((1.5*pAB),2) +pow((1.5*pAB) ,4)) );
      fD = sqrt( (9+3*pow((5.0*pD),2)+pow((5.0*pD),4))/(9+3*pow((5.0*pDAB),2)+pow((5.0*pDAB),4)) );
      power = 5;
      break;
    default:
      cout << "Incorrect spin in D0ToKLpipi2018::EvtResonance2.cc\n" <<endl;
  }

  double gammaAB= gammaR*pow(pAB/pR,power)*(mR/mAB)*fR*fR;
  switch (spin) {
    case 0:
      ampl=mag*complex<double>(cos(theta*pi180inv),sin(theta*pi180inv))*fR*fD/(mR*mR-mAB*mAB-complex<double>(0.0,mR*gammaAB));
      break;
    case 1:
      ampl=mag*complex<double>(cos(theta*pi180inv),sin(theta*pi180inv))*
           (fR*fD*(mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mAB*mAB)))/(mR*mR-mAB*mAB-complex<double>(0.0,mR*gammaAB)));
      break;
    case 2:
      ampl=mag*complex<double>(cos(theta*pi180inv),sin(theta*pi180inv))*
        (fR*fD/(mR*mR-mAB*mAB-complex<double>(0.0,mR*gammaAB)))*
        (pow((mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB)),2)-
        (1.0/3.0)*(mAB*mAB-2*mD*mD-2*mC*mC+pow((mD*mD- mC*mC)/mAB, 2))*
        (mAB*mAB-2*mA*mA-2*mB*mB+pow((mA*mA-mB*mB)/mAB,2)));
      break;
    default:
     cout << "Incorrect spin in D0ToKSpipi::Resonance2.cc\n" <<endl;
  }

  return ampl;
}

complex<double> D0ToKLpipi2018::K_matrix(vector<double> p_pip, vector<double> p_pim) { 

	//double pi180inv = 1.0/EvtConst::radToDegrees;
	bool reject=false;

	const double mD0 = 1.86483;
	const double mKl = 0.49761;
	const double mPi = 0.13957;

	TLorentzVector _p_pip(p_pip[0],p_pip[1],p_pip[2],p_pip[3]); 
	TLorentzVector _p_pim(p_pim[0],p_pim[1],p_pim[2],p_pim[3]); 

  double mAB = (_p_pip + _p_pim).Mag() ;

	double s = mAB*mAB;

	// Define the complex coupling constants
	//Double_t g[5][5]; // g[Physical pole]Decay channel]

	// pi+pi- channel
	//g[0][0]=0.22889;
	//g[1][0]=0.94128;
	//g[2][0]=0.36856;
	//g[3][0]=0.33650;
	//g[4][0]=0.18171;

	//// K+K- channel
	//g[0][1]=-0.55377;
	//g[1][1]=0.55095;
	//g[2][1]=0.23888;
	//g[3][1]=0.40907;
	//g[4][1]=-0.17558;

	//// 4pi channel
	//g[0][2]=0;
	//g[1][2]=0;
	//g[2][2]=0.55639;
	//g[3][2]=0.85679;
	//g[4][2]=-0.79658;

	//// eta eta channel
	//g[0][3]=-0.39899;
	//g[1][3]=0.39065;
	//g[2][3]=0.18340;
	//g[3][3]=0.19906;
	//g[4][3]=-0.00355;

	////eta eta' channel
	//g[0][4]=-0.34639;
	//g[1][4]=0.31503;
	//g[2][4]=0.18681;
	//g[3][4]=-0.00984;
	//g[4][4]=0.22358;

	// Define masses of the physical poles (in GeV)
	//Double_t ma[5];

	//ma[0]=0.651;
	//ma[1]=1.20360;
	//ma[2]=1.55817;
	//ma[3]=1.21000;
	//ma[4]=1.82206;

	// Define variables
	complex<double> n11,n12,n13,n14,n15,n21,n22,n23,n24,n25,n31,n32,n33,n34,n35,n41,n42,n43,n44,n45,n51,n52,n53,n54,n55;
	double rho1sq,rho2sq,rho4sq,rho5sq;
	complex<double> rho1,rho2,rho3,rho4,rho5;
	complex<double> rho[5];
	complex<double> pole,SVT,Adler;
	complex<double> det;
	complex<double> i[5][5];
	double f[5][5];

	// pi+, K+, eta, and eta' PDG masses
	double mpi=0.13957;
	double mK=0.493677;
	double meta=0.54775;
	double metap=0.95778;

	// Init matrices and vectors with zeros
	complex<double> K[5][5];
	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			i[k][l]=complex<double>(0.,0.);
			K[k][l]=complex<double>(0.,0.);
			f[k][l]=0.;
		}
		rho[k]=0.;
	}

	// Fill scattering data values
	Double_t s_scatt=-3.92637;
	Double_t sa=1.0;
	Double_t sa_0=-0.15;

	// f_scattering
	f[0][0]=0.23399;
	f[0][1]=0.15044;
	f[0][2]=-0.20545;
	f[0][3]=0.32825;
	f[0][4]=0.35412;

	f[1][0]=f[0][1];
	f[2][0]=f[0][2];
	f[3][0]=f[0][3];
	f[4][0]=f[0][4];

	// Compute phase space factors
	rho1sq=(1.0-(pow((mpi+mpi),2)/s));
	if(rho1sq >=0.) {
		rho1=complex<double>(sqrt(rho1sq),0.);
	}
	else{
		rho1=complex<double>(0.,sqrt(-rho1sq));
	}
	rho[0]=rho1;

	rho2sq=(1.0-(pow((mK+mK),2)/s));
	if(rho2sq >=0.) {
		rho2=complex<double>(sqrt(rho2sq),0.);
	}
	else{
		rho2=complex<double>(0.,sqrt(-rho2sq));
	}

	rho[1]=rho2;

	rho3=complex<double>(0.,0.);

	if(s<=1) {
		Double_t real = 1.2274+0.00370909/(s*s) - (0.111203)/(s) - 6.39017*s +16.8358*s*s - 21.8845*s*s*s + 11.3153*s*s*s*s;
		Double_t cont32=sqrt(1.0-(16.0*mpi*mpi));
		rho3=complex<double>(cont32*real,0.);
	}
	else{
		rho3=complex<double>(sqrt(1.0-(16.0*mpi*mpi/s)),0.);
	}
	rho[2]=rho3;

	rho4sq=(1.0-(pow((meta+meta),2)/s));
	if(rho4sq>=0.) {
		rho4=complex<double>(sqrt(rho4sq),0.);
	}
	else{
		rho4=complex<double>(0.,sqrt(-rho4sq));
	}
	rho[3]=rho4;

	rho5sq=(1.0-(pow((meta+metap),2)/s));
	if(rho5sq >=0.) {
		rho5=complex<double>(sqrt(rho5sq),0.);
	}
	else{
		rho5=complex<double>(0.,sqrt(-rho5sq));
	}
	rho[4]=rho5;

	// Sum over the poles
	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			for (Int_t pole_index=0;pole_index<5;pole_index++) {
				Double_t A=g[pole_index][k]*g[pole_index][l];
				Double_t B=ma[pole_index]*ma[pole_index]-s;
				K[k][l]=K[k][l]+complex<double>(A/B,0.);
			}
		}
	}

	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			Double_t C=f[k][l]*(1.0-s_scatt);
			Double_t D=(s-s_scatt);
			K[k][l]=K[k][l]+complex<double>(C/D,0.);
		}
	}

	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			Double_t E=(s-(sa*mpi*mpi*0.5))*(1.0-sa_0);
			Double_t F=(s-sa_0);
			K[k][l]=K[k][l]*complex<double>(E/F,0.);
		}
	}

	n11=complex<double>(1.,0.)-complex<double>(0.,1.)*K[0][0]*rho[0];
	n12=complex<double>(0.,0.)-complex<double>(0.,1.)*K[0][1]*rho[1];
	n13=complex<double>(0.,0.)-complex<double>(0.,1.)*K[0][2]*rho[2];
	n14=complex<double>(0.,0.)-complex<double>(0.,1.)*K[0][3]*rho[3];
	n15=complex<double>(0.,0.)-complex<double>(0.,1.)*K[0][4]*rho[4];

	n21=complex<double>(0.,0.)-complex<double>(0.,1.)*K[1][0]*rho[0];
	n22=complex<double>(1.,0.)-complex<double>(0.,1.)*K[1][1]*rho[1];
	n23=complex<double>(0.,0.)-complex<double>(0.,1.)*K[1][2]*rho[2];
	n24=complex<double>(0.,0.)-complex<double>(0.,1.)*K[1][3]*rho[3];
	n25=complex<double>(0.,0.)-complex<double>(0.,1.)*K[1][4]*rho[4];

	n31=complex<double>(0.,0.)-complex<double>(0.,1.)*K[2][0]*rho[0];
	n32=complex<double>(0.,0.)-complex<double>(0.,1.)*K[2][1]*rho[1];
	n33=complex<double>(1.,0.)-complex<double>(0.,1.)*K[2][2]*rho[2];
	n34=complex<double>(0.,0.)-complex<double>(0.,1.)*K[2][3]*rho[3];
	n35=complex<double>(0.,0.)-complex<double>(0.,1.)*K[2][4]*rho[4];

	n41=complex<double>(0.,0.)-complex<double>(0.,1.)*K[3][0]*rho[0];
	n42=complex<double>(0.,0.)-complex<double>(0.,1.)*K[3][1]*rho[1];
	n43=complex<double>(0.,0.)-complex<double>(0.,1.)*K[3][2]*rho[2];
	n44=complex<double>(1.,0.)-complex<double>(0.,1.)*K[3][3]*rho[3];
	n45=complex<double>(0.,0.)-complex<double>(0.,1.)*K[3][4]*rho[4];

	n51=complex<double>(0.,0.)-complex<double>(0.,1.)*K[4][0]*rho[0];
	n52=complex<double>(0.,0.)-complex<double>(0.,1.)*K[4][1]*rho[1];
	n53=complex<double>(0.,0.)-complex<double>(0.,1.)*K[4][2]*rho[2];
	n54=complex<double>(0.,0.)-complex<double>(0.,1.)*K[4][3]*rho[3];
	n55=complex<double>(1.,0.)-complex<double>(0.,1.)*K[4][4]*rho[4];

	// Compute the determinant
	det = (n15*n24*n33*n42*n51 - n14*n25*n33*n42*n51 - n15*n23*n34*n42*n51 +
			n13*n25*n34*n42*n51 + n14*n23*n35*n42*n51 - n13*n24*n35*n42*n51 -
			n15*n24*n32*n43*n51 + n14*n25*n32*n43*n51 + n15*n22*n34*n43*n51 -
			n12*n25*n34*n43*n51 - n14*n22*n35*n43*n51 + n12*n24*n35*n43*n51 +
			n15*n23*n32*n44*n51 - n13*n25*n32*n44*n51 - n15*n22*n33*n44*n51 +
			n12*n25*n33*n44*n51 + n13*n22*n35*n44*n51 - n12*n23*n35*n44*n51 -
			n14*n23*n32*n45*n51 + n13*n24*n32*n45*n51 + n14*n22*n33*n45*n51 -
			n12*n24*n33*n45*n51 - n13*n22*n34*n45*n51 + n12*n23*n34*n45*n51 -
			n15*n24*n33*n41*n52 + n14*n25*n33*n41*n52 + n15*n23*n34*n41*n52 -
			n13*n25*n34*n41*n52 - n14*n23*n35*n41*n52 + n13*n24*n35*n41*n52 +
			n15*n24*n31*n43*n52 - n14*n25*n31*n43*n52 - n15*n21*n34*n43*n52 +
			n11*n25*n34*n43*n52 + n14*n21*n35*n43*n52 - n11*n24*n35*n43*n52 -
			n15*n23*n31*n44*n52 + n13*n25*n31*n44*n52 + n15*n21*n33*n44*n52 -
			n11*n25*n33*n44*n52 - n13*n21*n35*n44*n52 + n11*n23*n35*n44*n52 +
			n14*n23*n31*n45*n52 - n13*n24*n31*n45*n52 - n14*n21*n33*n45*n52 +
			n11*n24*n33*n45*n52 + n13*n21*n34*n45*n52 - n11*n23*n34*n45*n52 +
			n15*n24*n32*n41*n53 - n14*n25*n32*n41*n53 - n15*n22*n34*n41*n53 +
			n12*n25*n34*n41*n53 + n14*n22*n35*n41*n53 - n12*n24*n35*n41*n53 -
			n15*n24*n31*n42*n53 + n14*n25*n31*n42*n53 + n15*n21*n34*n42*n53 -
			n11*n25*n34*n42*n53 - n14*n21*n35*n42*n53 + n11*n24*n35*n42*n53 +
			n15*n22*n31*n44*n53 - n12*n25*n31*n44*n53 - n15*n21*n32*n44*n53 +
			n11*n25*n32*n44*n53 + n12*n21*n35*n44*n53 - n11*n22*n35*n44*n53 -
			n14*n22*n31*n45*n53 + n12*n24*n31*n45*n53 + n14*n21*n32*n45*n53 -
			n11*n24*n32*n45*n53 - n12*n21*n34*n45*n53 + n11*n22*n34*n45*n53 -
			n15*n23*n32*n41*n54 + n13*n25*n32*n41*n54 + n15*n22*n33*n41*n54 -
			n12*n25*n33*n41*n54 - n13*n22*n35*n41*n54 + n12*n23*n35*n41*n54 +
			n15*n23*n31*n42*n54 - n13*n25*n31*n42*n54 - n15*n21*n33*n42*n54 +
			n11*n25*n33*n42*n54 + n13*n21*n35*n42*n54 - n11*n23*n35*n42*n54 -
			n15*n22*n31*n43*n54 + n12*n25*n31*n43*n54 + n15*n21*n32*n43*n54 -
			n11*n25*n32*n43*n54 - n12*n21*n35*n43*n54 + n11*n22*n35*n43*n54 +
			n13*n22*n31*n45*n54 - n12*n23*n31*n45*n54 - n13*n21*n32*n45*n54 +
			n11*n23*n32*n45*n54 + n12*n21*n33*n45*n54 - n11*n22*n33*n45*n54 +
			n14*n23*n32*n41*n55 - n13*n24*n32*n41*n55 - n14*n22*n33*n41*n55 +
			n12*n24*n33*n41*n55 + n13*n22*n34*n41*n55 - n12*n23*n34*n41*n55 -
			n14*n23*n31*n42*n55 + n13*n24*n31*n42*n55 + n14*n21*n33*n42*n55 -
			n11*n24*n33*n42*n55 - n13*n21*n34*n42*n55 + n11*n23*n34*n42*n55 +
			n14*n22*n31*n43*n55 - n12*n24*n31*n43*n55 - n14*n21*n32*n43*n55 +
			n11*n24*n32*n43*n55 + n12*n21*n34*n43*n55 - n11*n22*n34*n43*n55 -
			n13*n22*n31*n44*n55 + n12*n23*n31*n44*n55 + n13*n21*n32*n44*n55 -
			n11*n23*n32*n44*n55 - n12*n21*n33*n44*n55 + n11*n22*n33*n44*n55);

	if(det == complex<double>(0., 0.)) reject=true;

	// The 1st row of the inverse matrix {(I-iKp)^-1}_0j
	i[0][0] = (n25*n34*n43*n52 -
			n24*n35*n43*n52 - n25*n33*n44*n52 + n23*n35*n44*n52 +
			n24*n33*n45*n52 - n23*n34*n45*n52 - n25*n34*n42*n53 +
			n24*n35*n42*n53 + n25*n32*n44*n53 - n22*n35*n44*n53 -
			n24*n32*n45*n53 + n22*n34*n45*n53 + n25*n33*n42*n54 -
			n23*n35*n42*n54 - n25*n32*n43*n54 + n22*n35*n43*n54 +
			n23*n32*n45*n54 - n22*n33*n45*n54 - n24*n33*n42*n55 +
			n23*n34*n42*n55 + n24*n32*n43*n55 - n22*n34*n43*n55 -
			n23*n32*n44*n55 + n22*n33*n44*n55)/det;

	i[0][1] = (-n15*n34*n43*n52 +
			n14*n35*n43*n52 + n15*n33*n44*n52 - n13*n35*n44*n52 -
			n14*n33*n45*n52 + n13*n34*n45*n52 + n15*n34*n42*n53 -
			n14*n35*n42*n53 - n15*n32*n44*n53 + n12*n35*n44*n53 +
			n14*n32*n45*n53 - n12*n34*n45*n53 - n15*n33*n42*n54 +
			n13*n35*n42*n54 + n15*n32*n43*n54 - n12*n35*n43*n54 -
			n13*n32*n45*n54 + n12*n33*n45*n54 + n14*n33*n42*n55 -
			n13*n34*n42*n55 - n14*n32*n43*n55 + n12*n34*n43*n55 +
			n13*n32*n44*n55 - n12*n33*n44*n55)/det;

	i[0][2] = (n15*n24*n43*n52 -
			n14*n25*n43*n52 - n15*n23*n44*n52 + n13*n25*n44*n52 +
			n14*n23*n45*n52 - n13*n24*n45*n52 - n15*n24*n42*n53 +
			n14*n25*n42*n53 + n15*n22*n44*n53 - n12*n25*n44*n53 -
			n14*n22*n45*n53 + n12*n24*n45*n53 + n15*n23*n42*n54 -
			n13*n25*n42*n54 - n15*n22*n43*n54 + n12*n25*n43*n54 +
			n13*n22*n45*n54 - n12*n23*n45*n54 - n14*n23*n42*n55 +
			n13*n24*n42*n55 + n14*n22*n43*n55 - n12*n24*n43*n55 -
			n13*n22*n44*n55 + n12*n23*n44*n55)/det;

	i[0][3] = (-n15*n24*n33*n52 +
			n14*n25*n33*n52 + n15*n23*n34*n52 - n13*n25*n34*n52 -
			n14*n23*n35*n52 + n13*n24*n35*n52 + n15*n24*n32*n53 -
			n14*n25*n32*n53 - n15*n22*n34*n53 + n12*n25*n34*n53 +
			n14*n22*n35*n53 - n12*n24*n35*n53 - n15*n23*n32*n54 +
			n13*n25*n32*n54 + n15*n22*n33*n54 - n12*n25*n33*n54 -
			n13*n22*n35*n54 + n12*n23*n35*n54 + n14*n23*n32*n55 -
			n13*n24*n32*n55 - n14*n22*n33*n55 + n12*n24*n33*n55 +
			n13*n22*n34*n55 - n12*n23*n34*n55)/det;

	i[0][4] = (n15*n24*n33*n42 -
			n14*n25*n33*n42 - n15*n23*n34*n42 + n13*n25*n34*n42 +
			n14*n23*n35*n42 - n13*n24*n35*n42 - n15*n24*n32*n43 +
			n14*n25*n32*n43 + n15*n22*n34*n43 - n12*n25*n34*n43 -
			n14*n22*n35*n43 + n12*n24*n35*n43 + n15*n23*n32*n44 -
			n13*n25*n32*n44 - n15*n22*n33*n44 + n12*n25*n33*n44 +
			n13*n22*n35*n44 - n12*n23*n35*n44 - n14*n23*n32*n45 +
			n13*n24*n32*n45 + n14*n22*n33*n45 - n12*n24*n33*n45 -
			n13*n22*n34*n45 + n12*n23*n34*n45)/det;


	//--------------------------------------------------------------------------------------------------------------------------------
	double s0_prod = -0.07;

	double u1j_re_limit[5][2], u1j_im_limit[5][2] ;
	u1j_re_limit[0][0] =  0.   ; u1j_re_limit[0][1] = 1.    ;
	u1j_re_limit[1][0] = -0.29 ; u1j_re_limit[1][1] = 0.12  ;
	u1j_re_limit[2][0] = -0.17 ; u1j_re_limit[2][1] = 0.065 ;
	u1j_re_limit[3][0] = -0.66 ; u1j_re_limit[3][1] = 0.1   ;
	u1j_re_limit[4][0] = -1.36 ; u1j_re_limit[4][1] = 0.18  ;

	u1j_im_limit[0][0] = -0.58  ; u1j_im_limit[0][1] = 0.58 ;
	u1j_im_limit[1][0] =  0.00  ; u1j_im_limit[1][1] = 0.28 ;
	u1j_im_limit[2][0] = -0.135 ; u1j_im_limit[2][1] = 0.10 ;
	u1j_im_limit[3][0] = -0.13  ; u1j_im_limit[3][1] = 0.40 ;
	u1j_im_limit[4][0] = -0.36  ; u1j_im_limit[4][1] = 0.80 ;

	//for(int kk=0 ; kk<5 ; kk++)
	//{
	//        cout<<"beta["<<kk<<"]: "<<abs(beta[kk])<<", \t"<<arg(beta[kk])*EvtConst::radToDegrees<<endl;
	//        cout<<"fprod["<<kk<<"]: "<<abs(fprod[kk])<<", \t"<<arg(fprod[kk])*EvtConst::radToDegrees<<endl;
	//}

	complex<double> value0(0., 0.) ;
	complex<double> value1(0., 0.) ;

	for(int l=0;l<5;l++)
	{
		double u1j_re = real(i[0][l]) ;
		double u1j_im = imag(i[0][l]) ;
		if(u1j_re==0. || u1j_im==0.) reject=true;
		if(u1j_re < u1j_re_limit[l][0] || u1j_re > u1j_re_limit[l][1] || u1j_im < u1j_im_limit[l][0] || u1j_im > u1j_im_limit[l][1]) reject=true;

		for(int pole_index=0;pole_index<5;pole_index++)
		{
			complex<double> A = beta[pole_index]*g[pole_index][l];
			value0 += (i[0][l]*A)/(ma[pole_index]*ma[pole_index]-s) ;
			//cout<<"value0["<<l<<"]["<<pole_index<<"] = "<<value0<<endl;
		}
	}

	value1 += i[0][0]*fprod[0];
	value1 += i[0][1]*fprod[1];
	value1 += i[0][2]*fprod[2];
	value1 += i[0][3]*fprod[3];
	value1 += i[0][4]*fprod[4];
	value1 *= (1.-s0_prod)/(s-s0_prod) ;
	//cout<<"value1 = "<<value1<<endl;

	if(reject==true) return complex<double>(9999., 9999.) ;
	else return (value0+value1) ;
	//return (value0+value1) ;

}

complex<double> D0ToKLpipi2018::amplitude_LASS(vector<double> p_k0l, vector<double> p_pip, vector<double> p_pim, string reso, double A_r, double Phi_r) {

	double mR = 1.425 ;
	double gammaR = 0.27 ;
	double mab2=0.0;
	TLorentzVector _p_k0l(p_k0l[0],p_k0l[1],p_k0l[2],p_k0l[3]);
	TLorentzVector _p_pip(p_pip[0],p_pip[1],p_pip[2],p_pip[3]);
	TLorentzVector _p_pim(p_pim[0],p_pim[1],p_pim[2],p_pim[3]);
  if     (reso == "k0lpim") mab2 = pow((_p_k0l + _p_pim).Mag(),2);
  else if(reso == "k0lpip") mab2 = pow((_p_k0l + _p_pip).Mag(),2);

	double s = mab2;

	const double mD0 = 1.86483;
	const double mKl = 0.49761;
	const double mPi = 0.13957;

	double _a =     0.113   ;     
	double _r =     -33.8   ;
	double _R =     1.0     ;
	double _F =     0.96    ;
	double _phiR =  -1.9146 ;
	double _phiF =  0.0017  ;

	double fR=1.0; // K*0(1430) has spin zero
	int power=1; // Power is 1 for spin zero

	double mAB = sqrt(mab2) ; // (_p4_d1+_p4_d2).mass();

	double mA=mKl; // _p4_d1.mass();
	double mB=mPi; // _p4_d2.mass();
	double mC=mPi;
	double mD = mD0;

	double pAB=sqrt( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mAB*mAB));
	double q=pAB;

	double pR=sqrt( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mR*mR));
	double q0=pR;

	// Running width.
	double g = gammaR*pow(q/q0,power)*(mR/mAB)*fR*fR;

	complex<double> propagator_relativistic_BreitWigner = 1./(mR*mR - mAB*mAB - complex<double>(0.,mR*g));

	// Non-resonant phase shift
	double cot_deltaF = 1.0/(_a*q) + 0.5*_r*q;
	double qcot_deltaF = 1.0/_a + 0.5*_r*q*q;

	// Compute resonant part
	complex<double> expi2deltaF = complex<double>(qcot_deltaF, q)/ complex<double>(qcot_deltaF, -q);

	complex<double> resonant_term_T = _R * complex<double>(cos(_phiR + 2 * _phiF), sin(_phiR + 2 * _phiF)) * propagator_relativistic_BreitWigner * mR * gammaR * mR / q0 * expi2deltaF;

	// Compute non-resonant part
	complex<double> non_resonant_term_F = _F * complex<double>(cos(_phiF), sin(_phiF)) * (cos(_phiF) + cot_deltaF * sin(_phiF)) * sqrt(s) / complex<double>(qcot_deltaF, -q);

	// Add non-resonant and resonant terms
	complex<double> LASS_contribution = non_resonant_term_F + resonant_term_T;

	return complex<double>(A_r*cos(Phi_r), A_r*sin(Phi_r)) * LASS_contribution;

}
