#ifndef D0TO2KLPIPI2018_H
#define D0TO2KLPIPI2018_H

#include <vector>
#include <complex>

using namespace std;

class D0ToKLpipi2018{

public:
  
  D0ToKLpipi2018() {}
  virtual ~D0ToKLpipi2018();

  void init();

  complex<double> Amp_PFT(vector<double> k0l, vector<double> pip, vector<double> pim);

protected:


private:

  complex<double> K_matrix(vector<double> p_pip, vector<double> p_pim);
  complex<double> amplitude_LASS(vector<double> p_k0l, vector<double> p_pip, vector<double> p_pim, std::string reso, double A_r, double Phi_r);
  complex<double> Resonance2(vector<double> p4_p, vector<double> p4_d1, const vector<double> p4_d2, double mag, double theta, double gamma, double bwm, int spin);

  int _nd;

  float ar[12], phir[12];
  //vector < complex<double> > CP_mult, beta, fprod;
  complex<double> CP_mult[5], beta[5], fprod[5];
  double tan2thetaC;
  double pi180inv;
  double mass_R[12], width_R[12];
  int spin_R[12];
  double frac1[3], frac2[3], frac3[3];
  double rd[4], deltad[4], Rf[4];
  double ma[5], g[5][5];  // Kmatrix_couplings
};

#endif

