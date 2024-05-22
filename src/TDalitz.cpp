#include "AmpGen/TDalitz.h"
#include "TMath.h"
#include <iostream>

dcomplex pol(double a, double phi) {
  double radphi = gPi/180.*phi;
  return dcomplex(a*cos(radphi), a*sin(radphi));
}

dcomplex breitWigner(double s, double m, double gamma) {
  return 1./dcomplex(m*m-s,-m*gamma);
}

dcomplex gounarisSakurai(double s, double m, double gamma) {
  double m_pi = 139.57017e-3;

  double m2 = m*m;
  double m_pi2 = m_pi*m_pi;
  double ss = sqrt(s);

  double ppi2 = (s-4.*m_pi2)/4.;
  double p02 = (m2-4.*m_pi2)/4.;
  double p0 = sqrt(p02);
  double ppi = sqrt(ppi2);

  double hs = 2.*ppi/gPi/ss*log((ss+2.*ppi)/2./m_pi);
  double hm = 2.*p0/gPi/m*log((m+2.*ppi)/2./m_pi);

  double dhdq = hm*(1./8./p02 - 1./2./m2) + 1./2./gPi/m2;
  double f = gamma*m2/pow(p0,3)*(ppi2*(hs-hm) - p02*(s-m2)*dhdq);

  double gamma_s = gamma*m2*pow(ppi,3)/s/pow(p0,3);

  double dr = m2-s+f;
  double di = ss*gamma_s;

  double r = dr/(dr*dr+di*di);
  double i = di/(dr*dr+di*di);

//  printf("%f %f %f %f %f %f\n", r, i, dr, di, f, gamma_s);

  return dcomplex(r,i);
}

double log_poisson(double nu, double x) {
  if (x<0 || nu<0) return -1e40;
  double arg = x+1;
  return (x*TMath::Log(nu)-nu-TMath::LnGamma(arg));
}

TDalitz::TDalitz(double md, double ma, double mb, double mc) {
  fMa = ma;
  fMb = mb;
  fMc = mc;
  fMd = md;
  fM2a = ma*ma;
  fM2b = mb*mb;
  fM2c = mc*mc;
  fM2d = md*md;
  fM2sum = fM2a + fM2b + fM2c + fM2d;
  
  fWidthDep = 1;
  fMabForVectors = 1;
  fMabForTensors = 1;
  
  fR2r = 1.5*1.5;
  fR2d = 5.0*5.0;
}

dcomplex TDalitz::scalar_ampl(double x, double y, double m, double gamma, int chan) {
  double ma, mb, mc, m2ab;
  switch(chan) {
    case CH_AB : ma = fMa; mb = fMb; mc = fMc; m2ab = x; break;
    case CH_AC : ma = fMa; mb = fMc; mc = fMb; m2ab = y; break;
    case CH_BC : ma = fMb; mb = fMc; mc = fMa; m2ab = fM2sum-x-y; break;
    default: return(0.);
  }
  double m2r = pow(m,2);
  double m2sumab = pow(ma+mb,2);
  double m2difab = pow(ma-mb,2);
  double p2res1 = (m2r-m2sumab)*(m2r-m2difab);
  double p2ab1 = (m2ab-m2sumab)*(m2ab-m2difab);
  double gamma2;
  if (fWidthDep) {
    gamma2 = gamma*sqrt(p2ab1/p2res1)*(m2r/m2ab);
  } else {
    gamma2 = gamma;
  }
  return breitWigner(m2ab, m, gamma2);
}

dcomplex TDalitz::f0_ampl(double x, double y, double m, int chan) {
  double m2ab;
  switch(chan) {
    case CH_AB : m2ab = x; break;
    case CH_AC : m2ab = y; break;
    case CH_BC : m2ab = fM2sum-x-y; break;
    default: return(0.);
  }

  double m2pi = pow(139.57017e-3,2);
  double m2k = pow(497.61401e-3,2);
  double m2kch = pow(493.677e-3,2);


  double gamma = 0.09*sqrt(m2ab/4. - m2pi);
  
  if (m2ab/4. > m2kch) {
    gamma = gamma + 0.02/2.*sqrt(m2ab/4.-m2kch);
  }
  if (m2ab/4. > m2k) {
    gamma = gamma + 0.02/2.*sqrt(m2ab/4.-m2k);
  }

  return breitWigner(m2ab, m, gamma);
}

dcomplex TDalitz::vector_ampl(double x, double y, double m, double gamma, int chan) {
  double ma, mb, mc;
  double m2ab, m2ac, m2bc;
  double z = fM2sum-x-y;

  switch(chan) {
    case CH_AB : ma = fMa; mb = fMb; mc = fMc; m2ab = x; m2ac = y; m2bc = z; break;
    case CH_AC : ma = fMa; mb = fMc; mc = fMb; m2ab = y; m2ac = x; m2bc = z; break;
    case CH_BC : ma = fMb; mb = fMc; mc = fMa; m2ab = z; m2ac = x; m2bc = y; break;
    default: return(0.);
  }
  double mab = sqrt(m2ab);
  double m2a = ma*ma;
  double m2b = mb*mb;
  double m2c = mc*mc;
  double m2r = pow(m,2);

  double m2div = m2r;
  if (fMabForVectors) m2div = m2ab;

  double num = m2ac-m2bc+(fM2d-m2c)*(m2b-m2a)/m2div;
  double m2sumab = pow(ma+mb,2);
  double m2difab = pow(ma-mb,2);
  double p2res = (m2r-m2sumab)*(m2r-m2difab)/4./m2r;
  double p2ab = (m2ab-m2sumab)*(m2ab-m2difab)/4./m2ab;

  double p2c = (fM2d-pow(m+mc,2))*(fM2d-pow(m-mc,2))/4./fM2d;
  double p2d = (fM2d-pow(mc+mab,2))*(fM2d-pow(mc-mab,2))/4./fM2d;

  double fr = sqrt((1.+fR2r*p2res)/(1.+fR2r*p2ab));
  double fd = sqrt((1.+fR2d*p2c)/(1.+fR2d*p2d));

  double gamma2;
  if (fWidthDep) {
    gamma2 = gamma*pow(sqrt(p2ab/p2res),3)*(m/mab)*fr*fr;
  } else {
    gamma2 = gamma;
  }
  
  if (chan == CH_BC) num =-num;

  return num*breitWigner(m2ab, m, gamma2)*fr*fd;
}

dcomplex TDalitz::gs_ampl(double x, double y, double m, double gamma, int chan) {
  double ma, mb, mc;
  double m2ab, m2ac, m2bc;
  double z = fM2sum-x-y;

  switch(chan) {
    case CH_AB : ma = fMa; mb = fMb; mc = fMc; m2ab = x; m2ac = y; m2bc = z; break;
    case CH_AC : ma = fMa; mb = fMc; mc = fMb; m2ab = y; m2ac = x; m2bc = z; break;
    case CH_BC : ma = fMb; mb = fMc; mc = fMa; m2ab = z; m2ac = x; m2bc = y; break;
    default: return(0.);
  }
  double m2a = ma*ma;
  double m2b = mb*mb;
  double m2c = mc*mc;
  double m2r = m*m;

  double m2div = m2r;
  if (fMabForVectors) m2div = m2ab;

  double num = m2ac-m2bc+(fM2d-m2c)*(m2b-m2a)/m2div;

//  printf("gs: %f %f %f\n", fA.real(), fA.imag(), num);

  if (chan == CH_BC) num =-num;

  return num*gounarisSakurai(m2ab, m, gamma);
}

dcomplex TDalitz::tensor_ampl(double x, double y, double m, double gamma, int chan) {
  double ma, mb, mc;
  double m2ab, m2ac, m2bc;
  double z = fM2sum-x-y;

  switch(chan) {
    case CH_AB : ma = fMa; mb = fMb; mc = fMc; m2ab = x; m2ac = y; m2bc = z; break;
    case CH_AC : ma = fMa; mb = fMc; mc = fMb; m2ab = y; m2ac = x; m2bc = z; break;
    case CH_BC : ma = fMb; mb = fMc; mc = fMa; m2ab = z; m2ac = x; m2bc = y; break;
    default: return(0.);
  }
  double mab = sqrt(m2ab);
  double m2a = ma*ma;
  double m2b = mb*mb;
  double m2c = mc*mc;
  double m2r = pow(m,2);
  
  double m2div = m2r;
  if (fMabForTensors) m2div = m2ab;
  
  double num = pow(m2bc-m2ac+(fM2d-m2c)*(m2a-m2b)/m2div,2)-
               1./3.*(m2ab-2.*(fM2d+m2c)+pow(fM2d-m2c,2)/m2div)*
               (m2ab-2.*(m2a+m2b)+pow(m2a-m2b,2)/m2div);

  double m2sumab = pow(ma+mb,2);
  double m2difab = pow(ma-mb,2);
  double p2res = (m2r-m2sumab)*(m2r-m2difab)/4./m2r;
  double p2ab = (m2ab-m2sumab)*(m2ab-m2difab)/4./m2ab;

  double p2c = (fM2d-pow(m+mc,2))*(fM2d-pow(m-mc,2))/4./fM2d;
  double p2d = (fM2d-pow(mc+mab,2))*(fM2d-pow(mc-mab,2))/4./fM2d;

  double fr = sqrt((9.+3.*fR2r*p2res+pow(fR2r*p2res,2))/
                   (9.+3.*fR2r*p2ab +pow(fR2r*p2ab ,2)));
  double fd = sqrt((9.+3.*fR2d*p2c  +pow(fR2d*p2c  ,2))/
                   (9.+3.*fR2d*p2d  +pow(fR2d*p2d  ,2)));

  double gamma2;
  if (fWidthDep) {
    gamma2 = gamma*pow(sqrt(p2ab/p2res),5)*(m/mab)*fr*fr;
  } else {
    gamma2 = gamma;
  }
  return num*breitWigner(m2ab, m, gamma2)*fr*fd;
}

int TDalitz::kine_limits(double m2ab, double m2ac) {
  double m2bc = fM2sum - m2ab - m2ac;
  double mab = sqrt(m2ab);
  double mac = sqrt(m2ac);
  double mbc = sqrt(m2bc);
  double p2a = 0.25/fM2d*(fM2d-pow(mbc+fMa,2))*(fM2d-pow(mbc-fMa,2));
  double p2b = 0.25/fM2d*(fM2d-pow(mac+fMb,2))*(fM2d-pow(mac-fMb,2));
  double p2c = 0.25/fM2d*(fM2d-pow(mab+fMc,2))*(fM2d-pow(mab-fMc,2));
  double eb = (m2ab-fM2a+fM2b)/2./mab;
  double ec = (fM2d-m2ab-fM2c)/2./mab;
  if (eb<fMb || ec<fMc) return 0;
  double pb = sqrt(eb*eb-fM2b);
  double pc = sqrt(ec*ec-fM2c);
  double e2sum = pow(eb+ec,2);
  double m2bc_max = e2sum-pow(pb-pc,2);
  double m2bc_min = e2sum-pow(pb+pc,2);
  if (m2bc < m2bc_min || m2bc > m2bc_max) return 0;
  if (p2a>0 && p2b>0 && p2c>0) return 1;
  return 0;
}

