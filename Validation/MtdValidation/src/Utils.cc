#include "Validation/MtdValidation/interface/Utils.h"
#include "DataFormats/Math/interface/deltaR.h"
// #include "../interface/ReadParameters.h"
// #include "B0KstarMM/B0KstMuMu/interface/Utils.h"
// #include "B0KstarMM/B0KstMuMu/interface/ReadParameters.h"

#include <TMath.h>
#include <TFile.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

// ####################
// # Global constants #
// ####################
#define YvalueOutsideLimits 20.0 // Value given to bins with zero error in order not to show them

Utils::Utils (bool rightFlavorTag)
{
  muonMass     = 0.10565837;
  pionMass     = 0.13957018;
  kaonMass     = 0.493677;
  protonMass   = 0.93827;
  kstMass      = 0.896;
  B0Mass       = 5.27958;
  JPsiMass     = 3.096916;
  PsiPMass     = 3.686109;
  D0Mass       = 1.86483;

  JPsiBF       =  7.87e-5; // B0 --> J/psi(mu+mu-) K*0          (1.32+/-0.06 * 5.96+/-0.03)
  JPsiKpiBF    =  5.25e-5; // B0 --> J/psi(mu+mu-) K*0(K+pi-)   (1.32+/-0.06 * 5.96+/-0.03 * 2/3)

  KstMuMuBF    =  1.05e-6; // B0 --> K*0 mu+mu-
  KstKpiMuMuBF =  7e-7;    // B0 --> K*0(K+pi-) mu+mu-          (1.05+/-0.1 * 2/3)

  PsiPBF       = 47.72e-7; // B0 --> psi(2S)(mu+mu-) K*0        (6.04+/-0.4 * 7.9+/-0.9)
  PsiPKpiBF    = 31.81e-7; // B0 --> psi(2S)(mu+mu-) K*0(K+pi-) (6.04+/-0.4 * 7.9+/-0.9 * 2/3)

  muonMassErr  = 3.5e-9;
  pionMassErr  = 3.5e-7;
  kaonMassErr  = 1.6e-5;
  B0MassErr    = 1.7e-4;
  kstSigma     = 0.05;

  PI = 3.141592653589793;

//   KstMassShape = new TF1("KstMassShape",
// 			 "2*sqrt(2)*[0]*[1]* sqrt([0]*[0] * ([0]*[0] + [1]*[1])) / (TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0] * ([0]*[0] + [1]*[1])))) / ((x*x - [0]*[0]) * (x*x - [0]*[0]) + [0]*[0]*[1]*[1])",
// 			 0.0,kstMass*2.);
//   // Breit-Wigner distribution:
//   // [0]: mass of the resonance
//   // [1]: width of the resonance
//   KstMassShape->SetParameter(0,kstMass);
//   KstMassShape->SetParameter(1,kstSigma);
//   KstMassShape->SetParNames("Mean","Width");

  // Define whether to compute the efficiency with good-tagged or mis-tagged events
  RIGHTflavorTAG = rightFlavorTag;


  // ###############################################
  // # ===> Define codes to identify MC type <===  #
  // #                                             #
  // # 1 = B0    --> K*0(K+pi-) mu+mu-             #
  // # 2 = B0bar --> K*0bar(K-pi+) mu+mu-          #
  // #                                             #
  // # 3 = B0    --> K*0(K+pi-) J/psi(mu+mu-)      #
  // # 4 = B0bar --> K*0bar(K-pi+) J/psi(mu+mu-)   #
  // #                                             #
  // # 5 = B0    --> K*0(K-pi+) psi(2S)(mu+mu-)    #
  // # 6 = B0bar --> K*0bar(K-pi+) psi(2S)(mu+mu-) #
  // ###############################################
  B0ToKstMuMu  = 1;
  B0ToJPsiKst  = 3;
  B0ToPsi2SKst = 5;

  // ################################
  // # Print out internal variables #
  // ################################
  std::cout << "\n@@@@@@ Utils class settings : public  @@@@@@" << std::endl;
  std::cout << "RIGHTflavorTAG: "            << RIGHTflavorTAG << std::endl;
  std::cout << "B0ToKstMuMu: "               << B0ToKstMuMu << std::endl;
  std::cout << "B0ToJPsiKst: "               << B0ToJPsiKst << std::endl;
  std::cout << "B0ToPsi2SKst: "              << B0ToPsi2SKst << std::endl;
}

double Utils::computeInvMass (double Px1,
			      double Py1,
			      double Pz1,
			      double mass1,
			      double Px2,
			      double Py2,
			      double Pz2,
			      double mass2,
			      double Px3,
			      double Py3,
			      double Pz3,
			      double mass3)
{
  double Energy1 = sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass1*mass1);
  double Energy2 = sqrt(Px2*Px2 + Py2*Py2 + Pz2*Pz2 + mass2*mass2);
  double Energy3 = sqrt(Px3*Px3 + Py3*Py3 + Pz3*Pz3 + mass3*mass3);
  return sqrt((Energy1+Energy2+Energy3) * (Energy1+Energy2+Energy3) - ((Px1+Px2+Px3) * (Px1+Px2+Px3) + (Py1+Py2+Py3) * (Py1+Py2+Py3) + (Pz1+Pz2+Pz3) * (Pz1+Pz2+Pz3)));
}

double Utils::computeEta (double Px,
			  double Py,
			  double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double Utils::computePhi (double Px, double Py)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double Utils::computeEtaPhiDistance (double Px1,
				     double Py1,
				     double Pz1,
				     double Px2,
				     double Py2,
				     double Pz2)
{
  double phi1 = computePhi (Px1,Py1);
  double eta1 = computeEta (Px1,Py1,Pz1);
  double phi2 = computePhi (Px2,Py2);
  double eta2 = computeEta (Px2,Py2,Pz2);
  return deltaR(eta1,eta2,phi1,phi2);
//   return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}

void Utils::computeLS (double Vx,
		       double Vy,
		       double Vz,
		       double Wx,
		       double Wy,
		       double Wz,
		       double VxErr2,
		       double VyErr2,
		       double VzErr2,
		       double VxyCov,
		       double VxzCov,
		       double VyzCov,
		       double WxErr2,
		       double WyErr2,
		       double WzErr2,
		       double WxyCov,
		       double WxzCov,
		       double WyzCov,
		       double* deltaD,
		       double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}

void Utils::computeCosAlpha (double Vx,
			     double Vy,
			     double Vz,
			     double Wx,
			     double Wy,
			     double Wz,
			     double VxErr2,
			     double VyErr2,
			     double VzErr2,
			     double VxyCov,
			     double VxzCov,
			     double VyzCov,
			     double WxErr2,
			     double WyErr2,
			     double WzErr2,
			     double WxyCov,
			     double WxzCov,
			     double WyzCov,
			     double* cosAlpha,
			     double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  
  if ((Vnorm > 0.) && (Wnorm > 0.))
    {
      *cosAlpha = VdotW / (Vnorm * Wnorm);
      *cosAlphaErr = sqrt( ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			   
			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
    }
  else
    {
      *cosAlpha = 0.;
      *cosAlphaErr = 0.;
    }
}

