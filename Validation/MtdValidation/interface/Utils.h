#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#include "B0KstMuMuTreeContent.h"

class Utils {
public:
  Utils(bool rightFlavorTag = true);
  ~Utils(){};

  double computeInvMass(double Px1,
                        double Py1,
                        double Pz1,
                        double mass1,
                        double Px2,
                        double Py2,
                        double Pz2,
                        double mass2,
                        double Px3 = 0,
                        double Py3 = 0,
                        double Pz3 = 0,
                        double mass3 = 0);

  double computeEta(double Px, double Py, double Pz);

  double computePhi(double Px, double Py);

  double computeEtaPhiDistance(double Px1, double Py1, double Pz1, double Px2, double Py2, double Pz2);

  void computeLS(double Vx,
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
                 double* deltaDErr);

  void computeCosAlpha(double Vx,
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
                       double* cosAlphaErr);

  void Readq2Bins(std::string fileName, std::vector<double>* q2Bins);

  double muonMass;
  double pionMass;
  double kaonMass;
  double protonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;
  double D0Mass;

  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  bool RIGHTflavorTAG;

  int B0ToKstMuMu;
  int B0ToJPsiKst;
  int B0ToPsi2SKst;

  unsigned int nConfigParam;

private:
  std::vector<std::string> HLTpath;
  std::vector<std::string> TrigTable;
};

#endif
