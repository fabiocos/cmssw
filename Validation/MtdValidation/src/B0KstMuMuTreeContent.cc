#include "Validation/MtdValidation/interface/B0KstMuMuTreeContent.h"
#include <iostream>

B0KstMuMuTreeContent::B0KstMuMuTreeContent() {
  ClearScalars();

  // ### B0 Mass ###
  bMass = nullptr;
  bMassE = nullptr;
  bBarMass = nullptr;
  bBarMassE = nullptr;
  bPx = nullptr;
  bPy = nullptr;
  bPz = nullptr;

  // ### Pileup information in MC ###
  bunchXingMC = nullptr;
  numInteractionsMC = nullptr;
  trueNumInteractionsMC = nullptr;

  // ### B0 Vtx ###
  bVtxCL = nullptr;
  bVtxX = nullptr;
  bVtxY = nullptr;
  bVtxZ = nullptr;
  bCosAlphaBS = nullptr;
  bCosAlphaBSE = nullptr;
  bLBS = nullptr;
  bLBSE = nullptr;
  bDCABS = nullptr;
  bDCABSE = nullptr;

  // ### K*0 Mass ###
  kstMass = nullptr;
  kstMassE = nullptr;
  kstBarMass = nullptr;
  kstBarMassE = nullptr;
  kstPx = nullptr;
  kstPy = nullptr;
  kstPz = nullptr;
  kstPxxE = nullptr;
  kstPyyE = nullptr;
  kstPzzE = nullptr;
  kstPxyE = nullptr;
  kstPxzE = nullptr;
  kstPyzE = nullptr;

  // ### K*0 Vtx ###
  kstVtxCL = nullptr;
  kstVtxX = nullptr;
  kstVtxY = nullptr;
  kstVtxZ = nullptr;

  // ### mu+ mu- Mass ###
  mumuMass = nullptr;
  mumuMassE = nullptr;
  mumuPx = nullptr;
  mumuPy = nullptr;
  mumuPz = nullptr;

  // ### mu+ mu- Vtx ###
  mumuVtxCL = nullptr;
  mumuVtxX = nullptr;
  mumuVtxY = nullptr;
  mumuVtxZ = nullptr;
  mumuCosAlphaBS = nullptr;
  mumuCosAlphaBSE = nullptr;
  mumuLBS = nullptr;
  mumuLBSE = nullptr;
  mumuDCA = nullptr;

  // ### mu- ###
  mumHighPurity = nullptr;
  mumCL = nullptr;
  mumNormChi2 = nullptr;
  mumPx = nullptr;
  mumPy = nullptr;
  mumPz = nullptr;
  mumDCAVtx = nullptr;
  mumDCAVtxE = nullptr;
  mumDCABS = nullptr;
  mumDCABSE = nullptr;
  mumKinkChi2 = nullptr;
  mumFracHits = nullptr;
  mumdxyBS = nullptr;
  mumdzBS = nullptr;
  mumMinIP2D = nullptr;
  mumMinIP2DE = nullptr;
  mumMinIP = nullptr;
  mumMinIPS = nullptr;
  mumDeltaRwithMC = nullptr;
  mumCat = nullptr;
  mumNPixHits = nullptr;
  mumNPixLayers = nullptr;
  mumNTrkHits = nullptr;
  mumNTrkLayers = nullptr;
  mumNMuonHits = nullptr;
  mumNMatchStation = nullptr;
  mumIso = nullptr;
  mumIsoPt = nullptr;
  mumIsodR = nullptr;

  // ### mu+ ###
  mupHighPurity = nullptr;
  mupCL = nullptr;
  mupNormChi2 = nullptr;
  mupPx = nullptr;
  mupPy = nullptr;
  mupPz = nullptr;
  mupDCAVtx = nullptr;
  mupDCAVtxE = nullptr;
  mupDCABS = nullptr;
  mupDCABSE = nullptr;
  mupKinkChi2 = nullptr;
  mupFracHits = nullptr;
  mupdxyBS = nullptr;
  mupdzBS = nullptr;
  mupMinIP2D = nullptr;
  mupMinIP2DE = nullptr;
  mupMinIP = nullptr;
  mupMinIPS = nullptr;
  mupDeltaRwithMC = nullptr;
  mupCat = nullptr;
  mupNPixHits = nullptr;
  mupNPixLayers = nullptr;
  mupNTrkHits = nullptr;
  mupNTrkLayers = nullptr;
  mupNMuonHits = nullptr;
  mupNMatchStation = nullptr;
  mupIso = nullptr;
  mupIsoPt = nullptr;
  mupIsodR = nullptr;

  // ### K*0 track- ###
  kstTrkmHighPurity = nullptr;
  kstTrkmCL = nullptr;
  kstTrkmNormChi2 = nullptr;
  kstTrkmPx = nullptr;
  kstTrkmPy = nullptr;
  kstTrkmPz = nullptr;
  kstTrkmPxxE = nullptr;
  kstTrkmPyyE = nullptr;
  kstTrkmPzzE = nullptr;
  kstTrkmPxyE = nullptr;
  kstTrkmPxzE = nullptr;
  kstTrkmPyzE = nullptr;
  kstTrkmDCAVtx = nullptr;
  kstTrkmDCAVtxE = nullptr;
  kstTrkmDCABS = nullptr;
  kstTrkmDCABSE = nullptr;
  kstTrkmFracHits = nullptr;
  kstTrkmdxyVtx = nullptr;
  kstTrkmdzVtx = nullptr;
  kstTrkmMinIP2D = nullptr;
  kstTrkmMinIP2DE = nullptr;
  kstTrkmMinIP = nullptr;
  kstTrkmMinIPS = nullptr;
  kstTrkmDeltaRwithMC = nullptr;
  kstTrkmNPixHits = nullptr;
  kstTrkmNPixLayers = nullptr;
  kstTrkmNTrkHits = nullptr;
  kstTrkmNTrkLayers = nullptr;
  kstTrkmMuMatch = nullptr;
  kstTrkmIso = nullptr;
  kstTrkmIsoPt = nullptr;
  kstTrkmIsodR = nullptr;

  // ### K*0 track+ ###
  kstTrkpHighPurity = nullptr;
  kstTrkpCL = nullptr;
  kstTrkpNormChi2 = nullptr;
  kstTrkpPx = nullptr;
  kstTrkpPy = nullptr;
  kstTrkpPz = nullptr;
  kstTrkpPxxE = nullptr;
  kstTrkpPyyE = nullptr;
  kstTrkpPzzE = nullptr;
  kstTrkpPxyE = nullptr;
  kstTrkpPxzE = nullptr;
  kstTrkpPyzE = nullptr;
  kstTrkpDCAVtx = nullptr;
  kstTrkpDCAVtxE = nullptr;
  kstTrkpDCABS = nullptr;
  kstTrkpDCABSE = nullptr;
  kstTrkpFracHits = nullptr;
  kstTrkpdxyVtx = nullptr;
  kstTrkpdzVtx = nullptr;
  kstTrkpMinIP2D = nullptr;
  kstTrkpMinIP2DE = nullptr;
  kstTrkpMinIP = nullptr;
  kstTrkpMinIPS = nullptr;
  kstTrkpDeltaRwithMC = nullptr;
  kstTrkpNPixHits = nullptr;
  kstTrkpNPixLayers = nullptr;
  kstTrkpNTrkHits = nullptr;
  kstTrkpNTrkLayers = nullptr;
  kstTrkpMuMatch = nullptr;
  kstTrkpIso = nullptr;
  kstTrkpIsoPt = nullptr;
  kstTrkpIsodR = nullptr;

  rawmumPt = nullptr;
  rawmumPhi = nullptr;
  rawmumEta = nullptr;
  rawmupPt = nullptr;
  rawmupPhi = nullptr;
  rawmupEta = nullptr;
  rawkstTrkmPt = nullptr;
  rawkstTrkmPhi = nullptr;
  rawkstTrkmEta = nullptr;
  rawkstTrkpPt = nullptr;
  rawkstTrkpPhi = nullptr;
  rawkstTrkpEta = nullptr;

  bMinusVtxCL = nullptr;
  bMinusCosAlphaBS = nullptr;
  bPlusVtxCL = nullptr;
  bPlusCosAlphaBS = nullptr;

  // ### Matching Between Reconstructed and Generated ###
  truthMatchSignal = nullptr;
  truthMatchMum = nullptr;
  truthMatchMup = nullptr;
  truthMatchTrkm = nullptr;
  truthMatchTrkp = nullptr;
}

void B0KstMuMuTreeContent::Init() {
  // ### B0 Mass ###
  bMass = new std::vector<double>;
  bMassE = new std::vector<double>;
  bBarMass = new std::vector<double>;
  bBarMassE = new std::vector<double>;
  bPx = new std::vector<double>;
  bPy = new std::vector<double>;
  bPz = new std::vector<double>;

  // ### Pileup information in MC ###
  bunchXingMC = new std::vector<double>;
  numInteractionsMC = new std::vector<double>;
  trueNumInteractionsMC = new std::vector<double>;

  // ### B0 Vtx ###
  bVtxCL = new std::vector<double>;
  bVtxX = new std::vector<double>;
  bVtxY = new std::vector<double>;
  bVtxZ = new std::vector<double>;
  bCosAlphaBS = new std::vector<double>;
  bCosAlphaBSE = new std::vector<double>;
  bLBS = new std::vector<double>;
  bLBSE = new std::vector<double>;
  bDCABS = new std::vector<double>;
  bDCABSE = new std::vector<double>;

  // ### K*0 Mass ###
  kstMass = new std::vector<double>;
  kstMassE = new std::vector<double>;
  kstBarMass = new std::vector<double>;
  kstBarMassE = new std::vector<double>;
  kstPx = new std::vector<double>;
  kstPy = new std::vector<double>;
  kstPz = new std::vector<double>;
  kstPxxE = new std::vector<double>;
  kstPyyE = new std::vector<double>;
  kstPzzE = new std::vector<double>;
  kstPxyE = new std::vector<double>;
  kstPxzE = new std::vector<double>;
  kstPyzE = new std::vector<double>;

  // ### K*0 Vtx ###
  kstVtxCL = new std::vector<double>;
  kstVtxX = new std::vector<double>;
  kstVtxY = new std::vector<double>;
  kstVtxZ = new std::vector<double>;

  // ### mu+ mu- Mass ###
  mumuMass = new std::vector<double>;
  mumuMassE = new std::vector<double>;
  mumuPx = new std::vector<double>;
  mumuPy = new std::vector<double>;
  mumuPz = new std::vector<double>;

  // ### mu+ mu- Vtx ###
  mumuVtxCL = new std::vector<double>;
  mumuVtxX = new std::vector<double>;
  mumuVtxY = new std::vector<double>;
  mumuVtxZ = new std::vector<double>;
  mumuCosAlphaBS = new std::vector<double>;
  mumuCosAlphaBSE = new std::vector<double>;
  mumuLBS = new std::vector<double>;
  mumuLBSE = new std::vector<double>;
  mumuDCA = new std::vector<double>;

  // ### mu- ###
  //   mumHighPurity    = new std::vector<int>;
  mumHighPurity = new std::vector<bool>;
  mumCL = new std::vector<double>;
  mumNormChi2 = new std::vector<double>;
  mumPx = new std::vector<double>;
  mumPy = new std::vector<double>;
  mumPz = new std::vector<double>;
  mumDCAVtx = new std::vector<double>;
  mumDCAVtxE = new std::vector<double>;
  mumDCABS = new std::vector<double>;
  mumDCABSE = new std::vector<double>;
  mumKinkChi2 = new std::vector<double>;
  mumFracHits = new std::vector<double>;
  mumdxyBS = new std::vector<double>;
  mumdzBS = new std::vector<double>;
  mumMinIP2D = new std::vector<double>;
  mumMinIP2DE = new std::vector<double>;
  mumMinIP = new std::vector<double>;
  mumMinIPS = new std::vector<double>;
  mumDeltaRwithMC = new std::vector<double>;
  mumCat = new std::vector<std::string>;
  mumNPixHits = new std::vector<int>;
  mumNPixLayers = new std::vector<int>;
  mumNTrkHits = new std::vector<int>;
  mumNTrkLayers = new std::vector<int>;
  mumNMuonHits = new std::vector<int>;
  mumNMatchStation = new std::vector<int>;
  mumIso = new std::vector<std::vector<float> >;
  mumIsoPt = new std::vector<std::vector<float> >;
  mumIsodR = new std::vector<std::vector<float> >;
  matchMum = new std::vector<bool>;
  trkrefMum = new std::vector<reco::TrackRef>;
  mumPathToVtx = new std::vector<double>;

  // ### mu+ ###
  //   mupHighPurity    = new std::vector<int>;
  mupHighPurity = new std::vector<bool>;
  mupCL = new std::vector<double>;
  mupNormChi2 = new std::vector<double>;
  mupPx = new std::vector<double>;
  mupPy = new std::vector<double>;
  mupPz = new std::vector<double>;
  mupDCAVtx = new std::vector<double>;
  mupDCAVtxE = new std::vector<double>;
  mupDCABS = new std::vector<double>;
  mupDCABSE = new std::vector<double>;
  mupKinkChi2 = new std::vector<double>;
  mupFracHits = new std::vector<double>;
  mupdxyBS = new std::vector<double>;
  mupdzBS = new std::vector<double>;
  mupMinIP2D = new std::vector<double>;
  mupMinIP2DE = new std::vector<double>;
  mupMinIP = new std::vector<double>;
  mupMinIPS = new std::vector<double>;
  mupDeltaRwithMC = new std::vector<double>;
  mupCat = new std::vector<std::string>;
  mupNPixHits = new std::vector<int>;
  mupNPixLayers = new std::vector<int>;
  mupNTrkHits = new std::vector<int>;
  mupNTrkLayers = new std::vector<int>;
  mupNMuonHits = new std::vector<int>;
  mupNMatchStation = new std::vector<int>;
  mupIso = new std::vector<std::vector<float> >;
  mupIsoPt = new std::vector<std::vector<float> >;
  mupIsodR = new std::vector<std::vector<float> >;
  matchMup = new std::vector<bool>;
  trkrefMup = new std::vector<reco::TrackRef>;
  mupPathToVtx = new std::vector<double>;

  // ### K*0 track- ###
  //   kstTrkmHighPurity   = new std::vector<int>;
  kstTrkmHighPurity = new std::vector<bool>;
  kstTrkmCL = new std::vector<double>;
  kstTrkmNormChi2 = new std::vector<double>;
  kstTrkmPx = new std::vector<double>;
  kstTrkmPy = new std::vector<double>;
  kstTrkmPz = new std::vector<double>;
  kstTrkmPxxE = new std::vector<double>;
  kstTrkmPyyE = new std::vector<double>;
  kstTrkmPzzE = new std::vector<double>;
  kstTrkmPxyE = new std::vector<double>;
  kstTrkmPxzE = new std::vector<double>;
  kstTrkmPyzE = new std::vector<double>;
  kstTrkmDCAVtx = new std::vector<double>;
  kstTrkmDCAVtxE = new std::vector<double>;
  kstTrkmDCABS = new std::vector<double>;
  kstTrkmDCABSE = new std::vector<double>;
  kstTrkmFracHits = new std::vector<double>;
  kstTrkmdxyVtx = new std::vector<double>;
  kstTrkmdzVtx = new std::vector<double>;
  kstTrkmMinIP2D = new std::vector<double>;
  kstTrkmMinIP2DE = new std::vector<double>;
  kstTrkmMinIP = new std::vector<double>;
  kstTrkmMinIPS = new std::vector<double>;
  kstTrkmDeltaRwithMC = new std::vector<double>;
  kstTrkmNPixHits = new std::vector<int>;
  kstTrkmNPixLayers = new std::vector<int>;
  kstTrkmNTrkHits = new std::vector<int>;
  kstTrkmNTrkLayers = new std::vector<int>;
  kstTrkmMuMatch = new std::vector<std::string>;
  kstTrkmIso = new std::vector<std::vector<float> >;
  kstTrkmIsoPt = new std::vector<std::vector<float> >;
  kstTrkmIsodR = new std::vector<std::vector<float> >;
  matchTkm = new std::vector<bool>;
  trkrefTkm = new std::vector<reco::TrackRef>;
  trkmPathToVtx = new std::vector<double>;

  // ### K*0 track+ ###
  //   kstTrkpHighPurity   = new std::vector<int>;
  kstTrkpHighPurity = new std::vector<bool>;
  kstTrkpCL = new std::vector<double>;
  kstTrkpNormChi2 = new std::vector<double>;
  kstTrkpPx = new std::vector<double>;
  kstTrkpPy = new std::vector<double>;
  kstTrkpPz = new std::vector<double>;
  kstTrkpPxxE = new std::vector<double>;
  kstTrkpPyyE = new std::vector<double>;
  kstTrkpPzzE = new std::vector<double>;
  kstTrkpPxyE = new std::vector<double>;
  kstTrkpPxzE = new std::vector<double>;
  kstTrkpPyzE = new std::vector<double>;
  kstTrkpDCAVtx = new std::vector<double>;
  kstTrkpDCAVtxE = new std::vector<double>;
  kstTrkpDCABS = new std::vector<double>;
  kstTrkpDCABSE = new std::vector<double>;
  kstTrkpFracHits = new std::vector<double>;
  kstTrkpdxyVtx = new std::vector<double>;
  kstTrkpdzVtx = new std::vector<double>;
  kstTrkpMinIP2D = new std::vector<double>;
  kstTrkpMinIP2DE = new std::vector<double>;
  kstTrkpMinIP = new std::vector<double>;
  kstTrkpMinIPS = new std::vector<double>;
  kstTrkpDeltaRwithMC = new std::vector<double>;
  kstTrkpNPixHits = new std::vector<int>;
  kstTrkpNPixLayers = new std::vector<int>;
  kstTrkpNTrkHits = new std::vector<int>;
  kstTrkpNTrkLayers = new std::vector<int>;
  kstTrkpMuMatch = new std::vector<std::string>;
  kstTrkpIso = new std::vector<std::vector<float> >;
  kstTrkpIsoPt = new std::vector<std::vector<float> >;
  kstTrkpIsodR = new std::vector<std::vector<float> >;
  matchTkp = new std::vector<bool>;
  trkrefTkp = new std::vector<reco::TrackRef>;
  trkpPathToVtx = new std::vector<double>;

  rawmumPt = new std::vector<double>;
  rawmumPhi = new std::vector<double>;
  rawmumEta = new std::vector<double>;
  rawmupPt = new std::vector<double>;
  rawmupPhi = new std::vector<double>;
  rawmupEta = new std::vector<double>;
  rawkstTrkmPt = new std::vector<double>;
  rawkstTrkmPhi = new std::vector<double>;
  rawkstTrkmEta = new std::vector<double>;
  rawkstTrkpPt = new std::vector<double>;
  rawkstTrkpPhi = new std::vector<double>;
  rawkstTrkpEta = new std::vector<double>;

  bMinusVtxCL = new std::vector<double>;
  bMinusCosAlphaBS = new std::vector<double>;
  bPlusVtxCL = new std::vector<double>;
  bPlusCosAlphaBS = new std::vector<double>;

  // ### Matching Between Reconstructed and Generated ###
  //   truthMatchSignal = new std::vector<int>;
  //   truthMatchMum    = new std::vector<int>;
  //   truthMatchMup    = new std::vector<int>;
  //   truthMatchTrkm   = new std::vector<int>;
  //   truthMatchTrkp   = new std::vector<int>;
  truthMatchSignal = new std::vector<bool>;
  truthMatchMum = new std::vector<bool>;
  truthMatchMup = new std::vector<bool>;
  truthMatchTrkm = new std::vector<bool>;
  truthMatchTrkp = new std::vector<bool>;
}

B0KstMuMuTreeContent::~B0KstMuMuTreeContent() {
  // ### B0 Mass ###
  delete bMass;
  delete bMassE;
  delete bBarMass;
  delete bBarMassE;
  delete bPx;
  delete bPy;
  delete bPz;

  // ### Pileup information in MC ###
  delete bunchXingMC;
  delete numInteractionsMC;
  delete trueNumInteractionsMC;

  // ### B0 Vtx ###
  delete bVtxCL;
  delete bVtxX;
  delete bVtxY;
  delete bVtxZ;
  delete bCosAlphaBS;
  delete bCosAlphaBSE;
  delete bLBS;
  delete bLBSE;
  delete bDCABS;
  delete bDCABSE;

  // ### K*0 Mass ###
  delete kstMass;
  delete kstMassE;
  delete kstBarMass;
  delete kstBarMassE;
  delete kstPx;
  delete kstPy;
  delete kstPz;
  delete kstPxxE;
  delete kstPyyE;
  delete kstPzzE;
  delete kstPxyE;
  delete kstPxzE;
  delete kstPyzE;

  // ### K*0 Vtx ###
  delete kstVtxCL;
  delete kstVtxX;
  delete kstVtxY;
  delete kstVtxZ;

  // ### mu+ mu- Mass ###
  delete mumuMass;
  delete mumuMassE;
  delete mumuPx;
  delete mumuPy;
  delete mumuPz;

  // ### mu+ mu- Vtx ###
  delete mumuVtxCL;
  delete mumuVtxX;
  delete mumuVtxY;
  delete mumuVtxZ;
  delete mumuCosAlphaBS;
  delete mumuCosAlphaBSE;
  delete mumuLBS;
  delete mumuLBSE;
  delete mumuDCA;

  // ### mu- ###
  delete mumHighPurity;
  delete mumCL;
  delete mumNormChi2;
  delete mumPx;
  delete mumPy;
  delete mumPz;
  delete mumDCAVtx;
  delete mumDCAVtxE;
  delete mumDCABS;
  delete mumDCABSE;
  delete mumKinkChi2;
  delete mumFracHits;
  delete mumdxyBS;
  delete mumdzBS;
  delete mumMinIP2D;
  delete mumMinIP2DE;
  delete mumMinIP;
  delete mumMinIPS;
  delete mumDeltaRwithMC;
  delete mumCat;
  delete mumNPixHits;
  delete mumNPixLayers;
  delete mumNTrkHits;
  delete mumNTrkLayers;
  delete mumNMuonHits;
  delete mumNMatchStation;
  delete mumIso;
  delete mumIsoPt;
  delete mumIsodR;
  delete matchMum;
  delete trkrefMum;
  delete mumPathToVtx;

  // ### mu+ ###
  delete mupHighPurity;
  delete mupCL;
  delete mupNormChi2;
  delete mupPx;
  delete mupPy;
  delete mupPz;
  delete mupDCAVtx;
  delete mupDCAVtxE;
  delete mupDCABS;
  delete mupDCABSE;
  delete mupKinkChi2;
  delete mupFracHits;
  delete mupdxyBS;
  delete mupdzBS;
  delete mupMinIP2D;
  delete mupMinIP2DE;
  delete mupMinIP;
  delete mupMinIPS;
  delete mupDeltaRwithMC;
  delete mupCat;
  delete mupNPixHits;
  delete mupNPixLayers;
  delete mupNTrkHits;
  delete mupNTrkLayers;
  delete mupNMuonHits;
  delete mupNMatchStation;
  delete mupIso;
  delete mupIsoPt;
  delete mupIsodR;
  delete matchMup;
  delete trkrefMup;
  delete mupPathToVtx;

  // ### K*0 track- ###
  delete kstTrkmHighPurity;
  delete kstTrkmCL;
  delete kstTrkmNormChi2;
  delete kstTrkmPx;
  delete kstTrkmPy;
  delete kstTrkmPz;
  delete kstTrkmPxxE;
  delete kstTrkmPyyE;
  delete kstTrkmPzzE;
  delete kstTrkmPxyE;
  delete kstTrkmPxzE;
  delete kstTrkmPyzE;
  delete kstTrkmDCAVtx;
  delete kstTrkmDCAVtxE;
  delete kstTrkmDCABS;
  delete kstTrkmDCABSE;
  delete kstTrkmFracHits;
  delete kstTrkmdxyVtx;
  delete kstTrkmdzVtx;
  delete kstTrkmMinIP2D;
  delete kstTrkmMinIP2DE;
  delete kstTrkmMinIP;
  delete kstTrkmMinIPS;
  delete kstTrkmDeltaRwithMC;
  delete kstTrkmNPixHits;
  delete kstTrkmNPixLayers;
  delete kstTrkmNTrkHits;
  delete kstTrkmNTrkLayers;
  delete kstTrkmMuMatch;
  delete kstTrkmIso;
  delete kstTrkmIsoPt;
  delete kstTrkmIsodR;
  delete matchTkm;
  delete trkrefTkm;
  delete trkmPathToVtx;

  // ### K*0 track+ ###
  delete kstTrkpHighPurity;
  delete kstTrkpCL;
  delete kstTrkpNormChi2;
  delete kstTrkpPx;
  delete kstTrkpPy;
  delete kstTrkpPz;
  delete kstTrkpPxxE;
  delete kstTrkpPyyE;
  delete kstTrkpPzzE;
  delete kstTrkpPxyE;
  delete kstTrkpPxzE;
  delete kstTrkpPyzE;
  delete kstTrkpDCAVtx;
  delete kstTrkpDCAVtxE;
  delete kstTrkpDCABS;
  delete kstTrkpDCABSE;
  delete kstTrkpFracHits;
  delete kstTrkpdxyVtx;
  delete kstTrkpdzVtx;
  delete kstTrkpMinIP2D;
  delete kstTrkpMinIP2DE;
  delete kstTrkpMinIP;
  delete kstTrkpMinIPS;
  delete kstTrkpDeltaRwithMC;
  delete kstTrkpNPixHits;
  delete kstTrkpNPixLayers;
  delete kstTrkpNTrkHits;
  delete kstTrkpNTrkLayers;
  delete kstTrkpMuMatch;
  delete kstTrkpIso;
  delete kstTrkpIsoPt;
  delete kstTrkpIsodR;
  delete matchTkp;
  delete trkrefTkp;
  delete trkpPathToVtx;

  delete rawmumPt;
  delete rawmumPhi;
  delete rawmumEta;
  delete rawmupPt;
  delete rawmupPhi;
  delete rawmupEta;
  delete rawkstTrkmPt;
  delete rawkstTrkmPhi;
  delete rawkstTrkmEta;
  delete rawkstTrkpPt;
  delete rawkstTrkpPhi;
  delete rawkstTrkpEta;

  delete bMinusVtxCL;
  delete bMinusCosAlphaBS;
  delete bPlusVtxCL;
  delete bPlusCosAlphaBS;

  // ### Matching Between Reconstructed and Generated ###
  delete truthMatchSignal;
  delete truthMatchMum;
  delete truthMatchMup;
  delete truthMatchTrkm;
  delete truthMatchTrkp;
}

void B0KstMuMuTreeContent::ClearScalars() {
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  runN = 0;
  eventN = 0;
  recoVtxN = 0;
  evWeight = 1;
  evWeightE2 = 0;
  numEventsTried = 0;
  numEventsPassed = 0;

  nB = 0;

  // ### Primary Vertex and Beam Spot ###
  bsX = 0;
  bsY = 0;

  ClearScalarsMonteCarlo();
}

void B0KstMuMuTreeContent::ClearScalarsMonteCarlo() {
  // ### Generated Observables ###
  genSignal = 0;
  genMuMuBG = 0;
  genMuMuBGnTrksm = 0;
  genMuMuBGnTrksp = 0;
  genPsiPrompt = false;
  genSignHasFSR = false;
  genSignKstHasFSR = false;
  genSignPsiHasFSR = false;

  // # Generated Primary Vertex #
  genPriVtxX = 0;
  genPriVtxY = 0;
  genPriVtxZ = 0;

  // ### Generated B0 Mass ###
  genB0Mass = 0;
  genB0Px = 0;
  genB0Py = 0;
  genB0Pz = 0;

  // ### Generated B0 Vtx ###
  genB0VtxX = 0;
  genB0VtxY = 0;
  genB0VtxZ = 0;

  // ### Generated K*0 Mass ###
  genKstMass = 0;
  genKstPx = 0;
  genKstPy = 0;
  genKstPz = 0;

  // ### Generated K*0 Vtx ###
  genKstVtxX = 0;
  genKstVtxY = 0;
  genKstVtxZ = 0;

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  genPsiMass = 0;
  genPsiVtxX = 0;
  genPsiVtxY = 0;
  genPsiVtxZ = 0;

  // ### Generated mu- ###
  genMumMother = 0;
  genMumPx = 0;
  genMumPy = 0;
  genMumPz = 0;

  // ### Generated mu+ ###
  genMupMother = 0;
  genMupPx = 0;
  genMupPy = 0;
  genMupPz = 0;

  // ### Generated K*0 track- ###
  genKstTrkmMother = 0;
  genKstTrkmID = 0;
  genKstTrkmPx = 0;
  genKstTrkmPy = 0;
  genKstTrkmPz = 0;

  // ### Generated K*0 track+ ###
  genKstTrkpMother = 0;
  genKstTrkpID = 0;
  genKstTrkpPx = 0;
  genKstTrkpPy = 0;
  genKstTrkpPz = 0;
}

void B0KstMuMuTreeContent::ClearVectors() {
  // ### B0 Mass ###
  bMass->clear();
  bMassE->clear();
  bBarMass->clear();
  bBarMassE->clear();
  bPx->clear();
  bPy->clear();
  bPz->clear();

  // ### Pileup information in MC ###
  bunchXingMC->clear();
  numInteractionsMC->clear();
  trueNumInteractionsMC->clear();

  // ### B0 Vtx ###
  bVtxCL->clear();
  bVtxX->clear();
  bVtxY->clear();
  bVtxZ->clear();
  bCosAlphaBS->clear();
  bCosAlphaBSE->clear();
  bLBS->clear();
  bLBSE->clear();
  bDCABS->clear();
  bDCABSE->clear();

  // ### K*0 Mass ###
  kstMass->clear();
  kstMassE->clear();
  kstBarMass->clear();
  kstBarMassE->clear();
  kstPx->clear();
  kstPy->clear();
  kstPz->clear();
  kstPxxE->clear();
  kstPyyE->clear();
  kstPzzE->clear();
  kstPxyE->clear();
  kstPxzE->clear();
  kstPyzE->clear();

  // ### K*0 Vtx ###
  kstVtxCL->clear();
  kstVtxX->clear();
  kstVtxY->clear();
  kstVtxZ->clear();

  // ### mu+ mu- Mass ###
  mumuMass->clear();
  mumuMassE->clear();
  mumuPx->clear();
  mumuPy->clear();
  mumuPz->clear();

  // ### mu+ mu- Vtx ###
  mumuVtxCL->clear();
  mumuVtxX->clear();
  mumuVtxY->clear();
  mumuVtxZ->clear();
  mumuCosAlphaBS->clear();
  mumuCosAlphaBSE->clear();
  mumuLBS->clear();
  mumuLBSE->clear();
  mumuDCA->clear();

  // ### mu- ###
  mumHighPurity->clear();
  mumCL->clear();
  mumNormChi2->clear();
  mumPx->clear();
  mumPy->clear();
  mumPz->clear();
  mumDCAVtx->clear();
  mumDCAVtxE->clear();
  mumDCABS->clear();
  mumDCABSE->clear();
  mumKinkChi2->clear();
  mumFracHits->clear();
  mumdxyBS->clear();
  mumdzBS->clear();
  mumMinIP2D->clear();
  mumMinIP2DE->clear();
  mumMinIP->clear();
  mumMinIPS->clear();
  mumDeltaRwithMC->clear();
  mumCat->clear();
  mumNPixHits->clear();
  mumNPixLayers->clear();
  mumNTrkHits->clear();
  mumNTrkLayers->clear();
  mumNMuonHits->clear();
  mumNMatchStation->clear();
  mumIso->clear();
  mumIsoPt->clear();
  mumIsodR->clear();
  matchMum->clear();
  trkrefMum->clear();
  mumPathToVtx->clear();

  // ### mu+ ###
  mupHighPurity->clear();
  mupCL->clear();
  mupNormChi2->clear();
  mupPx->clear();
  mupPy->clear();
  mupPz->clear();
  mupDCAVtx->clear();
  mupDCAVtxE->clear();
  mupDCABS->clear();
  mupDCABSE->clear();
  mupKinkChi2->clear();
  mupFracHits->clear();
  mupdxyBS->clear();
  mupdzBS->clear();
  mupMinIP2D->clear();
  mupMinIP2DE->clear();
  mupMinIP->clear();
  mupMinIPS->clear();
  mupDeltaRwithMC->clear();
  mupCat->clear();
  mupNPixHits->clear();
  mupNPixLayers->clear();
  mupNTrkHits->clear();
  mupNTrkLayers->clear();
  mupNMuonHits->clear();
  mupNMatchStation->clear();
  mupIso->clear();
  mupIsoPt->clear();
  mupIsodR->clear();
  matchMup->clear();
  trkrefMup->clear();
  mupPathToVtx->clear();

  // ### K*0 track- ###
  kstTrkmHighPurity->clear();
  kstTrkmCL->clear();
  kstTrkmNormChi2->clear();
  kstTrkmPx->clear();
  kstTrkmPy->clear();
  kstTrkmPz->clear();
  kstTrkmPxxE->clear();
  kstTrkmPyyE->clear();
  kstTrkmPzzE->clear();
  kstTrkmPxyE->clear();
  kstTrkmPxzE->clear();
  kstTrkmPyzE->clear();
  kstTrkmDCAVtx->clear();
  kstTrkmDCAVtxE->clear();
  kstTrkmDCABS->clear();
  kstTrkmDCABSE->clear();
  kstTrkmFracHits->clear();
  kstTrkmdxyVtx->clear();
  kstTrkmdzVtx->clear();
  kstTrkmMinIP2D->clear();
  kstTrkmMinIP2DE->clear();
  kstTrkmMinIP->clear();
  kstTrkmMinIPS->clear();
  kstTrkmDeltaRwithMC->clear();
  kstTrkmNPixHits->clear();
  kstTrkmNPixLayers->clear();
  kstTrkmNTrkHits->clear();
  kstTrkmNTrkLayers->clear();
  kstTrkmMuMatch->clear();
  kstTrkmIso->clear();
  kstTrkmIsoPt->clear();
  kstTrkmIsodR->clear();
  matchTkm->clear();
  trkrefTkm->clear();

  // ### K*0 track+ ###
  kstTrkpHighPurity->clear();
  kstTrkpCL->clear();
  kstTrkpNormChi2->clear();
  kstTrkpPx->clear();
  kstTrkpPy->clear();
  kstTrkpPz->clear();
  kstTrkpPxxE->clear();
  kstTrkpPyyE->clear();
  kstTrkpPzzE->clear();
  kstTrkpPxyE->clear();
  kstTrkpPxzE->clear();
  kstTrkpPyzE->clear();
  kstTrkpDCAVtx->clear();
  kstTrkpDCAVtxE->clear();
  kstTrkpDCABS->clear();
  kstTrkpDCABSE->clear();
  kstTrkpFracHits->clear();
  kstTrkpdxyVtx->clear();
  kstTrkpdzVtx->clear();
  kstTrkpMinIP2D->clear();
  kstTrkpMinIP2DE->clear();
  kstTrkpMinIP->clear();
  kstTrkpMinIPS->clear();
  kstTrkpDeltaRwithMC->clear();
  kstTrkpNPixHits->clear();
  kstTrkpNPixLayers->clear();
  kstTrkpNTrkHits->clear();
  kstTrkpNTrkLayers->clear();
  kstTrkpMuMatch->clear();
  kstTrkpIso->clear();
  kstTrkpIsoPt->clear();
  kstTrkpIsodR->clear();
  matchTkp->clear();
  trkrefTkp->clear();

  rawmumPt->clear();
  rawmumPhi->clear();
  rawmumEta->clear();
  rawmupPt->clear();
  rawmupPhi->clear();
  rawmupEta->clear();
  rawkstTrkmPt->clear();
  rawkstTrkmPhi->clear();
  rawkstTrkmEta->clear();
  rawkstTrkpPt->clear();
  rawkstTrkpPhi->clear();
  rawkstTrkpEta->clear();

  bMinusVtxCL->clear();
  bMinusCosAlphaBS->clear();
  bPlusVtxCL->clear();
  bPlusCosAlphaBS->clear();

  ClearVectorsMonteCarlo();
}

void B0KstMuMuTreeContent::ClearVectorsMonteCarlo() {
  // ### Matching Between Reconstructed and Generated ###
  truthMatchSignal->clear();
  truthMatchMum->clear();
  truthMatchMup->clear();
  truthMatchTrkm->clear();
  truthMatchTrkp->clear();
}

void B0KstMuMuTreeContent::ClearNTuple() {
  ClearScalars();
  ClearVectors();
}

void B0KstMuMuTreeContent::ClearMonteCarlo() {
  ClearScalarsMonteCarlo();
  ClearVectorsMonteCarlo();
}

void B0KstMuMuTreeContent::MakeTreeBranches(TTree* theTree) {
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  theTree->Branch("runN", &runN, "runN/i");
  theTree->Branch("eventN", &eventN, "eventN/i");
  theTree->Branch("recoVtxN", &recoVtxN, "recoVtxN/i");
  theTree->Branch("evWeight", &evWeight, "evWeight/D");
  theTree->Branch("evWeightE2", &evWeightE2, "evWeightE2/D");
  theTree->Branch("numEventsTried", &numEventsTried, "numEventsTried/i");
  theTree->Branch("numEventsPassed", &numEventsPassed, "numEventsPassed/i");

  theTree->Branch("nB", &nB, "nB/i");

  // ### Primary Vertex and Beam Spot ###
  theTree->Branch("bsX", &bsX, "bsX/D");
  theTree->Branch("bsY", &bsY, "bsY/D");

  // ### B0 Mass ###
  theTree->Branch("bMass", &bMass);
  theTree->Branch("bMassE", &bMassE);
  theTree->Branch("bBarMass", &bBarMass);
  theTree->Branch("bBarMassE", &bBarMassE);
  theTree->Branch("bPx", &bPx);
  theTree->Branch("bPy", &bPy);
  theTree->Branch("bPz", &bPz);

  // ### Pileup information in MC ###
  theTree->Branch("bunchXingMC", &bunchXingMC);
  theTree->Branch("numInteractionsMC", &numInteractionsMC);
  theTree->Branch("trueNumInteractionsMC", &trueNumInteractionsMC);

  // ### B0 Vtx ###
  theTree->Branch("bVtxCL", &bVtxCL);
  theTree->Branch("bVtxX", &bVtxX);
  theTree->Branch("bVtxY", &bVtxY);
  theTree->Branch("bVtxZ", &bVtxZ);
  theTree->Branch("bCosAlphaBS", &bCosAlphaBS);
  theTree->Branch("bCosAlphaBSE", &bCosAlphaBSE);
  theTree->Branch("bLBS", &bLBS);
  theTree->Branch("bLBSE", &bLBSE);
  theTree->Branch("bDCABS", &bDCABS);
  theTree->Branch("bDCABSE", &bDCABSE);

  // ### K*0 Mass ###
  theTree->Branch("kstMass", &kstMass);
  theTree->Branch("kstMassE", &kstMassE);
  theTree->Branch("kstBarMass", &kstBarMass);
  theTree->Branch("kstBarMassE", &kstBarMassE);
  theTree->Branch("kstPx", &kstPx);
  theTree->Branch("kstPy", &kstPy);
  theTree->Branch("kstPz", &kstPz);
  theTree->Branch("kstPxxE", &kstPxxE);
  theTree->Branch("kstPyyE", &kstPyyE);
  theTree->Branch("kstPzzE", &kstPzzE);
  theTree->Branch("kstPxyE", &kstPxyE);
  theTree->Branch("kstPxzE", &kstPxzE);
  theTree->Branch("kstPyzE", &kstPyzE);

  // ### K*0 Vtx ###
  theTree->Branch("kstVtxCL", &kstVtxCL);
  theTree->Branch("kstVtxX", &kstVtxX);
  theTree->Branch("kstVtxY", &kstVtxY);
  theTree->Branch("kstVtxZ", &kstVtxZ);

  // ### mu+ mu- Mass ###
  theTree->Branch("mumuMass", &mumuMass);
  theTree->Branch("mumuMassE", &mumuMassE);
  theTree->Branch("mumuPx", &mumuPx);
  theTree->Branch("mumuPy", &mumuPy);
  theTree->Branch("mumuPz", &mumuPz);

  // ### mu+ mu- Vtx ###
  theTree->Branch("mumuVtxCL", &mumuVtxCL);
  theTree->Branch("mumuVtxX", &mumuVtxX);
  theTree->Branch("mumuVtxY", &mumuVtxY);
  theTree->Branch("mumuVtxZ", &mumuVtxZ);
  theTree->Branch("mumuCosAlphaBS", &mumuCosAlphaBS);
  theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE);
  theTree->Branch("mumuLBS", &mumuLBS);
  theTree->Branch("mumuLBSE", &mumuLBSE);
  theTree->Branch("mumuDCA", &mumuDCA);

  // ### mu- ###
  theTree->Branch("mumHighPurity", &mumHighPurity);
  theTree->Branch("mumCL", &mumCL);
  theTree->Branch("mumNormChi2", &mumNormChi2);
  theTree->Branch("mumPx", &mumPx);
  theTree->Branch("mumPy", &mumPy);
  theTree->Branch("mumPz", &mumPz);
  theTree->Branch("mumDCAVtx", &mumDCAVtx);
  theTree->Branch("mumDCAVtxE", &mumDCAVtxE);
  theTree->Branch("mumDCABS", &mumDCABS);
  theTree->Branch("mumDCABSE", &mumDCABSE);
  theTree->Branch("mumKinkChi2", &mumKinkChi2);
  theTree->Branch("mumFracHits", &mumFracHits);
  theTree->Branch("mumdxyBS", &mumdxyBS);
  theTree->Branch("mumdzBS", &mumdzBS);
  theTree->Branch("mumMinIP2D", &mumMinIP2D);
  theTree->Branch("mumMinIP2DE", &mumMinIP2DE);
  theTree->Branch("mumMinIP", &mumMinIP);
  theTree->Branch("mumMinIPS", &mumMinIPS);
  theTree->Branch("mumDeltaRwithMC", &mumDeltaRwithMC);
  theTree->Branch("mumCat", &mumCat);
  theTree->Branch("mumNPixHits", &mumNPixHits);
  theTree->Branch("mumNPixLayers", &mumNPixLayers);
  theTree->Branch("mumNTrkHits", &mumNTrkHits);
  theTree->Branch("mumNTrkLayers", &mumNTrkLayers);
  theTree->Branch("mumNMuonHits", &mumNMuonHits);
  theTree->Branch("mumNMatchStation", &mumNMatchStation);
  theTree->Branch("mumIso", &mumIso);
  theTree->Branch("mumIsoPt", &mumIsoPt);
  theTree->Branch("mumIsodR", &mumIsodR);

  // ### mu+ ###
  theTree->Branch("mupHighPurity", &mupHighPurity);
  theTree->Branch("mupCL", &mupCL);
  theTree->Branch("mupNormChi2", &mupNormChi2);
  theTree->Branch("mupPx", &mupPx);
  theTree->Branch("mupPy", &mupPy);
  theTree->Branch("mupPz", &mupPz);
  theTree->Branch("mupDCAVtx", &mupDCAVtx);
  theTree->Branch("mupDCAVtxE", &mupDCAVtxE);
  theTree->Branch("mupDCABS", &mupDCABS);
  theTree->Branch("mupDCABSE", &mupDCABSE);
  theTree->Branch("mupKinkChi2", &mupKinkChi2);
  theTree->Branch("mupFracHits", &mupFracHits);
  theTree->Branch("mupdxyBS", &mupdxyBS);
  theTree->Branch("mupdzBS", &mupdzBS);
  theTree->Branch("mupMinIP2D", &mupMinIP2D);
  theTree->Branch("mupMinIP2DE", &mupMinIP2DE);
  theTree->Branch("mupMinIP", &mupMinIP);
  theTree->Branch("mupMinIPS", &mupMinIPS);
  theTree->Branch("mupDeltaRwithMC", &mupDeltaRwithMC);
  theTree->Branch("mupCat", &mupCat);
  theTree->Branch("mupNPixHits", &mupNPixHits);
  theTree->Branch("mupNPixLayers", &mupNPixLayers);
  theTree->Branch("mupNTrkHits", &mupNTrkHits);
  theTree->Branch("mupNTrkLayers", &mupNTrkLayers);
  theTree->Branch("mupNMuonHits", &mupNMuonHits);
  theTree->Branch("mupNMatchStation", &mupNMatchStation);
  theTree->Branch("mupIso", &mupIso);
  theTree->Branch("mupIsoPt", &mupIsoPt);
  theTree->Branch("mupIsodR", &mupIsodR);

  // ### K*0 track- ###
  theTree->Branch("kstTrkmHighPurity", &kstTrkmHighPurity);
  theTree->Branch("kstTrkmCL", &kstTrkmCL);
  theTree->Branch("kstTrkmNormChi2", &kstTrkmNormChi2);
  theTree->Branch("kstTrkmPx", &kstTrkmPx);
  theTree->Branch("kstTrkmPy", &kstTrkmPy);
  theTree->Branch("kstTrkmPz", &kstTrkmPz);
  theTree->Branch("kstTrkmPxxE", &kstTrkmPxxE);
  theTree->Branch("kstTrkmPyyE", &kstTrkmPyyE);
  theTree->Branch("kstTrkmPzzE", &kstTrkmPzzE);
  theTree->Branch("kstTrkmPxyE", &kstTrkmPxyE);
  theTree->Branch("kstTrkmPxzE", &kstTrkmPxzE);
  theTree->Branch("kstTrkmPyzE", &kstTrkmPyzE);
  theTree->Branch("kstTrkmDCAVtx", &kstTrkmDCAVtx);
  theTree->Branch("kstTrkmDCAVtxE", &kstTrkmDCAVtxE);
  theTree->Branch("kstTrkmDCABS", &kstTrkmDCABS);
  theTree->Branch("kstTrkmDCABSE", &kstTrkmDCABSE);
  theTree->Branch("kstTrkmFracHits", &kstTrkmFracHits);
  theTree->Branch("kstTrkmdxyVtx", &kstTrkmdxyVtx);
  theTree->Branch("kstTrkmdzVtx", &kstTrkmdzVtx);
  theTree->Branch("kstTrkmMinIP2D", &kstTrkmMinIP2D);
  theTree->Branch("kstTrkmMinIP2DE", &kstTrkmMinIP2DE);
  theTree->Branch("kstTrkmMinIP", &kstTrkmMinIP);
  theTree->Branch("kstTrkmMinIPS", &kstTrkmMinIPS);
  theTree->Branch("kstTrkmDeltaRwithMC", &kstTrkmDeltaRwithMC);
  theTree->Branch("kstTrkmNPixHits", &kstTrkmNPixHits);
  theTree->Branch("kstTrkmNPixLayers", &kstTrkmNPixLayers);
  theTree->Branch("kstTrkmNTrkHits", &kstTrkmNTrkHits);
  theTree->Branch("kstTrkmNTrkLayers", &kstTrkmNTrkLayers);
  theTree->Branch("kstTrkmMuMatch", &kstTrkmMuMatch);
  theTree->Branch("kstTrkmIso", &kstTrkmIso);
  theTree->Branch("kstTrkmIsoPt", &kstTrkmIsoPt);
  theTree->Branch("kstTrkmIsodR", &kstTrkmIsodR);

  // ### K*0 track+ ###
  theTree->Branch("kstTrkpHighPurity", &kstTrkpHighPurity);
  theTree->Branch("kstTrkpCL", &kstTrkpCL);
  theTree->Branch("kstTrkpNormChi2", &kstTrkpNormChi2);
  theTree->Branch("kstTrkpPx", &kstTrkpPx);
  theTree->Branch("kstTrkpPy", &kstTrkpPy);
  theTree->Branch("kstTrkpPz", &kstTrkpPz);
  theTree->Branch("kstTrkpPxxE", &kstTrkpPxxE);
  theTree->Branch("kstTrkpPyyE", &kstTrkpPyyE);
  theTree->Branch("kstTrkpPzzE", &kstTrkpPzzE);
  theTree->Branch("kstTrkpPxyE", &kstTrkpPxyE);
  theTree->Branch("kstTrkpPxzE", &kstTrkpPxzE);
  theTree->Branch("kstTrkpPyzE", &kstTrkpPyzE);
  theTree->Branch("kstTrkpDCAVtx", &kstTrkpDCAVtx);
  theTree->Branch("kstTrkpDCAVtxE", &kstTrkpDCAVtxE);
  theTree->Branch("kstTrkpDCABS", &kstTrkpDCABS);
  theTree->Branch("kstTrkpDCABSE", &kstTrkpDCABSE);
  theTree->Branch("kstTrkpFracHits", &kstTrkpFracHits);
  theTree->Branch("kstTrkpdxyVtx", &kstTrkpdxyVtx);
  theTree->Branch("kstTrkpdzVtx", &kstTrkpdzVtx);
  theTree->Branch("kstTrkpMinIP2D", &kstTrkpMinIP2D);
  theTree->Branch("kstTrkpMinIP2DE", &kstTrkpMinIP2DE);
  theTree->Branch("kstTrkpMinIP", &kstTrkpMinIP);
  theTree->Branch("kstTrkpMinIPS", &kstTrkpMinIPS);
  theTree->Branch("kstTrkpDeltaRwithMC", &kstTrkpDeltaRwithMC);
  theTree->Branch("kstTrkpNPixHits", &kstTrkpNPixHits);
  theTree->Branch("kstTrkpNPixLayers", &kstTrkpNPixLayers);
  theTree->Branch("kstTrkpNTrkHits", &kstTrkpNTrkHits);
  theTree->Branch("kstTrkpNTrkLayers", &kstTrkpNTrkLayers);
  theTree->Branch("kstTrkpMuMatch", &kstTrkpMuMatch);
  theTree->Branch("kstTrkpIso", &kstTrkpIso);
  theTree->Branch("kstTrkpIsoPt", &kstTrkpIsoPt);
  theTree->Branch("kstTrkpIsodR", &kstTrkpIsodR);

  theTree->Branch("rawmumPt", &rawmumPt);
  theTree->Branch("rawmumPhi", &rawmumPhi);
  theTree->Branch("rawmumEta", &rawmumEta);
  theTree->Branch("rawmupPt", &rawmupPt);
  theTree->Branch("rawmupPhi", &rawmupPhi);
  theTree->Branch("rawmupEta", &rawmupEta);
  theTree->Branch("rawkstTrkmPt", &rawkstTrkmPt);
  theTree->Branch("rawkstTrkmPhi", &rawkstTrkmPhi);
  theTree->Branch("rawkstTrkmEta", &rawkstTrkmEta);
  theTree->Branch("rawkstTrkpPt", &rawkstTrkpPt);
  theTree->Branch("rawkstTrkpPhi", &rawkstTrkpPhi);
  theTree->Branch("rawkstTrkpEta", &rawkstTrkpEta);

  theTree->Branch("bMinusVtxCL", &bMinusVtxCL);
  theTree->Branch("bMinusCosAlphaBS", &bMinusCosAlphaBS);
  theTree->Branch("bPlusVtxCL", &bPlusVtxCL);
  theTree->Branch("bPlusCosAlphaBS", &bPlusCosAlphaBS);

  // ### Generated Observables ###
  theTree->Branch("genSignal", &genSignal, "genSignal/I");
  theTree->Branch("genMuMuBG", &genMuMuBG, "genMuMuBG/I");
  theTree->Branch("genMuMuBGnTrksm", &genMuMuBGnTrksm, "genMuMuBGnTrksm/I");
  theTree->Branch("genMuMuBGnTrksp", &genMuMuBGnTrksp, "genMuMuBGnTrksp/I");
  theTree->Branch("genPsiPrompt", &genPsiPrompt, "genPsiPrompt/O");
  theTree->Branch("genSignHasFSR", &genSignHasFSR, "genSignHasFSR/O");
  theTree->Branch("genSignKstHasFSR", &genSignKstHasFSR, "genSignKstHasFSR/O");
  theTree->Branch("genSignPsiHasFSR", &genSignPsiHasFSR, "genSignPsiHasFSR/O");

  // ### Generated Primary Vertex ###
  theTree->Branch("genPriVtxX", &genPriVtxX, "genPriVtxX/D");
  theTree->Branch("genPriVtxY", &genPriVtxY, "genPriVtxY/D");
  theTree->Branch("genPriVtxZ", &genPriVtxZ, "genPriVtxZ/D");

  // ### Generated B0 Mass ###
  theTree->Branch("genB0Mass", &genB0Mass, "genB0Mass/D");
  theTree->Branch("genB0Px", &genB0Px, "genB0Px/D");
  theTree->Branch("genB0Py", &genB0Py, "genB0Py/D");
  theTree->Branch("genB0Pz", &genB0Pz, "genB0Pz/D");

  // ### Generated B0 Vtx ###
  theTree->Branch("genB0VtxX", &genB0VtxX, "genB0VtxX/D");
  theTree->Branch("genB0VtxY", &genB0VtxY, "genB0VtxY/D");
  theTree->Branch("genB0VtxZ", &genB0VtxZ, "genB0VtxZ/D");

  // ### Generated K*0 Mass ###
  theTree->Branch("genKstMass", &genKstMass, "genKstMass/D");
  theTree->Branch("genKstPx", &genKstPx, "genKstPx/D");
  theTree->Branch("genKstPy", &genKstPy, "genKstPy/D");
  theTree->Branch("genKstPz", &genKstPz, "genKstPz/D");

  // ### Generated K*0 Vtx ###
  theTree->Branch("genKstVtxX", &genKstVtxX, "genKstVtxX/D");
  theTree->Branch("genKstVtxY", &genKstVtxY, "genKstVtxY/D");
  theTree->Branch("genKstVtxZ", &genKstVtxZ, "genKstVtxZ/D");

  // ### Generated J/psi or psi(2S) Mass and Vtx ###
  theTree->Branch("genPsiMass", &genPsiMass, "genPsiMass/D");
  theTree->Branch("genPsiVtxX", &genPsiVtxX, "genPsiVtxX/D");
  theTree->Branch("genPsiVtxY", &genPsiVtxY, "genPsiVtxY/D");
  theTree->Branch("genPsiVtxZ", &genPsiVtxZ, "genPsiVtxZ/D");

  // ### Generated mu- ###
  theTree->Branch("genMumMother", &genMumMother, "genMumMother/I");
  theTree->Branch("genMumPx", &genMumPx, "genMumPx/D");
  theTree->Branch("genMumPy", &genMumPy, "genMumPy/D");
  theTree->Branch("genMumPz", &genMumPz, "genMumPz/D");

  // ### Generated mu+ ###
  theTree->Branch("genMupMother", &genMupMother, "genMupMother/I");
  theTree->Branch("genMupPx", &genMupPx, "genMupPx/D");
  theTree->Branch("genMupPy", &genMupPy, "genMupPy/D");
  theTree->Branch("genMupPz", &genMupPz, "genMupPz/D");

  // ### Generated K*0 track- ###
  theTree->Branch("genKstTrkmMother", &genKstTrkmMother, "genKstTrkmMother/I");
  theTree->Branch("genKstTrkmID", &genKstTrkmID, "genKstTrkmID/I");
  theTree->Branch("genKstTrkmPx", &genKstTrkmPx, "genKstTrkmPx/D");
  theTree->Branch("genKstTrkmPy", &genKstTrkmPy, "genKstTrkmPy/D");
  theTree->Branch("genKstTrkmPz", &genKstTrkmPz, "genKstTrkmPz/D");

  // ### Generated K*0 track+ ###
  theTree->Branch("genKstTrkpMother", &genKstTrkpMother, "genKstTrkpMother/I");
  theTree->Branch("genKstTrkpID", &genKstTrkpID, "genKstTrkpID/I");
  theTree->Branch("genKstTrkpPx", &genKstTrkpPx, "genKstTrkpPx/D");
  theTree->Branch("genKstTrkpPy", &genKstTrkpPy, "genKstTrkpPy/D");
  theTree->Branch("genKstTrkpPz", &genKstTrkpPz, "genKstTrkpPz/D");

  // ### Matching Between Reconstructed and Generated ###
  theTree->Branch("truthMatchSignal", &truthMatchSignal);
  theTree->Branch("truthMatchMum", &truthMatchMum);
  theTree->Branch("truthMatchMup", &truthMatchMup);
  theTree->Branch("truthMatchTrkm", &truthMatchTrkm);
  theTree->Branch("truthMatchTrkp", &truthMatchTrkp);
}
