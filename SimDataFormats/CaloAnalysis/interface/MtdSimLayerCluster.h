// Author: Aurora Perego, Fabio Cossutti - aurora.perego@cern.ch, fabio.cossutti@ts.infn.it
// Date: 05/2023

#ifndef SimDataFormats_CaloAnalysis_MtdSimLayerCluster_h
#define SimDataFormats_CaloAnalysis_MtdSimLayerCluster_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include <vector>

//
// Forward declarations
//
class SimTrack;
class EncodedEventId;

class MtdSimLayerCluster : public SimCluster {
  friend std::ostream &operator<<(std::ostream &s, MtdSimLayerCluster const &tp);

public:
  MtdSimLayerCluster();

  MtdSimLayerCluster(const SimTrack &simtrk);
  MtdSimLayerCluster(EncodedEventId eventID, uint32_t particleID);  // for PU

  // destructor
  ~MtdSimLayerCluster();

  /** @brief add hit time */
  void addHitTime(float time) {
    times_.emplace_back(time);
    ++nsimhits_;
  }

  /** @brief computes the energy of the cluster */
  void addCluEnergy(float energy) { simLC_energy_ = energy; }

  /** @brief computes the position of the cluster */
  void addCluLocalPos(LocalPoint pos) { simLC_pos_ = pos; }

  /** @brief add the index of the simcluster */
  void addCluIndex(const uint32_t index) { seedId_ = index; }

  /** @brief returns the time of the cluster */
  float simTime() const { return simLC_time_; }

  /** @brief returns the local position of the cluster */
  LocalPoint simPos() const { return simLC_pos_; }

  /** @brief returns the accumulated sim energy in the cluster */
  float simEnergy() const { return simLC_energy_; }

  uint32_t seedId() const { return seedId_; }

private:
  // id of the simCluster it comes from
  uint32_t seedId_;

  uint64_t nsimhits_{0};
  EncodedEventId event_;

  uint32_t particleId_{0};
  float simLC_time_{0.f};

  float simLC_energy_{0.f};
  LocalPoint simLC_pos_;

  std::vector<uint32_t> hits_;
  std::vector<float> fractions_;
  std::vector<float> energies_;
  std::vector<float> times_;

  math::XYZTLorentzVectorF theMomentum_;

  // references to G4 and reco::GenParticle tracks
  std::vector<SimTrack> g4Tracks_;
  reco::GenParticleRefVector genParticles_;
};

#endif
