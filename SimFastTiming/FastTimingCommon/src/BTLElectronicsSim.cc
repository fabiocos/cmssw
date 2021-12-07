#include "SimFastTiming/FastTimingCommon/interface/BTLElectronicsSim.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"

using namespace mtd;

BTLElectronicsSim::BTLElectronicsSim(const edm::ParameterSet& pset, edm::ConsumesCollector iC)
    : debug_(pset.getUntrackedParameter<bool>("debug", false)),
      bxTime_(pset.getParameter<double>("bxTime")),
      testBeamMIPTimeRes_(pset.getParameter<double>("TestBeamMIPTimeRes")),
      scintillatorRiseTime_(pset.getParameter<double>("ScintillatorRiseTime")),
      scintillatorDecayTime_(pset.getParameter<double>("ScintillatorDecayTime")),
      channelTimeOffset_(pset.getParameter<double>("ChannelTimeOffset")),
      smearChannelTimeOffset_(pset.getParameter<double>("smearChannelTimeOffset")),
      energyThreshold_(pset.getParameter<double>("EnergyThreshold")),
      timeThreshold1_(pset.getParameter<double>("TimeThreshold1")),
      timeThreshold2_(pset.getParameter<double>("TimeThreshold2")),
      referencePulseNpe_(pset.getParameter<double>("ReferencePulseNpe")),
      sigmaDigitization_(pset.getParameter<double>("SigmaDigitization")),
      sigmaClock_(pset.getParameter<double>("SigmaClock")),
      paramDCR_(pset.getParameter<std::vector<double> >("DCRParam")),
      darkCountRate_(pset.getParameter<double>("DarkCountRate")),
      paramSR_(pset.getParameter<std::vector<double> >("SlewRateParam")),
      sigmaElectronicNoise_(pset.getParameter<double>("SigmaElectronicNoise")),
      sigmaElectronicNoiseConst_(pset.getParameter<double>("SigmaElectronicNoiseConst")),
      electronicGain_(pset.getParameter<double>("ElectronicGain")),
      smearTimeForOOTtails_(pset.getParameter<bool>("SmearTimeForOOTtails")),
      npe_to_pC_(pset.getParameter<double>("Npe_to_pC")),
      npe_to_V_(pset.getParameter<double>("Npe_to_V")),
      adcNbits_(pset.getParameter<uint32_t>("adcNbits")),
      tdcNbits_(pset.getParameter<uint32_t>("tdcNbits")),
      adcSaturation_MIP_(pset.getParameter<double>("adcSaturation_MIP")),
      adcBitSaturation_(std::pow(2, adcNbits_) - 1),
      adcLSB_MIP_(adcSaturation_MIP_ / adcBitSaturation_),
      adcThreshold_MIP_(pset.getParameter<double>("adcThreshold_MIP")),
      toaLSB_ns_(pset.getParameter<double>("toaLSB_ns")),
      tdcBitSaturation_(std::pow(2, tdcNbits_) - 1),
      corrCoeff_(pset.getParameter<double>("CorrelationCoefficient")),
      cosPhi_(0.5 * (sqrt(1. + corrCoeff_) + sqrt(1. - corrCoeff_))),
      sinPhi_(0.5 * corrCoeff_ / cosPhi_),
      scintillatorDecayTime2_(scintillatorDecayTime_ * scintillatorDecayTime_),
      scintillatorDecayTimeInv_(1. / scintillatorDecayTime_),
      sigmaElectronicNoiseConst2_(sigmaElectronicNoiseConst_ * sigmaElectronicNoiseConst_),
      sigmaConst2_(sigmaDigitization_ * sigmaDigitization_ + sigmaClock_ * sigmaClock_) {
#ifdef EDM_ML_DEBUG
  float lightOutput = 4.2f * pset.getParameter<double>("LightOutput");  // average npe for 4.2 MeV
  float s1 = sigma_stochastic(lightOutput);
  float s2 = sigma_DCR(lightOutput);
  float s3 = sigma_electronics(lightOutput);
  float s4 = sigmaDigitization_;
  float s5 = sigmaClock_;
  LogDebug("BTLElectronicsSim") << " BTL resolution model, for an average light output of " << std::fixed
                                << std::setw(14) << lightOutput << " :"
                                << "\n sigma stochastic   = " << std::setw(14) << sigma_stochastic(lightOutput)
                                << "\n sigma DCR          = " << std::setw(14) << sigma_DCR(lightOutput)
                                << "\n sigma electronics  = " << std::setw(14) << sigma_electronics(lightOutput)
                                << "\n sigma digitization = " << std::setw(14) << sigmaDigitization_
                                << "\n sigma clock        = " << std::setw(14) << sigmaClock_
                                << "\n ---------------------"
                                << "\n sigma total        = " << std::setw(14)
                                << std::sqrt(s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4 + s5 * s5);
#endif
}

void BTLElectronicsSim::run(const mtd::MTDSimHitDataAccumulator& input,
                            BTLDigiCollection& output,
                            CLHEP::HepRandomEngine* hre) const {
  MTDSimHitData chargeColl, toa1, toa2;

  for (MTDSimHitDataAccumulator::const_iterator it = input.begin(); it != input.end(); it++) {
    // --- Digitize only the in-time bucket:
    const unsigned int iBX = mtd_digitizer::kInTimeBX;

    chargeColl.fill(0.f);
    toa1.fill(0.f);
    toa2.fill(0.f);
    for (size_t iside = 0; iside < 2; iside++) {
      // --- Fluctuate the total number of photo-electrons
      float npe = CLHEP::RandPoissonQ::shoot(hre, (it->second).hit_info[2 * iside][iBX]);
      if (npe < energyThreshold_)
        continue;

      // --- Get the time of arrival and add a channel time offset
      float finalToA1 = (it->second).hit_info[1 + 2 * iside][iBX] + channelTimeOffset_;

      if (smearChannelTimeOffset_ > 0.) {
        float timeSmearing = CLHEP::RandGaussQ::shoot(hre, 0., smearChannelTimeOffset_);
        finalToA1 += timeSmearing;
      }

      // --- Calculate and add the time walk: the time of arrival is read in correspondence
      //                                      with two thresholds on the signal pulse
      std::array<float, 3> times =
          btlPulseShape_.timeAtThr(npe / referencePulseNpe_, timeThreshold1_ * npe_to_V_, timeThreshold2_ * npe_to_V_);

      // --- If the pulse amplitude is smaller than timeThreshold2, the trigger does not fire
      if (times[1] == 0.)
        continue;

      float finalToA2 = finalToA1 + times[1];
      finalToA1 += times[0];

      // --- Estimate the time uncertainty due to photons from earlier OOT hits in the current BTL cell
      if (smearTimeForOOTtails_) {
        float rate_oot = 0.;
        // Loop on earlier OOT hits
        for (int ibx = 0; ibx < mtd_digitizer::kInTimeBX; ++ibx) {
          if ((it->second).hit_info[2 * iside][ibx] > 0.) {
            float hit_time = (it->second).hit_info[1 + 2 * iside][ibx] + bxTime_ * (ibx - mtd_digitizer::kInTimeBX);
            float npe_oot = CLHEP::RandPoissonQ::shoot(hre, (it->second).hit_info[2 * iside][ibx]);
            rate_oot += npe_oot * exp(hit_time * scintillatorDecayTimeInv_) * scintillatorDecayTimeInv_;
          }
        }  // ibx loop

        if (rate_oot > 0.) {
          float sigma_oot = sqrt(rate_oot * scintillatorRiseTime_) * scintillatorDecayTime_ / npe;
          float smearing_oot = CLHEP::RandGaussQ::shoot(hre, 0., sigma_oot);
          finalToA1 += smearing_oot;
          finalToA2 += smearing_oot;
        }
      }  // if smearTimeForOOTtails_

      // --- Stochastich term, uncertainty due to the fluctuations of the n-th photon arrival time:
      if (testBeamMIPTimeRes_ > 0.) {
        // In this case the time resolution is parametrized from the testbeam.
        // The same parameterization is used for both thresholds.
        // the uncertainty is provided for the combination of two SiPMs
        float sigma = sqrt2_ * sigma_stochastic(npe);
        float smearing_stat_thr1 = CLHEP::RandGaussQ::shoot(hre, 0., sigma);
        float smearing_stat_thr2 = CLHEP::RandGaussQ::shoot(hre, 0., sigma);

        finalToA1 += smearing_stat_thr1;
        finalToA2 += smearing_stat_thr2;

      } else {
        // In this case the time resolution is taken from the literature.
        // The fluctuations due to the first timeThreshold1_ p.e. are common to both times
        float smearing_stat_thr1 =
            CLHEP::RandGaussQ::shoot(hre, 0., scintillatorDecayTime_ * sqrt(sigma2_pe(timeThreshold1_, npe)));
        float smearing_stat_thr2 = CLHEP::RandGaussQ::shoot(
            hre, 0., scintillatorDecayTime_ * sqrt(sigma2_pe(timeThreshold2_ - timeThreshold1_, npe)));
        finalToA1 += smearing_stat_thr1;
        finalToA2 += smearing_stat_thr1 + smearing_stat_thr2;
      }

      // --- Add in quadrature the uncertainties due to the SiPM timing resolution, the SiPM DCR,
      //     the electronic noise and the clock distribution:

      float sigmaDCR = sigma_DCR(npe);
      float sigmaElec = sigma_electronics(npe);
      float sigma2_tot_thr1 = sigmaDCR * sigmaDCR + sigmaElec * sigmaElec;

      // --- Add in quadrature uncertainties independent on npe: digitization and clock distribution

      sigma2_tot_thr1 += sigmaConst2_;
      sigma2_tot_thr1 *= 2.f;  // all uncertainties are provided for a combination of two SiPMs

      float sigma2_tot_thr2 = sigma2_tot_thr1;

      // --- Smear the arrival times using the correlated uncertainties:
      float smearing_thr1_uncorr = CLHEP::RandGaussQ::shoot(hre, 0., sqrt(sigma2_tot_thr1));
      float smearing_thr2_uncorr = CLHEP::RandGaussQ::shoot(hre, 0., sqrt(sigma2_tot_thr2));

      finalToA1 += cosPhi_ * smearing_thr1_uncorr + sinPhi_ * smearing_thr2_uncorr;
      finalToA2 += sinPhi_ * smearing_thr1_uncorr + cosPhi_ * smearing_thr2_uncorr;

      chargeColl[iside] = npe * npe_to_pC_;  // the p.e. number is here converted to pC

      toa1[iside] = finalToA1;
      toa2[iside] = finalToA2;

    }  // iside loop

    //run the shaper to create a new data frame
    BTLDataFrame rawDataFrame(it->first.detid_);
    runTrivialShaper(rawDataFrame, chargeColl, toa1, toa2, it->first.row_, it->first.column_);
    updateOutput(output, rawDataFrame);

  }  // MTDSimHitDataAccumulator loop
}

void BTLElectronicsSim::runTrivialShaper(BTLDataFrame& dataFrame,
                                         const mtd::MTDSimHitData& chargeColl,
                                         const mtd::MTDSimHitData& toa1,
                                         const mtd::MTDSimHitData& toa2,
                                         const uint8_t row,
                                         const uint8_t col) const {
  bool debug = debug_;
#ifdef EDM_ML_DEBUG
  for (int it = 0; it < (int)(chargeColl.size()); it++)
    debug |= (chargeColl[it] > adcThreshold_MIP_);
#endif

  if (debug)
    edm::LogVerbatim("BTLElectronicsSim") << "[runTrivialShaper]" << std::endl;

  //set new ADCs
  for (int it = 0; it < (int)(chargeColl.size()); it++) {
    BTLSample newSample;
    newSample.set(false, false, 0, 0, 0, row, col);

    //brute force saturation, maybe could to better with an exponential like saturation
    const uint32_t adc = std::min((uint32_t)std::floor(chargeColl[it] / adcLSB_MIP_), adcBitSaturation_);
    const uint32_t tdc_time1 = std::min((uint32_t)std::floor(toa1[it] / toaLSB_ns_), tdcBitSaturation_);
    const uint32_t tdc_time2 = std::min((uint32_t)std::floor(toa2[it] / toaLSB_ns_), tdcBitSaturation_);

    newSample.set(
        chargeColl[it] > adcThreshold_MIP_, tdc_time1 == tdcBitSaturation_, tdc_time2, tdc_time1, adc, row, col);
    dataFrame.setSample(it, newSample);

    if (debug)
      edm::LogVerbatim("BTLElectronicsSim") << adc << " (" << chargeColl[it] << "/" << adcLSB_MIP_ << ") ";
  }

  if (debug) {
    std::ostringstream msg;
    dataFrame.print(msg);
    edm::LogVerbatim("BTLElectronicsSim") << msg.str() << std::endl;
  }
}

void BTLElectronicsSim::updateOutput(BTLDigiCollection& coll, const BTLDataFrame& rawDataFrame) const {
  BTLDataFrame dataFrame(rawDataFrame.id());
  dataFrame.resize(dfSIZE);
  bool putInEvent(false);
  for (int it = 0; it < dfSIZE; ++it) {
    dataFrame.setSample(it, rawDataFrame[it]);
    if (it == 0)
      putInEvent = rawDataFrame[it].threshold();
  }

  if (putInEvent) {
    coll.push_back(dataFrame);
  }
}

float BTLElectronicsSim::sigma2_pe(const float& Q, const float& R) const {
  float OneOverR = 1. / R;
  float OneOverR2 = OneOverR * OneOverR;

  // --- This is Eq. (17) from Nucl. Instr. Meth. A 564 (2006) 185
  float sigma2 = Q * OneOverR2 *
                 (1. + 2. * (Q + 1.) * OneOverR + (Q + 1.) * (6. * Q + 11) * OneOverR2 +
                  (Q + 1.) * (Q + 2.) * (2. * Q + 5.) * OneOverR2 * OneOverR);

  return sigma2;
}

float BTLElectronicsSim::sigma_stochastic(const float& npe) const {
  return testBeamMIPTimeRes_ * std::sqrt(scintillatorDecayTime_ / npe);
}

float BTLElectronicsSim::sigma_DCR(const float& npe) const {
  // trick to safely switch off the electronics contribution for resolution studies

  if (darkCountRate_ == 0.) {
    return 0.;
  }
  return paramDCR_[0] * std::pow((darkCountRate_ / paramDCR_[1]), paramDCR_[2]) * scintillatorDecayTime_ / npe;
}

float BTLElectronicsSim::sigma_electronics(const float npe) const {
  // trick to safely switch off the electronics contribution for resolution studies

  if (electronicGain_ == 0.) {
    return 0.;
  }

  float gainXnpe = electronicGain_ * npe;
  float res = sigmaElectronicNoise_ / sqrt2_;
  if (gainXnpe < paramSR_[0]) {
    res /= (paramSR_[2] * gainXnpe + paramSR_[1]);
  } else {
    res /= (paramSR_[3] * std::log(gainXnpe) + paramSR_[2] * paramSR_[0] - paramSR_[3] * std::log(paramSR_[0]));
  }
  return std::sqrt(res * res + sigmaElectronicNoiseConst2_);
}
