import FWCore.ParameterSet.Config as cms

_common_BTLparameters = cms.PSet(
  bxTime                   = cms.double(25),      # [ns]
  LightOutput              = cms.double(2000.),  # [photons/MeV], including Light Yield, Light Collection Efficincy and Photon Detection Efficiency
)

_barrel_MTDDigitizer = cms.PSet(
    digitizerName     = cms.string("BTLDigitizer"),
    inputSimHits      = cms.InputTag("g4SimHits:FastTimerHitsBarrel"),
    digiCollectionTag = cms.string("FTLBarrel"),
    maxSimHitsAccTime = cms.uint32(100),
    premixStage1      = cms.bool(False),
    premixStage1MinCharge = cms.double(1e-4),
    premixStage1MaxCharge = cms.double(1e6),
    DeviceSimulation = cms.PSet(
        _common_BTLparameters,
        LightCollectionSlope     = cms.double(0.075),   # [ns/cm]
        ),
    ElectronicsSimulation = cms.PSet(
        _common_BTLparameters,
        TestBeamMIPTimeRes         = cms.double(0.2697), # = 0.020[ns]*sqrt(7000.[npe]/38.5[ps])
        ScintillatorRiseTime       = cms.double(1.1),   # [ns]
        ScintillatorDecayTime      = cms.double(40.),   # [ns]
        ChannelTimeOffset          = cms.double(0.),    # [ns]
        smearChannelTimeOffset     = cms.double(0.),    # [ns]
        EnergyThreshold            = cms.double(4.),    # [photo-electrons]
        TimeThreshold1             = cms.double(20.),   # [photo-electrons]
        TimeThreshold2             = cms.double(50.),   # [photo-electrons]
        ReferencePulseNpe          = cms.double(100.),  # [photo-electrons]
        SigmaDigitization          = cms.double(0.007), # [ns]
        SigmaClock                 = cms.double(0.015), # [ns], 0.015 ps uncertainty on the combination of SiPMs
        DCRParam                   = cms.vdouble(6.234,30.,0.41), # 0.040[ns]*6000[pe]/38.5[ns], 30 [GHz], optimal exponent from fit to labo measurements
        DarkCountRate              = cms.double(10.), # [GHz]
        SlewRateParam              = cms.vdouble(5.32470e-01,0.,2.92152e+01,7.79368e+00), # parameterization of slew rate vs Gain * npe
        SigmaElectronicNoise       = cms.double(0.335), # [ns]
        SigmaElectronicNoiseConst  = cms.double(0.0167), # 0.0167[ns]
        ElectronicGain             = cms.double(0.0001457), # best gain / gain(3.5 Vov) / 9500. [pe]
        CorrelationCoefficient     = cms.double(1.),
        SmearTimeForOOTtails       = cms.bool(True),
        Npe_to_pC                  = cms.double(0.016), # [pC]
        Npe_to_V                   = cms.double(0.0064),# [V]

        # n bits for the ADC 
        adcNbits          = cms.uint32(10),
        # n bits for the TDC
        tdcNbits          = cms.uint32(10),
        # ADC saturation
        adcSaturation_MIP = cms.double(600.),           # [pC]
        # for different thickness
        adcThreshold_MIP   = cms.double(0.064),         # [pC]
        # LSB for time of arrival estimate from TDC
        toaLSB_ns         = cms.double(0.020),          # [ns]
        )


)

_endcap_MTDDigitizer = cms.PSet(
    digitizerName     = cms.string("ETLDigitizer"),
    inputSimHits      = cms.InputTag("g4SimHits:FastTimerHitsEndcap"),
    digiCollectionTag = cms.string("FTLEndcap"),
    maxSimHitsAccTime = cms.uint32(100),
    premixStage1      = cms.bool(False),
    premixStage1MinCharge = cms.double(1e-4),
    premixStage1MaxCharge = cms.double(1e6),
    DeviceSimulation  = cms.PSet(
        bxTime            = cms.double(25),
        tofDelay          = cms.double(1),
        meVPerMIP         = cms.double(0.085), # from HGCal
        ),
    ElectronicsSimulation = cms.PSet(
        bxTime               = cms.double(25),
        IntegratedLuminosity = cms.double(1000.),      # [1/fb]
        FluenceVsRadius      = cms.string("1.937*TMath::Power(x,-1.706)"),
        LGADGainVsFluence    = cms.string("TMath::Min(15.,30.-x)"),
        TimeResolution2      = cms.string("0.0225/x"), # [ns^2]
        # n bits for the ADC 
        adcNbits             = cms.uint32(8),
        # n bits for the TDC
        tdcNbits             = cms.uint32(11),
        # ADC saturation
        adcSaturation_MIP  = cms.double(25),
        # for different thickness
        adcThreshold_MIP   = cms.double(0.025),
        # LSB for time of arrival estimate from TDC in ns
        toaLSB_ns          = cms.double(0.013),
        )
)

from Configuration.Eras.Modifier_phase2_etlV4_cff import phase2_etlV4
phase2_etlV4.toModify(_endcap_MTDDigitizer.DeviceSimulation, meVPerMIP = 0.015 )

from Configuration.ProcessModifiers.premix_stage1_cff import premix_stage1
for _m in [_barrel_MTDDigitizer, _endcap_MTDDigitizer]:
    premix_stage1.toModify(_m, premixStage1 = True)

# Fast Timing
mtdDigitizer = cms.PSet( 
    accumulatorType   = cms.string("MTDDigiProducer"),
    makeDigiSimLinks  = cms.bool(False),
    verbosity         = cms.untracked.uint32(0),

    barrelDigitizer = _barrel_MTDDigitizer,
    endcapDigitizer = _endcap_MTDDigitizer
)
