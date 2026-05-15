import FWCore.ParameterSet.Config as cms

from SimFastTiming.FastTimingCommon.mtdDigitizer_cfi import mtdDigitizer

_barrelAlgo = cms.PSet(
    algoName = cms.string("BTLUncalibRecHitAlgo"),
    invLightSpeedLYSO = mtdDigitizer.barrelDigitizer.DeviceSimulation.LightCollectionSlope, # [ns/cm]
    npeToADC = mtdDigitizer.barrelDigitizer.ElectronicsSimulation.PulseQParam, # Npe to ADC counts conversion
    npePerMeV = mtdDigitizer.barrelDigitizer.DeviceSimulation.LightOutput, # [Npe/MeV]
    npeSaturationCorrection = mtdDigitizer.barrelDigitizer.ElectronicsSimulation.SiPMSaturationParam, # effective Npe to true Npe conversion (SiPM saturation)
    tdcLSB_ns = cms.double(0.020), # [ns]
    timeResolutionInNs = cms.string("0.0593858*pow(x,-1.02826)+0.0156719"), # [ns]
    timeWalkCorrection = cms.string("1.16946*pow(x,-0.671018)+0.0443454") # [ns]
)

_endcapAlgo = cms.PSet(
    algoName      = cms.string("ETLUncalibRecHitAlgo"),
    adcNbits      = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcNbits,
    adcSaturation = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcSaturation_MIP,
    toaLSB_ns     = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.toaLSB_ns,
    timeResolutionInNs = cms.string("0.0370"), # [ns]
    timeCorr_p0 = cms.double(0.967683), # 0.974683 - 0.007, ad hoc correction for bias from global delay removal
    timeCorr_p1 = cms.double(-0.237274),
    timeCorr_p2 = cms.double(0.021455),
    timeCorr_p3 = cms.double(-0.000727429)
)

mtdUncalibratedRecHits = cms.EDProducer(
    "MTDUncalibratedRecHitProducer",
    barrel = _barrelAlgo,
    endcap = _endcapAlgo,
    barrelDigis = cms.InputTag('mix:FTLBarrel'),
    endcapDigis = cms.InputTag('mix:FTLEndcap'),
    BarrelHitsName = cms.string('FTLBarrel'),
    EndcapHitsName = cms.string('FTLEndcap')
)

btlBaseRecHitsSoA = cms.EDProducer(
    'btlrechit::BTLBaseRecHitSoAProducer@alpaka',
    digi = cms.InputTag("mix", "FTLBarrelSoA"),
    npeToADC0 = cms.double(-22.5),
    npeToADC1 = cms.double(0.0348), # Npe to ADC counts conversion
    npeSaturationCorr0 = cms.double(-8.54e-06),
    npeSaturationCorr1 = cms.double(1.034),
    npePerMeV = mtdDigitizer.barrelDigitizer.DeviceSimulation.LightOutput # [Npe/MeV]
)

etlBaseRecHitsSoA = cms.EDProducer(
    'etlrechit::ETLBaseRecHitSoAProducer@alpaka',
    digi = cms.InputTag("mix", "FTLEndcapSoA"),
    adcNbits      = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcNbits,
    adcSaturation = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcSaturation_MIP,
    toaLSB_ns     = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.toaLSB_ns,
    tdcWindowStart     = cms.double(9.375), # now set to 3 x ETROC_clock, phase can be adjusted to set the start at any value
    timeCorr_p0 = cms.double(0.967683), # 0.974683 - 0.007, ad hoc correction for bias from global delay removal
    timeCorr_p1 = cms.double(-0.237274),
    timeCorr_p2 = cms.double(0.021455),
    timeCorr_p3 = cms.double(-0.000727429)
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(mtdUncalibratedRecHits,
    barrelDigis = 'mixData:FTLBarrel',
    endcapDigis = 'mixData:FTLEndcap',
)
