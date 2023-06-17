import FWCore.ParameterSet.Config as cms


from Configuration.StandardSequences.Eras import eras
process = cms.Process('mtdHarvesting',eras.Phase2C4_timing_layer_bar)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2023D38Reco_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:step3_inDQM.root')
)

# Path and EndPath definitions

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")

process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor)

process.p = cms.Path( process.harvesting )

process.schedule = cms.Schedule( process.edmtome_step , process.p , process.dqmsave_step )
