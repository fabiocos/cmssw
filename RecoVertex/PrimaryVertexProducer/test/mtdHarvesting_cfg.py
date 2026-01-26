import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting', Phase2C17I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source - DQM file from mtdValidation_GNN_cfg.py
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:mtdValidation_DQM.root')
)

# Path and EndPath definitions
process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing for all MTD validators
process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

process.harvesting = cms.Sequence(
    process.btlSimHitsPostProcessor + 
    process.btlLocalRecoPostProcessor + 
    process.MtdTracksPostProcessor + 
    process.MtdEleIsoPostProcessor + 
    process.Primary4DVertexPostProcessor
)

process.p = cms.Path(process.harvesting)

process.schedule = cms.Schedule(process.edmtome_step, process.p, process.dqmsave_step)

print("=" * 60)
print("mtdHarvesting_cfg.py for CMSSW_16_1_0_pre1")
print("Input: mtdValidation_DQM.root")
print("Output: DQM_V0001_R00000001__Global__CMSSW_X_Y_Z__RECO.root")
print("=" * 60)
