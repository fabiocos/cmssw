import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdValidation', Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
)

# Input: use revtx_step3.root from vertexTask_cfg.py output
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:revtx_step3.root'
    )
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- BTL Validation
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlDigiHitsValid + process.btlLocalRecoValid)

# --- ETL Validation
process.load("Validation.MtdValidation.etlSimHitsValid_cfi")
process.load("Validation.MtdValidation.etlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.etlLocalRecoValid_cfi")
etlValidation = cms.Sequence(process.etlSimHitsValid + process.etlDigiHitsValid + process.etlLocalRecoValid)

# --- Global Validation
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.load("Validation.MtdValidation.mtdEleIsoValid_cfi")
process.load("Validation.MtdValidation.vertices4DValid_cfi")

# Enable optional plots for vertex validation
process.vertices4DValid.optionalPlots = True

# ============================================================================
# Configure validation to use GNN vertex outputs
# ============================================================================
process.mtdTracksValid.inputTagV = 'offlinePrimaryVerticesGNN'
process.mtdTracksValid.t0SafePID = 'tofPIDGNN:t0safe'
process.mtdTracksValid.sigmat0SafePID = 'tofPIDGNN:sigmat0safe'
process.mtdTracksValid.sigmat0PID = 'tofPIDGNN:sigmat0'
process.mtdTracksValid.t0PID = 'tofPIDGNN:t0'

process.vertices4DValid.offline4DPV = 'offlinePrimaryVerticesGNN'
process.vertices4DValid.t0PID = 'tofPIDGNN:t0'
process.vertices4DValid.t0SafePID = 'tofPIDGNN:t0safe'
process.vertices4DValid.sigmat0SafePID = 'tofPIDGNN:sigmat0safe'
process.vertices4DValid.probPi = 'tofPIDGNN:probPi'
process.vertices4DValid.probK = 'tofPIDGNN:probK'
process.vertices4DValid.probP = 'tofPIDGNN:probP'

process.validation = cms.Sequence(
    btlValidation + 
    etlValidation + 
    process.mtdTracksValid + 
    process.mtdEleIsoValid + 
    process.vertices4DValid
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:mtdValidation_DQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path(process.mix + process.mtdTrackingRecHits + process.validation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(process.p, process.endjob_step, process.DQMoutput_step)

print("=" * 60)
print("mtdValidation_GNN_cfg.py for CMSSW_16_1_0_pre1")
print("Input: revtx_step3.root (GNN vertexing output)")
print("Output: mtdValidation_DQM.root")
print("=" * 60)
