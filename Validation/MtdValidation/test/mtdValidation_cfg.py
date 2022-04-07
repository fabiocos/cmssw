import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('mtdValidation',Phase2C11I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.Geometry.GeometryExtended2026D77Reco_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    # reportEvery = cms.untracked.int32(100),
    reportEvery = cms.untracked.int32(1),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:step3_1.root',
        'file:step3_2.root',
        'file:step3_3.root',
        'file:step3_4.root',
        'file:step3_5.root',
        'file:step3_6.root',
        'file:step3_7.root',
        'file:step3_8.root',
        'file:step3_9.root',
        'file:step3_10.root',
        'file:step3_11.root',
        'file:step3_12.root',
        'file:step3_13.root',
        'file:step3_14.root',
        'file:step3_15.root',
        'file:step3_16.root',
        'file:step3_17.root',
        'file:step3_18.root',
        'file:step3_19.root',
        'file:step3_20.root'
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
process.load("Validation.MtdValidation.vertices4DValid_cfi")
process.load("Validation.MtdValidation.mtdSecondaryPvValid_cfi")

# --- Optional plots activation

# process.btlDigiHitsValid.LocalPositionDebug = True
# process.etlDigiHitsValid.LocalPositionDebug = True
# process.btlLocalRecoValid.LocalPositionDebug = True
# process.etlLocalRecoValid.LocalPositionDebug = True
process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True
process.mtdSecondaryPvValid.optionalPlots = True
process.mtdSecondaryPvValid.printMsg = cms.untracked.bool(True)

# process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.vertices4DValid)
process.validation = cms.Sequence( process.mtdTracksValid + process.vertices4DValid + process.mtdSecondaryPvValid)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path( process.mix + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
