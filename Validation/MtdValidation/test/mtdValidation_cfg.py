import FWCore.ParameterSet.Config as cms


from Configuration.StandardSequences.Eras import eras
process = cms.Process('mtdValidation',eras.Phase2C4_timing_layer_bar)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2023D38Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
# process.MessageLogger.cerr.MtdTracksValidation = cms.untracked.PSet(
    # limit = cms.untracked.int32(-1)
# )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
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
'file:step3_1.root',
'file:step3_20.root',
'file:step3_21.root',
'file:step3_22.root',
'file:step3_23.root',
'file:step3_24.root',
'file:step3_25.root',
'file:step3_26.root',
'file:step3_27.root',
'file:step3_28.root',
'file:step3_29.root',
'file:step3_2.root',
'file:step3_30.root',
'file:step3_31.root',
'file:step3_32.root',
'file:step3_33.root',
'file:step3_34.root',
'file:step3_35.root',
'file:step3_36.root',
'file:step3_37.root',
'file:step3_38.root',
'file:step3_39.root',
'file:step3_3.root',
'file:step3_40.root',
'file:step3_41.root',
'file:step3_42.root',
'file:step3_43.root',
'file:step3_44.root',
'file:step3_45.root',
'file:step3_46.root',
'file:step3_47.root',
'file:step3_48.root',
'file:step3_49.root',
'file:step3_4.root',
'file:step3_50.root',
'file:step3_5.root',
'file:step3_6.root',
'file:step3_7.root',
'file:step3_8.root',
'file:step3_9.root')
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- BTL Validation
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlLocalRecoValid)


# --- Global Validation
process.load("Validation.MtdValidation.mtdTracksValid_cfi")

# process.btlDigiHitsValid.optionalPlots = True
# process.etlDigiHitsValid.optionalPlots = True
# process.btlLocalRecoValid.optionalPlots = True
# process.etlLocalRecoValid.optionalPlots = True
# process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True

process.validation = cms.Sequence(btlValidation + process.mtdTracksValid)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path( process.mix + process.mtdTrackingRecHits + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
