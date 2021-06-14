## import skeleton process
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger

process = cms.Process("PAT", stage2L1Trigger, Phase2C11I13M9)

process.Tracer = cms.Service("Tracer")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3_inMINIAODSIM.root'),
    secondaryFileNames = cms.untracked.vstring()
)

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.load("Configuration.StandardSequences.MagneticField_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patTuple.root'),
    ## save only events passing the full path
    #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    ## save PAT output; you need a '*' to unpack the list of commands
    ## 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
    )

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)
process.outpath = cms.EndPath(process.out, patAlgosToolsTask)

## uncomment the following line to update different jet collections
## and add them to the event content
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

from CommonTools.PileupAlgos.Puppi_cff import *
process.puppi4D = puppi.clone(
    candName = cms.InputTag('packedPFCandidates'),
    vertexName = cms.InputTag('offlineSlimmedPrimaryVertices4D'))
patAlgosToolsTask.add(process.puppi4D)

from RecoJets.Configuration.RecoPFJets_cff import *
process.ak4PFJetsPuppi4D = ak4PFJetsPuppi.clone(
    src = cms.InputTag("puppi4D"),
    srcPVs = cms.InputTag("offlineSlimmedPrimaryVertices4D"),
    applyWeight = False,
    # srcWeights = cms.InputTag("puppi4D")
)
patAlgosToolsTask.add(process.ak4PFJetsPuppi4D)

labelAK4PFPUPPI4D = "Puppi4D"
addJetCollection(
    process,
    labelName = labelAK4PFPUPPI4D,
    jetSource = cms.InputTag('ak4PFJetsPuppi4D'),
    algo = 'AK',
    rParam = 0.4,
    pfCandidates=cms.InputTag('puppi4D'),
    pvSource=cms.InputTag('offlineSlimmedPrimaryVertices4D'),
    genJetCollection = cms.InputTag('slimmedGenJets'),
    jetCorrections = ('AK4PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = ['None'], # turn-off b tagging
    getJetMCFlavour = False # jet flavor needs to be disabled for groomed fat jets
)

process.patJetPartonMatchPuppi4D.matched = cms.InputTag("prunedGenParticles")

from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent
process.out.outputCommands = MINIAODSIMEventContent.outputCommands
process.out.outputCommands.append('keep *_selectedPatJetsPuppi4D_*_*')
process.out.outputCommands.append('drop CaloTowers_selectedPatJetsPuppi4D_*_*')
#                                         ##
process.out.fileName = 'puppi4D_fromMiniAOD.root'
#                                         ##
# process.options.wantSummary = False   ##  (to suppress the long output at the end of the job)
