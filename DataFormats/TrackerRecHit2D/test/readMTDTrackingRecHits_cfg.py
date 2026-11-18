import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.filePrepend = ""
options.parseArguments()

process = cms.Process("MTDREFREAD")
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring("file:"+options.filePrepend+"step3_mtdRereco.root"),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.reader = cms.EDAnalyzer(
    "TestReadMTDTrackingRecHits",
    src=cms.InputTag("mtdTrackingRecHits"),
)
process.p = cms.Path(process.reader)
