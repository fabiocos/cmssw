import FWCore.ParameterSet.Config as cms

##____________________________________________________________________________||
from RecoMET.METProducers.pfMet_cfi import *
from RecoMET.METProducers.pfChMet_cfi import *
from CommonTools.PileupAlgos.Puppi_cff import puppiNoLep
from RecoMET.METProducers.pfMetPuppi_cfi import *

##____________________________________________________________________________||
recoPFMETTask = cms.Task(pfMet , particleFlowForChargedMET , pfChMet, puppiNoLep, pfMetPuppi)
recoPFMET = cms.Sequence(recoPFMETTask)


from CommonTools.PileupAlgos.Puppi_cff import puppi4DNoLep
from Configuration.Eras.Modifier_phase2_timing_layer_cff import phase2_timing_layer
pfMetPuppi4D = pfMetPuppi.clone()
phase2_timing_layer.toModify(pfMetPuppi4D, srcWeights = cms.InputTag("puppi4DNoLep") )

recoPFMET4DOnlyTask = cms.Task(puppi4DNoLep,
                           pfMetPuppi4D
   )
recoPFMET4DOnly = cms.Sequence(recoPFMET4DOnlyTask)
recoPFMET4DTask = recoPFMETTask.copy()
recoPFMET4DTask.add(puppi4DNoLep, pfMetPuppi4D)
phase2_timing_layer.toReplaceWith(recoPFMETTask, recoPFMET4DTask)

##____________________________________________________________________________||
from Configuration.ProcessModifiers.pp_on_AA_cff import pp_on_AA
pp_on_AA.toModify(pfMet,  globalThreshold = 999.)
pp_on_AA.toModify(pfChMet, globalThreshold = 999.)
pp_on_AA.toModify(pfMetPuppi,  globalThreshold = 999.)
