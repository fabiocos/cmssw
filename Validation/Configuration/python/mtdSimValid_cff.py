import FWCore.ParameterSet.Config as cms

# MTD validation sequences
from Validation.MtdValidation.btlSimHitsValid_cfi import btlSimHitsValid
from Validation.MtdValidation.btlLocalRecoValid_cfi import btlLocalRecoValid
from Validation.MtdValidation.mtdTracksValid_cfi import mtdTracksValid

mtdSimValid  = cms.Sequence(btlSimHitsValid)
mtdRecoValid = cms.Sequence(btlLocalRecoValid  + mtdTracksValid)
