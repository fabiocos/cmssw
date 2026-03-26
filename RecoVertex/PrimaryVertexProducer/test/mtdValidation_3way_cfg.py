"""
mtdValidation_3way_cfg.py - MTD Validation with 3D, 4D and GNN vertex comparisons

Runs Primary4DVertexValidation and MtdTracksValidation for all three vertex collections
in a single job, writing to separate DQM folders for comparison.

Usage: cmsRun mtdValidation_3way_cfg.py backend=alpaka
"""

import FWCore.ParameterSet.Config as cms
import sys

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdValidation', Phase2C17I13M9)

# ============================================================================
# Parse command line for backend selection
# ============================================================================
backend = "alpaka"  # default
for arg in sys.argv[1:]:
    if arg.startswith("backend="):
        backend = arg.split("=")[1].lower()

if backend not in ["onnx", "alpaka"]:
    raise ValueError(f"Unknown backend: {backend}. Use 'onnx' or 'alpaka'")

input_file = f"file:revtx_step3_{backend}_v17s7.root"
output_file = f"file:mtdValidation_3way_{backend}_v17s7.root"

print("=" * 60)
print(f"MTD 3-Way Validation (3D + 4D + GNN) - Backend: {backend.upper()}")
print(f"Input:  {input_file}")
print(f"Output: {output_file}")
print("=" * 60)

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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

# Multi-threaded
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
)

# Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(input_file)
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# ============================================================================
# BTL/ETL Validation (shared, run once)
# ============================================================================
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlDigiHitsValid + process.btlLocalRecoValid)

process.load("Validation.MtdValidation.etlSimHitsValid_cfi")
process.load("Validation.MtdValidation.etlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.etlLocalRecoValid_cfi")
etlValidation = cms.Sequence(process.etlSimHitsValid + process.etlDigiHitsValid + process.etlLocalRecoValid)

process.load("Validation.MtdValidation.mtdEleIsoValid_cfi")

# ============================================================================
# Load base modules (we'll clone from these)
# ============================================================================
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.load("Validation.MtdValidation.vertices4DValid_cfi")

# ============================================================================
# 3D Vertex Validation
# Uses: offlinePrimaryVertices (3D, no timing)
# tofPID: tofPID3D (3D PID)
# ============================================================================
process.vertices3DValid = process.vertices4DValid.clone(
    folder = cms.string('MTD/Vertices/3D/'),
    offline4DPV = 'offlinePrimaryVertices',
    t0PID = 'tofPID3D:t0',
    t0SafePID = 'tofPID3D:t0safe',
    sigmat0SafePID = 'tofPID3D:sigmat0safe',
    probPi = 'tofPID3D:probPi',
    probK = 'tofPID3D:probK',
    probP = 'tofPID3D:probP',
    optionalPlots = True,
)

process.mtdTracks3DValid = process.mtdTracksValid.clone(
    folder = cms.string('MTD/Tracks/3D/'),
    inputTagV = 'offlinePrimaryVertices',
    t0PID = 'tofPID3D:t0',
    t0SafePID = 'tofPID3D:t0safe',
    sigmat0SafePID = 'tofPID3D:sigmat0safe',
    sigmat0PID = 'tofPID3D:sigmat0',
)

# ============================================================================
# 4D Vertex Validation
# Uses: offlinePrimaryVertices4D (4D, with timing)
# tofPID: tofPID (default 4D PID)
# ============================================================================
process.vertices4DValid.folder = cms.string('MTD/Vertices/4D/')
process.vertices4DValid.offline4DPV = 'offlinePrimaryVertices4D'
process.vertices4DValid.t0PID = 'tofPID:t0'
process.vertices4DValid.t0SafePID = 'tofPID:t0safe'
process.vertices4DValid.sigmat0SafePID = 'tofPID:sigmat0safe'
process.vertices4DValid.probPi = 'tofPID:probPi'
process.vertices4DValid.probK = 'tofPID:probK'
process.vertices4DValid.probP = 'tofPID:probP'
process.vertices4DValid.optionalPlots = True

process.mtdTracksValid.folder = cms.string('MTD/Tracks/4D/')
process.mtdTracksValid.inputTagV = 'offlinePrimaryVertices4D'
process.mtdTracksValid.t0PID = 'tofPID:t0'
process.mtdTracksValid.t0SafePID = 'tofPID:t0safe'
process.mtdTracksValid.sigmat0SafePID = 'tofPID:sigmat0safe'
process.mtdTracksValid.sigmat0PID = 'tofPID:sigmat0'

# ============================================================================
# GNN Vertex Validation
# Uses: offlinePrimaryVerticesGNN (GNN)
# tofPID: tofPIDGNN
# ============================================================================
process.verticesGNNValid = process.vertices4DValid.clone(
    folder = cms.string('MTD/Vertices/GNN/'),
    offline4DPV = 'offlinePrimaryVerticesGNN',
    t0PID = 'tofPIDGNN:t0',
    t0SafePID = 'tofPIDGNN:t0safe',
    sigmat0SafePID = 'tofPIDGNN:sigmat0safe',
    probPi = 'unsortedOfflinePrimaryVerticesGNN:gnnPiWeight0',
    probK = 'unsortedOfflinePrimaryVerticesGNN:gnnPiWeight1',
    probP = 'unsortedOfflinePrimaryVerticesGNN:gnnPiWeight2',
    trackweightTh = 0.0,  # match producer's trackAssignmentThreshold = 0.0
    optionalPlots = True,
)

process.mtdTracksGNNValid = process.mtdTracksValid.clone(
    folder = cms.string('MTD/Tracks/GNN/'),
    inputTagV = 'offlinePrimaryVerticesGNN',
    t0PID = 'tofPIDGNN:t0',
    t0SafePID = 'tofPIDGNN:t0safe',
    sigmat0SafePID = 'tofPIDGNN:sigmat0safe',
    sigmat0PID = 'tofPIDGNN:sigmat0',
)

# ============================================================================
# Sequences
# ============================================================================
process.validation = cms.Sequence(
    btlValidation +
    etlValidation +
    # 3D
    process.vertices3DValid +
    process.mtdTracks3DValid +
    # 4D
    process.mtdTracksValid +
    process.vertices4DValid +
    # GNN
    process.verticesGNNValid +
    process.mtdTracksGNNValid +
    # Shared
    process.mtdEleIsoValid
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(output_file),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path(process.mix + process.mtdTrackingRecHits + process.validation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(process.p, process.endjob_step, process.DQMoutput_step)

print("\nValidation modules:")
print("  - vertices3DValid    -> MTD/Vertices/3D/  (offlinePrimaryVertices)")
print("  - vertices4DValid    -> MTD/Vertices/4D/  (offlinePrimaryVertices4D)")
print("  - verticesGNNValid   -> MTD/Vertices/GNN/ (offlinePrimaryVerticesGNN)")
print("  - mtdTracks3DValid   -> MTD/Tracks/3D/")
print("  - mtdTracksValid     -> MTD/Tracks/4D/")
print("  - mtdTracksGNNValid  -> MTD/Tracks/GNN/")
