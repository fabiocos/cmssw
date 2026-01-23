"""
ONNX-only GPU Test (single thread)
Extracted from working test_parity_onnx_alpaka_cfg.py
"""

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('ONNXGPU', Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root'
    ),
)

# Logging
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.PrimaryVertexProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizer = cms.untracked.PSet(limit = cms.untracked.int32(-1))

# Single thread for isolation
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# TOF PID (required for timing products)
from RecoVertex.Configuration.RecoVertex_phase2_timing_cff import tofPID4DnoPID
process.tofPID4DnoPID = tofPID4DnoPID.clone()

# TRACK FILTER PARAMETERS (exactly as in working parity config)
TkFilterParams = cms.PSet(
    algorithm = cms.string('filter'),
    maxNormalizedChi2 = cms.double(10.0),
    minPixelLayersWithHits = cms.int32(2),
    minSiliconLayersWithHits = cms.int32(5),
    maxD0Significance = cms.double(4.0),
    maxD0Error = cms.double(1.0),
    maxDzError = cms.double(1.0),
    minPt = cms.double(0.0),
    maxEta = cms.double(4.0),
    trackQuality = cms.string('any'),
    minValidStripHits = cms.int32(0),
)

# ONNX PATH: GNN2D_vect with dummy ONNX model (exactly as parity config)
process.unsortedOfflinePrimaryVerticesONNX = process.unsortedOfflinePrimaryVertices4D.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_vect"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.5),
            trackAssignmentThreshold = cms.double(0.0),
            numSlots = cms.int32(200),
            nnVersion = cms.string("dummy_parity_test"),
            onnxBackend = cms.string("CUDA"),
            onnxModelPath = cms.FileInPath('RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx'),
            Tmin = cms.double(4.0),
            Tpurge = cms.double(4.0),
            Tstop = cms.double(2.0),
            vertexSize = cms.double(0.01),
            d0CutOff = cms.double(3.0),
            verbose = cms.untracked.bool(True),
        )
    ),
    TkFilterParameters = TkFilterParams,
    TrackTimesLabel = cms.InputTag("tofPID4DnoPID:t0safe"),
    TrackTimeResosLabel = cms.InputTag("tofPID4DnoPID:sigmat0safe"),
    verbose = cms.untracked.bool(True),
)

# Inspector
process.gnnInspectorONNX = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesONNX", "", "ONNXGPU"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(10),
    dropNaNs = cms.bool(True),
)

# TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("gnn_onnx_gpu_only.root")
)

# Remove default GNN producer
process.vertexreco.remove(process.unsortedOfflinePrimaryVerticesGNN)

# PATHS
process.common_path = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.vertexreco *
    process.tofPID4DnoPID
)

process.onnx_path = cms.Path(
    process.unsortedOfflinePrimaryVerticesONNX
)

process.inspect_path = cms.EndPath(process.gnnInspectorONNX)

# SCHEDULE
process.schedule = cms.Schedule(
    process.common_path,
    process.onnx_path,
    process.inspect_path
)

print("=" * 70)
print("ONNX GPU-only Test (single thread)")
print("=" * 70)
