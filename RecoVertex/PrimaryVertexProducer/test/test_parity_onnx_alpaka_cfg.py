"""
Parity Test: ONNX vs Alpaka GNN Paths

This config runs BOTH GNN paths using the SAME dummy model weights:
1. ONNX path: GNN2D_vect with dummy_vertex_slot.onnx
2. Alpaka path: GNN2D_alpaka with dummy_vertex_slot.pt

Outputs:
- revtx_parity_test.root - vertex collections from both paths
- gnn_parity_test.root   - inspector histograms (subdirs: gnnInspectorONNX, gnnInspectorAlpaka)
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('PARITY', Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Alpaka + CUDA backend for GPU inference
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')

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
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.PrimaryVertexProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizerFromAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.TrackFeatureProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(10),
    numberOfStreams = cms.untracked.uint32(10),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# =============================================================================
# TOF PID (required for timing products)
# =============================================================================
from RecoVertex.Configuration.RecoVertex_phase2_timing_cff import tofPID4DnoPID
process.tofPID4DnoPID = tofPID4DnoPID.clone()

# =============================================================================
# TRACK FILTER PARAMETERS (Phase2 production defaults, same for both paths)
# =============================================================================
TkFilterParams = cms.PSet(
    algorithm = cms.string('filter'),
    maxNormalizedChi2 = cms.double(10.0),
    minPixelLayersWithHits = cms.int32(2),
    minSiliconLayersWithHits = cms.int32(5),
    maxD0Significance = cms.double(4.0),
    maxD0Error = cms.double(1.0),
    maxDzError = cms.double(1.0),
    minPt = cms.double(0.0),   # Phase2 production default
    maxEta = cms.double(4.0),  # Phase2 production default
    trackQuality = cms.string('any'),
    minValidStripHits = cms.int32(0),
)

# =============================================================================
# 1. ONNX PATH: GNN2D_vect with dummy ONNX model
# =============================================================================
process.unsortedOfflinePrimaryVerticesONNX = process.unsortedOfflinePrimaryVertices4D.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_vect"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.5),
            trackAssignmentThreshold = cms.double(0.4),
            numSlots = cms.int32(180),  # Must match model's num_slots
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

# =============================================================================
# 2. ALPAKA PATH: TrackFeatureProducer + GNNVertexProducerAlpaka + PVP
# =============================================================================

# Step 2a: TrackFeatureProducer (uses same filter as ONNX path)
process.trackFeatureProducer = cms.EDProducer("vertexgnn::TrackFeatureProducer",
    tracks = cms.InputTag("generalTracks"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    trackTimesLabel = cms.InputTag("tofPID4DnoPID:t0safe"),
    trackTimeResosLabel = cms.InputTag("tofPID4DnoPID:sigmat0safe"),
    trackMTDTimeQualityVMapTag = cms.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    trackAssocSrc = cms.InputTag("trackExtenderWithMTD:generalTrackassoc"),
    tmtdSrc = cms.InputTag("trackExtenderWithMTD:generalTracktmtd"),
    sigmatmtdSrc = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd"),
    pathmtd = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength"),
    btlMatchChi2Src = cms.InputTag("trackExtenderWithMTD", "btlMatchChi2"),
    btlMatchTimeChi2Src = cms.InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"),
    etlMatchChi2Src = cms.InputTag("trackExtenderWithMTD", "etlMatchChi2"),
    etlMatchTimeChi2Src = cms.InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"),
    tofPi = cms.InputTag("trackExtenderWithMTD:generalTrackTofPi"),
    tofK = cms.InputTag("trackExtenderWithMTD:generalTrackTofK"),
    tofP = cms.InputTag("trackExtenderWithMTD:generalTrackTofP"),
    sigmatofpiSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"),
    sigmatofkSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"),
    sigmatofpSrc = cms.InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"),
    npixBarrelSrc = cms.InputTag("trackExtenderWithMTD", "npixBarrel"),
    npixEndcapSrc = cms.InputTag("trackExtenderWithMTD", "npixEndcap"),
    minTrackTimeQuality = cms.double(0.8),
    useMVACut = cms.bool(False),
    useTrackFilter = cms.bool(True),
    verbose = cms.untracked.bool(True),
    TkFilterParameters = TkFilterParams,
)

# Step 2b: GNNVertexProducerAlpaka (runs inference)
process.gnnVertexProducer = cms.EDProducer("vertexgnn::GNNVertexProducerAlpaka@alpaka",
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    verbose = cms.untracked.bool(True),
)

# Step 2c: PrimaryVertexProducer with GNN2D_alpaka
process.unsortedOfflinePrimaryVerticesAlpaka = process.unsortedOfflinePrimaryVertices4D.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_alpaka"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.5),
            trackAssignmentThreshold = cms.double(0.4),
            gnnOutput = cms.InputTag("gnnVertexProducer"),
            verbose = cms.untracked.bool(True),
        )
    ),
    TkFilterParameters = TkFilterParams,
    TrackTimesLabel = cms.InputTag("tofPID4DnoPID:t0safe"),
    TrackTimeResosLabel = cms.InputTag("tofPID4DnoPID:sigmat0safe"),
    verbose = cms.untracked.bool(True),
)

# =============================================================================
# GNN TRACK INSPECTORS (one per path)
# =============================================================================
process.gnnInspectorONNX = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesONNX", "", "PARITY"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(10),
    dropNaNs = cms.bool(True),
)

process.gnnInspectorAlpaka = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesAlpaka", "", "PARITY"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(10),
    dropNaNs = cms.bool(True),
)

# =============================================================================
# TFileService - separate for each inspector
# =============================================================================
# Note: TFileService can only write one file, so we'll use a combined file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("gnn_parity_test.root")
)

# =============================================================================
# OUTPUT
# =============================================================================
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('revtx_parity_test.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_unsortedOfflinePrimaryVerticesONNX_*_*',
        'keep *_unsortedOfflinePrimaryVerticesAlpaka_*_*',
        'keep *_gnnVertexProducer_*_*',
    )
)

# =============================================================================
# PATHS - Use vertexreco sequence to provide dependencies
# =============================================================================
# vertexreco provides unsortedOfflinePrimaryVertices which tofPID4DnoPID needs
# BUT it also includes unsortedOfflinePrimaryVerticesGNN which we don't want
# So we remove it to only test ONNX vs Alpaka

# Remove the default GNN producer from vertexreco to avoid running 3 producers
process.vertexreco.remove(process.unsortedOfflinePrimaryVerticesGNN)

# Common base path with vertex reconstruction
process.common_path = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.vertexreco *
    process.tofPID4DnoPID
)

# ONNX path - runs after common
process.onnx_path = cms.Path(
    process.unsortedOfflinePrimaryVerticesONNX
)

# Alpaka path - runs after common
process.alpaka_path = cms.Path(
    process.trackFeatureProducer *
    process.gnnVertexProducer *
    process.unsortedOfflinePrimaryVerticesAlpaka
)

# Inspectors
process.inspect_path = cms.EndPath(
    process.gnnInspectorONNX *
    process.gnnInspectorAlpaka
)

process.output_step = cms.EndPath(process.output)

# =============================================================================
# SCHEDULE
# =============================================================================
process.schedule = cms.Schedule(
    process.common_path,
    process.onnx_path,
    process.alpaka_path,
    process.inspect_path,
    process.output_step
)

print("=" * 70)
print("GNN Parity Test: ONNX vs Alpaka")
print("=" * 70)
print("Both paths use IDENTICAL model weights (seed=42)")
print("Outputs:")
print("  - revtx_parity_test.root: vertex collections")
print("  - gnn_parity_test.root:   inspector histograms")
print("=" * 70)
