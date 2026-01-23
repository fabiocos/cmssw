"""
Production Test Configuration for Alpaka GNN Vertex Pipeline with Real Tracks

Pipeline:
  generalTracks + MTD timing 
            ↓
  TrackFeatureProducer (standard) → TrackFeaturesHostCollection
            ↓
  GNNVertexProducerAlpaka (Alpaka) → SlotPredictions + Assignments
            ↓
  GNNVertexBuilderFromAlpaka (standard) → reco::VertexCollection

Usage:
  cmsRun RecoVertex/PrimaryVertexProducer/test/test_production_alpaka_cfg.py
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('AlpakaGNNProduction', Phase2C17I13M9)

# Standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Accelerators_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root'
    ),
)

# Enable logging
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.TrackFeatureProducer = cms.untracked.PSet()
process.MessageLogger.GNNVertexProducerAlpaka = cms.untracked.PSet()
process.MessageLogger.GNNVertexBuilderFromAlpaka = cms.untracked.PSet()
process.MessageLogger.GNNClusterizerFromAlpaka = cms.untracked.PSet()
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# =============================================================================
# LOAD TOF PID PRODUCERS (required for timing products)
# =============================================================================
from RecoVertex.Configuration.RecoVertex_phase2_timing_cff import tofPID4DnoPID
process.tofPID4DnoPID = tofPID4DnoPID.clone()

# =============================================================================
# STEP 1: TrackFeatureProducer (Standard EDProducer)
# Consumes generalTracks + MTD timing → TrackFeaturesHostCollection
# =============================================================================
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
    verbose = cms.untracked.bool(True)
)

# =============================================================================
# STEP 2: GNNVertexProducerAlpaka (Alpaka Producer)
# Consumes TrackFeaturesHostCollection → SlotPredictions + Assignments
# NOTE: Currently consumes DeviceCollection, but HostCollection should 
# auto-transfer via framework. If not, we need a Host→Device copier.
# =============================================================================
process.gnnVertexProducer = cms.EDProducer("alpaka_serial_sync::vertexgnn::GNNVertexProducerAlpaka",
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    numSlots = cms.int32(200),
    verbose = cms.untracked.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("")
    )
)

# =============================================================================
# STEP 3: GNNVertexBuilderFromAlpaka (Standard EDProducer)
# Consumes SlotPredictions + Assignments → reco::VertexCollection
# =============================================================================
process.gnnVertexBuilderAlpaka = cms.EDProducer("GNNVertexBuilderFromAlpaka",
    tracks = cms.InputTag("generalTracks"),
    slotPredictions = cms.InputTag("gnnVertexProducer"),
    assignments = cms.InputTag("gnnVertexProducer"),
    numSlots = cms.int32(200),
    existenceThreshold = cms.double(0.5),
    trackAssignmentThreshold = cms.double(0.0),
    verbose = cms.untracked.bool(True)
)

# =============================================================================
# PATH - need vertexreco first for unsortedOfflinePrimaryVertices
# =============================================================================
process.alpaka_gnn_path = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.vertexreco *
    process.tofPID4DnoPID +
    process.trackFeatureProducer +
    process.gnnVertexProducer +
    process.gnnVertexBuilderAlpaka
)

# Output
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('production_alpaka_gnn_vertices.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_gnnVertexBuilderAlpaka_*_*',
        'keep *_trackFeatureProducer_*_*',
    )
)
process.output_step = cms.EndPath(process.output)

print("=" * 70)
print("Production Alpaka GNN Full Pipeline Test")
print("=" * 70)
print("Pipeline:")
print("  TrackFeatureProducer (real tracks) → GNNVertexProducerAlpaka → Vertices")
print("=" * 70)
