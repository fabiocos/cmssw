"""
vertexTask_alpaka_cfg.py - Alpaka GNN Vertex Reconstruction (Drop-in replacement)

Same structure as vertexTask_cfg.py but uses Alpaka backend for GPU inference.
Replaces unsortedOfflinePrimaryVerticesGNN in vertexreco with Alpaka pipeline.
"""

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('REVTX', Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Alpaka + CUDA for GPU inference
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32, cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0220b466-b1b6-47f6-ad50-1ee4609c6de6.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/06190f6c-a9eb-4ed2-b764-36bcd6ae5147.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0721da40-e1ea-493d-86fa-e505b1c93b53.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/08d81f57-3c82-442f-877f-28be46da3441.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/09696683-5b79-4cfa-8ab2-64b3ad94ab80.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0a9ad4a9-9a7a-43f6-86ca-873dec4a26a6.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/1017a25b-dfa2-4888-aa25-617a42e3deea.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/1049f975-3edd-4edb-a472-720cf974bda7.root',
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/16741a4f-db51-4bb3-ae11-76188d9b8a19.root',
    ),
    secondaryFileNames = cms.untracked.vstring()
)

# Debug logging for GNN modules
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.cerr.DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.PrimaryVertexProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNClusterizerFromAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.SequentialPrimaryVertexFitterAdapter = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.TrackFeatureProducer = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.MessageLogger.cerr.GNNVertexProducerAlpaka = cms.untracked.PSet(limit = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# ============================================================================
# ALPAKA GNN Vertex Producer - GPU Inference
# ============================================================================

# Step 1: TrackFeatureProducer (CPU, produces HostCollection)
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
    verbose = cms.untracked.bool(False),
    # Use same TkFilterParameters as unsortedOfflinePrimaryVerticesGNN
    TkFilterParameters = process.unsortedOfflinePrimaryVerticesGNN.TkFilterParameters.clone(),
)

# Step 2: GNNVertexProducerAlpaka (GPU inference)
process.gnnVertexProducer = cms.EDProducer("vertexgnn::GNNVertexProducerAlpaka@alpaka",
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/test/vertex_slot_v17s7_alpaka.pt"),
    verbose = cms.untracked.bool(False),
)

# Step 3: Replace unsortedOfflinePrimaryVerticesGNN with Alpaka version
process.unsortedOfflinePrimaryVerticesGNN = process.unsortedOfflinePrimaryVerticesGNN.clone(
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_alpaka"),
        TkDAClusParameters = cms.PSet(
            existenceThreshold = cms.double(0.01),  # Low threshold for untrained model debugging
            trackAssignmentThreshold = cms.double(0.5),
            gnnOutput = cms.InputTag("gnnVertexProducer"),
            verbose = cms.untracked.bool(False),
        )
    ),
)

# Use fitter geometric weights instead of GNN assignment probabilities for comparison
# Set to False to use AdaptiveVertexFitter's chi2-based weights
# Set to True (default) to use GNN assignment probabilities
for vc in process.unsortedOfflinePrimaryVerticesGNN.vertexCollections:
    vc.useClusterWeights = cms.untracked.bool(False)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Alpaka GNN Vertex Reconstruction (GPU)'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
outputCommand = process.FEVTDEBUGHLTEventContent.outputCommands
outputCommand.append('drop *_offlinePrimaryVertices_*_RECO')
outputCommand.append('drop *_offlinePrimaryVertices4D_*_RECO')
outputCommand.append('drop *_offlinePrimaryVertices4DWithBS_*_RECO')
outputCommand.append('drop *_offlinePrimaryVerticesWithBS_*_RECO')
outputCommand.append('drop *_TriggerResults_*_RECO')
outputCommand.append('drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModel_RECO')
outputCommand.append('drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModelResolution_RECO')
outputCommand.append('drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModel_RECO')
outputCommand.append('drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModelResolution_RECO')
outputCommand.append('drop *_tofPID_probK_RECO')
outputCommand.append('drop *_tofPID_probP_RECO')
outputCommand.append('drop *_tofPID_probPi_RECO')
outputCommand.append('drop *_tofPID_sigmat0_RECO')
outputCommand.append('drop *_tofPID_sigmat0safe_RECO')
outputCommand.append('drop *_tofPID_t0_RECO')
outputCommand.append('drop *_tofPID_t0safe_RECO')
outputCommand.append('drop *_ak4CaloJetsForTrk_*_RECO')
outputCommand.append('drop *_inclusiveSecondaryVertices_*_RECO')
outputCommand.append('drop *_generalV0Candidates_Kshort_RECO')
outputCommand.append('drop *_generalV0Candidates_Lambda_RECO')

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:revtx_step3_alpaka_v17s7.root'),
    outputCommands = outputCommand,
    splitLevel = cms.untracked.int32(0)
)

# TFileService for GNN track inspector output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("gnn_inspector_v17s7.root")
)

# ============================================================================
# GNN Track Inspector (same as vertexTask)
# ============================================================================
process.gnnInspector = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesGNN", "", "REVTX"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(15),
    dropNaNs = cms.bool(True),
)

# ============================================================================
# Path definitions (same structure as vertexTask)
# ============================================================================
# Add Alpaka producers before vertexreco since we're replacing the GNN in it
process.alpaka_producers = cms.Sequence(
    process.trackFeatureProducer *
    process.gnnVertexProducer
)

process.exe = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.alpaka_producers *
    process.vertexreco
)
process.inspect_path = cms.EndPath(process.gnnInspector)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Performance monitoring
from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)

print("=" * 70)
print("vertexTask_alpaka_cfg.py - Alpaka GNN VertexSlotModel (GPU)")
print("=" * 70)
print("Mode: GPU Inference via Alpaka/PyTorchAlpaka")
print("Parameters:")
print("  - existenceThreshold: 0.01")
print("  - trackAssignmentThreshold: 0.0")
print("  - useClusterWeights: False (using adaptive fitter weights, not GNN A)")
print("  - model: vertex_slot_v17s7_alpaka.pt (v17s7)")
print("  - backend: CUDA (GPU)")
