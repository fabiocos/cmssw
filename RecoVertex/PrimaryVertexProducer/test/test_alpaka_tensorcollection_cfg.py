#!/usr/bin/env python3
"""
Test the new proper Alpaka TensorCollection-based GNNVertexProducerAlpaka.

This tests:
- TrackFeatureSource: produces TrackFeaturesDeviceCollection [N, 13]
- GNNVertexProducerAlpaka: uses model_.forward(queue, inputs, outputs)
- Output: GNNOutputDeviceCollection [N, K=200] with Eigen columns

The model outputs batch=N tensors (slot predictions replicated per track)
for proper TensorCollection integration.
"""

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = "INFO"

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(2))

# Input file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/relval/CMSSW_16_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/02098e56-4e94-4d1d-afea-d901f9842456.root'
    )
)

# Standard track building
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.load("Configuration.StandardSequences.MagneticField_cff")

# First step vertices (for initial track collection)
process.firstStepPrimaryVerticesUnsorted = cms.EDProducer("PrimaryVertexProducer",
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    verbose = cms.untracked.bool(False),
    TrackLabel = cms.InputTag("generalTracks"),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Significance = cms.double(4.0),
        maxD0Error = cms.double(1.0),
        maxDzError = cms.double(1.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any'),
        maxNormalizedChi2 = cms.double(10.0),
        maxEta = cms.double(4.0),
    ),
    TkClusParameters = cms.PSet(
        algorithm = cms.string('DA_vect'),
        TkDAClusParameters = cms.PSet(
            coolingFactor = cms.double(0.6),
            zmerge = cms.double(0.01),
            Tmin = cms.double(2.0),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            d0CutOff = cms.double(3.0),
            dzCutOff = cms.double(3.0),
            vertexSize = cms.double(0.006),
            uniquetrkweight = cms.double(0.9),
            uniquetrkminp = cms.double(0.0),
        )
    ),
    vertexCollections = cms.VPSet(
        cms.PSet(
            label = cms.string(''),
            algorithm = cms.string('AdaptiveVertexFitter'),
            chi2cutoff = cms.double(2.5),
            minNdof = cms.double(0.0),
            useBeamConstraint = cms.bool(False),
            maxDistanceToBeam = cms.double(1.0),
        )
    ),
    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag(""),
    useMVACut = cms.bool(False),
)

# Track feature producer
process.trackFeatureProducer = cms.EDProducer("vertexgnn::TrackFeatureSource@alpaka",
    vertices = cms.InputTag("firstStepPrimaryVerticesUnsorted"),
    tracks = cms.InputTag("generalTracks"),
    tofPID = cms.InputTag("tofPID4DnoPID"),
    verbose = cms.untracked.bool(True),
)

# GNN producer using proper TensorCollection pattern
process.gnnVertexProducer = cms.EDProducer("vertexgnn::GNNVertexProducerAlpaka@alpaka",
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    verbose = cms.untracked.bool(True),
)

# Sequence
process.p = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.trackFeatureProducer *
    process.gnnVertexProducer
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(True),
)
