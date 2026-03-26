"""
mtdHarvesting_3way_cfg.py - MTD Harvesting for 3D, 4D and GNN vertex comparison

Runs post-processors for all three vertex types, reading from the DQM folders
written by mtdValidation_3way_cfg.py.

Usage: cmsRun mtdHarvesting_3way_cfg.py backend=alpaka
"""

import FWCore.ParameterSet.Config as cms
import sys

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting', Phase2C17I13M9)

# ============================================================================
# Parse command line for backend selection
# ============================================================================
backend = "alpaka"  # default
for arg in sys.argv[1:]:
    if arg.startswith("backend="):
        backend = arg.split("=")[1].lower()

if backend not in ["onnx", "alpaka"]:
    raise ValueError(f"Unknown backend: {backend}. Use 'onnx' or 'alpaka'")

input_file = f"file:mtdValidation_3way_{backend}_v17s7.root"

print("=" * 60)
print(f"MTD 3-Way Harvesting (3D + 4D + GNN) - Backend: {backend.upper()}")
print(f"Input: {input_file}")
print("=" * 60)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring(input_file)
)

# Path and EndPath definitions
process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# ============================================================================
# Load base post-processors
# ============================================================================
process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

# ============================================================================
# Vertex Post-Processors (3D, 4D, GNN)
# ============================================================================

# 3D
process.Primary4DVertexPostProcessor3D = process.Primary4DVertexPostProcessor.clone(
    folder = cms.string('MTD/Vertices/3D/')
)
process.MtdTracksPostProcessor3D = process.MtdTracksPostProcessor.clone(
    folder = cms.string('MTD/Tracks/3D/')
)

# 4D (reuse base module, just change folder)
process.Primary4DVertexPostProcessor.folder = cms.string('MTD/Vertices/4D/')
process.MtdTracksPostProcessor.folder = cms.string('MTD/Tracks/4D/')

# GNN
process.Primary4DVertexPostProcessorGNN = process.Primary4DVertexPostProcessor.clone(
    folder = cms.string('MTD/Vertices/GNN/')
)
process.MtdTracksPostProcessorGNN = process.MtdTracksPostProcessor.clone(
    folder = cms.string('MTD/Tracks/GNN/')
)

# ============================================================================
# Harvesting Sequence
# ============================================================================
process.harvesting = cms.Sequence(
    process.btlSimHitsPostProcessor +
    process.btlLocalRecoPostProcessor +
    process.MtdEleIsoPostProcessor +
    # 3D
    process.Primary4DVertexPostProcessor3D +
    process.MtdTracksPostProcessor3D +
    # 4D
    process.Primary4DVertexPostProcessor +
    process.MtdTracksPostProcessor +
    # GNN
    process.Primary4DVertexPostProcessorGNN +
    process.MtdTracksPostProcessorGNN
)

process.p = cms.Path(process.harvesting)

process.schedule = cms.Schedule(process.edmtome_step, process.p, process.dqmsave_step)

print("\nHarvesting post-processors:")
print("  - Primary4DVertexPostProcessor3D  -> MTD/Vertices/3D/")
print("  - Primary4DVertexPostProcessor    -> MTD/Vertices/4D/")
print("  - Primary4DVertexPostProcessorGNN -> MTD/Vertices/GNN/")
print("  - MtdTracksPostProcessor3D        -> MTD/Tracks/3D/")
print("  - MtdTracksPostProcessor          -> MTD/Tracks/4D/")
print("  - MtdTracksPostProcessorGNN       -> MTD/Tracks/GNN/")
