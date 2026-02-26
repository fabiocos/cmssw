# GNN Vertex Reconstruction - Alpaka Backend Implementation

This document provides comprehensive documentation of the GNN-based primary vertex reconstruction system, covering both the **ONNX backend** (original) and the **Alpaka backend** (newly developed).

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture Comparison](#architecture-comparison)
3. [ONNX Backend (Original)](#onnx-backend-original)
4. [Alpaka Backend (New)](#alpaka-backend-new)
5. [SoA Data Structures](#soa-data-structures)
6. [Model Specification](#model-specification)
7. [Configuration Guide](#configuration-guide)
8. [Testing & Parity Verification](#testing--parity-verification)
9. [File Reference](#file-reference)
10. [Troubleshooting](#troubleshooting)

---

## Overview

The GNN Vertex Producer uses a **VertexSlotModel** neural network to cluster tracks into primary vertices. The model outputs:
- **Assignment probabilities** A[N, K]: Which slot each track belongs to
- **Slot positions** z_hat[K], t_hat[K]: Predicted z and t for each slot
- **Existence probabilities** p[K]: Whether each slot is a real vertex
- **PID weights** pi[N, 3]: Particle ID weights per track (pion, kaon, proton)

Two inference backends are supported:
- **ONNX Runtime** (`GNN2D_vect`): CPU/CUDA via ONNX Runtime
- **PyTorchAlpaka** (`GNN2D_alpaka`): CPU/CUDA/ROCm via TorchScript + Alpaka

---

## Architecture Comparison

### ONNX Path (Single Module)

```
┌─────────────────────────────────────────────────────────────┐
│                   PrimaryVertexProducer                      │
│                   (algorithm = "GNN2D_vect")                 │
├─────────────────────────────────────────────────────────────┤
│  1. TrackFilterForPVFinding::select()                       │
│  2. GNNClusterizer::vertices()                              │
│     - Extract 13 features per track                         │
│     - ONNX inference: model([N,13]) → A, z, t, p, pi        │
│     - Build TransientVertex from outputs                    │
│  3. AdaptiveVertexFitter                                     │
│  4. Produce reco::VertexCollection                          │
└─────────────────────────────────────────────────────────────┘
```

### Alpaka Path (Multi-Module Pipeline)

```
┌──────────────────────┐     ┌──────────────────────────┐     ┌─────────────────────────┐
│  TrackFeatureProducer │ ──> │  GNNVertexProducerAlpaka │ ──> │  PrimaryVertexProducer  │
│  (Standard EDProducer)│     │  (Alpaka stream::EDProducer)   │  (algorithm = "GNN2D_alpaka")
├──────────────────────┤     ├──────────────────────────┤     ├─────────────────────────┤
│  • Consume tracks     │     │  • TensorCollection input│     │  • Consume SoA outputs  │
│  • Apply track filter │     │  • TorchScript inference │     │  • GNNClusterizerFromAlpaka
│  • Extract 13 features│     │  • TensorCollection output     │  • Build vertices       │
│  • Produce SoA [N,13] │     │  • Produce SoA [N,K]+[N,3]     │  • Produce ValueMaps    │
└──────────────────────┘     └──────────────────────────┘     └─────────────────────────┘
```

---

## ONNX Backend (Original)

### Key Files

| File | Purpose |
|------|---------|
| `src/GNNClusterizer.cc` | Main clustering logic + ONNX inference |
| `interface/GNNClusterizer.h` | Header with ValueMap accessors |

### Data Flow

1. **Track Selection**: `TrackFilterForPVFinding::select()` filters `generalTracks`
2. **Feature Extraction**: 13 features per track extracted from TransientTrack
3. **ONNX Inference**: 
   ```cpp
   input_values.push_back(features_flat);  // [1, N, 13]
   auto outputs = onnxRuntime_->run(input_names, input_values, ...);
   ```
4. **Output Parsing**:
   - `A_flat[N*K]` - assignment probabilities
   - `z_hat_flat[K]` - z positions
   - `t_hat_flat[K]` - times
   - `p_flat[K]` - existence probabilities
   - `pi_flat[N*3]` - PID weights
5. **Vertex Building**: For each slot with `p[k] > threshold`, collect assigned tracks

### 13 Input Features

| Index | Name | Description |
|-------|------|-------------|
| 0 | vz | Track z at closest approach |
| 1 | dz | Track z error |
| 2 | pt | Transverse momentum |
| 3 | eta | Pseudorapidity |
| 4 | mva | MVA quality |
| 5 | pl | Path length |
| 6 | t_pi | Time of flight (pion hypothesis) |
| 7 | t_k | Time of flight (kaon hypothesis) |
| 8 | t_p | Time of flight (proton hypothesis) |
| 9 | s_pi | Time uncertainty (pion) |
| 10 | s_k | Time uncertainty (kaon) |
| 11 | s_p | Time uncertainty (proton) |
| 12 | has_time | 1.0 if MTD timing available, else 0.0 |

---

## Alpaka Backend (New)

### Key Files

| File | Purpose |
|------|---------|
| `interface/VertexGNNSoA.h` | SoA layout definitions |
| `interface/VertexGNNHostCollection.h` | Host collection typedefs |
| `interface/alpaka/VertexGNNDeviceCollection.h` | Device collection typedefs |
| `plugins/TrackFeatureProducer.cc` | Feature extraction → SoA |
| `plugins/alpaka/GNNVertexProducerAlpaka.cc` | TorchScript inference |
| `plugins/GNNVertexBuilderFromAlpaka.cc` | Standalone vertex builder |
| `src/GNNClusterizerFromAlpaka.cc` | Vertex building logic |
| `plugins/PrimaryVertexProducer.cc` | Integrated GNN2D_alpaka path |

### Module 1: TrackFeatureProducer

**Purpose**: Standard EDProducer that extracts features from tracks.

```cpp
produces<TrackFeaturesHostCollection>();
```

**Critical**: Uses same `TrackFilterForPVFinding` as ONNX path to ensure identical track selection:
```cpp
trackFilter_ = std::make_unique<TrackFilterForPVFinding>(
    params.getParameter<edm::ParameterSet>("TkFilterParameters"));
```

### Module 2: GNNVertexProducerAlpaka

**Purpose**: Alpaka stream::EDProducer that runs TorchScript inference on GPU.

**GPU Data Flow**:
```cpp
// 1. Consume HostCollection from standard EDProducer
const auto& hostInput = event.get(trackFeaturesToken_);

// 2. Copy to device for GPU inference
TrackFeaturesDeviceCollection deviceInput(N, event.queue());
alpaka::memcpy(event.queue(), deviceInput.buffer(), hostInput.buffer());

// 3. Allocate output on device
GNNOutputDeviceCollection deviceOutput(N, event.queue());

// 4. Run inference on device (GPU)
model_.forward(event.queue(), inputs, outputs);

// 5. Put DeviceCollection - framework handles automatic D2H for downstream consumers
event.emplace(gnnOutputToken_, std::move(deviceOutput));
```

**TensorCollection Input Setup** (13 features → [N, 13]):
```cpp
inputs.add<::vertexgnn::TrackFeaturesSoA>("features",
    inputRecords.vz(), inputRecords.dz(), inputRecords.pt(), inputRecords.eta(),
    inputRecords.mva(), inputRecords.pl(), inputRecords.t_pi(), inputRecords.t_k(),
    inputRecords.t_p(), inputRecords.s_pi(), inputRecords.s_k(), inputRecords.s_p(),
    inputRecords.has_time());
```

**TensorCollection Output Setup** (5 tensors):
```cpp
outputs.add<::vertexgnn::GNNOutputSoA>("A", outputRecords.A());        // [N, K]
outputs.add<::vertexgnn::GNNOutputSoA>("z_hat", outputRecords.z_hat());  // [N, K]
outputs.add<::vertexgnn::GNNOutputSoA>("t_hat", outputRecords.t_hat());  // [N, K]
outputs.add<::vertexgnn::GNNOutputSoA>("p", outputRecords.p());          // [N, K]
outputs.add<::vertexgnn::GNNOutputSoA>("pi", outputRecords.pi());        // [N, 3]
```

**Order matters!** Must match TorchScript model output tuple order.

### Module 3: PrimaryVertexProducer (GNN2D_alpaka)

**Purpose**: Consumes Alpaka output, builds vertices, produces ValueMaps.

**Reading from SoA**:
```cpp
auto gnnView = gnnOutput.const_view();
for (int i = 0; i < gnn_N; ++i) {
    auto elem = gnnView[i];
    // Eigen columns accessed as: elem.A()[k], elem.pi()[j]
    gnn_pi_0[i] = elem.pi()[0];  // Per-track PID
    gnn_pi_1[i] = elem.pi()[1];
    gnn_pi_2[i] = elem.pi()[2];
}
```

---

## SoA Data Structures

### VertexGNNSoA.h

```cpp
namespace vertexgnn {
    constexpr int kNumSlots = 180;  // K (matches v17p1 model)
    
    using SlotVector = Eigen::Vector<float, kNumSlots>;  // [K]
    using PIDVector = Eigen::Vector<float, 3>;           // [3]
    
    // Input: [N] tracks × 13 features
    GENERATE_SOA_LAYOUT(TrackFeaturesLayout,
        SOA_COLUMN(float, vz), SOA_COLUMN(float, dz), ...);  // 13 columns
    
    // Output: [N] tracks × (K slots + 3 PID)
    GENERATE_SOA_LAYOUT(GNNOutputLayout,
        SOA_EIGEN_COLUMN(SlotVector, A),      // [N, K]
        SOA_EIGEN_COLUMN(SlotVector, z_hat),  // [N, K]
        SOA_EIGEN_COLUMN(SlotVector, t_hat),  // [N, K]
        SOA_EIGEN_COLUMN(SlotVector, p),      // [N, K]
        SOA_EIGEN_COLUMN(PIDVector, pi));     // [N, 3]
}
```

### Why Slot Predictions are Replicated

TensorCollection requires all outputs to have the same batch dimension (N tracks). The TorchScript model replicates global slot predictions [K] to [N, K]:
```python
z_hat_expanded = self.z_hat.unsqueeze(0).expand(N, K)  # [K] → [N, K]
```

When consuming, we only read from track 0 since all tracks have identical values:
```cpp
z_hat[k] = gnnView[0].z_hat()[k];  // Same for all tracks
```

---

## Model Specification

### TorchScript Model (Alpaka)

**Output**: Tuple of 5 tensors
```python
return (A, z_hat_expanded, t_hat_expanded, p_expanded, pi)
# Shapes: ([N,K], [N,K], [N,K], [N,K], [N,3])
```

### ONNX Model (ONNX Path)

**Output**: 5 named tensors
```
A: [N, K], z_hat: [K], t_hat: [K], p: [K], pi: [N, 3]
```

Note: ONNX path uses non-replicated [K] for slot predictions.

### Dummy Model for Parity Testing

`test/create_dummy_model.py` generates both formats with identical weights (seed=42):
- `data/dummy_vertex_slot.pt` - TorchScript
- `data/dummy_vertex_slot.onnx` - ONNX

---

## Configuration Guide

### ONNX Path Configuration

```python
process.producer = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_vect"),
        TkDAClusParameters = cms.PSet(
            onnxModelPath = cms.FileInPath('RecoVertex/PrimaryVertexProducer/data/vertex_slot_model.onnx'),
            onnxBackend = cms.string("CPU"),  # or "CUDA"
            numSlots = cms.int32(180),  # Must match model
            existenceThreshold = cms.double(0.5),
            ...
        )
    ),
    TkFilterParameters = TkFilterParams,  # Track filter
)
```

### Alpaka Path Configuration

```python
# Step 1: Feature extraction
process.trackFeatureProducer = cms.EDProducer("vertexgnn::TrackFeatureProducer",
    tracks = cms.InputTag("generalTracks"),
    TkFilterParameters = TkFilterParams,  # SAME filter as PrimaryVertexProducer
    ...
)

# Step 2: Inference
process.gnnVertexProducer = cms.EDProducer("vertexgnn::GNNVertexProducerAlpaka@alpaka",
    trackFeatures = cms.InputTag("trackFeatureProducer"),
    model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"),
)

# Step 3: Vertex building
process.producer = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        algorithm = cms.string("GNN2D_alpaka"),
        TkDAClusParameters = cms.PSet(
            gnnOutput = cms.InputTag("gnnVertexProducer"),
            existenceThreshold = cms.double(0.5),
        )
    ),
    TkFilterParameters = TkFilterParams,
)
```

### Track Filter Parameters (Phase 2 Production Defaults)

```python
TkFilterParams = cms.PSet(
    algorithm = cms.string('filter'),
    maxNormalizedChi2 = cms.double(10.0),
    minPixelLayersWithHits = cms.int32(2),
    minSiliconLayersWithHits = cms.int32(5),
    maxD0Significance = cms.double(4.0),
    maxD0Error = cms.double(1.0),
    maxDzError = cms.double(1.0),
    minPt = cms.double(0.0),    # Phase2 default
    maxEta = cms.double(4.0),   # Phase2 default
    trackQuality = cms.string('any'),
)
```

---

## Testing & Parity Verification

### Configuration File

`test/test_parity_onnx_alpaka_cfg.py` runs both paths on identical events.

### Running Parity Test

```bash
cd $CMSSW_BASE/src

# Generate dummy models
python3 RecoVertex/PrimaryVertexProducer/test/create_dummy_model.py

# Build
scram b -j8

# Run parity test (50 events)
cmsRun RecoVertex/PrimaryVertexProducer/test/test_parity_onnx_alpaka_cfg.py

# Compare outputs
python3 compare_gnn_inspectors.py
```

### Expected Output

```
slotAssign:  ✓ MATCH (entries identical)
maxProb:     ✓ MATCH (entries identical)
piWeight0:   ✓ MATCH (mean identical)
piWeight1:   ✓ MATCH (mean identical)
piWeight2:   ✓ MATCH (mean identical)
```

---

## File Reference

### Interface Headers

| File | Description |
|------|-------------|
| `interface/VertexGNNSoA.h` | SoA layouts for input/output |
| `interface/VertexGNNHostCollection.h` | Host collection typedefs |
| `interface/alpaka/VertexGNNDeviceCollection.h` | Device collection typedefs |
| `interface/GNNClusterizer.h` | ONNX clusterizer interface |
| `interface/GNNClusterizerFromAlpaka.h` | Alpaka clusterizer interface |

### Plugin Sources

| File | Description |
|------|-------------|
| `plugins/PrimaryVertexProducer.cc` | Main producer (both paths) |
| `plugins/TrackFeatureProducer.cc` | Feature extraction → SoA |
| `plugins/alpaka/GNNVertexProducerAlpaka.cc` | TorchScript inference |
| `plugins/GNNVertexBuilderFromAlpaka.cc` | Standalone builder (unused in integrated path) |

### Implementation Sources

| File | Description |
|------|-------------|
| `src/GNNClusterizer.cc` | ONNX inference + vertex building |
| `src/GNNClusterizerFromAlpaka.cc` | Vertex building from SoA |

### Alpaka Dictionary Files (CUDA SoA Serialization)

| File | Description |
|------|-------------|
| `src/alpaka/classes_cuda.h` | CUDA device dictionary headers |
| `src/alpaka/classes_cuda_def.xml` | CUDA collection wrappers for EDM |
| `src/alpaka/BuildFile.xml` | Build rules for CUDA dictionaries |

### Test Files

| File | Description |
|------|-------------|
| `test/create_dummy_model.py` | Generate TorchScript + ONNX models |
| `test/test_parity_onnx_alpaka_cfg.py` | Parity test configuration |
| `test/test_onnx_gpu_only_cfg.py` | ONNX GPU-only isolated test |
| `test/test_alpaka_gpu_only_cfg.py` | Alpaka GPU-only isolated test |
| `test/compare_gnn_inspectors.py` | Compare histogram outputs |

### Production Configs (src/)

| File | Description |
|------|-------------|
| `vertexTask_cfg.py` | Original ONNX-based GNN vertex task |
| `vertexTask_alpaka_cfg.py` | **Drop-in Alpaka replacement** (GPU)|

---

## Drop-in Replacement: vertexTask_alpaka_cfg.py

The `vertexTask_alpaka_cfg.py` is a **production-ready drop-in replacement** for `vertexTask_cfg.py` that uses GPU inference via the Alpaka backend.

### Key Differences from ONNX vertexTask

| Aspect | vertexTask_cfg.py | vertexTask_alpaka_cfg.py |
|--------|------------------|--------------------------|
| Algorithm | `GNN2D_vect` (ONNX) | `GNN2D_alpaka` (TorchScript) |
| Backend | CPU or CUDA (ONNX Runtime) | CUDA (PyTorchAlpaka) |
| Config Lines | Single module | 3 modules in sequence |
| Model Format | `.onnx` | `.pt` (TorchScript) |

### Path Structure

```python
# Alpaka producers run before vertexreco
process.alpaka_producers = cms.Sequence(
    process.trackFeatureProducer *
    process.gnnVertexProducer
)

process.exe = cms.Path(
    process.firstStepPrimaryVerticesUnsorted *
    process.alpaka_producers *
    process.vertexreco
)
```

### Running

```bash
# GPU inference with Alpaka
cmsRun vertexTask_alpaka_cfg.py

# Compare with ONNX version
cmsRun vertexTask_cfg.py
```

---

## Troubleshooting

### Shape Mismatch Errors

**Symptom**: `output with shape [X, 1] doesn't match broadcast shape [X, X]`

**Cause**: TensorCollection expects Eigen columns for outputs. Scalar columns (`SOA_COLUMN(float, ...)`) expect `[N, 1]` tensors, but model outputs `[N]`.

**Solution**: Use `SOA_EIGEN_COLUMN` for all outputs, even scalars:
```cpp
using PIDVector = Eigen::Vector<float, 3>;
SOA_EIGEN_COLUMN(PIDVector, pi)  // [N, 3]
```

### Track Count Mismatch

**Symptom**: Different track counts between ONNX and Alpaka paths

**Cause**: Filter parameters not synchronized

**Solution**: Ensure both paths use identical `TkFilterParameters`:
```python
process.trackFeatureProducer.TkFilterParameters = TkFilterParams
process.producer.TkFilterParameters = TkFilterParams
```

### ONNX Runtime Initialization for Alpaka Path

**Symptom**: ONNX runtime loading when using Alpaka path

**Cause**: `initializeGlobalCache` must return `nullptr` for GNN2D_alpaka

**Solution**: Check `PrimaryVertexProducer.cc`:
```cpp
if (algo == "GNN2D_alpaka" || algo == "DA2D_vect" || ...) {
    return nullptr;
}
```

---

## Author Notes

- **Parity Status**: 100% verified between ONNX and Alpaka backends
- **Phase 2 Ready**: Uses MTD timing features
- **GPU Support**: Alpaka path supports CUDA and ROCm backends

---

*Last updated: January 2026*
