#!/usr/bin/env python3
"""
Compare ONNX vs TorchScript GNN inference in a producer-like context.
Simulates what GNNClusterizer does: prepare features, run inference, extract vertices.

Usage:
    source /afs/cern.ch/work/p/prsolank/public/fastgraphenv/bin/activate
    python3 RecoVertex/PrimaryVertexProducer/test/compare_producers.py
"""

import numpy as np
import torch
import onnxruntime as ort
from pathlib import Path
import argparse


def generate_track_features(num_tracks: int = 100, seed: int = 42) -> np.ndarray:
    """
    Generate dummy track features similar to what GNNClusterizer prepares.
    
    13 features per track:
    - vz: track z vertex
    - dz: dz error  
    - pt: transverse momentum
    - eta: pseudorapidity
    - mva: track quality MVA
    - pl: significance (dz/dz_error)
    - t_pi, t_k, t_p: timing PID hypotheses
    - s_pi, s_k, s_p: timing uncertainties
    - has_time: MTD flag
    """
    np.random.seed(seed)
    
    # Generate realistic-ish track features
    vz = np.random.normal(0, 5, num_tracks)  # z vertex in cm
    dz = np.abs(np.random.exponential(0.1, num_tracks))  # dz error
    pt = np.abs(np.random.exponential(2, num_tracks) + 0.5)  # pt in GeV
    eta = np.random.uniform(-2.5, 2.5, num_tracks)  # pseudorapidity
    mva = np.random.uniform(0.5, 1.0, num_tracks)  # quality MVA
    pl = np.abs(vz / (dz + 1e-6))  # significance
    
    # Timing features - some tracks have MTD info
    has_time = (np.random.random(num_tracks) > 0.3).astype(np.float32)
    t_pi = np.where(has_time > 0.5, np.random.normal(0, 0.05, num_tracks), 0)
    t_k = np.where(has_time > 0.5, np.random.normal(0, 0.05, num_tracks), 0)
    t_p = np.where(has_time > 0.5, np.random.normal(0, 0.05, num_tracks), 0)
    s_pi = np.where(has_time > 0.5, np.abs(np.random.exponential(0.03, num_tracks)), 0.2)
    s_k = np.where(has_time > 0.5, np.abs(np.random.exponential(0.03, num_tracks)), 0.2)
    s_p = np.where(has_time > 0.5, np.abs(np.random.exponential(0.03, num_tracks)), 0.2)
    
    # Stack features: [N, 13]
    features = np.stack([
        vz, dz, pt, eta, mva, pl,
        t_pi, t_k, t_p, s_pi, s_k, s_p, has_time
    ], axis=1).astype(np.float32)
    
    return features


def run_onnx_inference(model_path: str, features: np.ndarray):
    """Run ONNX inference like GNNClusterizer does."""
    sess = ort.InferenceSession(model_path)
    
    # Input shape: [1, N, 13]
    x = features[np.newaxis, :, :]
    
    outputs = sess.run(None, {"x": x})
    
    # Parse outputs
    A = outputs[0]      # [N, K] assignment
    z_hat = outputs[1]  # [K] z positions
    t_hat = outputs[2]  # [K] t positions
    p = outputs[3]      # [K] existence prob
    pi = outputs[4]     # [K, 4] PID weights
    
    return {
        "A": A,
        "z_hat": z_hat,
        "t_hat": t_hat,
        "p": p,
        "pi": pi
    }


def run_torchscript_inference(model_path: str, features: np.ndarray):
    """Run TorchScript inference like GNNVertexProducerAlpaka would."""
    model = torch.jit.load(model_path)
    model.eval()
    
    # Input shape: [1, N, 13]
    x = torch.from_numpy(features[np.newaxis, :, :])
    
    with torch.no_grad():
        outputs = model(x)
    
    # Parse tuple outputs
    A, z_hat, t_hat, p, pi = outputs
    
    return {
        "A": A.numpy(),
        "z_hat": z_hat.numpy(),
        "t_hat": t_hat.numpy(),
        "p": p.numpy(),
        "pi": pi.numpy()
    }


def extract_vertices(outputs: dict, existence_threshold: float = 0.5):
    """
    Extract vertices from model outputs - simulating what the producers do.
    
    Returns list of vertex dictionaries with:
    - z: vertex z position
    - t: vertex t position
    - tracks: list of (track_idx, weight) assigned to this vertex
    - pid_weights: [4] array of PID weights
    """
    A = outputs["A"]  # [N, K]
    z_hat = outputs["z_hat"]  # [K]
    t_hat = outputs["t_hat"]  # [K]
    p = outputs["p"]  # [K]
    pi = outputs["pi"]  # [K, 4]
    
    N, K = A.shape
    vertices = []
    
    for k in range(K):
        if p[k] < existence_threshold:
            continue
        
        # Find tracks assigned to this slot
        track_to_slot = A.argmax(axis=1)  # [N]
        assigned_mask = track_to_slot == k
        
        if not assigned_mask.any():
            continue
        
        # Get track indices and weights
        tracks = []
        for i in range(N):
            if assigned_mask[i]:
                weight = A[i, k]
                tracks.append((i, float(weight)))
        
        vertices.append({
            "slot": k,
            "z": float(z_hat[k]),
            "t": float(t_hat[k]),
            "existence": float(p[k]),
            "tracks": tracks,
            "pid_weights": pi[k].tolist() if hasattr(pi[k], 'tolist') else list(pi[k])
        })
    
    return vertices


def compare_vertices(onnx_vtx: list, ts_vtx: list, tolerance: float = 1e-5):
    """Compare vertex lists from ONNX and TorchScript."""
    print(f"\n{'='*60}")
    print("VERTEX COMPARISON")
    print(f"{'='*60}")
    
    print(f"ONNX vertices: {len(onnx_vtx)}")
    print(f"TorchScript vertices: {len(ts_vtx)}")
    
    if len(onnx_vtx) != len(ts_vtx):
        print("⚠️  Different number of vertices!")
        # Still compare what we can
    
    # Match by slot
    onnx_by_slot = {v["slot"]: v for v in onnx_vtx}
    ts_by_slot = {v["slot"]: v for v in ts_vtx}
    
    all_slots = set(onnx_by_slot.keys()) | set(ts_by_slot.keys())
    
    all_match = True
    for slot in sorted(all_slots):
        if slot not in onnx_by_slot:
            print(f"  Slot {slot}: Only in TorchScript")
            all_match = False
            continue
        if slot not in ts_by_slot:
            print(f"  Slot {slot}: Only in ONNX")
            all_match = False
            continue
        
        ov = onnx_by_slot[slot]
        tv = ts_by_slot[slot]
        
        z_diff = abs(ov["z"] - tv["z"])
        t_diff = abs(ov["t"] - tv["t"])
        p_diff = abs(ov["existence"] - tv["existence"])
        
        match = z_diff < tolerance and t_diff < tolerance and p_diff < tolerance
        status = "✓" if match else "✗"
        
        print(f"  Slot {slot}: z_diff={z_diff:.2e}, t_diff={t_diff:.2e}, p_diff={p_diff:.2e} {status}")
        print(f"    Tracks: ONNX={len(ov['tracks'])}, TS={len(tv['tracks'])}")
        
        if not match:
            all_match = False
    
    return all_match


def main():
    parser = argparse.ArgumentParser(description="Compare ONNX vs TorchScript producer outputs")
    parser.add_argument("--onnx-model", type=str, 
                        default="RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx")
    parser.add_argument("--ts-model", type=str,
                        default="RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt")
    parser.add_argument("--num-tracks", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.5)
    args = parser.parse_args()
    
    print(f"{'='*60}")
    print("GNN VERTEX PRODUCER COMPARISON: ONNX vs TorchScript")
    print(f"{'='*60}")
    print(f"Num tracks: {args.num_tracks}")
    print(f"Seed: {args.seed}")
    print(f"Existence threshold: {args.threshold}")
    
    # Generate features
    print(f"\n1. Generating {args.num_tracks} track features...")
    features = generate_track_features(args.num_tracks, args.seed)
    print(f"   Feature shape: {features.shape}")
    
    # ONNX inference
    print(f"\n2. Running ONNX inference...")
    onnx_outputs = run_onnx_inference(args.onnx_model, features)
    print(f"   A shape: {onnx_outputs['A'].shape}")
    print(f"   z_hat shape: {onnx_outputs['z_hat'].shape}")
    
    # TorchScript inference
    print(f"\n3. Running TorchScript inference...")
    ts_outputs = run_torchscript_inference(args.ts_model, features)
    print(f"   A shape: {ts_outputs['A'].shape}")
    print(f"   z_hat shape: {ts_outputs['z_hat'].shape}")
    
    # Compare raw outputs
    print(f"\n{'='*60}")
    print("RAW OUTPUT COMPARISON")
    print(f"{'='*60}")
    
    for name in ["A", "z_hat", "t_hat", "p", "pi"]:
        diff = np.abs(onnx_outputs[name] - ts_outputs[name]).max()
        status = "✓" if diff < 1e-5 else "✗"
        print(f"  {name:8s}: max_diff = {diff:.2e} {status}")
    
    # Extract vertices
    print(f"\n4. Extracting vertices (threshold={args.threshold})...")
    onnx_vertices = extract_vertices(onnx_outputs, args.threshold)
    ts_vertices = extract_vertices(ts_outputs, args.threshold)
    
    # Compare vertices
    all_match = compare_vertices(onnx_vertices, ts_vertices)
    
    print(f"\n{'='*60}")
    if all_match:
        print("✅ All vertices match between ONNX and TorchScript!")
    else:
        print("⚠️  Some differences detected (may be numerical precision)")
    print(f"{'='*60}")
    
    return 0 if all_match else 1


if __name__ == "__main__":
    exit(main())
