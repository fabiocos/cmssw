#!/usr/bin/env python3
"""
Create a dummy model with the same I/O signature as VertexSlotModel.
Exports to both ONNX and TorchScript formats for backend comparison testing.

Input:  x [1, N, 13] - batch of N tracks with 13 features each
Output: A [N, K] - assignment matrix
        z_hat [K] - vertex z positions
        t_hat [K] - vertex t positions
        p [K] - existence probabilities
        pi [K, 4] - PID weights
"""

import torch
import torch.nn as nn
import numpy as np
import argparse
from pathlib import Path


class DummyVertexSlotModel(nn.Module):
    """
    Dummy model with same I/O as VertexSlotModel.
    Uses simple linear layers - not a real vertex model, just for testing infrastructure.
    """
    
    def __init__(self, num_features: int = 13, num_slots: int = 200, num_pid: int = 4):
        super().__init__()
        self.num_features = num_features
        self.num_slots = num_slots
        self.num_pid = num_pid
        
        hidden = 64
        
        # Track encoder
        self.track_encoder = nn.Sequential(
            nn.Linear(num_features, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU()
        )
        
        # Assignment head: [N, hidden] -> [N, K]
        self.assignment_head = nn.Linear(hidden, num_slots)
        
        # Slot predictions (global pooling then per-slot)
        self.slot_z = nn.Linear(hidden, num_slots)
        self.slot_t = nn.Linear(hidden, num_slots)
        self.slot_p = nn.Linear(hidden, num_slots)
        self.slot_pi = nn.Linear(hidden, num_slots * num_pid)
        
    def forward(self, x: torch.Tensor):
        """
        Args:
            x: [1, N, 13] batch of track features
            
        Returns:
            A: [N, K] assignment probabilities (softmax)
            z_hat: [K] predicted z positions
            t_hat: [K] predicted t positions
            p: [K] existence probabilities (sigmoid)
            pi: [K, 4] PID weights (softmax per slot)
        """
        # Remove batch dim: [1, N, 13] -> [N, 13]
        x = x.squeeze(0)
        N = x.shape[0]
        
        # Encode tracks: [N, 13] -> [N, hidden]
        h = self.track_encoder(x)
        
        # Assignment: [N, K] with softmax over slots
        A = torch.softmax(self.assignment_head(h), dim=-1)
        
        # Global pooled features for slot predictions
        h_global = h.mean(dim=0, keepdim=True)  # [1, hidden]
        
        # Slot predictions
        z_hat = self.slot_z(h_global).squeeze(0)  # [K]
        t_hat = self.slot_t(h_global).squeeze(0)  # [K]
        p = torch.sigmoid(self.slot_p(h_global).squeeze(0))  # [K]
        
        # PID weights: [K, 4] with softmax over PID classes
        pi_raw = self.slot_pi(h_global).view(self.num_slots, self.num_pid)  # [K, 4]
        pi = torch.softmax(pi_raw, dim=-1)
        
        return A, z_hat, t_hat, p, pi


def create_dummy_input(num_tracks: int = 50, num_features: int = 13):
    """Create dummy input tensor for testing."""
    return torch.randn(1, num_tracks, num_features)


def export_onnx(model: nn.Module, output_path: Path, num_tracks: int = 50):
    """Export model to ONNX format."""
    model.eval()
    dummy_input = create_dummy_input(num_tracks, model.num_features)
    
    torch.onnx.export(
        model,
        dummy_input,
        str(output_path),
        input_names=["x"],
        output_names=["A", "z_hat", "t_hat", "p", "pi"],
        dynamic_axes={
            "x": {1: "num_tracks"},
            "A": {0: "num_tracks"}
        },
        opset_version=17,
        do_constant_folding=True
    )
    print(f"Exported ONNX model to: {output_path}")


def export_torchscript(model: nn.Module, output_path: Path, num_tracks: int = 50):
    """Export model to TorchScript format."""
    model.eval()
    
    # Use torch.jit.script for full Python support
    try:
        scripted = torch.jit.script(model)
    except Exception as e:
        print(f"torch.jit.script failed, falling back to trace: {e}")
        dummy_input = create_dummy_input(num_tracks, model.num_features)
        scripted = torch.jit.trace(model, dummy_input)
    
    scripted.save(str(output_path))
    print(f"Exported TorchScript model to: {output_path}")


def test_inference(model: nn.Module, num_tracks: int = 50):
    """Run inference and print output shapes."""
    model.eval()
    x = create_dummy_input(num_tracks, model.num_features)
    
    with torch.no_grad():
        A, z_hat, t_hat, p, pi = model(x)
    
    print(f"\nInput shape: {x.shape}")
    print(f"Output shapes:")
    print(f"  A (assignment):  {A.shape}")
    print(f"  z_hat (z pos):   {z_hat.shape}")
    print(f"  t_hat (t pos):   {t_hat.shape}")
    print(f"  p (existence):   {p.shape}")
    print(f"  pi (PID):        {pi.shape}")
    
    return A, z_hat, t_hat, p, pi


def compare_onnx_torchscript(onnx_path: Path, ts_path: Path, num_tracks: int = 50):
    """Compare outputs from ONNX and TorchScript models."""
    import onnxruntime as ort
    
    # Create identical input
    np.random.seed(42)
    x_np = np.random.randn(1, num_tracks, 13).astype(np.float32)
    x_torch = torch.from_numpy(x_np)
    
    # ONNX inference
    sess = ort.InferenceSession(str(onnx_path))
    onnx_outputs = sess.run(None, {"x": x_np})
    A_onnx, z_onnx, t_onnx, p_onnx, pi_onnx = onnx_outputs
    
    # TorchScript inference
    ts_model = torch.jit.load(str(ts_path))
    ts_model.eval()
    with torch.no_grad():
        A_ts, z_ts, t_ts, p_ts, pi_ts = ts_model(x_torch)
    
    # Compare
    print("\n=== ONNX vs TorchScript Comparison ===")
    outputs = [
        ("A", A_onnx, A_ts.numpy()),
        ("z_hat", z_onnx, z_ts.numpy()),
        ("t_hat", t_onnx, t_ts.numpy()),
        ("p", p_onnx, p_ts.numpy()),
        ("pi", pi_onnx, pi_ts.numpy())
    ]
    
    all_match = True
    for name, onnx_out, ts_out in outputs:
        max_diff = np.abs(onnx_out - ts_out).max()
        match = max_diff < 1e-5
        status = "✓" if match else "✗"
        print(f"  {name:8s}: max_diff = {max_diff:.2e} {status}")
        if not match:
            all_match = False
    
    if all_match:
        print("\n✅ All outputs match between ONNX and TorchScript!")
    else:
        print("\n⚠️  Some outputs differ (may be due to numerical precision)")
    
    return all_match


def main():
    parser = argparse.ArgumentParser(description="Create and test dummy VertexSlotModel")
    parser.add_argument("--output-dir", type=Path, default=Path("RecoVertex/PrimaryVertexProducer/data"),
                        help="Output directory for model files")
    parser.add_argument("--num-tracks", type=int, default=50, help="Number of tracks for testing")
    parser.add_argument("--num-slots", type=int, default=200, help="Number of slots (K)")
    parser.add_argument("--compare", action="store_true", help="Compare ONNX and TorchScript outputs")
    args = parser.parse_args()
    
    # Create model
    model = DummyVertexSlotModel(num_features=13, num_slots=args.num_slots)
    print(f"Created DummyVertexSlotModel with {args.num_slots} slots")
    
    # Test inference
    test_inference(model, args.num_tracks)
    
    # Export both formats
    args.output_dir.mkdir(parents=True, exist_ok=True)
    onnx_path = args.output_dir / "dummy_vertex_slot.onnx"
    ts_path = args.output_dir / "dummy_vertex_slot.pt"
    
    export_onnx(model, onnx_path, args.num_tracks)
    export_torchscript(model, ts_path, args.num_tracks)
    
    # Compare if requested
    if args.compare:
        compare_onnx_torchscript(onnx_path, ts_path, args.num_tracks)
    
    print(f"\nModel files saved to {args.output_dir}/")
    print(f"  ONNX:        {onnx_path.name}")
    print(f"  TorchScript: {ts_path.name}")


if __name__ == "__main__":
    main()
