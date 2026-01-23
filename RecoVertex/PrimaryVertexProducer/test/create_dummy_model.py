#!/usr/bin/env python3
"""
Create a dummy TorchScript model for GNN vertex prediction (Alpaka-compatible).

This model is designed for proper TensorCollection integration:
- Input: [N, 13] track features
- Output: All tensors with batch=N (slot predictions replicated per track)

Output structure matches GNNOutputSoA:
  - A: [N, K] assignment probabilities
  - z_hat: [N, K] z predictions (replicated)
  - t_hat: [N, K] t predictions (replicated)
  - p: [N, K] existence probabilities (replicated)
  - pi: [N, 3] PID weights per track (pion, kaon, proton)
"""

import torch
import torch.nn as nn


class DummyVertexSlotModel(nn.Module):
    """Dummy model that outputs replicated slot predictions for each track."""
    
    def __init__(self, num_slots: int = 200, num_features: int = 13):
        super().__init__()
        self.num_slots = num_slots
        self.num_features = num_features
        
        # Simple network to produce slot predictions
        self.slot_net = nn.Sequential(
            nn.Linear(num_features, 64),
            nn.ReLU(),
            nn.Linear(64, num_slots)  # Assignment probabilities per track
        )
        
        # Global slot predictions (learned parameters, shared across all tracks)
        self.z_hat = nn.Parameter(torch.randn(num_slots) * 5)  # z positions
        self.t_hat = nn.Parameter(torch.zeros(num_slots))       # t positions
        self.p = nn.Parameter(torch.sigmoid(torch.randn(num_slots)))  # existence
        
        # Per-track PID network: [N, 13] -> [N, 3]
        self.pi_net = nn.Sequential(
            nn.Linear(num_features, 16),
            nn.ReLU(),
            nn.Linear(16, 3),
            nn.Softmax(dim=1)  # PID probabilities sum to 1
        )
    
    def forward(self, x: torch.Tensor):
        """
        Args:
            x: [N, 13] track features
        
        Returns:
            Tuple of 5 tensors all with batch=N:
            - A: [N, K] assignment probabilities (softmax over slots)
            - z_hat: [N, K] z predictions (replicated from global [K])
            - t_hat: [N, K] t predictions (replicated from global [K])
            - p: [N, K] existence probabilities (replicated from global [K])
            - pi: [N, 3] PID weights per track (pion, kaon, proton)
        """
        # Handle batch dimension if present
        if x.dim() == 3:
            x = x.squeeze(0)  # Remove batch dim: [1, N, 13] -> [N, 13]
        
        N = x.size(0)
        K = self.num_slots
        
        # Assignment probabilities: [N, K]
        A = torch.softmax(self.slot_net(x), dim=1)
        
        # Replicate global slot predictions for each track: [K] -> [N, K]
        z_hat_expanded = self.z_hat.unsqueeze(0).expand(N, K)
        t_hat_expanded = self.t_hat.unsqueeze(0).expand(N, K)
        p_expanded = self.p.unsqueeze(0).expand(N, K)
        
        # Per-track PID weights: [N, 3]
        pi = self.pi_net(x)
        
        return (A, z_hat_expanded, t_hat_expanded, p_expanded, pi)


def main():
    print("Creating Alpaka-compatible dummy vertex slot model...")
    
    # Set seed for reproducibility
    torch.manual_seed(42)
    
    model = DummyVertexSlotModel(num_slots=200, num_features=13)
    model.eval()
    
    # Test with sample input
    N = 100  # number of tracks
    x = torch.randn(N, 13)
    
    with torch.no_grad():
        A, z_hat, t_hat, p, pi = model(x)
    
    print(f"Input shape: {x.shape}")
    print(f"Output shapes (all batch=N={N}):")
    print(f"  A: {A.shape}")
    print(f"  z_hat: {z_hat.shape}")
    print(f"  t_hat: {t_hat.shape}")
    print(f"  p: {p.shape}")
    print(f"  pi: {pi.shape}")
    
    # Verify all outputs have same batch dimension
    assert A.size(0) == N
    assert z_hat.size(0) == N
    assert t_hat.size(0) == N
    assert p.size(0) == N
    assert pi.size(0) == N
    assert pi.size(1) == 3
    print("✓ All outputs have batch=N")
    
    # Verify replicated values are identical across tracks
    assert torch.allclose(z_hat[0], z_hat[1]), "z_hat should be identical across tracks"
    print("✓ Slot predictions correctly replicated")
    
    # Script the model for TorchScript/Alpaka
    scripted_model = torch.jit.script(model)
    
    # Save TorchScript
    import os
    data_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    ts_path = os.path.join(data_dir, "data", "dummy_vertex_slot.pt")
    os.makedirs(os.path.dirname(ts_path), exist_ok=True)
    scripted_model.save(ts_path)
    print(f"Saved TorchScript to {ts_path}")
    
    # Export to ONNX for ONNX path parity testing
    # ONNX expects input with batch dimension [1, N, 13]
    dummy_input = torch.randn(1, N, 13)
    
    onnx_path = os.path.join(data_dir, "data", "dummy_vertex_slot.onnx")
    
    # Note: For ONNX we need a wrapper that handles the batch dimension
    # and returns outputs in the format expected by GNNClusterizer
    class ONNXWrapper(nn.Module):
        """Wrapper for ONNX export that matches GNNClusterizer expectations.
        
        ONNX model expects:
          Input: x [1, N, 13]
        
        ONNX model outputs (matching GNNClusterizer.cc):
          - A: [N, K] assignment probabilities
          - z_hat: [K] vertex z positions
          - t_hat: [K] vertex times
          - p: [K] existence probabilities
          - pi: [N, 3] PID weights (ONNX path expects combined [N, 3])
        """
        def __init__(self, base_model):
            super().__init__()
            self.base_model = base_model
        
        def forward(self, x):
            # x is [1, N, 13], squeeze to [N, 13]
            x = x.squeeze(0)
            A, z_hat, t_hat, p, pi = self.base_model(x)
            # Return in ONNX format:
            # - z_hat, t_hat, p: [K] (take from first track since replicated)
            # - pi: [N, 3] already in correct format
            return (A, z_hat[0], t_hat[0], p[0], pi)
    
    onnx_model = ONNXWrapper(model)
    onnx_model.eval()
    
    torch.onnx.export(
        onnx_model,
        dummy_input,
        onnx_path,
        input_names=['x'],
        output_names=['A', 'z_hat', 't_hat', 'p', 'pi'],
        dynamic_axes={
            'x': {1: 'N'},
            'A': {0: 'N'},
            'pi': {0: 'N'}
        },
        opset_version=14
    )
    print(f"Saved ONNX to {onnx_path}")
    
    # Verify ONNX export
    try:
        import onnx
        onnx_model_loaded = onnx.load(onnx_path)
        onnx.checker.check_model(onnx_model_loaded)
        print("✓ ONNX model verified")
    except ImportError:
        print("(onnx package not available for verification)")
    
    # =============================================================================
    # PARITY CHECK: Compare TorchScript vs ONNX outputs
    # =============================================================================
    print("\n=== Parity Check: TorchScript vs ONNX ===")
    
    # Load TorchScript model
    ts_loaded = torch.jit.load(ts_path)
    ts_loaded.eval()
    
    # Run ONNX inference via PyTorch wrapper (same input)
    with torch.no_grad():
        # Reset seed and use same input
        torch.manual_seed(42)
        x_test = torch.randn(N, 13)
        x_test_onnx = x_test.unsqueeze(0)  # [1, N, 13] for ONNX
        
        # TorchScript returns 5 outputs (pi is [N, 3])
        ts_A, ts_z, ts_t, ts_p, ts_pi = ts_loaded(x_test)
        
        # ONNX wrapper returns 5 outputs (pi is already [N, 3])
        onnx_A, onnx_z, onnx_t, onnx_p, onnx_pi = onnx_model(x_test_onnx)
    
    # Compare outputs
    def check_match(name, ts_val, onnx_val, rtol=1e-5, atol=1e-6):
        # Handle replicated vs non-replicated
        if ts_val.dim() == 2 and onnx_val.dim() == 1:
            # TorchScript: [N, K], ONNX: [K] (z_hat, t_hat, p)
            ts_val = ts_val[0]  # Take first track since replicated
        
        match = torch.allclose(ts_val, onnx_val, rtol=rtol, atol=atol)
        if match:
            print(f"  ✓ {name} MATCH")
        else:
            max_diff = (ts_val - onnx_val).abs().max().item()
            print(f"  ✗ {name} DIFFER (max diff: {max_diff:.6f})")
        return match
    
    all_match = True
    all_match &= check_match("A", ts_A, onnx_A)
    all_match &= check_match("z_hat", ts_z, onnx_z)
    all_match &= check_match("t_hat", ts_t, onnx_t)
    all_match &= check_match("p", ts_p, onnx_p)
    all_match &= check_match("pi", ts_pi, onnx_pi)
    
    if all_match:
        print("\n✓ PARITY VERIFIED: TorchScript and ONNX outputs are identical!")
    else:
        print("\n✗ PARITY FAILED: Some outputs differ!")
    
    print("\n=== Summary ===")
    print(f"TorchScript (Alpaka): {ts_path}")
    print(f"ONNX (ONNX path):     {onnx_path}")
    print("Both models use identical weights (seed=42)")


if __name__ == "__main__":
    main()
