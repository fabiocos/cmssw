#!/usr/bin/env python3
"""
export_production_model.py - Export VertexSlotModel to TorchScript

Creates a fresh model with random weights, loads scaler stats, and exports
to TorchScript format compatible with the Alpaka C++ SoA.

Output shapes (for C++ compatibility):
- A:     [N, K]  - assignment probabilities
- z_hat: [N, K]  - replicated from [K]
- t_hat: [N, K]  - replicated from [K]
- p:     [N, K]  - replicated from [K]
- pi:    [N, 3]  - PID weights

Usage:
    python export_production_model.py
"""

import torch
import torch.nn as nn
from typing import Tuple

# Import from local model.py
from model import VertexSlotModel


class VertexSlotModelCMSSW(nn.Module):
    """
    CMSSW-compatible wrapper for VertexSlotModel.
    
    Replicates slot outputs [K] -> [N, K] for SoA compatibility.
    Output order must match GNNOutputSoA: A, z_hat, t_hat, p, pi
    """
    
    def __init__(self, model: VertexSlotModel):
        super().__init__()
        self.model = model
    
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward with replicated outputs for SoA compatibility.
        
        Args:
            x: Input features [1, N, 13] or [N, 13]
            
        Returns:
            A:     [N, K] assignment probabilities
            z_hat: [N, K] replicated z predictions
            t_hat: [N, K] replicated t predictions
            p:     [N, K] replicated existence probabilities
            pi:    [N, 3] PID weights
        """
        # Get native outputs
        A, z_hat, t_hat, p, pi, _, _ = self.model(x)
        
        # A is already [N, K], pi is already [N, 3]
        N = A.size(0)
        
        # Replicate slot outputs [K] -> [N, K]
        z_hat_rep = z_hat.unsqueeze(0).expand(N, -1)  # [K] -> [N, K]
        t_hat_rep = t_hat.unsqueeze(0).expand(N, -1)  # [K] -> [N, K]
        p_rep = p.unsqueeze(0).expand(N, -1)          # [K] -> [N, K]
        
        return A, z_hat_rep, t_hat_rep, p_rep, pi


def main():
    print("=" * 60)
    print("Export Production VertexSlotModel to TorchScript")
    print("=" * 60)
    
    # Configuration matching C++ kNumSlots
    NUM_SLOTS = 180  # Must match kNumSlots in VertexGNNSoA.h
    INPUT_DIM = 13
    SLOT_DIM = 48
    
    # Create fresh model
    print(f"\n1. Creating VertexSlotModel (num_slots={NUM_SLOTS})...")
    model = VertexSlotModel(
        input_dim=INPUT_DIM,
        slot_dim=SLOT_DIM,
        num_slots=NUM_SLOTS,
        num_iterations=4,
        hidden_dim=32,
        k=30,
        dropout=0.0,
        n_heads_encoder=2,
        n_heads_slot=4,
    )
    
    # Load scaler stats
    print("\n2. Loading scaler statistics from scaler_stats.pt...")
    scaler_stats = torch.load("scaler_stats.pt", map_location="cpu")
    model.featnorm.mean.copy_(scaler_stats["mean"])
    model.featnorm.std.copy_(scaler_stats["std"])
    print(f"   mean: {model.featnorm.mean[:3].tolist()}...")
    print(f"   std:  {model.featnorm.std[:3].tolist()}...")
    
    # Wrap for CMSSW compatibility
    print("\n3. Wrapping model for CMSSW output shapes...")
    cmssw_model = VertexSlotModelCMSSW(model)
    cmssw_model.eval()
    
    # Test forward pass
    print("\n4. Testing forward pass...")
    N_test = 100
    x_test = torch.randn(1, N_test, INPUT_DIM)
    x_test[0, :, 1] = x_test[0, :, 1].abs() + 0.01  # sigma_z > 0
    x_test[0, :, 9:12] = x_test[0, :, 9:12].abs() + 0.01  # sigma_t > 0
    x_test[0, :, 12] = (torch.rand(N_test) > 0.3).float()  # has_time
    
    with torch.no_grad():
        A, z_hat, t_hat, p, pi = cmssw_model(x_test)
    
    print(f"   A:     {A.shape}")        # [N, K]
    print(f"   z_hat: {z_hat.shape}")    # [N, K]
    print(f"   t_hat: {t_hat.shape}")    # [N, K]
    print(f"   p:     {p.shape}")        # [N, K]
    print(f"   pi:    {pi.shape}")       # [N, 3]
    
    # Verify wrapper replication preserves original values
    print("\n4b. Verifying wrapper replication preserves values...")
    with torch.no_grad():
        A_orig, z_orig, t_orig, p_orig, pi_orig, _, _ = model(x_test)
    
    # z_hat[i,:] should equal z_orig[:] for all i
    z_diff = (z_hat[0, :] - z_orig).abs().max().item()
    t_diff = (t_hat[0, :] - t_orig).abs().max().item()
    p_diff = (p[0, :] - p_orig).abs().max().item()
    
    # All rows should be identical (replicated)
    z_row_diff = (z_hat - z_hat[0:1, :]).abs().max().item()
    t_row_diff = (t_hat - t_hat[0:1, :]).abs().max().item()
    p_row_diff = (p - p[0:1, :]).abs().max().item()
    
    print(f"   z_hat vs original: {z_diff:.2e}, row consistency: {z_row_diff:.2e}")
    print(f"   t_hat vs original: {t_diff:.2e}, row consistency: {t_row_diff:.2e}")
    print(f"   p vs original:     {p_diff:.2e}, row consistency: {p_row_diff:.2e}")
    
    if max(z_diff, t_diff, p_diff, z_row_diff, t_row_diff, p_row_diff) < 1e-6:
        print("   ✓ Replication verified - wrapper is correct!")
    
    # Export to TorchScript
    print("\n5. Exporting to TorchScript...")
    output_path = "vertex_slot_production.pt"
    scripted = torch.jit.script(cmssw_model)
    scripted.save(output_path)
    print(f"   ✓ Saved to {output_path}")
    
    # Verify loaded model
    print("\n6. Verifying loaded TorchScript...")
    loaded = torch.jit.load(output_path)
    loaded.eval()
    
    with torch.no_grad():
        A2, z2, t2, p2, pi2 = loaded(x_test)
    
    max_diff = max(
        (A - A2).abs().max().item(),
        (z_hat - z2).abs().max().item(),
        (t_hat - t2).abs().max().item(),
        (p - p2).abs().max().item(),
        (pi - pi2).abs().max().item(),
    )
    print(f"   Max diff: {max_diff:.2e}")
    if max_diff < 1e-5:
        print("   ✓ Verification passed!")
    else:
        print("   ✗ Outputs differ - investigate!")
    
    # Copy to CMSSW data directory
    print("\n7. Copying to CMSSW data directory...")
    import shutil
    dest = "RecoVertex/PrimaryVertexProducer/data/vertex_slot_production.pt"
    shutil.copy(output_path, dest)
    print(f"   ✓ Copied to {dest}")
    
    print("\n" + "=" * 60)
    print("Export complete!")
    print("=" * 60)
    print(f"\nTo use in configs:")
    print(f'  model = cms.FileInPath("RecoVertex/PrimaryVertexProducer/data/vertex_slot_production.pt")')


if __name__ == "__main__":
    main()
