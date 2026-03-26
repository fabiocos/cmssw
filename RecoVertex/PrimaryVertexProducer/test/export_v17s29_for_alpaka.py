#!/usr/bin/env python3
"""
export_v17s29_for_alpaka.py - Wrap pre-exported v17s29 TorchScript for Alpaka SoA

The input vertex_slot_model_v17s29.pt is already a TorchScript model with 5 outputs:
  A:     [N, K=220]
  z_hat: [K]
  t_hat: [K]
  p:     [K]
  pi:    [N, 3]

This script wraps it to replicate [K] → [N,K] for Alpaka SoA compatibility.

Usage:
    source /afs/cern.ch/work/p/prsolank/public/fastgraphenv/bin/activate
    cd /eos/user/p/prsolank/www/CMSSW_16_1_0_pre1/src/RecoVertex/PrimaryVertexProducer/test
    python export_v17s29_for_alpaka.py
"""

import os
from typing import Tuple

import torch
import torch.nn as nn


class AlpakaWrapper(nn.Module):
    """Wraps a pre-exported TorchScript model: replicates [K] → [N,K]."""
    
    def __init__(self, inner: torch.jit.ScriptModule):
        super().__init__()
        self.inner = inner
    
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        A, z_hat, t_hat, p, pi = self.inner(x)
        N = A.size(0)
        z_hat_rep = z_hat.unsqueeze(0).expand(N, -1)
        t_hat_rep = t_hat.unsqueeze(0).expand(N, -1)
        p_rep = p.unsqueeze(0).expand(N, -1)
        return A, z_hat_rep, t_hat_rep, p_rep, pi


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(script_dir, "vertex_slot_model_v17s29.pt")
    output_path = os.path.join(script_dir, "vertex_slot_v17s29_alpaka.pt")
    
    print("=" * 70)
    print("Wrap v17s29 TorchScript for Alpaka SoA ([K] → [N,K])")
    print("=" * 70)
    
    # Load original TorchScript
    print(f"\n1. Loading original TorchScript: {input_path}")
    original = torch.jit.load(input_path, map_location='cpu')
    original.eval()
    
    # Test original outputs
    print(f"\n2. Testing original model outputs...")
    N_test = 100
    x_test = torch.randn(1, N_test, 13)
    x_test[0, :, 1] = x_test[0, :, 1].abs() + 0.01
    x_test[0, :, 9:12] = x_test[0, :, 9:12].abs() + 0.01
    x_test[0, :, 12] = (torch.rand(N_test) > 0.3).float()
    
    with torch.no_grad():
        A_orig, z_orig, t_orig, p_orig, pi_orig = original(x_test)
    print(f"   A:     {A_orig.shape}")
    print(f"   z_hat: {z_orig.shape}")
    print(f"   t_hat: {t_orig.shape}")
    print(f"   p:     {p_orig.shape}")
    print(f"   pi:    {pi_orig.shape}")
    
    K = z_orig.shape[0]
    print(f"   → K={K} slots detected")
    
    if K != 220:
        print(f"   ⚠ WARNING: Expected K=220 (kNumSlots in VertexGNNSoA.h), got K={K}!")
    
    # Wrap for Alpaka SoA
    print(f"\n3. Wrapping for Alpaka SoA compatibility [K] → [N,K]")
    wrapped = AlpakaWrapper(original)
    wrapped.eval()
    
    with torch.no_grad():
        A_w, z_w, t_w, p_w, pi_w = wrapped(x_test)
    print(f"   A:     {A_w.shape} (expected: [N, K={K}])")
    print(f"   z_hat: {z_w.shape} (expected: [N, K={K}])")
    print(f"   t_hat: {t_w.shape} (expected: [N, K={K}])")
    print(f"   p:     {p_w.shape} (expected: [N, K={K}])")
    print(f"   pi:    {pi_w.shape} (expected: [N, 3])")
    
    # Parity: original vs wrapped
    print(f"\n4. Parity check: original vs wrapped")
    diffs = {
        'A': (A_orig - A_w).abs().max().item(),
        'z_hat': (z_orig - z_w[0, :]).abs().max().item(),
        't_hat': (t_orig - t_w[0, :]).abs().max().item(),
        'p': (p_orig - p_w[0, :]).abs().max().item(),
        'pi': (pi_orig - pi_w).abs().max().item(),
    }
    for k, v in diffs.items():
        status = "✓" if v < 1e-6 else "✗"
        print(f"   {k:8s}: {v:.2e} {status}")
    if all(v < 1e-6 for v in diffs.values()):
        print("   ✓ ALL MATCH — wrapper is transparent!")
    else:
        print("   ✗ MISMATCH — investigate!")
    
    # Replication consistency
    print(f"\n5. Replication consistency (all rows identical):")
    z_row_diff = (z_w - z_w[0:1, :]).abs().max().item()
    t_row_diff = (t_w - t_w[0:1, :]).abs().max().item()
    p_row_diff = (p_w - p_w[0:1, :]).abs().max().item()
    print(f"   z_hat rows: {z_row_diff:.2e} {'✓' if z_row_diff < 1e-6 else '✗'}")
    print(f"   t_hat rows: {t_row_diff:.2e} {'✓' if t_row_diff < 1e-6 else '✗'}")
    print(f"   p rows:     {p_row_diff:.2e} {'✓' if p_row_diff < 1e-6 else '✗'}")
    
    # Export to TorchScript
    print(f"\n6. Exporting wrapped model to TorchScript: {output_path}")
    scripted = torch.jit.script(wrapped)
    scripted.save(output_path)
    print(f"   ✓ Saved!")
    
    # Verify loaded model
    print(f"\n7. Verifying loaded TorchScript...")
    loaded = torch.jit.load(output_path, map_location='cpu')
    with torch.no_grad():
        A_l, z_l, t_l, p_l, pi_l = loaded(x_test)
    max_diff = max(
        (A_w - A_l).abs().max().item(),
        (z_w - z_l).abs().max().item(),
        (t_w - t_l).abs().max().item(),
        (p_w - p_l).abs().max().item(),
        (pi_w - pi_l).abs().max().item(),
    )
    print(f"   Max diff (wrapped vs loaded): {max_diff:.2e}")
    if max_diff < 1e-5:
        print("   ✓ Verification passed!")
    else:
        print("   ✗ Verification FAILED!")
    
    # Input shape flexibility
    print(f"\n8. Input shape flexibility ([1,N,13] vs [N,13]):")
    x_2d = x_test.squeeze(0)
    with torch.no_grad():
        A_2d, z_2d, t_2d, p_2d, pi_2d = loaded(x_2d)
    diff_2d = max(
        (A_w - A_2d).abs().max().item(),
        (z_w - z_2d).abs().max().item(),
        (t_w - t_2d).abs().max().item(),
        (p_w - p_2d).abs().max().item(),
        (pi_w - pi_2d).abs().max().item(),
    )
    print(f"   Max diff: {diff_2d:.2e}")
    if diff_2d < 1e-6:
        print(f"   ✓ Both input shapes produce identical outputs!")
    else:
        print(f"   ⚠ Outputs differ — investigate!")

    print("\n" + "=" * 70)
    print("SUCCESS!")
    print(f"Output: {output_path}")
    print(f"Copy to CMSSW test dir (already there if run from test/)")
    print("=" * 70)


if __name__ == "__main__":
    main()
