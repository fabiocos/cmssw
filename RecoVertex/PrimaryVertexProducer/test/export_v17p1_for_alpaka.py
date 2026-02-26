#!/usr/bin/env python3
"""
export_v17p1_for_alpaka.py - Export v17p1 model with Alpaka SoA-compatible shapes

Follows the exact export logic from vertex_slot_train_v17p1.py:
1. Load model weights from ckpt["model"]
2. Load scaler stats from scaler_stats.pt
3. Use export_to_torchscript() but wrap in CMSSW wrapper for [K] -> [N,K] replication

Usage:
    source /afs/cern.ch/work/p/prsolank/public/fastgraphenv/bin/activate
    cd /eos/user/p/prsolank/work/Develop/CheckReco/CMSSW_16_1_0_pre1/src/RecoVertex/PrimaryVertexProducer/test
    python export_v17p1_for_alpaka.py
"""

import sys
import os
import argparse
import importlib.util
from typing import Tuple

import torch
import torch.nn as nn


# =============================================================================
# CMSSW-Compatible Wrapper (replicates slot outputs [K] → [N, K])
# =============================================================================

class VertexSlotModelCMSSW(nn.Module):
    """Replicates slot outputs [K] → [N, K] for Alpaka SoA compatibility."""
    
    def __init__(self, model: nn.Module):
        super().__init__()
        self.model = model
    
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        # v17p1 VertexSlotModel returns 7 outputs: A, z_hat, t_hat, p, pi, pi_logits, s
        A, z_hat, t_hat, p, pi, _, _ = self.model(x)
        N = A.size(0)
        z_hat_rep = z_hat.unsqueeze(0).expand(N, -1)
        t_hat_rep = t_hat.unsqueeze(0).expand(N, -1)
        p_rep = p.unsqueeze(0).expand(N, -1)
        return A, z_hat_rep, t_hat_rep, p_rep, pi


def load_module_from_path(module_name: str, file_path: str):
    """Dynamically load a Python module from a file path."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def main():
    # Paths - exactly matching training script
    ckpt_dir = "/afs/cern.ch/work/p/prsolank/private/OC_FastGraph/checkpoints_vertex_slot_v17p1"
    model_src = "/afs/cern.ch/work/p/prsolank/private/OC_FastGraph/vertex_slot_model_v17p1.py"
    ckpt_path = os.path.join(ckpt_dir, "best_vertex_slot.pt")
    scaler_path = os.path.join(ckpt_dir, "scaler_stats.pt")
    output_path = "vertex_slot_v17p1_alpaka.pt"
    
    print("=" * 70)
    print("Export v17p1 for Alpaka SoA (following train script logic)")
    print("=" * 70)
    
    # Load model definition
    print(f"\n1. Loading model definition from {model_src}")
    model_module = load_module_from_path("vertex_slot_model_v17p1", model_src)
    VertexSlotModel = model_module.VertexSlotModel
    
    # Load checkpoint - exactly as in train script line 1284
    print(f"\n2. Loading checkpoint from {ckpt_path}")
    ckpt = torch.load(ckpt_path, map_location='cpu', weights_only=False)
    print(f"   Checkpoint keys: {list(ckpt.keys())}")
    print(f"   Epoch: {ckpt.get('epoch', 'N/A')}, best_val: {ckpt.get('best_val', 'N/A')}")
    
    # Create model with EXACT same args as train script lines 1285-1300
    print(f"\n3. Creating model (matching train script)")
    export_model = VertexSlotModel(
        input_dim=13,
        slot_dim=48,
        num_slots=180,
        num_iterations=4,
        hidden_dim=64,
        k=24,
        num_tokens=128,
        heads_global=4,
        heads_local=4,
        edge_dim=32,
        knn_space_dim=16,
        local_layers=1,
        dropout=0.0,
        n_heads_slot=4,
    )
    
    # Load weights - exactly as train script line 1301
    export_model.load_state_dict(ckpt["model"])
    print("   ✓ Loaded model weights")
    
    # Load scaler stats - exactly as train script lines 1304-1308
    print(f"\n4. Loading scaler stats from {scaler_path}")
    scaler_stats = torch.load(scaler_path, map_location='cpu', weights_only=False)
    export_model.featnorm.mean.copy_(scaler_stats['mean'])
    export_model.featnorm.std.copy_(scaler_stats['std'])
    print(f"   mean[:3]: {export_model.featnorm.mean[:3].tolist()}")
    print(f"   std[:3]:  {export_model.featnorm.std[:3].tolist()}")
    
    # Wrap for Alpaka SoA
    print(f"\n5. Wrapping for Alpaka SoA compatibility [K] -> [N,K]")
    cmssw_model = VertexSlotModelCMSSW(export_model)
    cmssw_model.eval()
    
    # Test shapes
    print(f"\n6. Testing shapes...")
    N_test = 100
    x_test = torch.randn(1, N_test, 13)
    x_test[0, :, 1] = x_test[0, :, 1].abs() + 0.01
    x_test[0, :, 9:12] = x_test[0, :, 9:12].abs() + 0.01
    x_test[0, :, 12] = (torch.rand(N_test) > 0.3).float()
    
    with torch.no_grad():
        A, z, t, p, pi = cmssw_model(x_test)
    print(f"   A:     {A.shape} (expected: [N, K])")
    print(f"   z_hat: {z.shape} (expected: [N, K])")
    print(f"   t_hat: {t.shape} (expected: [N, K])")
    print(f"   p:     {p.shape} (expected: [N, K])")
    print(f"   pi:    {pi.shape} (expected: [N, 3])")
    
    # Test input shape flexibility: [1,N,13] vs [N,13]
    print(f"\n   --- Input Shape Flexibility ---")
    x_3d = x_test                      # [1, N, 13]
    x_2d = x_test.squeeze(0)           # [N, 13]
    with torch.no_grad():
        A_3d, z_3d, t_3d, p_3d, pi_3d = cmssw_model(x_3d)
        A_2d, z_2d, t_2d, p_2d, pi_2d = cmssw_model(x_2d)
    diff_max = max(
        (A_3d - A_2d).abs().max().item(),
        (z_3d - z_2d).abs().max().item(),
        (t_3d - t_2d).abs().max().item(),
        (p_3d - p_2d).abs().max().item(),
        (pi_3d - pi_2d).abs().max().item(),
    )
    print(f"   [1,N,13] vs [N,13]: max_diff = {diff_max:.2e}")
    if diff_max < 1e-6:
        print(f"   ✓ Both input shapes produce identical outputs!")
    else:
        print(f"   ⚠ Outputs differ - investigate!")
    
    # Export to TorchScript
    print(f"\n7. Exporting to TorchScript: {output_path}")
    scripted = torch.jit.script(cmssw_model)
    scripted.save(output_path)
    print(f"   ✓ Saved!")
    
    # Verify
    print(f"\n8. Verifying loaded TorchScript...")
    loaded = torch.jit.load(output_path)
    with torch.no_grad():
        A2, z2, t2, p2, pi2 = loaded(x_test)
    max_diff = max(
        (A - A2).abs().max().item(),
        (z - z2).abs().max().item(),
        (t - t2).abs().max().item(),
        (p - p2).abs().max().item(),
        (pi - pi2).abs().max().item(),
    )
    print(f"   Max diff: {max_diff:.2e}")
    if max_diff < 1e-5:
        print("   ✓ Verification passed!")
    
    # =============================================================================
    # THREE-WAY PARITY TEST
    # =============================================================================
    print(f"\n" + "=" * 70)
    print("THREE-WAY PARITY TEST")
    print("=" * 70)
    
    # 1. Native PyTorch (7 outputs)
    print(f"\n1. Native PyTorch model (7 outputs)")
    with torch.no_grad():
        A_native, z_native, t_native, p_native, pi_native, _, _ = export_model(x_test)
    print(f"   A:     {A_native.shape}")
    print(f"   z_hat: {z_native.shape}")
    print(f"   t_hat: {t_native.shape}")
    print(f"   p:     {p_native.shape}")
    print(f"   pi:    {pi_native.shape}")
    
    # 2. Old TorchScript (5 outputs, [K] shapes)
    old_ts_path = os.path.join(ckpt_dir, "vertex_slot_model.pt")
    print(f"\n2. Old TorchScript: {old_ts_path}")
    if os.path.exists(old_ts_path):
        old_ts = torch.jit.load(old_ts_path)
        old_ts.eval()
        with torch.no_grad():
            A_old, z_old, t_old, p_old, pi_old = old_ts(x_test)
        print(f"   A:     {A_old.shape}")
        print(f"   z_hat: {z_old.shape}")
        print(f"   t_hat: {t_old.shape}")
        print(f"   p:     {p_old.shape}")
        print(f"   pi:    {pi_old.shape}")
    else:
        print(f"   ⚠ Not found, skipping")
        A_old = z_old = t_old = p_old = pi_old = None
    
    # 3. New CMSSW TorchScript (5 outputs, [N,K] shapes)
    print(f"\n3. New CMSSW TorchScript: {output_path}")
    print(f"   A:     {A2.shape}")
    print(f"   z_hat: {z2.shape}")
    print(f"   t_hat: {t2.shape}")
    print(f"   p:     {p2.shape}")
    print(f"   pi:    {pi2.shape}")
    
    # Parity checks
    print(f"\n--- PARITY RESULTS ---")
    
    # Native vs Old TorchScript
    if A_old is not None:
        print(f"\nNative PyTorch vs Old TorchScript:")
        diffs = {
            'A': (A_native - A_old).abs().max().item(),
            'z_hat': (z_native - z_old).abs().max().item(),
            't_hat': (t_native - t_old).abs().max().item(),
            'p': (p_native - p_old).abs().max().item(),
            'pi': (pi_native - pi_old).abs().max().item(),
        }
        for k, v in diffs.items():
            status = "✓" if v < 1e-5 else "✗"
            print(f"   {k:8s}: {v:.2e} {status}")
        if all(v < 1e-5 for v in diffs.values()):
            print("   ✓ ALL MATCH!")
    
    # Native vs New CMSSW TorchScript
    print(f"\nNative PyTorch vs New CMSSW TorchScript:")
    # A and pi should match directly
    # z/t/p: new has [N,K] replicated, native has [K] - compare row 0
    diffs = {
        'A': (A_native - A2).abs().max().item(),
        'z_hat': (z_native - z2[0, :]).abs().max().item(),
        't_hat': (t_native - t2[0, :]).abs().max().item(),
        'p': (p_native - p2[0, :]).abs().max().item(),
        'pi': (pi_native - pi2).abs().max().item(),
    }
    for k, v in diffs.items():
        status = "✓" if v < 1e-5 else "✗"
        print(f"   {k:8s}: {v:.2e} {status}")
    if all(v < 1e-5 for v in diffs.values()):
        print("   ✓ ALL MATCH!")
    
    # Replication consistency
    print(f"\nReplication consistency (all rows identical):")
    z_row_diff = (z2 - z2[0:1, :]).abs().max().item()
    t_row_diff = (t2 - t2[0:1, :]).abs().max().item()
    p_row_diff = (p2 - p2[0:1, :]).abs().max().item()
    print(f"   z_hat rows: {z_row_diff:.2e} {'✓' if z_row_diff < 1e-6 else '✗'}")
    print(f"   t_hat rows: {t_row_diff:.2e} {'✓' if t_row_diff < 1e-6 else '✗'}")
    print(f"   p rows:     {p_row_diff:.2e} {'✓' if p_row_diff < 1e-6 else '✗'}")
    
    print("\n" + "=" * 70)
    print("SUCCESS!")
    print(f"Output: {output_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
