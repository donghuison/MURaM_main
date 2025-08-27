#!/usr/bin/env python3
"""
Validation script for MURaM Test_3D output
Adapted from verify_lite_3D.py for local testing
"""

import numpy as np
from numpy import linalg 
import sys
import os

print(f"NumPy version: {np.version.version}")

# Path to our test output (current test run)
path_test = '/Users/donghuison/workspace/myGit/MURaM-study/MURaM_main/TEST/Test_3D/3D/'

# Check what iterations are available
available_iters = []
for file in os.listdir(path_test):
    if file.startswith('Header.'):
        iter_num = int(file.split('.')[1])
        available_iters.append(iter_num)

available_iters.sort()
print(f"Available iterations: {available_iters}")

# Use the latest available iteration for validation
if available_iters:
    iter_to_check = available_iters[-1]  # Use the last available iteration
else:
    print("No output files found!")
    sys.exit(1)

print(f"\nValidating iteration: {iter_to_check}")

# Threshold for relative error
thresh = 1e-4

# Transpose for data layout
trans = [0, 1, 2]

def inttostring(ii, ts_size=6):
    """Convert integer to zero-padded string"""
    str_num = str(ii)
    for bb in range(len(str_num), ts_size, 1):
        str_num = '0' + str_num
    return str_num

def read_var_3d(dir, var, iter, layout=None):
    """Read 3D variable from MURaM output files"""
    h = np.loadtxt(dir + 'Header.' + inttostring(iter, ts_size=6))
    
    size = h[0:3].astype(int)
    dx = h[3:6]
    time = h[6]
    
    tmp = np.fromfile(dir + var + '.' + inttostring(iter, ts_size=6), dtype=np.float32)
    tmp = tmp.reshape([size[2], size[1], size[0]])
    
    if layout is not None:
        tmp = tmp.transpose(layout)
    
    return tmp, dx, size, time

# Read various output variables
print("\nReading output files...")

try:
    rho, dx, size, time = read_var_3d(path_test, 'result_prim_0', iter_to_check, trans)
    print(f"✓ Density (rho): shape={rho.shape}, min={rho.min():.6e}, max={rho.max():.6e}")
    
    T, _, _, _ = read_var_3d(path_test, 'eosT', iter_to_check, trans)
    print(f"✓ Temperature (T): shape={T.shape}, min={T.min():.6e}, max={T.max():.6e}")
    
    V, _, _, _ = read_var_3d(path_test, 'result_prim_1', iter_to_check, trans)
    print(f"✓ Velocity (V): shape={V.shape}, min={V.min():.6e}, max={V.max():.6e}")
    
    E, _, _, _ = read_var_3d(path_test, 'result_prim_4', iter_to_check, trans)
    print(f"✓ Energy (E): shape={E.shape}, min={E.min():.6e}, max={E.max():.6e}")
    
    B, _, _, _ = read_var_3d(path_test, 'result_prim_5', iter_to_check, trans)
    print(f"✓ Magnetic Field (B): shape={B.shape}, min={B.min():.6e}, max={B.max():.6e}")
    
    Qrad, _, _, _ = read_var_3d(path_test, 'Qtot', iter_to_check, trans)
    print(f"✓ Radiation (Qrad): shape={Qrad.shape}, min={Qrad.min():.6e}, max={Qrad.max():.6e}")
    
    Jtot, _, _, _ = read_var_3d(path_test, 'Jtot', iter_to_check, trans)
    print(f"✓ Current (Jtot): shape={Jtot.shape}, min={Jtot.min():.6e}, max={Jtot.max():.6e}")
    
    Qcor, _, _, _ = read_var_3d(path_test, 'QxCor', iter_to_check, trans)
    print(f"✓ Corona heating (Qcor): shape={Qcor.shape}, min={Qcor.min():.6e}, max={Qcor.max():.6e}")
    
    cond, _, _, _ = read_var_3d(path_test, 'result_prim_8', iter_to_check, trans)
    print(f"✓ Conduction (cond): shape={cond.shape}, min={cond.min():.6e}, max={cond.max():.6e}")
    
except Exception as e:
    print(f"Error reading files: {e}")
    sys.exit(1)

print(f"\nSimulation time: {time:.6f}")
print(f"Grid size: {size}")
print(f"Grid spacing (dx): {dx}")

# Basic validation checks
print("\n" + "="*50)
print("VALIDATION CHECKS")
print("="*50)

# Check for NaN or Inf values
variables = {
    'Density': rho,
    'Temperature': T,
    'Velocity': V,
    'Energy': E,
    'Magnetic Field': B,
    'Radiation': Qrad,
    'Current': Jtot,
    'Corona heating': Qcor,
    'Conduction': cond
}

all_valid = True
for name, var in variables.items():
    has_nan = np.any(np.isnan(var))
    has_inf = np.any(np.isinf(var))
    
    if has_nan or has_inf:
        print(f"❌ {name}: Contains {'NaN' if has_nan else ''} {'Inf' if has_inf else ''}")
        all_valid = False
    else:
        print(f"✅ {name}: No NaN or Inf values")

# Physical checks
print("\n" + "-"*50)
print("PHYSICAL VALIDITY CHECKS")
print("-"*50)

# Density should be positive
if np.all(rho > 0):
    print(f"✅ Density is positive everywhere")
else:
    print(f"❌ Negative density found!")
    all_valid = False

# Temperature should be positive
if np.all(T > 0):
    print(f"✅ Temperature is positive everywhere")
else:
    print(f"❌ Negative temperature found!")
    all_valid = False

# Energy should be positive
if np.all(E > 0):
    print(f"✅ Energy is positive everywhere")
else:
    print(f"❌ Negative energy found!")
    all_valid = False

# Calculate some statistics
print("\n" + "-"*50)
print("STATISTICAL SUMMARY")
print("-"*50)

for name, var in variables.items():
    mean_val = np.mean(var)
    std_val = np.std(var)
    print(f"{name:20s}: mean={mean_val:12.6e}, std={std_val:12.6e}")

# Conservation checks (if we had reference data)
print("\n" + "-"*50)
print("CONSERVATION CHECKS")
print("-"*50)

# Total mass
total_mass = np.sum(rho) * np.prod(dx)
print(f"Total mass: {total_mass:.6e}")

# Total energy
total_energy = np.sum(E) * np.prod(dx)
print(f"Total energy: {total_energy:.6e}")

# Magnetic energy
B_energy = np.sum(B**2) * np.prod(dx) / (8 * np.pi)
print(f"Magnetic energy: {B_energy:.6e}")

# Final verdict
print("\n" + "="*50)
if all_valid:
    print("✅ VALIDATION PASSED: All basic checks passed")
    print("   Output files appear physically valid")
else:
    print("❌ VALIDATION FAILED: Some checks failed")
    print("   Please review the output above")
print("="*50)

# Compare with initial conditions if available
ini_path = '/Users/donghuison/workspace/myGit/MURaM-study/MURaM_main/TEST/Test_3D/ini/'
if os.path.exists(ini_path):
    print("\n" + "="*50)
    print("COMPARISON WITH INITIAL CONDITIONS")
    print("="*50)
    
    try:
        rho_ini, _, _, _ = read_var_3d(ini_path, 'result_prim_0', 0, trans)
        
        # Calculate relative change from initial conditions
        rho_change = np.mean(np.abs(rho - rho_ini)) / np.mean(np.abs(rho_ini))
        print(f"Mean relative density change: {rho_change:.6e}")
        
        if rho_change < 1.0:
            print(f"✅ Density evolution appears reasonable (< 100% change)")
        else:
            print(f"⚠️  Large density changes detected (> 100% change)")
            
    except Exception as e:
        print(f"Could not read initial conditions: {e}")

print("\n✅ Validation script completed successfully")