#!/usr/bin/env python3
"""
Simple visualization script for MURaM Test_3D output
Creates plots of key variables at different heights
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Path to test output
path_test = '/Users/donghuison/workspace/myGit/MURaM-study/MURaM_main/TEST/Test_3D/3D/'

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

# Use the latest iteration
iter_to_plot = 20
trans = [0, 1, 2]

print(f"Creating visualizations for iteration {iter_to_plot}...")

# Read data
rho, dx, size, time = read_var_3d(path_test, 'result_prim_0', iter_to_plot, trans)
T, _, _, _ = read_var_3d(path_test, 'eosT', iter_to_plot, trans)
V, _, _, _ = read_var_3d(path_test, 'result_prim_1', iter_to_plot, trans)
B, _, _, _ = read_var_3d(path_test, 'result_prim_5', iter_to_plot, trans)

# Create figure with subplots
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
fig.suptitle(f'MURaM Test_3D Output - Iteration {iter_to_plot} (t={time:.3f}s)', fontsize=14)

# Select heights to plot (bottom, middle, top)
z_indices = [10, size[2]//2, size[2]-10]
z_labels = ['Near Bottom', 'Middle', 'Near Top']

# Density plots
for i, (z_idx, label) in enumerate(zip(z_indices[:2], z_labels[:2])):
    ax = axes[0, i]
    im = ax.imshow(np.log10(rho[:, :, z_idx]), cmap='viridis', aspect='auto')
    ax.set_title(f'Log Density - {label}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='log10(ρ)')

# Temperature plots
for i, (z_idx, label) in enumerate(zip(z_indices[:2], z_labels[:2])):
    ax = axes[0, i+2]
    im = ax.imshow(T[:, :, z_idx]/1000, cmap='hot', aspect='auto')
    ax.set_title(f'Temperature - {label}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='T (kK)')

# Velocity plots
for i, (z_idx, label) in enumerate(zip(z_indices[:2], z_labels[:2])):
    ax = axes[1, i]
    im = ax.imshow(V[:, :, z_idx]/1e5, cmap='RdBu_r', aspect='auto')
    ax.set_title(f'Velocity - {label}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='V (km/s)')

# Magnetic field plots
for i, (z_idx, label) in enumerate(zip(z_indices[:2], z_labels[:2])):
    ax = axes[1, i+2]
    im = ax.imshow(B[:, :, z_idx], cmap='PuOr', aspect='auto')
    ax.set_title(f'Magnetic Field - {label}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.colorbar(im, ax=ax, label='B (G)')

plt.tight_layout()

# Save figure
output_file = 'muram_test3d_visualization.png'
plt.savefig(output_file, dpi=100, bbox_inches='tight')
print(f"✅ Visualization saved to: {output_file}")

# Create a vertical profile plot
fig2, axes2 = plt.subplots(1, 4, figsize=(16, 4))
fig2.suptitle(f'Vertical Profiles (horizontally averaged) - Iteration {iter_to_plot}', fontsize=14)

# Calculate horizontal averages
rho_avg = np.mean(rho, axis=(0, 1))
T_avg = np.mean(T, axis=(0, 1))
V_avg = np.mean(V, axis=(0, 1))
B_rms = np.sqrt(np.mean(B**2, axis=(0, 1)))

# Heights in Mm - note size is [nx, ny, nz] but data is [nx, ny, nz] after transpose
# The z-dimension is the first one after transpose
z = np.arange(rho.shape[2]) * dx[2] / 1e8

# Density profile
axes2[0].semilogy(z, rho_avg)
axes2[0].set_xlabel('Height (Mm)')
axes2[0].set_ylabel('Density (g/cm³)')
axes2[0].set_title('Density Profile')
axes2[0].grid(True, alpha=0.3)

# Temperature profile
axes2[1].plot(z, T_avg/1000)
axes2[1].set_xlabel('Height (Mm)')
axes2[1].set_ylabel('Temperature (kK)')
axes2[1].set_title('Temperature Profile')
axes2[1].grid(True, alpha=0.3)

# Velocity profile
axes2[2].plot(z, V_avg/1e5)
axes2[2].set_xlabel('Height (Mm)')
axes2[2].set_ylabel('Velocity (km/s)')
axes2[2].set_title('Vertical Velocity Profile')
axes2[2].axhline(y=0, color='k', linestyle='--', alpha=0.5)
axes2[2].grid(True, alpha=0.3)

# Magnetic field profile
axes2[3].plot(z, B_rms)
axes2[3].set_xlabel('Height (Mm)')
axes2[3].set_ylabel('B_rms (G)')
axes2[3].set_title('RMS Magnetic Field')
axes2[3].grid(True, alpha=0.3)

plt.tight_layout()

# Save profile figure
profile_file = 'muram_test3d_profiles.png'
plt.savefig(profile_file, dpi=100, bbox_inches='tight')
print(f"✅ Profile plot saved to: {profile_file}")

print("\n✨ Visualization complete!")