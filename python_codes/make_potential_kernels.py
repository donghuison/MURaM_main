#!/usr/bin/env python3
import numpy as np

# Calculate the potential field kernels for use with the MURaM code.
# Based on the IDL routines of M.C. 2009.

# Vertical is z.

# Simulation Setup
out_dir = "../RUNDIR/"

# Set Ny = 1 for 2D
# NB dx must = dy for this to work
Nx = 128
Ny = 1
Nz = 128
n_kernel_levels = 3

# level 0: boundary-reference level
# level 1,2: first/second ghost-layer-corresponding levels
# therefore n_kernel_levels is not a ghost-only count
XY_LEVELS = (0, 1, 2)
Z_LEVELS = (1, 2)

dx = 8.0e8 / Nx
dy = dx
dz = 4.0e8 / Nz

res = 16

required_levels = max(max(XY_LEVELS), max(Z_LEVELS)) + 1
if n_kernel_levels < required_levels:
    raise ValueError(
        "n_kernel_levels is too small for requested XY_LEVELS/Z_LEVELS."
    )
if dx != dy:
    raise ValueError("dx must equal dy for this kernel construction.")

if not out_dir.endswith("/"):
    out_dir = out_dir + "/"

heightscale = dz / dx

Nxs = Nx * res
if Ny > 1:
    Nys = Ny * res
else:
    Nys = 1

# Define Fourier grid
kx = np.zeros([Nxs])
ky = np.zeros([Nys])

kx[0 : Nxs // 2] = np.arange(Nxs // 2) / Nxs
kx[Nxs // 2 : Nxs] = np.arange(Nxs // 2) / Nxs - 0.5

kx = np.array([kx[:]] * Nys).transpose()
kx *= 2 * np.pi

if Ny > 1:
    ky[0 : Nys // 2] = np.arange(Nys // 2) / Nys
    ky[Nys // 2 :] = np.arange(Nys // 2) / Nys - 0.5

ky = np.array([ky[:]] * Nxs)
ky *= 2 * np.pi

k2 = kx * kx + ky * ky
kabs = np.sqrt(k2)

# Complex Arrays
HxB = np.zeros([Nxs, Nys], dtype=complex)
HyB = np.zeros([Nxs, Nys], dtype=complex)
HzB = np.ones([Nxs, Nys], dtype=complex)

HxB[np.where(k2 != 0)] = -1j * kx[np.where(k2 != 0)] / kabs[np.where(k2 != 0)]
HxB[0, 0] = -1j

if Ny == 1:
    HxB[0, :] = -1j

HyB[np.where(k2 != 0)] = -1j * ky[np.where(k2 != 0)] / kabs[np.where(k2 != 0)]
HyB[0, 0] = -1j

if Ny == 1:
    HyB[0, :] = -1j

# Delta Function
delta = np.zeros([Nxs, Nys])

for i in range(-4 * res, 4 * res + 1):
    i0 = i
    while i0 < 0:
        i0 = i0 + Nxs
    for j in range(-4 * res, 4 * res + 1):
        j0 = j
        while j0 < 0:
            j0 = j0 + Nys
        if Ny == 1:
            delta[i0, :] = np.exp(-np.double(np.fmod(i, Nxs)) ** 2 / res**2)
        else:
            delta[i0, j0] = np.exp(
                -(
                    np.double(np.fmod(i, Nxs)) ** 2
                    + np.double(np.fmod(j, Nys)) ** 2
                )
                / res**2
            )

delta /= delta.sum()

# FFT delta function
FFTdelta = np.fft.fft2(delta) / Nxs / Nys


# define function symmetric rebin
def symmetric_rebin(a, oNx, oNy):
    Nx = a.shape[0]
    Ny = a.shape[1]
    dx = Nx / oNx
    dy = Ny / oNy

    dxi = int(dx)
    dyi = int(dy)

    out = np.zeros([oNx, oNy])
    b = np.zeros([Nx + 2 * dxi, Ny + 2 * dyi])

    b[dxi : dxi + Nx, dyi : dyi + Ny] = a

    # Fill y ghost cells
    b[0:dxi, dyi : Ny + dyi] = a[Nx - dxi : Nx, 0:Ny]
    b[Nx + dxi : Nx + 2 * dxi, dyi : Ny + dyi] = a[0:dxi, 0:Ny]

    # Fill x ghost cells
    b[:, 0:dyi] = b[:, Ny : Ny + dyi]
    b[:, Ny + dyi : Ny + 2 * dyi] = b[:, dyi : 2 * dyi]

    cx = int(dx / 2)
    cy = int(dy / 2)

    for i in range(dxi, Nx + dxi, dxi):
        for j in range(dyi, Ny + dyi, dyi):
            out[i // dxi - 1, j // dyi - 1] = b[
                i - cx : i + cx + 1, j - cy : j + cy + 1
            ].sum()

    return out


# Make and save kernels
xbkernel = np.zeros([Nx, Ny])
ybkernel = np.zeros([Nx, Ny])
zbkernel = np.zeros([Nx, Ny])

a = np.zeros(Nx * Ny + 3, dtype=np.single)
a[0:3] = [Nx, Ny, heightscale]

saved_files = []

for level_idx in range(n_kernel_levels):
    level_height = level_idx * heightscale * res

    xbkernel[:, :] = symmetric_rebin(
        np.fft.ifft2(FFTdelta * HxB * np.exp(-kabs * level_height)).real, Nx, Ny
    )
    ybkernel[:, :] = symmetric_rebin(
        np.fft.ifft2(FFTdelta * HyB * np.exp(-kabs * level_height)).real, Nx, Ny
    )
    zbkernel[:, :] = symmetric_rebin(
        np.fft.ifft2(FFTdelta * HzB * np.exp(-kabs * level_height)).real, Nx, Ny
    )

    flux = zbkernel[:, :].sum()

    xbkernel /= flux
    ybkernel /= flux
    zbkernel /= flux

    print(
        level_idx,
        flux,
        np.abs(xbkernel[:, :]).sum() / Nx / Ny,
        np.abs(ybkernel[:, :]).sum() / Nx / Ny,
        np.abs(zbkernel[:, :]).sum() / Nx / Ny,
    )

    if level_idx in XY_LEVELS:
        a[3:] = xbkernel[:, :].transpose().ravel()
        a.tofile(out_dir + "PSF-kernel-x-" + str(level_idx) + ".dat")
        saved_files.append("PSF-kernel-x-" + str(level_idx) + ".dat")

        a[3:] = ybkernel[:, :].transpose().ravel()
        a.tofile(out_dir + "PSF-kernel-y-" + str(level_idx) + ".dat")
        saved_files.append("PSF-kernel-y-" + str(level_idx) + ".dat")

    if level_idx in Z_LEVELS:
        a[3:] = zbkernel[:, :].transpose().ravel()
        a.tofile(out_dir + "PSF-kernel-z-" + str(level_idx) + ".dat")
        saved_files.append("PSF-kernel-z-" + str(level_idx) + ".dat")

print("Saved files:")
for name in saved_files:
    print("  " + name)
