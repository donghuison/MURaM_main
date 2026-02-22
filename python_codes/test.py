import read_muram as rmu
import dp_plot_tools as dplt
import muram_eos as eos
import numpy as np
import matplotlib.pyplot as plt
import os

G_SOLAR_REF = 2.74e4  # cm/s^2
R_SUN_MM = 695.7


def hydrostatic_gravity(z_mm, p_cgs, rho_cgs):
    mask = (
        np.isfinite(z_mm)
        & np.isfinite(p_cgs)
        & np.isfinite(rho_cgs)
        & (p_cgs > 0.0)
        & (rho_cgs > 0.0)
    )
    if np.count_nonzero(mask) < 3:
        raise ValueError("Need at least 3 valid points to compute hydrostatic gravity.")

    z = z_mm[mask]
    p = p_cgs[mask]
    rho = rho_cgs[mask]

    order = np.argsort(z)
    z_sorted = z[order]
    p_sorted = p[order]
    rho_sorted = rho[order]

    z_cm = z_sorted * 1.0e8
    dp_dz = np.gradient(p_sorted, z_cm, edge_order=2)
    g_hse = -dp_dz / rho_sorted

    return z_sorted, g_hse


def finite_stats(arr):
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.nan, np.nan, np.nan
    vals = arr[finite]
    return float(np.nanmin(vals)), float(np.nanmedian(vals)), float(np.nanmax(vals))


def invert_energy_from_rho_T(mu_eos, rho_arr, T_arr):
    rho_arr = np.asarray(rho_arr, dtype=float)
    T_arr = np.asarray(T_arr, dtype=float)
    out = np.full(rho_arr.shape, np.nan, dtype=float)
    eps_axis = np.exp(mu_eos.xeps) - mu_eos.eps_off
    for i in range(rho_arr.size):
        T_axis = mu_eos.interp_T(eps_axis, np.full(mu_eos.n_eps, rho_arr[i], dtype=float))
        valid = np.isfinite(T_axis) & np.isfinite(mu_eos.xeps)
        if np.count_nonzero(valid) < 2:
            continue
        out[i] = np.exp(np.interp(T_arr[i], T_axis[valid], mu_eos.xeps[valid])) - mu_eos.eps_off
    return out


def cp_from_gamma(gamma_arr, p_arr, rho_arr, T_arr, min_gamma=1.001):
    gamma_arr = np.asarray(gamma_arr, dtype=float)
    p_arr = np.asarray(p_arr, dtype=float)
    rho_arr = np.asarray(rho_arr, dtype=float)
    T_arr = np.asarray(T_arr, dtype=float)

    cp = np.full(p_arr.shape, np.nan, dtype=float)
    valid = (
        np.isfinite(gamma_arr)
        & np.isfinite(p_arr)
        & np.isfinite(rho_arr)
        & np.isfinite(T_arr)
        & (gamma_arr > min_gamma)
        & (p_arr > 0.0)
        & (rho_arr > 0.0)
        & (T_arr > 0.0)
    )
    if not np.any(valid):
        return cp

    cp[valid] = gamma_arr[valid] / (gamma_arr[valid] - 1.0) * p_arr[valid] / (rho_arr[valid] * T_arr[valid])
    return cp

## Load the EOS file
eos_file = "../RUNDIR/Uppsala_mergedeos_PI_A.dat"
mu_eos=eos.mu_eos(eos_file)

## Set the simulation size and range
## Iteration number

iter = 0

## Nz = vertical direction (dimension 0 in MURaM)
## Nx = first horizontal direction (set = 1 for 1D)
## Ny = second horizontal direection (set = 1 for 2D or 1D)

Nz = 4000
Nx = 2400
Ny = 1

## z_box_top = the top height above the surface in Mm
## z_box_range = the vertical extend in Mm

z_box_top = 1.0
# z_box_range = 4.0
z_box_range =4.0

x_box_range = 24.0
y_box_range = 24.0

dz = z_box_range/Nz
dx = x_box_range/Nx
dy = y_box_range/Ny

## Load the SSM Background
N_SSM = 2480
z_SSM = -1.0e-6*np.loadtxt('tables/R_SSM.dat',dtype=float).ravel()[::-1]
r_SSM = np.loadtxt('tables/rho_SSM.dat',dtype=float).ravel()[::-1]
t_SSM = np.loadtxt('tables/T_SSM.dat',dtype=float).ravel()[::-1]
p_SSM = np.loadtxt('tables/P_SSM.dat',dtype=float).ravel()[::-1]
gamma_SSM = np.loadtxt('tables/gamma_SSM.dat',dtype=float).ravel()[::-1]

## bkg values
z_bkg = z_box_top-((Nz-1)-np.arange(Nz,dtype=float))*float(z_box_range)/float(Nz)

## Interpolate for rho, pressure and density
r_bkg = np.exp(np.interp(z_bkg,z_SSM,np.log(r_SSM)))
p_bkg = np.exp(np.interp(z_bkg,z_SSM,np.log(p_SSM)))
t_bkg = np.exp(np.interp(z_bkg,z_SSM,np.log(t_SSM)))

## Above the SSM use an isothermal, constant pressure scale height extrapolation

r0 = r_SSM[-1]
p0 = p_SSM[-1]
t0 = t_SSM[-1]
z0 = z_SSM[-1]
H0 = p0/r0/G_SOLAR_REF

for kk in range(Nz):
    if(z_bkg[kk]> z_SSM[-1]):
      r_bkg[kk] = r0*np.exp((z0-z_bkg[kk])*1.0e8/H0)
      p_bkg[kk] = p0*np.exp((z0-z_bkg[kk])*1.0e8/H0)
      t_bkg[kk] = t0

## Invert temperature table for internal energy
e_bkg = invert_energy_from_rho_T(mu_eos, r_bkg, t_bkg)

s_bkg = mu_eos.interp_s(e_bkg,r_bkg)

## Bottom boundary conditions
zstag= z_bkg[0]-0.5*z_box_range/float(Nz)

rbc  = np.exp(np.interp(zstag,z_SSM,np.log(r_SSM)))
pbc  = np.exp(np.interp(zstag,z_SSM,np.log(p_SSM)))
Tbc  = np.exp(np.interp(zstag,z_SSM,np.log(t_SSM)))

## Plot computed physical quantities
zmin = z_bkg.min()
zmax = z_bkg.max()
ssm_mask = (z_SSM >= zmin) & (z_SSM <= zmax)

if np.count_nonzero(ssm_mask) < 2:
    raise ValueError("SSM and background z-ranges do not overlap sufficiently for plotting.")

z_SSM_overlap = z_SSM[ssm_mask]
r_SSM_overlap = r_SSM[ssm_mask]
p_SSM_overlap = p_SSM[ssm_mask]
t_SSM_overlap = t_SSM[ssm_mask]
gamma_SSM_overlap = gamma_SSM[ssm_mask]

n_plot = max(z_bkg.size, z_SSM_overlap.size)
z_plot = np.linspace(zmin, zmax, n_plot)

r_bkg_plot = np.exp(np.interp(z_plot, z_bkg, np.log(r_bkg)))
p_bkg_plot = np.exp(np.interp(z_plot, z_bkg, np.log(p_bkg)))
t_bkg_plot = np.exp(np.interp(z_plot, z_bkg, np.log(t_bkg)))
e_bkg_plot = np.exp(np.interp(z_plot, z_bkg, np.log(e_bkg)))
s_bkg_plot = np.interp(z_plot, z_bkg, s_bkg)

r_ssm_plot = np.exp(np.interp(z_plot, z_SSM_overlap, np.log(r_SSM_overlap)))
p_ssm_plot = np.exp(np.interp(z_plot, z_SSM_overlap, np.log(p_SSM_overlap)))
t_ssm_plot = np.exp(np.interp(z_plot, z_SSM_overlap, np.log(t_SSM_overlap)))

# Hydrostatic gravity estimate from SSM/background:
#   dP/dz = -rho * g  ->  g = -(1/rho) * dP/dz
z_ssm_g, g_ssm_hse = hydrostatic_gravity(z_SSM_overlap, p_SSM_overlap, r_SSM_overlap)
z_bkg_g, g_bkg_hse = hydrostatic_gravity(z_bkg, p_bkg, r_bkg)

g_ssm_min, g_ssm_med, g_ssm_max = finite_stats(g_ssm_hse)
g_bkg_min, g_bkg_med, g_bkg_max = finite_stats(g_bkg_hse)
g_ssm_rel_med = (g_ssm_med - G_SOLAR_REF) / G_SOLAR_REF
g_bkg_rel_med = (g_bkg_med - G_SOLAR_REF) / G_SOLAR_REF

# Additional thermodynamic diagnostics for figure-style plot
H_bkg_km = p_bkg / (r_bkg * g_bkg_hse) / 1.0e5
H_ssm_km = p_SSM_overlap / (r_SSM_overlap * g_ssm_hse) / 1.0e5

gamma_bkg = np.interp(z_bkg, z_SSM, gamma_SSM)
cp_bkg = cp_from_gamma(gamma_bkg, p_bkg, r_bkg, t_bkg)
cp_ssm = cp_from_gamma(gamma_SSM_overlap, p_SSM_overlap, r_SSM_overlap, t_SSM_overlap)

rnorm_bkg = 1.0 + z_bkg / R_SUN_MM
rnorm_ssm = 1.0 + z_SSM_overlap / R_SUN_MM
rnorm_bkg_g = 1.0 + z_bkg_g / R_SUN_MM
rnorm_ssm_g = 1.0 + z_ssm_g / R_SUN_MM

fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
axes = axes.ravel()

axes[0].semilogy(z_plot, r_bkg_plot, lw=2.0, label="Background", zorder=2, alpha=0.9)
axes[0].semilogy(z_plot, r_ssm_plot, "--", lw=1.8, label="SSM (overlap)", color="k", zorder=3)
axes[0].set_ylabel(r"$\rho$")
axes[0].set_title("Density")
axes[0].grid(alpha=0.3)
axes[0].legend()

axes[1].semilogy(z_plot, p_bkg_plot, lw=2.0, label="Background", zorder=2, alpha=0.9)
axes[1].semilogy(z_plot, p_ssm_plot, "--", lw=1.8, label="SSM (overlap)", color="k", zorder=3)
axes[1].set_ylabel("Pressure")
axes[1].set_title("Pressure")
axes[1].grid(alpha=0.3)
axes[1].legend()

axes[2].plot(z_plot, t_bkg_plot, lw=2.0, label="Background", zorder=2, alpha=0.9)
axes[2].plot(z_plot, t_ssm_plot, "--", lw=1.8, label="SSM (overlap)", color="k", zorder=3)
axes[2].set_ylabel("Temperature")
axes[2].set_title("Temperature")
axes[2].grid(alpha=0.3)
axes[2].legend()

axes[3].semilogy(z_plot, e_bkg_plot, lw=2.0, color="tab:purple")
axes[3].set_xlabel("z [Mm]")
axes[3].set_ylabel("Internal Energy")
axes[3].set_title("Internal Energy")
axes[3].grid(alpha=0.3)

axes[4].plot(z_plot, s_bkg_plot, lw=2.0, color="tab:brown")
axes[4].set_xlabel("z [Mm]")
axes[4].set_ylabel("Entropy")
axes[4].set_title("Entropy")
axes[4].grid(alpha=0.3)

for ax in axes[:5]:
    ax.set_xlim(zmin, zmax)

axes[5].axis("off")
axes[5].text(
    0.02,
    0.95,
    "\n".join(
        [
            f"N_SSM(total)   = {z_SSM.size}",
            f"N_SSM(overlap) = {z_SSM_overlap.size}",
            f"N_bkg          = {z_bkg.size}",
            f"N_plot         = {n_plot}",
            "",
            f"zstag = {zstag:.4f} Mm",
            f"rbc   = {rbc:.4e}",
            f"pbc   = {pbc:.4e}",
            f"Tbc   = {Tbc:.4e}",
            "",
            f"g_ref        = {G_SOLAR_REF:.4e} cm/s^2",
            f"g_ssm(med)   = {g_ssm_med:.4e} ({g_ssm_rel_med:+.2%})",
            f"g_bkg(med)   = {g_bkg_med:.4e} ({g_bkg_rel_med:+.2%})",
            f"g_ssm(range) = [{g_ssm_min:.4e}, {g_ssm_max:.4e}]",
            f"g_bkg(range) = [{g_bkg_min:.4e}, {g_bkg_max:.4e}]",
        ]
    ),
    va="top",
    ha="left",
    fontsize=11,
    family="monospace",
)
axes[5].set_title("Bottom Boundary")

fig.suptitle("Initial Condition Profiles", y=1.02)
fig.tight_layout()
plot_file = "initial_condition_profiles.png"
fig.savefig(plot_file, dpi=180, bbox_inches="tight")
print(f"Saved plot: {plot_file}")

# Additional diagnostic plot for hydrostatic gravity estimate
fig_g, ax_g = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax_g[0].plot(z_ssm_g, g_ssm_hse, lw=1.8, label="SSM: -(1/rho) dP/dz")
ax_g[0].plot(z_bkg_g, g_bkg_hse, lw=1.8, label="Background: -(1/rho) dP/dz")
ax_g[0].axhline(G_SOLAR_REF, color="k", ls="--", lw=1.2, label="Reference g = 2.74e4")
ax_g[0].set_ylabel("g [cm/s^2]")
ax_g[0].set_title("Hydrostatic Gravity from SSM / Background")
ax_g[0].grid(alpha=0.3)
ax_g[0].legend()

g_ssm_interp = np.interp(z_plot, z_ssm_g, g_ssm_hse)
g_bkg_interp = np.interp(z_plot, z_bkg_g, g_bkg_hse)
ax_g[1].plot(z_plot, (g_ssm_interp - G_SOLAR_REF) / G_SOLAR_REF, lw=1.6, label="SSM relative error")
ax_g[1].plot(z_plot, (g_bkg_interp - G_SOLAR_REF) / G_SOLAR_REF, lw=1.6, label="Background relative error")
ax_g[1].axhline(0.0, color="k", ls="--", lw=1.0)
ax_g[1].set_xlabel("z [Mm]")
ax_g[1].set_ylabel("(g_est - g_ref) / g_ref")
ax_g[1].grid(alpha=0.3)
ax_g[1].legend()

fig_g.tight_layout()
gravity_plot_file = "ssm_gravity_profile.png"
fig_g.savefig(gravity_plot_file, dpi=180, bbox_inches="tight")
print(f"Saved plot: {gravity_plot_file}")

# Figure style similar to the attached 2x3 reference plot
fig_ref, ax_ref = plt.subplots(2, 3, figsize=(8.2, 6.2), sharex=True)
ax_ref = ax_ref.ravel()

for ax in ax_ref:
    ax.set_xlim(0.7, 1.0)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0])
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.set_xlabel(r"r/$R_{\odot}$")

m = np.isfinite(rnorm_bkg) & np.isfinite(r_bkg) & (r_bkg > 0.0)
ax_ref[0].semilogy(rnorm_bkg[m], r_bkg[m], color="k", lw=1.2)
m = np.isfinite(rnorm_ssm) & np.isfinite(r_SSM_overlap) & (r_SSM_overlap > 0.0)
ax_ref[0].semilogy(rnorm_ssm[m], r_SSM_overlap[m], color="r", ls="--", lw=1.1)
ax_ref[0].set_ylim(1.0e-6, 1.0e0)
ax_ref[0].set_ylabel(r"$\rho_0$ [g cm$^{-3}$]")
ax_ref[0].set_title("(a)", pad=4)

m = np.isfinite(rnorm_bkg) & np.isfinite(p_bkg) & (p_bkg > 0.0)
ax_ref[1].semilogy(rnorm_bkg[m], p_bkg[m], color="k", lw=1.2)
m = np.isfinite(rnorm_ssm) & np.isfinite(p_SSM_overlap) & (p_SSM_overlap > 0.0)
ax_ref[1].semilogy(rnorm_ssm[m], p_SSM_overlap[m], color="r", ls="--", lw=1.1)
ax_ref[1].set_ylim(1.0e6, 1.0e14)
ax_ref[1].set_ylabel(r"$p_0$ [dyn cm$^{-2}$]")
ax_ref[1].set_title("(b)", pad=4)

m = np.isfinite(rnorm_bkg) & np.isfinite(t_bkg) & (t_bkg > 0.0)
ax_ref[2].semilogy(rnorm_bkg[m], t_bkg[m], color="k", lw=1.2)
m = np.isfinite(rnorm_ssm) & np.isfinite(t_SSM_overlap) & (t_SSM_overlap > 0.0)
ax_ref[2].semilogy(rnorm_ssm[m], t_SSM_overlap[m], color="r", ls="--", lw=1.1)
ax_ref[2].set_ylim(1.0e4, 1.0e7)
ax_ref[2].set_ylabel(r"$T_0$ [K]")
ax_ref[2].set_title("(c)", pad=4)

m = np.isfinite(rnorm_bkg_g) & np.isfinite(g_bkg_hse)
ax_ref[3].plot(rnorm_bkg_g[m], g_bkg_hse[m] / 1.0e4, color="k", lw=1.2)
m = np.isfinite(rnorm_ssm_g) & np.isfinite(g_ssm_hse)
ax_ref[3].plot(rnorm_ssm_g[m], g_ssm_hse[m] / 1.0e4, color="r", ls="--", lw=1.1)
ax_ref[3].set_ylim(0.0, 6.0)
ax_ref[3].set_ylabel(r"$g$ [$10^4$ cm s$^{-2}$]")
ax_ref[3].set_title("(d)", pad=4)

m = np.isfinite(rnorm_bkg_g) & np.isfinite(H_bkg_km) & (H_bkg_km > 0.0)
ax_ref[4].semilogy(rnorm_bkg_g[m], H_bkg_km[m], color="k", lw=1.2)
m = np.isfinite(rnorm_ssm_g) & np.isfinite(H_ssm_km) & (H_ssm_km > 0.0)
ax_ref[4].semilogy(rnorm_ssm_g[m], H_ssm_km[m], color="r", ls="--", lw=1.1)
ax_ref[4].set_ylim(1.0e2, 1.0e5)
ax_ref[4].set_ylabel(r"$H_{p0}$ [km]")
ax_ref[4].set_title("(e)", pad=4)

m = np.isfinite(rnorm_bkg) & np.isfinite(cp_bkg)
ax_ref[5].plot(rnorm_bkg[m], cp_bkg[m] / 1.0e9, color="k", lw=1.2)
m = np.isfinite(rnorm_ssm) & np.isfinite(cp_ssm)
ax_ref[5].plot(rnorm_ssm[m], cp_ssm[m] / 1.0e9, color="r", ls="--", lw=1.1)
ax_ref[5].set_ylim(0.0, 2.0)
ax_ref[5].set_ylabel(r"$C_{p0}$ [$10^9$ erg g$^{-1}$ K$^{-1}$]")
ax_ref[5].set_title("(f)", pad=4)

fig_ref.tight_layout()
reference_style_plot_file = "reference_style_profiles.png"
fig_ref.savefig(reference_style_plot_file, dpi=220, bbox_inches="tight")
print(f"Saved plot: {reference_style_plot_file}")

if os.environ.get("DISPLAY"):
    plt.show()
else:
    plt.close(fig)
    plt.close(fig_g)
    plt.close(fig_ref)