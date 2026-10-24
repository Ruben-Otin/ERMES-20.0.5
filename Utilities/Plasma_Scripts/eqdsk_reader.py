#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eqdsk_reader.py
===============
Reads a G-EQDSK equilibrium file and the kinetic profiles and plots them, to
check the inputs before generating Scotty or ERMES files:

    1. poloidal flux psi and normalised flux psi_N (with the LCFS),
    2. poloidal field |B_p| (with field lines) and toroidal field B_t,
    3. poloidal current density |J_p| (with current lines) and J_phi,
    4. electron density and temperature on the (R, Z) plane and vs rho_pol.

The EQDSK is read and B is computed with eqdsk_to_scotty.py (read_geqdsk,
compute_fields), i.e. exactly as in the Scotty and ERMES files:
    B_R = -(1/R) dpsi/dZ,  B_Z = (1/R) dpsi/dR,  B_t = F/R (vacuum F outside the LCFS)
    J_R = -(1/(mu0 R)) dF/dZ, J_Z = (1/(mu0 R)) dF/dR       (F = R B_t)
    J_phi = R p'(psi) + F F'(psi) / (mu0 R)                  (inside the LCFS)

Run:  python eqdsk_reader.py   (edit the file names below)
"""

###############################################################################
## File names
###############################################################################

EQDSK_FILE   = "./EQDSK_Samples/mast-u-sample.eqdsk"
NORM_RP_FILE = "./EQDSK_Samples/rho_flux.txt"   # rho_pol = sqrt(psi_N), 1 column
NORM_NE_FILE = "./EQDSK_Samples/rho_ne.txt"     # n_e [m^-3], 1 column
NORM_TE_FILE = "./EQDSK_Samples/rho_te.txt"     # T_e [eV],   1 column

###############################################################################

import numpy as np
import matplotlib.pyplot as plt

from scipy.constants   import mu_0
from scipy.interpolate import RectBivariateSpline, interp1d

import eqdsk_to_scotty as e2s


def plot_map(ax, R, Z, field, title, eq, levels=40, cmap="turbo", lines=None):
    """Filled contour of field[iR, iZ] with the LCFS; optional (U, V) stream lines."""
    cs = ax.contourf(R, Z, field.T, levels=levels, cmap=cmap)
    plt.colorbar(cs, ax=ax, shrink=0.8)
    if lines is not None:
        ax.streamplot(R, Z, lines[0].T, lines[1].T, color="white", density=1.2,
                      linewidth=0.6)
    if eq.rbbbs.size:
        ax.plot(eq.rbbbs, eq.zbbbs, "w--", lw=1)
    ax.set_aspect("equal")
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_title(title)


def main():
    eq = e2s.read_geqdsk(EQDSK_FILE)
    R, Z = eq.R, eq.Z
    B_R, B_t, B_Z, psi_n = e2s.compute_fields(eq)
    inside = e2s.inside_lcfs_mask(eq, psi_n)
    Rcol = R[:, None]

    print(f"{EQDSK_FILE}: {eq.nw} x {eq.nh} grid, axis (R, Z) = "
          f"({eq.rmaxis:.4f}, {eq.zmaxis:.4f}) m, B0 = {eq.bcentr:.3f} T "
          f"at R0 = {eq.rcentr:.3f} m, Ip (header) = {eq.current:.4e} A")

    # --- Current density -----------------------------------------------------
    F = B_t * Rcol                                     # F = R B_t [T m]
    F_spl = RectBivariateSpline(R, Z, F)
    J_R = -F_spl(R, Z, dx=0, dy=1) / (mu_0 * Rcol)
    J_Z =  F_spl(R, Z, dx=1, dy=0) / (mu_0 * Rcol)

    pprime = interp1d(eq.psi_n_1d, eq.pprime, kind="cubic")
    ffprim = interp1d(eq.psi_n_1d, eq.ffprim, kind="cubic")
    psi_c = np.clip(psi_n, 0.0, 1.0)
    J_phi = np.where(inside, Rcol * pprime(psi_c) + ffprim(psi_c) / (mu_0 * Rcol), 0.0)

    # --- Kinetic profiles on the grid (last value beyond the profile) -------
    rho = np.loadtxt(NORM_RP_FILE)
    ne = np.loadtxt(NORM_NE_FILE)
    te = np.loadtxt(NORM_TE_FILE)
    rho_grid = np.sqrt(np.clip(psi_n, 0.0, None))
    ne_2D = np.interp(rho_grid, rho, ne)
    te_2D = np.interp(rho_grid, rho, te)

    # --- Plots ---------------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(10, 7), constrained_layout=True)
    plot_map(axs[0], R, Z, eq.psirz, r"$\Psi$ (Wb/rad)", eq)
    plot_map(axs[1], R, Z, psi_n, r"$\Psi_N$", eq)
    axs[1].contour(R, Z, psi_n.T, levels=[1.0], colors="k")

    fig, axs = plt.subplots(1, 2, figsize=(10, 7), constrained_layout=True)
    plot_map(axs[0], R, Z, np.hypot(B_R, B_Z), r"$|B_p|$ (T)", eq, cmap="plasma",
             lines=(B_R, B_Z))
    plot_map(axs[1], R, Z, B_t, r"$B_t$ (T)", eq, cmap="RdBu_r")

    fig, axs = plt.subplots(1, 2, figsize=(10, 7), constrained_layout=True)
    plot_map(axs[0], R, Z, np.hypot(J_R, J_Z), r"$|J_p|$ (A/m$^2$)", eq, cmap="plasma",
             lines=(J_R, J_Z))
    plot_map(axs[1], R, Z, J_phi, r"$J_\phi$ (A/m$^2$)", eq, cmap="RdBu_r")

    fig, axs = plt.subplots(1, 2, figsize=(10, 7), constrained_layout=True)
    plot_map(axs[0], R, Z, ne_2D, r"$n_e$ (m$^{-3}$)", eq)
    plot_map(axs[1], R, Z, te_2D, r"$T_e$ (eV)", eq)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, y, lab in zip(axs, (ne, te), (r"$n_e$ (m$^{-3}$)", r"$T_e$ (eV)")):
        ax.plot(rho, y)
        ax.axvline(1.0, color="k", ls=":")
        ax.set_xlabel(r"$\rho_{pol} = \sqrt{\Psi_N}$")
        ax.set_ylabel(lab)
        ax.grid()

    plt.show()


if __name__ == "__main__":
    main()
