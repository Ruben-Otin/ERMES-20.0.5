#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eqdsk_to_scotty.py
==================

Convert an EFIT G-EQDSK equilibrium file (plus kinetic profiles) into the
input files read by the Scotty beam tracer.

HOW TO USE
----------
1. Edit the "USER PARAMETERS" block just below this docstring.
2. Run the file (python eqdsk_to_scotty.py, or Run in Spyder/VS Code).
3. In Scotty's beam_me_up(), use:
       find_B_method        = "torbeam"      ("omfit" if you use the JSON)
       density_fit_method   = "smoothing-spline-file"
       input_filename_suffix = SUFFIX         (the same value as below)
       magnetic_data_path   = OUTPUT_DIR
       ne_data_path         = OUTPUT_DIR
       Te_data_path         = OUTPUT_DIR      (only if relativistic_flag=True)

FILES PRODUCED  (<sfx> = SUFFIX)
--------------------------------
    topfile<sfx>        B_R, B_t, B_Z and psi_N on the (R, Z) grid.
                        Comes 100 % from the EQDSK.
    topfile<sfx>.json   Same data in JSON (only for find_B_method="omfit").
    ne<sfx>.dat         Electron density [1e19 m^-3] vs sqrt(psi_N).
                        NOT contained in the EQDSK: needs an external profile.
    Te<sfx>.dat         Electron temperature [keV] vs sqrt(psi_N).
                        NOT contained in the EQDSK. Scotty only reads it when
                        relativistic_flag=True.

All formats were checked against the Scotty source (torbeam.py,
beam_me_up.py, profile_fit.py).

CONVENTIONS
-----------
* psi in the EQDSK is assumed to be in Wb/rad (standard EFIT). If the
  plasma-current check printed at the end is off by ~2*pi, set
  PSI_IN_WEBER = True.
* B_R = -(1/R) dpsi/dZ,  B_Z = (1/R) dpsi/dR,  B_t = F(psi)/R.
  If your EQDSK uses the opposite poloidal sign (COCOS), set FLIP_BP = True.
* Outside the last closed flux surface F is set to its vacuum value.

Requirements: numpy, scipy, matplotlib (plots and LCFS masking).
"""

from __future__ import annotations

###############################################################################
###############################  USER PARAMETERS  #############################
###############################################################################
# Change only this block. Paths can be absolute (e.g. "C:/Users/me/Plasma/")
# or relative to the folder you run the script from.
###############################################################################

# ---- Input / output ---------------------------------------------------------

INPUT_DIR   = "Plasma/"                 # folder with the EQDSK and profile files
EQDSK_FILE  = "mast-u-sample.eqdsk"     # G-EQDSK file name (inside INPUT_DIR)
OUTPUT_DIR  = "Scotty_input/"           # where the Scotty files are written
SUFFIX      = "_sample"                 # -> topfile_test_1, ne_test_1.dat, ...
                                        # must equal Scotty's input_filename_suffix

# ---- Magnetic geometry options ----------------------------------------------

WRITE_JSON    = False    # also write topfile<SUFFIX>.json (find_B_method="omfit")
PSI_IN_WEBER  = False    # True if psi in the EQDSK is in Wb instead of Wb/rad
FLIP_BP       = False    # True to flip the sign of B_R and B_Z (other COCOS)
INTERP_ORDER  = 3        # spline order used to differentiate psi(R, Z)

# ---- Kinetic profiles (NOT in the EQDSK) ------------------------------------

# Each profile file can be:
#   * two columns  (radial coordinate, value), or
#   * one column   (values only) -> the radial coordinate is read from
#                  PROFILE_X_FILE (one column, same number of points).

# Set a file to None to skip it.
PROFILE_X_FILE = "flux_sample.dat"      # radial coordinate for 1-column files
PROFILE_X_TYPE = "psi_n"                # "psi_n" (normalised flux) or
                                        # "rho"   (rho_pol = sqrt(psi_N))

NE_FILE  = "ne_sample.dat"              # electron density file (or None)
NE_UNITS = "m-3"                        # "m-3" or "1e19" (units inside NE_FILE)

TE_FILE  = None                         # electron temperature file (or None)
TE_UNITS = "eV"                         # "eV" or "keV" (units inside TE_FILE)

# Rough fallback when no density is available: n_e = p / (e T_e (1 + Ti/Te)),
# using the EQDSK pressure, assuming n_i = n_e and no fast ions.
# Only used if NE_FILE is None; requires TE_FILE.
NE_FROM_PRESSURE = False
TI_OVER_TE       = 1.0

# ---- Diagnostics ------------------------------------------------------------

MAKE_PLOTS = True       # save check plots (PNG) in OUTPUT_DIR
SHOW_PLOTS = True       # also open the plot windows

###############################################################################
########################  END OF USER PARAMETERS  #############################
###############################################################################

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.interpolate import RectBivariateSpline, interp1d

MU0 = constants.mu_0
E_CHARGE = constants.e

# Regex that splits Fortran-formatted numbers even when they are glued
# together, e.g. "0.285254274E-02-0.599517148E+00".
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


# =============================================================================
# 1. Reading the G-EQDSK file
# =============================================================================
@dataclass
class GEqdsk:
    """Container for the contents of a G-EQDSK file (names follow the spec)."""

    header: str
    nw: int                 # number of R grid points
    nh: int                 # number of Z grid points
    rdim: float             # width of the R grid [m]
    zdim: float             # height of the Z grid [m]
    rcentr: float           # reference major radius for bcentr [m]
    rleft: float            # R of the left edge of the grid [m]
    zmid: float             # Z of the centre of the grid [m]
    rmaxis: float           # R of the magnetic axis [m]
    zmaxis: float           # Z of the magnetic axis [m]
    simag: float            # psi at the magnetic axis [Wb/rad]
    sibry: float            # psi at the plasma boundary (LCFS) [Wb/rad]
    bcentr: float           # vacuum toroidal field at rcentr [T]
    current: float          # plasma current [A]
    fpol: np.ndarray        # F = R*B_t on the uniform 1D psi grid [T m]
    pres: np.ndarray        # pressure [Pa]
    ffprim: np.ndarray      # F dF/dpsi
    pprime: np.ndarray      # dp/dpsi
    psirz: np.ndarray       # psi on the (R,Z) grid, shape (nw, nh) -> [iR, iZ]
    qpsi: np.ndarray        # safety factor
    rbbbs: np.ndarray       # R of the LCFS points [m]
    zbbbs: np.ndarray       # Z of the LCFS points [m]
    rlim: np.ndarray        # R of the limiter points [m]
    zlim: np.ndarray        # Z of the limiter points [m]

    @property
    def R(self) -> np.ndarray:
        """1D major-radius grid [m]."""
        return np.linspace(self.rleft, self.rleft + self.rdim, self.nw)

    @property
    def Z(self) -> np.ndarray:
        """1D vertical grid [m]."""
        return np.linspace(self.zmid - self.zdim / 2, self.zmid + self.zdim / 2, self.nh)

    @property
    def psi_n_1d(self) -> np.ndarray:
        """Normalised flux grid of the 1D profiles (fpol, pres, ...): 0 -> 1."""
        return np.linspace(0.0, 1.0, self.nw)


def read_geqdsk(path: str | Path) -> GEqdsk:
    """
    Read a G-EQDSK file.

    Line 1 is the header (format A48, 3I4): a free text string followed by
    three integers, the last two being nw and nh. Everything afterwards is a
    stream of numbers (5 per line, 16 characters each), so instead of reading
    line by line we tokenise the whole remainder and consume it in order.
    This is robust against glued negative numbers and short last lines.
    """
    with open(path, "r") as f:
        first_line = f.readline()
        body = f.read()

    # --- header: the last two integers on the first line are nw, nh --------
    ints = re.findall(r"-?\d+", first_line[48:]) or re.findall(r"-?\d+", first_line)
    nw, nh = int(ints[-2]), int(ints[-1])
    header = first_line[:48].strip()

    # --- numeric body as one flat stream -------------------------------------
    tokens = iter(float(t.replace("D", "E").replace("d", "e"))
                  for t in _FLOAT_RE.findall(body))

    def take(n: int) -> np.ndarray:
        return np.array([next(tokens) for _ in range(n)])

    # 20 scalars (4 lines of 5). Some entries are duplicates or dummies.
    (rdim, zdim, rcentr, rleft, zmid,
     rmaxis, zmaxis, simag, sibry, bcentr,
     current, _simag, _xdum, _rmaxis, _xdum,
     _zmaxis, _xdum, _sibry, _xdum, _xdum) = take(20)

    fpol = take(nw)
    pres = take(nw)
    ffprim = take(nw)
    pprime = take(nw)
    # psirz is stored with R running fastest: nh rows of nw values -> [iZ, iR]
    psirz = take(nw * nh).reshape(nh, nw).T          # -> [iR, iZ]
    qpsi = take(nw)

    # Boundary and limiter (may be absent in some files)
    try:
        nbbbs, limitr = int(next(tokens)), int(next(tokens))
        bbbs = take(2 * nbbbs).reshape(nbbbs, 2)      # pairs (R, Z)
        lim = take(2 * limitr).reshape(limitr, 2)
    except StopIteration:
        bbbs = np.zeros((0, 2))
        lim = np.zeros((0, 2))

    return GEqdsk(header, nw, nh, rdim, zdim, rcentr, rleft, zmid,
                  rmaxis, zmaxis, simag, sibry, bcentr, current,
                  fpol, pres, ffprim, pprime, psirz, qpsi,
                  bbbs[:, 0], bbbs[:, 1], lim[:, 0], lim[:, 1])


# =============================================================================
# 2. Magnetic field on the (R, Z) grid
# =============================================================================
def normalised_flux(eq: GEqdsk) -> np.ndarray:
    """psi_N = (psi - psi_axis) / (psi_LCFS - psi_axis): 0 on axis, 1 on LCFS.

    Works for either sign convention of psi (no assumption on which of
    simag/sibry is larger, unlike the min/max rescaling of the old script).
    """
    return (eq.psirz - eq.simag) / (eq.sibry - eq.simag)


def inside_lcfs_mask(eq: GEqdsk, psi_n: np.ndarray) -> np.ndarray:
    """Boolean mask [iR, iZ] of grid points inside the plasma.

    Uses the LCFS polygon when available (this excludes the private-flux
    region below an X-point, where psi_N < 1 but there is no plasma).
    Falls back to psi_N <= 1 otherwise.
    """
    RR, ZZ = np.meshgrid(eq.R, eq.Z, indexing="ij")
    mask = psi_n <= 1.0
    if eq.rbbbs.size >= 3:
        try:
            from matplotlib.path import Path as MplPath
            poly = MplPath(np.column_stack([eq.rbbbs, eq.zbbbs]))
            in_poly = poly.contains_points(np.column_stack([RR.ravel(), ZZ.ravel()]))
            mask &= in_poly.reshape(RR.shape)
        except ImportError:
            print("[warn] matplotlib not found: using psi_N <= 1 as plasma mask")
    return mask


def compute_fields(eq: GEqdsk, psi_in_weber: bool = False, flip_bp: bool = False,
                   interp_order: int = 3):
    """Return B_R, B_t, B_Z and psi_N on the EQDSK grid, all shaped [iR, iZ]."""
    R, Z = eq.R, eq.Z
    psi = eq.psirz / (2 * np.pi) if psi_in_weber else eq.psirz

    # Spline of psi(R, Z) to get smooth derivatives
    psi_spline = RectBivariateSpline(R, Z, psi, kx=interp_order, ky=interp_order, s=0)
    dpsi_dR = psi_spline(R, Z, dx=1, dy=0)
    dpsi_dZ = psi_spline(R, Z, dx=0, dy=1)

    Rcol = R[:, None]                       # broadcast R along Z
    sign = -1.0 if flip_bp else 1.0
    B_R = sign * (-dpsi_dZ / Rcol)
    B_Z = sign * (dpsi_dR / Rcol)

    # Toroidal field: B_t = F(psi_N) / R, with F = fpol inside the plasma
    # and F = F_vacuum = fpol[-1] outside (no spline extrapolation).
    psi_n = normalised_flux(eq)
    F_of_psi = interp1d(eq.psi_n_1d, eq.fpol, kind="cubic")
    F = np.full_like(psi_n, eq.fpol[-1])
    inside = inside_lcfs_mask(eq, psi_n)
    F[inside] = F_of_psi(np.clip(psi_n[inside], 0.0, 1.0))
    B_t = F / Rcol

    return B_R, B_t, B_Z, psi_n


# =============================================================================
# 3. Sanity checks (printed, not fatal)
# =============================================================================
def run_checks(eq: GEqdsk, B_R, B_t, B_Z, psi_n) -> None:
    print("\n--- Sanity checks -------------------------------------------------")
    print(f"Grid: {eq.nw} x {eq.nh}, R = [{eq.R[0]:.3f}, {eq.R[-1]:.3f}] m, "
          f"Z = [{eq.Z[0]:.3f}, {eq.Z[-1]:.3f}] m")
    print(f"Magnetic axis: R = {eq.rmaxis:.4f} m, Z = {eq.zmaxis:.4f} m")

    # (a) Vacuum toroidal field: F_vac should equal rcentr * bcentr
    print(f"F_vac = fpol[-1] = {eq.fpol[-1]:.5f} T m   |   "
          f"rcentr*bcentr = {eq.rcentr * eq.bcentr:.5f} T m")

    # (b) psi_N at the grid point nearest the axis should be ~0
    iR = np.argmin(abs(eq.R - eq.rmaxis))
    iZ = np.argmin(abs(eq.Z - eq.zmaxis))
    print(f"psi_N at grid point nearest the axis = {psi_n[iR, iZ]:.4f} (expect ~0)")

    # (c) Ampere's law around the LCFS: |closed integral B_p.dl| = mu0*|Ip|.
    #     Checks units of psi (Wb/rad vs Wb) and the B_p construction.
    if eq.rbbbs.size >= 3:
        sR = RectBivariateSpline(eq.R, eq.Z, B_R)
        sZ = RectBivariateSpline(eq.R, eq.Z, B_Z)
        Rb = np.append(eq.rbbbs, eq.rbbbs[0])
        Zb = np.append(eq.zbbbs, eq.zbbbs[0])
        Rm, Zm = 0.5 * (Rb[1:] + Rb[:-1]), 0.5 * (Zb[1:] + Zb[:-1])
        circ = np.sum(sR(Rm, Zm, grid=False) * np.diff(Rb)
                      + sZ(Rm, Zm, grid=False) * np.diff(Zb))
        Ip_calc = circ / MU0
        ratio = abs(Ip_calc / eq.current) if eq.current else np.nan
        print(f"Ip from Ampere's law = {Ip_calc:.4e} A  |  Ip in header = "
              f"{eq.current:.4e} A  |  |ratio| = {ratio:.3f}")
        if np.isfinite(ratio):
            if abs(ratio - 1) < 0.1:
                print("  OK: magnitudes agree.")
            elif abs(ratio * 2 * np.pi - 1) < 0.1:
                print("  ratio ~ 1/(2 pi): psi looks like Wb/rad was already applied "
                      "twice? Check the file units.")
            elif abs(ratio / (2 * np.pi) - 1) < 0.1:
                print("  ratio ~ 2 pi: psi is probably in Wb -> rerun with --psi-in-weber")
            else:
                print("  WARNING: no simple factor explains the mismatch; check the "
                      "header current and the boundary data.")
    print("-------------------------------------------------------------------\n")


# =============================================================================
# 4. Writers
# =============================================================================
def write_topfile(path: Path, eq: GEqdsk, B_R, B_t, B_Z, psi_n) -> None:
    """
    Write a TORBEAM 'topfile'.

    Scotty's reader (scotty.torbeam.Torbeam.from_file) skips everything until
    a line containing 'X-coordinates', then reads blocks separated by lines
    containing 'Z-coordinates', 'B_R', 'B_t', 'B_Z' and 'psi'. Each 2D block
    is read as a flat list reshaped to (nZ, nR), i.e. R runs fastest, which
    is why the [iR, iZ] arrays are transposed before writing.

    The four header lines follow the original TORBEAM layout (grid size and
    psi at the separatrix = 1 since psi is normalised), so the same file can
    also be fed to TORBEAM itself.
    """
    def block(f, label, values):
        f.write(label + "\n")
        np.savetxt(f, np.ravel(values), fmt="%.10e")

    with open(path, "w") as f:
        f.write("Number of radial and vertical grid points\n")
        f.write(f"{eq.nw} {eq.nh}\n")
        f.write("Inside and outside radius and psi_sep\n")
        f.write(f"{eq.R[0]:.10e} {eq.R[-1]:.10e} 1.0\n")
        block(f, "Grid: X-coordinates", eq.R)
        block(f, "Grid: Z-coordinates", eq.Z)
        block(f, "Magnetic field: B_R", B_R.T)
        block(f, "Magnetic field: B_t", B_t.T)
        block(f, "Magnetic field: B_Z", B_Z.T)
        block(f, "Poloidal flux: psi", psi_n.T)
    print(f"Wrote {path}")


def write_topfile_json(path: Path, eq: GEqdsk, B_R, B_t, B_Z, psi_n) -> None:
    """OMFIT-style JSON topfile (find_B_method='omfit'). 2D data flattened
    with R running fastest, same ordering as the text topfile."""
    data = {
        "R": eq.R.tolist(),
        "Z": eq.Z.tolist(),
        "Br": B_R.T.ravel().tolist(),
        "Bt": B_t.T.ravel().tolist(),
        "Bz": B_Z.T.ravel().tolist(),
        "pol_flux": psi_n.T.ravel().tolist(),
    }
    with open(path, "w") as f:
        json.dump(data, f)
    print(f"Wrote {path}")


def write_profile(path: Path, psi_n: np.ndarray, values: np.ndarray, what: str) -> None:
    """Write a Scotty profile file: N, then rows 'sqrt(psi_N)  value'."""
    rho = np.sqrt(psi_n)
    with open(path, "w") as f:
        f.write(f"{len(rho)}\n")
        np.savetxt(f, np.column_stack([rho, values]), fmt="%.8e")
    print(f"Wrote {path} ({what}, {len(rho)} points, "
          f"rho_pol = {rho[0]:.3f} -> {rho[-1]:.3f})")


# =============================================================================
# 5. Kinetic profiles (external data)
# =============================================================================
def load_profile(value_file: str, x_file: str | None, x_type: str):
    """
    Load a 1D profile and return (psi_N, values), sorted and de-duplicated.

    value_file : either one column (values) -> x comes from x_file,
                 or two columns (x, value).
    x_type     : 'psi_n' if x is normalised poloidal flux, 'rho' if x is
                 sqrt(psi_N).
    """
    data = np.loadtxt(value_file)
    if data.ndim == 2 and data.shape[1] >= 2:
        x, y = data[:, 0], data[:, 1]
    else:
        if x_file is None:
            raise ValueError(f"{value_file} has a single column: give --profile-x")
        x, y = np.loadtxt(x_file), np.ravel(data)
        if len(x) != len(y):
            raise ValueError(f"{x_file} and {value_file} have different lengths")

    psi_n = x**2 if x_type == "rho" else x
    psi_n, idx = np.unique(psi_n, return_index=True)   # sorts + removes duplicates
    return psi_n, y[idx]


def ne_from_pressure(eq: GEqdsk, psi_te: np.ndarray, te_ev: np.ndarray,
                     ti_over_te: float):
    """
    Rough n_e estimate from the EQDSK pressure: p = n_e e T_e + n_i e T_i with
    n_i = n_e (pure hydrogenic plasma, Z_eff = 1, no fast ions):
        n_e = p / (e T_e (1 + T_i/T_e)).
    Only valid inside the LCFS (psi_N <= 1). Note that in many EFITs p
    includes fast-ion pressure, which makes this an overestimate.
    """
    psi = eq.psi_n_1d
    te_on_grid = np.interp(psi, psi_te, te_ev)
    ne_m3 = eq.pres / (E_CHARGE * te_on_grid * (1.0 + ti_over_te))
    return psi, ne_m3


# =============================================================================
# 6. Optional plots
# =============================================================================
def make_plots(eq: GEqdsk, B_R, B_t, B_Z, psi_n, out_dir: Path, sfx: str,
               profiles: dict) -> None:
    import matplotlib.pyplot as plt

    RR, ZZ = np.meshgrid(eq.R, eq.Z, indexing="ij")
    fig, axs = plt.subplots(1, 4, figsize=(16, 6), constrained_layout=True)
    for ax, field, title in zip(axs, [psi_n, B_R, B_t, B_Z],
                                [r"$\psi_N$", r"$B_R$ [T]", r"$B_t$ [T]", r"$B_Z$ [T]"]):
        cs = ax.contourf(RR, ZZ, field, 60, cmap="RdBu_r")
        ax.contour(RR, ZZ, psi_n, levels=[1.0], colors="k", linewidths=1.5)
        if eq.rbbbs.size:
            ax.plot(eq.rbbbs, eq.zbbbs, "g--", lw=1)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("R [m]")
        fig.colorbar(cs, ax=ax, shrink=0.7)
    axs[0].set_ylabel("Z [m]")
    fig.savefig(out_dir / f"topfile{sfx}_check.png", dpi=150)

    if profiles:
        fig, axs = plt.subplots(1, len(profiles), figsize=(5 * len(profiles), 4),
                                constrained_layout=True, squeeze=False)
        for ax, (name, (psi, val, unit)) in zip(axs[0], profiles.items()):
            ax.plot(np.sqrt(psi), val)
            ax.axvline(1.0, color="k", ls=":")
            ax.set_xlabel(r"$\rho_{pol} = \sqrt{\psi_N}$")
            ax.set_ylabel(f"{name} [{unit}]")
        fig.savefig(out_dir / f"profiles{sfx}_check.png", dpi=150)
    print(f"Saved check plots in {out_dir}")


# =============================================================================
# 7. Main: runs with the USER PARAMETERS defined at the top of the file
# =============================================================================
def main():
    in_dir = Path(INPUT_DIR)
    out_dir = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    def in_path(name):
        """Full path of an input file, or None if the parameter is None."""
        return None if name is None else str(in_dir / name)

    # --- 1. Equilibrium -> topfile ------------------------------------------
    eq = read_geqdsk(in_dir / EQDSK_FILE)
    print(f"Read '{EQDSK_FILE}' ({eq.header}), grid {eq.nw} x {eq.nh}, "
          f"{eq.rbbbs.size} LCFS points")
    B_R, B_t, B_Z, psi_n = compute_fields(eq, PSI_IN_WEBER, FLIP_BP, INTERP_ORDER)
    run_checks(eq, B_R, B_t, B_Z, psi_n)
    write_topfile(out_dir / f"topfile{SUFFIX}", eq, B_R, B_t, B_Z, psi_n)
    if WRITE_JSON:
        write_topfile_json(out_dir / f"topfile{SUFFIX}.json", eq, B_R, B_t, B_Z, psi_n)

    # --- 2. Electron temperature -> Te<SUFFIX>.dat (keV) ---------------------
    plotted = {}
    te_data = None
    if TE_FILE is not None:
        psi_te, te = load_profile(in_path(TE_FILE), in_path(PROFILE_X_FILE), PROFILE_X_TYPE)
        te_ev = te if TE_UNITS == "eV" else te * 1e3
        te_data = (psi_te, te_ev)
        write_profile(out_dir / f"Te{SUFFIX}.dat", psi_te, te_ev / 1e3, "T_e [keV]")
        plotted["T_e"] = (psi_te, te_ev / 1e3, "keV")

    # --- 3. Electron density -> ne<SUFFIX>.dat (1e19 m^-3) -------------------
    if NE_FILE is not None:
        psi_ne, ne = load_profile(in_path(NE_FILE), in_path(PROFILE_X_FILE), PROFILE_X_TYPE)
        ne_1e19 = ne / 1e19 if NE_UNITS == "m-3" else ne
        write_profile(out_dir / f"ne{SUFFIX}.dat", psi_ne, ne_1e19, "n_e [1e19 m^-3]")
        plotted["n_e"] = (psi_ne, ne_1e19, r"$10^{19}$ m$^{-3}$")
    elif NE_FROM_PRESSURE:
        if te_data is None:
            raise SystemExit("NE_FROM_PRESSURE = True needs a temperature file (TE_FILE)")
        psi_ne, ne_m3 = ne_from_pressure(eq, *te_data, TI_OVER_TE)
        write_profile(out_dir / f"ne{SUFFIX}.dat", psi_ne, ne_m3 / 1e19,
                      "n_e [1e19 m^-3], ESTIMATED from pressure")
        plotted["n_e (from p)"] = (psi_ne, ne_m3 / 1e19, r"$10^{19}$ m$^{-3}$")
    else:
        print("[info] No density given: ne.dat not written (it cannot be "
              "derived from the EQDSK alone).")

    # --- 4. Check plots -------------------------------------------------------
    if MAKE_PLOTS:
        make_plots(eq, B_R, B_t, B_Z, psi_n, out_dir, SUFFIX, plotted)
        if SHOW_PLOTS:
            import matplotlib.pyplot as plt
            plt.show()


if __name__ == "__main__":
    main()
