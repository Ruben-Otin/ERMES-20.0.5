#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eqdsk_to_ermes_turb.py
======================
Convert an EFIT G-EQDSK equilibrium file (plus an electron density profile)
into the plasma files read by ERMES: the electron density and the external
magnetic field at every node of the GiD mesh, with an optional synthetic
turbulent density field and blobs/holes. Companion of eqdsk_to_scotty.py:
both use the SAME equilibrium calculation, so the Scotty beam and the ERMES
plasma agree.

This is a thin, problem_generator-friendly front end of torbeam_to_ermes.py
(the engine: mesh reading, fields, turbulence, blobs, check plot). It also
provides the equilibrium functions used by problem_generator.py:
    Generate_Plasma_Files       plasma files at the mesh nodes
    Get_Point_Equilibrium_Data  ne and B at a point
    Get_Pcrtc_Equilibrium_Data  ne and B where the beam reaches the critical
                                density (used for the launch polarisation)
(same signatures as the old functions of utilities.py, plus keyword options).

HOW TO USE
----------
1. Mesh the problem in GiD so that <PROBLEM_NAME>.gid/<PROBLEM_NAME>-1.dat exists.
2a. Standalone: edit the USER PARAMETERS below and run
       python eqdsk_to_ermes_turb.py
2b. From problem_generator.py (what it does):
       import eqdsk_to_ermes_turb as e2e
       e2e.Generate_Plasma_Files(EQDSK_File, Norm_ne_File, Norm_rp_File,
                                 Absolute_Path, Problem_Name, Plasma_Folder,
                                 Plasma_NameId, **Plasma_Options)
    Plasma_Options may contain toroidal_convention="right-handed"/"legacy"
    and any parameter of torbeam_to_ermes.py, e.g.
       dict(turbulence=True, turb_Lr=0.01, blobs=True, blob_amplitude=1.5)
    The USER PARAMETERS below are the defaults for anything not given.

FILES PRODUCED  (<name> = PLASMA_NAME, sorted by node id)
--------------------------------------------------------
    ne_<name>.dat          node  n_e [m^-3]
    Be_<name>.dat          node  Bx By Bz [T] (Cartesian, mesh frame)
    dn_<name>.dat          node  delta n / n        (turbulence / blobs only)
    check_<name>.png       cut through the mesh + profiles

CONVENTIONS
-----------
* B_R = -(1/R) dpsi/dZ, B_Z = (1/R) dpsi/dR, B_t = F(psi)/R, F = vacuum value
  outside the LCFS; psi_N = (psi - psi_axis)/(psi_LCFS - psi_axis)
  (eqdsk_to_scotty.compute_fields; PSI_IN_WEBER / FLIP_BP as in that script).
* Mesh frame: the tokamak axis is TOKAMAK_AXIS = "z" (R = sqrt(x^2 + y^2),
  phi = atan2(y, x)), as utilities.py and problem_generator.py.
* TOROIDAL_CONVENTION:
    "right-handed": (R, phi, Z) right-handed (phi counter-clockwise seen from
                    +Z), as the EQDSK, Scotty and utilities.Transform_FRZT_2_FXYZ.
    "legacy"      : toroidal field reversed (mirrored field-line helicity), as
                    the old y-axis version of utilities.py produced. Only to
                    reproduce old runs (see README.md).
* Nodes with psi_N < 1 outside the LCFS polygon (private-flux region below an
  X-point) are treated as outside the plasma.
* The density profile is used where it is defined (rho_flux.txt extends to
  rho = 1.3, into the SOL); beyond its last point NE_OUTSIDE decides:
  "constant" (its last value) or "zero".
* Point coordinates (P0 below) and the mesh are in MESH_UNITS.

TURBULENCE AND BLOBS
--------------------
Field-aligned random density field (k_par ~ 0, radial / binormal correlation
lengths of ~cm, rms delta n / n growing towards the edge, lognormal or Gaussian
statistics, reproducible seed) and optional field-aligned blobs / holes.
Model and parameters: docstring of torbeam_to_ermes.py.

Requirements: numpy, scipy, matplotlib; eqdsk_to_scotty.py,
torbeam_to_ermes.py and utilities.py in the same folder.
"""
from __future__ import annotations

###############################################################################
###############################  USER PARAMETERS  #############################
###############################################################################

# ---- Equilibrium and profile ------------------------------------------------
EQDSK_FILE    = "./EQDSK_Samples/mast-u-sample.eqdsk"   # G-EQDSK file
PSI_IN_WEBER  = False      # True if psi in the EQDSK is in Wb instead of Wb/rad
FLIP_BP       = False      # True to flip the sign of B_R and B_Z (other COCOS)

NE_FILE       = "./EQDSK_Samples/rho_ne.txt"    # density: 1 column (with NE_X_FILE)
                                                # or 2 columns (x, ne)
NE_X_FILE     = "./EQDSK_Samples/rho_flux.txt"  # radial coordinate of a 1-column NE_FILE
NE_X_TYPE     = "rho"      # "rho" (sqrt(psi_N)) or "psi_n"
NE_UNITS      = "m-3"      # "m-3" or "1e19" (units inside NE_FILE)
NE_OUTSIDE    = "constant" # beyond the last profile point: "constant" (last value)
                           # or "zero"

# ---- ERMES problem (mesh) ---------------------------------------------------
PROBLEM_NAME   = "MASTU-2D-t20_p20"   # GiD problem: <PROBLEM_FOLDER><NAME>.gid/<NAME>-1.dat
PROBLEM_FOLDER = "./"
MESH_UNITS     = "m"       # "m", "cm" or "mm": units of the GiD mesh

TOKAMAK_AXIS        = "z"             # Cartesian axis along the tokamak axis
TOROIDAL_CONVENTION = "right-handed"  # "right-handed" | "legacy" (see CONVENTIONS)

# ---- Output -----------------------------------------------------------------
PLASMA_NAME   = "2D-t20_p20"          # -> ne_<NAME>.dat, Be_<NAME>.dat
PLASMA_FOLDER = "./Plasma/"
MAKE_PLOT     = True                  # check_<NAME>.png in PLASMA_FOLDER

# ---- Turbulence -------------------------------------------------------------
TURBULENCE       = False   # True: add the synthetic turbulent density field
TURB_AMPLITUDE   = [(0.0, 0.003), (0.7, 0.005), (0.9, 0.02), (0.97, 0.05),
                    (1.0, 0.10), (1.2, 0.30)]     # rms delta n / n vs rho_pol
TURB_LR          = 0.015   # radial correlation length [m]
TURB_LB          = 0.030   # binormal correlation length [m]
TURB_LPAR        = 5.0     # parallel correlation length [m] (0 = k_par = 0)
TURB_KB_PEAK     = 0.0     # binormal spectral peak [1/m] (0 = none)
TURB_MODES       = 256     # random Fourier modes
TURB_STATISTICS  = "lognormal"   # "lognormal" | "gaussian"
TURB_SEED        = 12345
WRITE_DN_FILE    = True    # dn_<NAME>.dat with delta n / n

# ---- Optional exponential SOL (if the profile ends at rho = 1) -------------
SOL_NE_SEPARATRIX = 0.0    # ne at the separatrix [m^-3]; 0 = off
SOL_DECAY_LENGTH  = 0.01   # e-folding length outside the separatrix [m]

###############################################################################
########################  END OF USER PARAMETERS  #############################
###############################################################################

import io
import os
import sys
import contextlib

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import torbeam_to_ermes as P           # noqa: E402  (engine: mesh, fields, turbulence)
from scipy.constants import e as QE, m_e as ME, epsilon_0 as EPS0  # noqa: E402


def build_config(eqdsk_file=None, ne_file=None, rp_file=None, mesh_file=None,
                 out_dir=None, plasma_name=None, **overrides):
    """torbeam_to_ermes configuration from the USER PARAMETERS (+ overrides).

    overrides: any torbeam_to_ermes.py parameter, or
    toroidal_convention = "right-handed" | "legacy"."""
    cfg = P.default_parameters()
    name = plasma_name or PLASMA_NAME
    settings = dict(
        equilibrium="eqdsk",
        eqdsk_file=eqdsk_file or EQDSK_FILE,
        eqdsk_psi_in_weber=PSI_IN_WEBER, eqdsk_flip_bp=FLIP_BP, eqdsk_use_lcfs=True,
        ne_file=ne_file or NE_FILE,
        ne_x_file=rp_file if rp_file is not None else (NE_X_FILE or ""),
        ne_radial_coordinate="rho" if NE_X_TYPE == "rho" else "psi",
        ne_units=1.0 if NE_UNITS == "m-3" else 1e19,
        ne_outside=NE_OUTSIDE,
        mesh_file=mesh_file or os.path.join(PROBLEM_FOLDER, PROBLEM_NAME + ".gid",
                                            PROBLEM_NAME + "-1.dat"),
        mesh_units=MESH_UNITS,
        geometry="toroidal", tokamak_axis=TOKAMAK_AXIS,
        flip_bt=(TOROIDAL_CONVENTION == "legacy"),
        out_dir=out_dir or PLASMA_FOLDER,
        ne_out=f"ne_{name}.dat", mag_out=f"Be_{name}.dat", dn_out=f"dn_{name}.dat",
        write_turbulence_file=WRITE_DN_FILE, plot=MAKE_PLOT,
        turbulence=TURBULENCE, turb_amplitude_profile=list(TURB_AMPLITUDE),
        turb_Lr=TURB_LR, turb_Lb=TURB_LB, turb_Lpar=TURB_LPAR, turb_kb_peak=TURB_KB_PEAK,
        turb_modes=TURB_MODES, turb_statistics=TURB_STATISTICS, turb_seed=TURB_SEED,
        sol_ne_separatrix=SOL_NE_SEPARATRIX, sol_decay_length=SOL_DECAY_LENGTH,
    )
    if "toroidal_convention" in overrides:
        conv = overrides.pop("toroidal_convention")
        if conv not in ("right-handed", "legacy"):
            raise ValueError(f"toroidal_convention = '{conv}': use 'right-handed' or 'legacy'")
        settings["flip_bt"] = conv == "legacy"
    settings.update(overrides)
    for k, v in settings.items():
        if not hasattr(cfg, k):
            raise TypeError(f"unknown parameter '{k}'")
        setattr(cfg, k, v)
    cfg.plot_file = f"check_{name}.png"
    return cfg


def _run(cfg):
    """torbeam_to_ermes.main with the check plot renamed."""
    P.main(cfg)
    src = os.path.join(cfg.out_dir, "plasma_check.png")
    if cfg.plot and os.path.isfile(src):
        os.replace(src, os.path.join(cfg.out_dir, cfg.plot_file))
        print(f"Check plot: {os.path.join(cfg.out_dir, cfg.plot_file)}")


# ##########################################################################
# Functions used by problem_generator.py
# ##########################################################################

def Generate_Plasma_Files(EQDSK_File, Norm_ne_File, Norm_rp_File,
                          Problem_Folder, Problem_Name,
                          Plasma_Folder, Plasma_Name, **options):
    """Write ne_<Plasma_Name>.dat, Be_<Plasma_Name>.dat (and dn_...) in
    Plasma_Folder for the nodes of <Problem_Folder><Problem_Name>.gid/
    <Problem_Name>-1.dat.
    options: torbeam_to_ermes.py parameters (e.g. turbulence=True) or
    toroidal_convention="legacy" / "right-handed"."""
    mesh = os.path.join(Problem_Folder, Problem_Name + ".gid", Problem_Name + "-1.dat")
    if not os.path.exists(mesh):
        raise FileNotFoundError(f"Mesh nodes file {mesh} not found: mesh the problem "
                                "in GiD first (Files > WriteCalcFile).")
    os.makedirs(Plasma_Folder, exist_ok=True)
    cfg = build_config(EQDSK_File, Norm_ne_File, Norm_rp_File, mesh, Plasma_Folder,
                       Plasma_Name, **options)
    _run(cfg)


_CACHE = {}


def _equilibrium(EQDSK_File, Norm_ne_File, Norm_rp_File, options):
    """Equilibrium read once per set of files/options (cached)."""
    key = (EQDSK_File, Norm_ne_File, Norm_rp_File, repr(sorted(options.items())))
    if key not in _CACHE:
        cfg = build_config(EQDSK_File, Norm_ne_File, Norm_rp_File, mesh_file="-",
                           out_dir="-", **{**options, "turbulence": False, "plot": False})
        with contextlib.redirect_stdout(io.StringIO()):
            _CACHE[key] = (cfg, P.load_equilibrium(cfg))
    return _CACHE[key]


def Get_Point_Equilibrium_Data(P0, EQDSK_File, Norm_ne_File, Norm_rp_File, **options):
    """ne [m^-3] and B [T] (Cartesian, mesh frame) at the point P0 (in
    MESH_UNITS; an (n, 3) array gives n points), with the same equilibrium
    and convention as the plasma files (no turbulence).
    Returns (ne array, B (3,)) for one point, (ne (n,), B (n, 3)) for n."""
    cfg, loaded = _equilibrium(EQDSK_File, Norm_ne_File, Norm_rp_File, options)
    xyz = np.atleast_2d(np.asarray(P0, float)) * P.UNIT_SCALE[cfg.mesh_units]
    with contextlib.redirect_stdout(io.StringIO()):
        ne, B, _, _, _ = P.evaluate(np.arange(1, len(xyz) + 1), xyz, cfg, loaded)
    return (ne, B[0]) if len(xyz) == 1 else (ne, B)


def Get_Pcrtc_Equilibrium_Data(P0, k0, Freq, Rpl, Cpl,
                               EQDSK_File, Norm_ne_File, Norm_rp_File,
                               step=1e-3, tol=1e-3, max_iter=100000, **options):
    """Density and B where the beam reaches the critical density.

    Marches from the point where the ray P0 + t k0 enters the plasma sphere
    (radius Rpl, centre Cpl) in steps of `step` (forwards if ne < nc there,
    backwards otherwise) until ne - nc changes sign, then interpolates
    linearly. nc = eps0 m_e (2 pi Freq / e)^2.
    Returns (ne at that point ~ nc, B (3,) there)."""
    import utilities as gu
    wo = 2.0 * np.pi * Freq
    nc = EPS0 * ME * (wo / QE) ** 2
    get = lambda p: Get_Point_Equilibrium_Data(p, EQDSK_File, Norm_ne_File, Norm_rp_File,
                                               **options)
    k0 = np.asarray(k0, float)
    Pn = gu.LineSphere_Intersect(P0, k0, Rpl, Cpl)
    if Pn is None:
        raise ValueError("Get_Pcrtc_Equilibrium_Data: the beam does not hit the plasma sphere")
    ne_prev, B_prev = get(Pn)
    if np.isclose(ne_prev[0], nc, rtol=tol):
        return nc, B_prev
    direction = k0 if ne_prev[0] < nc else -k0
    for _ in range(max_iter):
        Pn_next = Pn + direction * step
        ne_next, B_next = get(Pn_next)
        if (ne_prev[0] - nc) * (ne_next[0] - nc) <= 0:
            t = (nc - ne_prev[0]) / (ne_next[0] - ne_prev[0])
            ni, Bcr = get(Pn + t * (Pn_next - Pn))
            return ni[0], Bcr
        Pn, ne_prev = Pn_next, ne_next
    raise RuntimeError("Critical density not found within max_iter")


if __name__ == "__main__":
    print(f"EQDSK file: {EQDSK_FILE}")
    print(f"Toroidal convention: {TOROIDAL_CONVENTION}, tokamak axis = {TOKAMAK_AXIS}, "
          f"turbulence = {TURBULENCE}")
    Generate_Plasma_Files(EQDSK_FILE, NE_FILE, NE_X_FILE, PROBLEM_FOLDER, PROBLEM_NAME,
                          PLASMA_FOLDER, PLASMA_NAME)
