#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eqdsk_to_ermes_simp.py
======================
Minimal EQDSK -> ERMES plasma files (no turbulence): the electron density and
the external magnetic field at every node of a meshed GiD/ERMES problem.

HOW TO USE
----------
1. Mesh the problem in GiD so that <PROBLEM_FOLDER><PROBLEM_NAME>.gid/
   <PROBLEM_NAME>-1.dat exists (problem_generator.py does it).
2. Edit the USER PARAMETERS and run:  python eqdsk_to_ermes_simp.py

FILES PRODUCED  (<name> = PLASMA_NAME, sorted by node id)
---------------------------------------------------------
    ne_<name>.dat    node  n_e [m^-3]
    Be_<name>.dat    node  Bx By Bz [T]   (Cartesian, mesh frame)

METHOD
------
* EQDSK read and B, psi_N computed with eqdsk_to_scotty.py (read_geqdsk,
  compute_fields): B_R = -(1/R) dpsi/dZ, B_Z = (1/R) dpsi/dR, B_t = F/R with
  the vacuum F outside the LCFS, psi_N = 0 on axis and 1 on the LCFS.
* ne(psi_N) from the 1-column files NE_FILE (values) and NE_X_FILE
  (rho_pol = sqrt(psi_N)); beyond the last profile point ne keeps its last
  value.
* Frame: z = tokamak axis, R = sqrt(x^2 + y^2), (R, phi, Z) right-handed
  (utilities.Transform_FRZT_2_FXYZ).

eqdsk_to_ermes_turb.py does the same (plus turbulence, blobs, private-flux
region handling) and is the version used by problem_generator.py. Both agree
to ~1e-3 in ne and ~5e-3 T in B on the MAST-U sample (the difference is the
private-flux region and the spline orders).

Requirements: numpy, scipy, matplotlib; eqdsk_to_scotty.py and utilities.py
in the same folder.
"""

###############################################################################
###############################  USER PARAMETERS  #############################
###############################################################################

EQDSK_FILE   = "./EQDSK_Samples/mast-u-sample.eqdsk"   # G-EQDSK file
NE_FILE      = "./EQDSK_Samples/rho_ne.txt"            # ne [m^-3], 1 column
NE_X_FILE    = "./EQDSK_Samples/rho_flux.txt"          # rho_pol, 1 column

PROBLEM_NAME   = "MASTU-2D-t20_p20"    # GiD problem name
PROBLEM_FOLDER = "./"                  # folder containing <PROBLEM_NAME>.gid

PLASMA_NAME   = "2D-t20_p20"           # -> ne_<NAME>.dat, Be_<NAME>.dat
PLASMA_FOLDER = "./Plasma/"

B_CARTESIAN  = True    # True : Bx, By, Bz in the mesh frame (ERMES 3D models)
                       # False: B_R, B_t, B_Z (cylindrical components)
SPLINE_ORDER = 5       # spline order used to evaluate psi_N, B and ne at the nodes

###############################################################################
########################  END OF USER PARAMETERS  #############################
###############################################################################

import os

import numpy as np
from scipy.interpolate import RectBivariateSpline, UnivariateSpline

import eqdsk_to_scotty as e2s
import utilities as gu


def main():
    mesh_file = os.path.join(PROBLEM_FOLDER, PROBLEM_NAME + ".gid", PROBLEM_NAME + "-1.dat")
    ne_out = os.path.join(PLASMA_FOLDER, "ne_" + PLASMA_NAME + ".dat")
    Be_out = os.path.join(PLASMA_FOLDER, "Be_" + PLASMA_NAME + ".dat")

    print("EQDSK file: " + EQDSK_FILE)
    print("Mesh  file: " + mesh_file)

    # --- Equilibrium on the EQDSK (R, Z) grid --------------------------------
    eq = e2s.read_geqdsk(EQDSK_FILE)
    B_R, B_t, B_Z, psi_n = e2s.compute_fields(eq)

    k = SPLINE_ORDER
    spline = lambda a: RectBivariateSpline(eq.R, eq.Z, a, kx=k, ky=k, s=0)
    Ps_spline, Br_spline, Bt_spline, Bz_spline = map(spline, (psi_n, B_R, B_t, B_Z))

    # --- Density profile ne(psi_N); ext=3: last value beyond the profile -----
    psi_ne, ne_prof = e2s.load_profile(NE_FILE, NE_X_FILE, "rho")
    ne_spline = UnivariateSpline(psi_ne, ne_prof, k=k, s=0, ext=3)

    # --- Mesh nodes -----------------------------------------------------------
    nId, nX, nY, nZ = gu.Read_GiD_Nodes(mesh_file)
    Rnod = np.sqrt(nX**2 + nY**2)          # tokamak axis = z
    Znod = nZ

    outside = ((Rnod < eq.R[0]) | (Rnod > eq.R[-1]) |
               (Znod < eq.Z[0]) | (Znod > eq.Z[-1]))
    if outside.any():
        print(f"WARNING: {outside.sum()} nodes are outside the EQDSK grid "
              "(B and ne there are extrapolated).")

    # --- Fields at the nodes --------------------------------------------------
    ne = np.clip(ne_spline(Ps_spline.ev(Rnod, Znod)), 0.0, None)
    Br = Br_spline.ev(Rnod, Znod)
    Bt = Bt_spline.ev(Rnod, Znod)
    Bz = Bz_spline.ev(Rnod, Znod)

    if B_CARTESIAN:
        B1, B2, B3 = gu.Transform_FRZT_2_FXYZ(nX, nY, nZ, Br, Bz, Bt)
    else:
        B1, B2, B3 = Br, Bt, Bz

    # --- Write, sorted by node id ---------------------------------------------
    os.makedirs(PLASMA_FOLDER, exist_ok=True)
    o = np.argsort(nId)
    np.savetxt(ne_out, np.column_stack([nId[o], ne[o]]), fmt=["%d", "%.8e"])
    np.savetxt(Be_out, np.column_stack([nId[o], B1[o], B2[o], B3[o]]),
               fmt=["%d", "%.8e", "%.8e", "%.8e"])
    print(f"Wrote {ne_out} and {Be_out} ({nId.size} nodes, "
          f"ne max = {ne.max():.4e} m^-3)")


if __name__ == "__main__":
    main()
