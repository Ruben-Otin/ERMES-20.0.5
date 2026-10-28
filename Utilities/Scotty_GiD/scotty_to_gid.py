#!/usr/bin/env python3
"""
scotty_to_gid.py
===================

Run a Scotty beam-tracing simulation (https://github.com/beam-tracing/Scotty)
and turn the resulting Gaussian beam into geometry that GiD can import:

  * STL   : triangulated "beam envelope" (tube of elliptical cross-sections,
            capped at both ends).  GiD: Files > Import > STL.
  * NURBS : the same envelope written as exact rational B-spline surfaces
            (IGES entity 128) in an .igs file.  GiD: Files > Import > IGES.
            Unlike the STL this is real CAD geometry, so GiD can re-mesh it
            at whatever element size the ERMES model needs.
  * FEM ENVELOPE : a closed polyhedron with PLANAR faces that encloses the
            beam at a minimum distance ``envelope_separation`` from it, to be
            used as the ERMES computational volume.  Its first and last faces
            lie in the planes of the beam's end caps.  Written as
            <stem>_envelope.igs (one trimmed planar surface per face) and/or
            <stem>_envelope.stl.  See section 10.

The script has two stages, which can be used together or separately:

  1. RUN     -- builds the Scotty keyword dictionary from the USER
                CONFIGURATION block below and calls ``beam_me_up``.
  2. EXPORT  -- reads a Scotty ``scotty_output*.h5`` file, (optionally)
                extends the beam through vacuum past its first/last point,
                and writes the STL / IGES / text outputs.

Geometry of the envelope
------------------------
At every point of the central ray, Scotty gives the local transverse basis
(x_hat, y_hat) and the complex beam tensor Psi projected on it.  The 1/e
amplitude contour of the beam is the ellipse

    w^T Im(Psi) w = 2      (w = transverse displacement)

so the semi-axes are W_i = sqrt(2 / lambda_i), lambda_i the eigenvalues of
Im(Psi).  Each cross-section is stored as: 
             P(theta) = C + A cos(theta) + B sin(theta)
(C = ray point, A/B = two conjugate semi-diameters).  A and B are re-phased
from one point to the next so the parametrisation never jumps (eigenvectors
returned by numpy have arbitrary sign -- in the old script this produced
twisted, self-intersecting triangles wherever a sign flipped).

The NURBS surface is exact in the circumferential direction (an ellipse is an
affine image of a circle, and a circle is an exact rational quadratic NURBS);
along the beam it is a cubic B-spline skinned through a set of cross-sections
that is chosen adaptively until the surface reproduces EVERY Scotty
cross-section to within ``--nurbs-tol``.

ERMES input / output
--------------------
With ``ermes_input = True`` the launch is defined as in ERMES: launch point
(Cartesian X, Y, Z), beam direction k0, waist w0 and waist position. The
script converts it to Scotty's (R, zeta, Z), TORBEAM angles, launch beam
width and curvature (section 11). After Scotty runs, the polarisation at the
plasma entry is written back in ERMES form (Epol, |Epar|, |Eper|, phases) in
<stem>_ermes_launch.txt and compared with the ERMES polarisation, if given.

O-X mode conversion diagnostics
-------------------------------
``ox_diagnostics = True`` evaluates, along the ray, N_par = N.b and
N_y = N.(b x grad psi / |grad psi|), the optimal N_par,opt = sqrt(Y/(1+Y)),
L_n and the Mjolhus transmission
    T = exp{-pi k0 L_n sqrt(Y/2) [2(1+Y)(|N_par|-N_par,opt)^2 + N_y^2]}
plus a beam-averaged T for the angular spectrum of the launched Gaussian
beam. B and psi come from topfile<suffix>.json, ne from ne<suffix>.dat.
Beam tracing is not valid AT the cut-off: the values there are estimates.

ERMES FEM cost
--------------
``fem_cost_estimate = True`` estimates the tetrahedra, nodes, unknowns,
matrix non-zeros, memory and solve time if the FEM envelope is meshed for
ERMES (section 11). Order-of-magnitude only: calibrate the constants with
one of your own runs.

Vacuum extension
----------------
``--extend-first`` / ``--extend-last`` prolong the beam by a given distance
beyond its first/last computed point as a Gaussian beam in free space: the ray
continues straight along the local wavevector and Psi evolves with the
free-space law   Psi^-1(s) = Psi^-1(0) + (s/k0) I.

Axis convention (Scotty's Cartesian outputs)
--------------------------------------------
X, Z -> tokamak R, Z ; Y -> toroidal, "into the page" (zeta = 0 plane is XZ)

Parameters
----------
All parameters you normally change are plain assignments in section 1
("USER PARAMETERS") right below the imports, each with its possible values
in a comment. Edit them there, or change any of them for one run from the
command line with  --name value  (same name, "-" instead of "_").
    python scotty_to_gid.py --show-params     # values in use
    python scotty_to_gid.py --help            # all flags with defaults

Usage examples
--------------
    # run Scotty with the parameters of section 1, then write the GiD files
    python scotty_to_gid.py

    # same, STL + NURBS, no plot window
    python scotty_to_gid.py --output-format both --no-show-plot

    # different frequency / angles / antenna position for this run only
    python scotty_to_gid.py --frequency 32e9 --poloidal-angle 10 \\
        --launch-position 1.0 0 1.4

    # do not run Scotty: only export an existing output file
    python scotty_to_gid.py --h5                              \\
        ./ScottyResults/scotty_output_freq28.0_pol0.0_tor0.0000.h5 \\
        --output-format nurbs                                      \\
        --extend-first 0.5 --extend-last 0.5 --no-show-plot

    # only run Scotty
    python scotty_to_gid.py --no-export-geometry

    # see what is inside a Scotty file
    python scotty_to_gid.py --h5 file.h5 --list-datasets

Requires Python >= 3.9, numpy and h5py (+ matplotlib for the plot, + scotty
for stage 1).
"""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

###############################################################################
# 1. USER PARAMETERS
###############################################################################
# Change the values here.
# Any of them can ALSO be changed for a single run from the command line:
#     --name value          (the same name, with "-" instead of "_")
#     --name / --no-name    for True/False parameters
# e.g.  python scotty_to_gid.py 
#              --frequency 32e9 --launch-position 1.0 0 1.4 \
#              --output-format both --no-show-plot
# "python scotty_to_gid.py --show-params" prints the values in use.
###############################################################################

_names_before_parameters = set(globals())

# ---------------------------- What to do -------------------------------------

run_scotty = True              # True : run Scotty, then export its output
                               # False: do NOT run Scotty; export scotty_h5_file
                               
scotty_h5_file = None          # Scotty output (.h5) to export when run_scotty = False
                               # None : newest scotty_output*.h5 in scotty_output_path
                               # (giving a file on the command line with --h5 FILE
                               #  implies run_scotty = False)
                               
export_geometry = True         # True : write the GiD files after Scotty
                               # False: stop after running Scotty
                               
show_plot = True               # True / False: 3D plot of the beam surface
                               # (the script waits until the window is closed)
                               
print_scotty_keywords = False  # True / False: print every keyword passed to Scotty

# ---------------------------- Scotty: device ---------------------------------

device = "DBS_UCLA_MAST-U"     # Scotty device preset (sets launch beam width,
                               # curvature and default antenna position)
                               
find_B_method = "torbeam"      # "torbeam" | "EFITpp" | "UDA_saved" | "UDA"

user = "rotin"                 # user name passed to the device preset

# ---------------------------- Scotty: launch ---------------------------------

frequency = 34.8e9             # launch frequency [Hz]

poloidal_angle = -9.276669     # poloidal launch angle [deg], TORBEAM convention
                               # (= ERMES k0 below, used when ermes_input = False)
                               
toroidal_angle = -15.045861    # toroidal launch angle [deg], TORBEAM convention
mode_flag = 1                  # -1 | +1  (device dependent; MAST-U: -1 = X, +1 = O)

launch_position = (2.073813918, 0.035598947, 0.03249)  # antenna position (R [m], zeta [rad], Z [m])
                                       # (= ERMES launch point X,Y,Z = 2.0725, 0.07381, 0.03249)
                                       # None: use the device preset position
                                       # zeta is normally 0 (it only rotates the picture)
                               
launch_beam_width = 17.92e-3   # launch beam width W (1/e field radius) [m]
                               # None: from the device preset
                               
launch_beam_curvature = 0.0    # launch wavefront curvature [1/m]: < 0 converging,
                               # 0 = waist at launch, > 0 diverging
                               # None: from the device preset

# ---------------------------- ERMES launch input -----------------------------

ermes_input = True             # True : define the launch with the ERMES data below;
                               # (the Scotty launch above is the SAME beam, so
                               #  True / False give identical runs)
                               # REPLACES launch_position, poloidal_angle,
                               # toroidal_angle, launch_beam_width/curvature
                               # (frequency and mode_flag are still used)
                               # False: use the Scotty launch parameters above
                               
ermes_launch_point = (2.0725, 0.07381, 0.03249)  # X, Y, Z [m]
                                                 # Cartesian, Z = tokamak axis, origin on the axis
                               
ermes_k0 = (-0.961603, 0.222113, 0.161202) # beam direction k0 (any length)

ermes_beam_width = 17.92e-3    # beam waist w0 (1/e field radius) [m]
ermes_waist_distance = 0.0     # launch point -> waist, along k0 [m]
                               # > 0: waist ahead (converging) | 0: at launch
                               # < 0: waist behind (diverging)
                               
ermes_epol = (0.025167, -0.513531, 0.857702)  # ERMES Epol (perpendicular to k0), to
                                              # compare polarisations; None: no comparison
                               
ermes_epar = (0.614842, 0.0)         # ERMES |Epar|, phase [rad]  (along Epol)
ermes_eper = (-0.78865, -1.5708)     # ERMES |Eper|, phase [rad]  (along k0 x Epol)
                              


# ---------------------------- Scotty: plasma input ---------------------------

magnetic_data_path = "./Scotty_input/"   # folder of topfile<suffix>.json
ne_data_path = "./Scotty_input/"         # folder of ne<suffix>.dat

input_filename_suffix = "_sample"  # Scotty reads ne<suffix>.dat, topfile<suffix>.json

poloidal_flux_zero_density = 1.134 # normalised poloidal flux where ne = 0
                                   # (ideally 1.0; depends on your ne file)
                                   
poloidal_flux_enter = 1.134        # where the beam enters the plasma;
                                   # keep equal to poloidal_flux_zero_density

# ---------------------------- Scotty: run options ----------------------------

scotty_output_path = "./Scotty_results/"  # folder where Scotty writes the .h5 file

vacuum_propagation = True      # True : launch point is in vacuum | False: in plasma

psi_bc_flag = "continuous"     # vacuum/plasma boundary condition for Psi
                               # (Scotty Psi_BC_flag, e.g. "continuous")
                               
scotty_figures = True          # True / False: Scotty's own diagnostic figures

# ---------------------------- Scotty: ODE solver (rarely changed) ------------

delta_R = -1e-5                # finite-difference step in R [m]
delta_Z = -1e-5                # finite-difference step in Z [m]
delta_K_R = 0.1                # finite-difference step in K_R [1/m]
delta_K_zeta = 0.1             # finite-difference step in K_zeta
delta_K_Z = 0.1                # finite-difference step in K_Z [1/m]
atol = 1e-7                    # ODE absolute tolerance
rtol = 1e-4                    # ODE relative tolerance
len_tau = 1002                 # number of points along the ray (= beam cross-sections)

# ---------------------------- O-X mode conversion diagnostics ----------------

ox_diagnostics = True            # True : N_par, N_y, L_n and Mjolhus O-X transmission
                                 # along the ray (reads topfile<suffix>.json and
                                 # ne<suffix>.dat) | False: skip
                               
ox_ne_radial_coordinate = "rho"  # 1st column of ne<suffix>.dat:
                                 # "rho" = sqrt(psi_n) (Scotty) | "psi" = psi_n
                               
ox_ne_units = 1e19               # factor from the ne file densities to m^-3

# ---------------------------- Vacuum extension of the beam -------------------

extend_first_from = "launch"   # how the beam IN FRONT of Scotty's first point is made:
                               # "launch": vacuum Gaussian beam from the LAUNCH
                               #   definition (antenna point, launch_beam_width and
                               #   curvature, read from the .h5 when present), from
                               #   extend_first BEHIND the antenna to the plasma entry.
                               #   Correct: Scotty's Psi at its first point is on the
                               #   PLASMA side of the boundary (Psi_BC jump), so it
                               #   must not be propagated back through vacuum.
                               # "entry": old behaviour, extend_first back from Scotty's
                               #   first point using its (plasma-side) Psi
                               
extend_first = -0.3             # "launch": where the beam starts, from the antenna [m]
                               #   > 0: this far BEHIND the antenna
                               #   0  : at the antenna / port
                               #   < 0: this far AHEAD of the antenna, towards the
                               #        plasma (smaller FEM volume); the beam there
                               #        is the launched beam after |extend_first| m
                               #        of vacuum (wider, curved, still circular);
                               #        <stem>_ermes_launch.txt gives the ERMES port
                               #        data (point, width, waist position) there
                               # "entry": distance back from Scotty's first point [m]
                               
extend_last = 0.3              # extend the beam FORWARDS from its last point [m]; 0 = none
                               # (0 for O-X: past the turning point it is the reflected beam)
                               
extend_points = 20             # number of cross-sections in each extension

# ---------------------------- Geometry output (GiD) --------------------------

output_format = "both"         # "stl"   : triangle mesh   (GiD: Import > STL)
                               # "nurbs" : NURBS surfaces  (GiD: Import > IGES)
                               # "both"  | "none"
                               
units = "m"                    # "m" | "cm" | "mm"  unit of all written geometry
caps = True                    # True : close both ends of the beam | False: open tube

out_dir = "./GiD_input/"       # folder for ALL files written after Scotty
                               # (.stl, .igs, central ray, E-field/k files);
                               # created automatically if it does not exist
                               
file_stem = None               # start of the output file names
                               # None: from the .h5 name ('.' replaced by '_')
                               
stl_points_per_section = 64    # STL: points around each cross-section
stl_stride = 1                 # STL: use every N-th cross-section (1 = all)
stl_ascii = False              # STL: True = ASCII file | False = binary file
nurbs_tolerance = 1e-4         # NURBS: max deviation from any Scotty cross-section [m]
nurbs_max_sections = 600       # NURBS: max cross-sections the surface goes through
nurbs_patches = 4              # NURBS: 4 = tube in 4 quarter surfaces (best for GiD)
                               #        1 = one closed surface
                               
nurbs_include_ray = False      # NURBS: True = add the central ray as a curve

# ---------------------------- FEM envelope (polyhedron around the beam) ------

envelope = True                # True : also build a closed polyhedron with planar
                               # faces that encloses the beam (FEM volume for ERMES)
                               # and write <stem>_envelope.igs / .stl (which of the
                               # two follows output_format)
                               
envelope_separation = 0.1    # minimum distance from the beam surface (1/e) to
                               # the side faces of the envelope [m]. The two end
                               # faces lie in the planes of the beam's end caps.
                               
envelope_sides = 10            # side faces around the beam in each cell (>= 3)

envelope_max_cell_length = 0.0 # longest cell along the beam [m]: shorter cells
                               # follow changes of the beam width more closely
                               # 0 = automatic (~ beam width + separation)
                               # -1 = no subdivision (fewest faces)
                               
envelope_ray_tolerance = 0.0   # max deviation of the ray from a straight line
                               # inside one cell [m]; 0 = separation / 2

# ---------------------------- ERMES FEM cost estimate ------------------------

fem_cost_estimate = True       # True : mesh size and solver cost if the FEM envelope
                               # is the ERMES domain (needs envelope = True) | False
                               
fem_elements_per_wavelength = 10.0  # mean tetrahedron edge h = lambda_local / this

fem_refractive_index = 1.0     # |N| in lambda_local = lambda0/|N| (use > 1 for
                               # short-wavelength regions, e.g. slow X / UHR);
                               # 0 = max |K|/k0 along the Scotty ray (>= 1)
                               
fem_element_order = 1          # 1 | 2: nodal Lagrange tetrahedra (4 / 10 nodes)
fem_dofs_per_node = 3          # unknowns per node (nodal E field: 3)

fem_vertices_per_tet = 0.164   # mesh vertices per tetrahedron (unstructured ~1/6.1)
fem_solver = "iterative"       # "direct" (MUMPS-like) | "iterative" (Krylov)
fem_fill_constant = 15.0       # direct: factor entries ~ C * Ndof^(4/3)  (calibrate)
fem_flop_constant = 60.0       # direct: factorisation flops ~ C * Ndof^2 (calibrate)

fem_cores = 128                # cores used by the solver
fem_gflops_per_core = 10.0     # sustained GFLOP/s per core
fem_iterations = 1000          # iterative: Krylov iterations assumed for the time

# ---------------------------- Other output files -----------------------------

write_central_ray = True       # True / False: <stem>_central_ray.dat (x,y,z per line)
write_efield = True            # True / False: E-field and k at both beam ends

efield_file = "E_final.dat"    # name of that end-point file

write_field_arrays = True      # True / False: E-field and k at every ray point
                               # (<stem>_Efield_along_ray.dat, <stem>_Kvector_along_ray.dat)
                               
plot_points_around = 36        # plot: points around each cross-section (the surface
                               # always uses ALL cross-sections)

###############################################################################
# (end of user parameters -- nothing below needs to be edited)
###############################################################################

USER_PARAMETERS = [n for n in globals()
                   if n not in _names_before_parameters and not n.startswith("_")]
_PARAMETER_DEFAULTS = {n: globals()[n] for n in USER_PARAMETERS}   # snapshot

# Values that are checked when given on the command line
PARAMETER_CHOICES = {
    "find_B_method": ["torbeam", "EFITpp", "UDA_saved", "UDA"],
    "mode_flag": [-1, 1],
    "output_format": ["stl", "nurbs", "both", "none"],
    "units": ["m", "cm", "mm"],
    "nurbs_patches": [1, 4],
    "ox_ne_radial_coordinate": ["rho", "psi"],
    "fem_element_order": [1, 2],
    "fem_solver": ["direct", "iterative"],
    "extend_first_from": ["launch", "entry"],
}
# Extra / shorter command-line names (also keep old command lines working)
PARAMETER_ALIASES = {
    "scotty_h5_file": ["--h5"],
    "frequency": ["--freq"],
    "poloidal_angle": ["--pol"],
    "toroidal_angle": ["--tor"],
    "output_format": ["--format"],
    "show_plot": ["--plot"],
    "write_efield": ["--efield"],
    "write_field_arrays": ["--field-arrays"],
    "write_central_ray": ["--central-ray", "--txt"],
    "nurbs_tolerance": ["--nurbs-tol"],
    "print_scotty_keywords": ["--verbose"],
}
PARAMETER_TYPES = {"scotty_h5_file": str, "file_stem": str, "launch_position": float,
                   "launch_beam_width": float, "launch_beam_curvature": float}
# Parameters that take several numbers on the command line
PARAMETER_NARGS = {"launch_position": 3, "ermes_launch_point": 3, "ermes_k0": 3,
                   "ermes_epol": 3, "ermes_epar": 2, "ermes_eper": 2}

SPEED_OF_LIGHT = 299_792_458.0  # m/s
E_CHARGE = 1.602176634e-19      # C
E_MASS = 9.1093837015e-31       # kg
EPS0 = 8.8541878128e-12         # F/m
UNIT_SCALE = {"m": 1.0, "cm": 100.0, "mm": 1000.0}


def default_parameters():
    """The values of section 1 as a namespace. From another script:
        import scotty_beam_to_gid as s
        cfg = s.default_parameters(); cfg.output_format = "nurbs"; s.main(cfg)"""
    return argparse.Namespace(**dict(_PARAMETER_DEFAULTS),
                              show_params=False, list_datasets=False)


# ##########################################################################
# 2. RUNNING SCOTTY
# ##########################################################################

def scotty_output_suffix(freq_GHz, pol_deg, tor_deg) -> str:
    """Suffix Scotty appends to its output files (same as the old run
    script, so existing result files keep matching)."""
    return f"_freq{freq_GHz:.1f}_pol{pol_deg:.1f}_tor{tor_deg:.4f}"


def run_scotty_simulation(cfg) -> str:
    """Run Scotty with the parameters in ``cfg``; return the .h5 path."""
    # Imported here so the export stage works on machines without Scotty.
    from scotty.beam_me_up import beam_me_up
    from scotty.init_bruv import get_parameters_for_Scotty

    freq_GHz = cfg.frequency / 1e9
    # Device preset (launch beam width/curvature, default antenna, ...)
    kw = get_parameters_for_Scotty(cfg.device, launch_freq_GHz=freq_GHz,
                                   find_B_method=cfg.find_B_method, user=cfg.user)

    # Translate the section-1 parameters into Scotty's keyword names
    kw.update(
        poloidal_launch_angle_Torbeam=cfg.poloidal_angle,
        toroidal_launch_angle_Torbeam=cfg.toroidal_angle,
        mode_flag=cfg.mode_flag,
        poloidal_flux_enter=cfg.poloidal_flux_enter,
        poloidal_flux_zero_density=cfg.poloidal_flux_zero_density,
        input_filename_suffix=cfg.input_filename_suffix,
        magnetic_data_path=cfg.magnetic_data_path,
        ne_data_path=cfg.ne_data_path,
        output_path=cfg.scotty_output_path,
        vacuum_propagation_flag=cfg.vacuum_propagation,
        Psi_BC_flag=cfg.psi_bc_flag,
        figure_flag=cfg.scotty_figures,
        delta_R=cfg.delta_R, delta_Z=cfg.delta_Z,
        delta_K_R=cfg.delta_K_R, delta_K_zeta=cfg.delta_K_zeta, delta_K_Z=cfg.delta_K_Z,
        atol=cfg.atol, rtol=cfg.rtol, len_tau=cfg.len_tau,
    )
    if cfg.launch_position is not None:
        # Scotty only handles launches at zeta = 0 (it treats (R, zeta, Z) like
        # Cartesian (X, Y, Z) with the angles in the zeta = 0 frame): launch at
        # zeta = 0 with the same LOCAL angles and rotate the output by zeta
        # afterwards (rotate_scotty_output). Exact: the plasma is axisymmetric.
        R0, zeta0, Z0 = map(float, cfg.launch_position)
        kw["launch_position"] = np.array([R0, 0.0, Z0])
        if zeta0 != 0.0:
            print(f"Launching Scotty at zeta = 0; its output is rotated by zeta = "
                  f"{zeta0:.6f} rad ({np.degrees(zeta0):.4f} deg) back to your frame.")
    if cfg.launch_beam_width is not None:
        kw["launch_beam_width"] = float(cfg.launch_beam_width)
    if cfg.launch_beam_curvature is not None:
        kw["launch_beam_curvature"] = float(cfg.launch_beam_curvature)
    kw["output_filename_suffix"] = scotty_output_suffix(
        kw["launch_freq_GHz"], cfg.poloidal_angle, cfg.toroidal_angle)

    cfg.used_launch = dict(position=(float(kw["launch_position"][0]), scotty_zeta_rotation(cfg),
                                     float(kw["launch_position"][2])),
                           width=float(kw["launch_beam_width"]),
                           curvature=float(kw["launch_beam_curvature"]))
    print("Scotty launch position (R, zeta, Z):", kw["launch_position"])
    print("Scotty launch beam width  [m]      :", kw["launch_beam_width"])
    print("Scotty launch beam curvature [1/m] :", kw["launch_beam_curvature"])
    if cfg.print_scotty_keywords:
        print(kw)

    os.makedirs(kw["output_path"], exist_ok=True)
    try:
        beam_me_up(**kw)
    except ValueError as err:
        if "t_eval" not in str(err):
            raise
        # Scotty first traces a simple ray to find where the beam leaves the
        # plasma (tau_leave) and then solves on linspace(0, tau_leave). If that
        # ray ends at once, tau_leave ~ 0 and scipy raises this error.
        R, zeta, Z = kw["launch_position"]
        sys.exit(
            "\nScotty stopped: the beam ray ended immediately (tau_leave ~ 0).\n"
            f"Launch position R = {R} m, Z = {Z} m, angles pol = {cfg.poloidal_angle} deg, "
            f"tor = {cfg.toroidal_angle} deg.\n"
            "Usual causes:\n"
            "  * the launch point is INSIDE the plasma (normalised poloidal flux below\n"
            "    poloidal_flux_enter) while vacuum_propagation = True. Move the antenna\n"
            "    outside the last closed flux surface;\n"
            "  * the launch direction points away from the plasma, so the beam never\n"
            "    enters it (check the sign of the TORBEAM angles);\n"
            "  * the wave is cut off right at the plasma edge for this mode/frequency.\n")

    h5 = os.path.join(kw["output_path"], f"scotty_output{kw['output_filename_suffix']}.h5")
    if not os.path.isfile(h5):
        h5 = newest_scotty_output(kw["output_path"])
    print(f"Scotty output: {h5}")
    return h5


def newest_scotty_output(folder) -> str:
    found = sorted(glob.glob(os.path.join(folder, "scotty_output*.h5")), key=os.path.getmtime)
    if not found:
        sys.exit(f"No scotty_output*.h5 file found in {folder}")
    return found[-1]


# ##########################################################################
# 3. READING THE SCOTTY HDF5 OUTPUT
# ##########################################################################

def find_dataset(h5file, candidate_names):
    """Search an HDF5 file recursively for a dataset whose leaf name matches
    one of ``candidate_names`` (case-insensitive, in order of preference).

    Exact leaf-name matches win; otherwise a looser match is accepted (leaf
    starting/ending with ``<name>_`` / ``_<name>``) and reported, because the
    dataset layout differs between Scotty versions. Returns None if absent.
    """
    import h5py
    wanted = [n.lower() for n in candidate_names]
    exact, fuzzy = {}, {}

    def visitor(name, obj):
        if not isinstance(obj, h5py.Dataset):
            return
        leaf = name.split("/")[-1].lower()
        for w in wanted:
            if leaf == w:
                exact.setdefault(w, (name, obj))
            elif leaf.startswith(w + "_") or leaf.endswith("_" + w):
                fuzzy.setdefault(w, (name, obj))

    h5file.visititems(visitor)
    for w in wanted:
        if w in exact:
            return exact[w][1]
    for w in wanted:
        if w in fuzzy:
            print(f"NOTE: using '{fuzzy[w][0]}' as a fuzzy match for '{w}'.")
            return fuzzy[w][1]
    return None


def list_all_datasets(h5file, with_shapes=False):
    import h5py
    out = []

    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset):
            out.append(f"{name}  {obj.shape}  {obj.dtype}" if with_shapes else name)

    h5file.visititems(visitor)
    return out


def to_complex(arr) -> np.ndarray:
    """Return a numpy array, converting h5netcdf/xarray-style compound
    (r, i) / (real, imag) records into genuine complex numbers."""
    arr = np.asarray(arr)
    if np.iscomplexobj(arr):
        return arr
    names = arr.dtype.names
    if names:
        low = {n.lower(): n for n in names}
        for re_key, im_key in (("r", "i"), ("real", "imag")):
            if re_key in low and im_key in low:
                return arr[low[re_key]].astype(float) + 1j * arr[low[im_key]].astype(float)
    return arr.astype(float)


def as_real(ds):
    return np.real(to_complex(ds[()]))


def ensure_n_by_3(arr) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 2 and arr.shape[1] == 3:
        return arr
    if arr.ndim == 2 and arr.shape[0] == 3:
        return arr.T
    raise ValueError(f"Expected an (N, 3) array, got shape {arr.shape}")


def broadcast_to_n(value, n):
    """Scalars / size-1 arrays -> constant length-n array; else unchanged."""
    arr = np.asarray(value)
    return np.full(n, arr.reshape(-1)[0]) if arr.size == 1 else arr


def load_scotty_beam(path: str) -> dict:
    """Read what the export needs from a Scotty .h5 file.

    Keys: beam, xhat, yhat (N,3 Cartesian); psi_xx/psi_xy/psi_yy (complex,
    local frame); distance; K_R, K_zeta, K_Z; e_hat; q_R, q_zeta; omega.
    Optional entries that are missing are None.
    """
    import h5py
    d = {}
    with h5py.File(path, "r") as f:
        beam_ds = find_dataset(f, ["beam_cartesian"])
        xhat_ds = find_dataset(f, ["x_hat_cartesian", "xhat", "x_hat"])
        yhat_ds = find_dataset(f, ["y_hat_cartesian", "yhat", "y_hat"])
        if beam_ds is None or xhat_ds is None or yhat_ds is None:
            raise KeyError("Could not find beam_cartesian / x_hat / y_hat. Datasets present:\n  "
                           + "\n  ".join(list_all_datasets(f)))
        d["beam"] = ensure_n_by_3(beam_ds[()])
        d["xhat"] = ensure_n_by_3(xhat_ds[()])
        d["yhat"] = ensure_n_by_3(yhat_ds[()])
        n = d["beam"].shape[0]

        # Beam tensor in the local (xhat, yhat) frame. Prefer Psi_xx/xy/yy;
        # otherwise project the full 3x3 Cartesian tensor.
        ds = [find_dataset(f, [k]) for k in ("psi_xx", "psi_xy", "psi_yy")]
        if all(x is not None for x in ds):
            d["psi_xx"], d["psi_xy"], d["psi_yy"] = (to_complex(x[()]) for x in ds)
        else:
            p3 = find_dataset(f, ["psi_3d_cartesian"])
            if p3 is not None:
                p3 = to_complex(p3[()])
                proj = lambda a, b: np.einsum("ni,nij,nj->n", a, p3, b)
                d["psi_xx"] = proj(d["xhat"], d["xhat"])
                d["psi_xy"] = proj(d["xhat"], d["yhat"])
                d["psi_yy"] = proj(d["yhat"], d["yhat"])
            else:
                print("WARNING: no Psi data found -- using a constant circular "
                      "cross-section (NOT the real beam width); vacuum "
                      "extension disabled.")
                d["psi_xx"] = d["psi_xy"] = d["psi_yy"] = None

        dist = find_dataset(f, ["distance_along_line", "l_lc"])
        d["distance"] = as_real(dist) if dist is not None else np.arange(n, dtype=float)

        # Wavevector in cylindrical components. K_zeta is the canonical
        # (dimensionless) component = R * physical toroidal wavenumber, and
        # is conserved (toroidal symmetry), so it may be stored once.
        k_r, k_z = find_dataset(f, ["k_r", "kr"]), find_dataset(f, ["k_z", "kz"])
        k_zeta = find_dataset(f, ["k_zeta", "kzeta"])
        if k_zeta is None:
            k_zeta = find_dataset(f, ["k_zeta_initial", "kzeta_initial"])
            if k_zeta is not None:
                print("NOTE: using the constant K_zeta_initial for every point.")
        if k_r is not None and k_zeta is not None and k_z is not None:
            d["K_R"], d["K_Z"] = as_real(k_r), as_real(k_z)
            d["K_zeta"] = broadcast_to_n(as_real(k_zeta), n)
        else:
            d["K_R"] = d["K_zeta"] = d["K_Z"] = None

        # Polarisation vector, physical cylindrical (R, zeta, Z) components.
        e_hat = find_dataset(f, ["e_hat", "ehat"])
        d["e_hat"] = ensure_n_by_3(to_complex(e_hat[()])) if e_hat is not None else None

        # Ray position (R, zeta): needed to rotate cylindrical -> Cartesian.
        q_r, q_zeta = find_dataset(f, ["q_r", "qr"]), find_dataset(f, ["q_zeta", "qzeta"])
        d["q_R"] = as_real(q_r) if q_r is not None else None
        d["q_zeta"] = as_real(q_zeta) if q_zeta is not None else None
        if d["q_zeta"] is not None and np.max(np.abs(d["q_zeta"])) > 3 * np.pi:
            print("NOTE: q_zeta looks like degrees -- converting to radians.")
            d["q_zeta"] = np.deg2rad(d["q_zeta"])

        # launch definition as Scotty used it (inputs), if stored
        for key, names in (("launch_position", ["launch_position"]),
                           ("launch_beam_width", ["launch_beam_width"]),
                           ("launch_beam_curvature", ["launch_beam_curvature"])):
            ds = find_dataset(f, names)
            d[key] = np.real(to_complex(ds[()])).ravel() if ds is not None else None

        omega = find_dataset(f, ["launch_angular_frequency"])
        f_ghz = find_dataset(f, ["launch_freq_ghz"])
        if omega is not None:
            d["omega"] = float(as_real(omega))
        elif f_ghz is not None:
            d["omega"] = 2 * np.pi * 1e9 * float(as_real(f_ghz))
        else:
            d["omega"] = None
    return d


# ##########################################################################
# 4. BEAM PHYSICS HELPERS
# ##########################################################################

@dataclass
class Beam:
    """Everything known at each point of the (possibly extended) ray."""
    points: np.ndarray                 # (N,3) central ray, Cartesian
    xhat: np.ndarray                   # (N,3) local transverse basis
    yhat: np.ndarray                   # (N,3)
    psi_xx: Optional[np.ndarray]       # (N,) complex, or None
    psi_xy: Optional[np.ndarray]
    psi_yy: Optional[np.ndarray]
    distance: np.ndarray               # (N,) distance along the ray
    E: Optional[np.ndarray] = None     # (N,3) complex polarisation, Cartesian
    K: Optional[np.ndarray] = None     # (N,3) wavevector, Cartesian [1/m]

    @property
    def n(self):
        return self.points.shape[0]


def cyl_to_cart_vector(v_r, v_zeta, v_z, zeta, canonical=False, R=None):
    """Rotate cylindrical (R, zeta, Z) vector components to Cartesian, with
    X = R cos(zeta), Y = R sin(zeta). With ``canonical`` the zeta component
    is first divided by R (Scotty's K_zeta; Hall-Chen, Parra & Hillesheim
    2022, eqs 56-58). e_hat is a physical (non-canonical) unit vector."""
    v_phi = v_zeta / R if canonical else v_zeta
    c, s = np.cos(zeta), np.sin(zeta)
    return np.stack([v_r * c - v_phi * s, v_r * s + v_phi * c, v_z], axis=1)


def cartesian_fields(d):
    """(E_cart, K_cart) for every Scotty point, or (None, None)."""
    need = ("K_R", "e_hat", "q_R", "q_zeta")
    if any(d[k] is None for k in need):
        return None, None
    K = cyl_to_cart_vector(d["K_R"], d["K_zeta"], d["K_Z"], d["q_zeta"],
                           canonical=True, R=d["q_R"])
    e = d["e_hat"]
    E = cyl_to_cart_vector(e[:, 0], e[:, 1], e[:, 2], d["q_zeta"])
    return E, K


def propagate_psi_vacuum(psi_xx, psi_xy, psi_yy, s, k0):
    """Free-space Gaussian beam: Psi^-1(s) = Psi^-1(0) + (s/k0) I.
    ``s`` may be an array of signed distances; returns three arrays."""
    M = np.array([[psi_xx, psi_xy], [psi_xy, psi_yy]], dtype=complex)
    s = np.atleast_1d(s)
    Minv = np.linalg.inv(M)[None] + (s / k0)[:, None, None] * np.eye(2)[None]
    Mn = np.linalg.inv(Minv)
    return Mn[:, 0, 0], Mn[:, 0, 1], Mn[:, 1, 1]


def extend_beam(beam: Beam, k0, dist_first, dist_last, n_ext) -> Beam:
    """Prolong the beam through vacuum by ``dist_first`` before its first
    point and ``dist_last`` after its last point (0 = no extension).

    The ray is continued in a straight line along the local wavevector;
    the transverse basis and the polarisation are kept constant, and Psi is
    propagated analytically. Needs Psi and the wavevector."""
    if beam.psi_xx is None or beam.K is None or k0 is None:
        print("WARNING: vacuum extension needs Psi, the wavevector and the "
              "launch frequency -- extension skipped.")
        return beam

    parts = {f: [getattr(beam, f)] for f in
             ("points", "xhat", "yhat", "psi_xx", "psi_xy", "psi_yy", "distance", "E", "K")}

    for which, dist in (("first", dist_first), ("last", dist_last)):
        if dist <= 0:
            continue
        i = 0 if which == "first" else -1
        sign = -1.0 if which == "first" else +1.0
        t = np.real(beam.K[i]) / np.linalg.norm(np.real(beam.K[i]))
        s = np.linspace(dist / n_ext, dist, n_ext)          # near -> far
        if which == "first":
            s = s[::-1]                                     # far -> near (ray order)
        new = {
            "points": beam.points[i] + sign * s[:, None] * t,
            "xhat": np.tile(beam.xhat[i], (n_ext, 1)),
            "yhat": np.tile(beam.yhat[i], (n_ext, 1)),
            "distance": beam.distance[i] + sign * s,
            "E": np.tile(beam.E[i], (n_ext, 1)) if beam.E is not None else None,
            "K": np.tile(beam.K[i], (n_ext, 1)),
        }
        new["psi_xx"], new["psi_xy"], new["psi_yy"] = propagate_psi_vacuum(
            beam.psi_xx[i], beam.psi_xy[i], beam.psi_yy[i], sign * s, k0)
        for f, arr in new.items():
            if which == "first":
                parts[f].insert(0, arr)
            else:
                parts[f].append(arr)
        print(f"Extended the beam {dist:.4g} m {'before the first' if which == 'first' else 'past the last'} "
              f"point ({n_ext} cross-sections).")

    cat = lambda lst: None if any(a is None for a in lst) else np.concatenate(lst, axis=0)
    return Beam(**{f: cat(v) for f, v in parts.items()})


def scotty_zeta_rotation(cfg):
    """Toroidal angle of the launch point: Scotty is run at zeta = 0 and its
    output is rotated by this angle (see run_scotty_simulation)."""
    if getattr(cfg, "launch_position", None) is None:
        return 0.0
    return float(cfg.launch_position[1])


def rotate_scotty_output(d, dzeta):
    """Rotate a loaded Scotty output (launched at zeta = 0) by dzeta about the
    tokamak axis Z, in place: Cartesian beam / xhat / yhat, q_zeta and the
    stored launch position."""
    if not dzeta:
        return d
    c, s = np.cos(dzeta), np.sin(dzeta)
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    for key in ("beam", "xhat", "yhat"):
        if d.get(key) is not None:
            d[key] = d[key] @ Rz.T
    if d.get("q_zeta") is not None:
        d["q_zeta"] = d["q_zeta"] + dzeta
    lp = d.get("launch_position")
    if lp is not None and np.size(lp) == 3 and abs(float(lp[1])) < 1e-12:
        d["launch_position"] = np.array([lp[0], dzeta, lp[2]], float)
    print(f"Scotty output rotated by zeta = {dzeta:.6f} rad ({np.degrees(dzeta):.4f} deg) "
          "about the tokamak axis.")
    return d


def launch_definition(d, cfg):
    """(launch point Cartesian, width, curvature) Scotty used: from the .h5
    inputs, else from the run just made, else from cfg. None if unknown."""
    used = getattr(cfg, "used_launch", None) or {}
    pick = lambda h5key, ukey, cfgval: (
        d.get(h5key) if d.get(h5key) is not None else
        (np.atleast_1d(used[ukey]) if ukey in used else
         (np.atleast_1d(cfgval) if cfgval is not None else None)))
    pos = pick("launch_position", "position", cfg.launch_position)
    W = pick("launch_beam_width", "width", cfg.launch_beam_width)
    c = pick("launch_beam_curvature", "curvature", cfg.launch_beam_curvature)
    if pos is None or W is None or c is None or np.size(pos) != 3:
        return None
    R, zeta, Z = map(float, pos)
    return np.array([R * np.cos(zeta), R * np.sin(zeta), Z]), float(W[0]), float(c[0])


def _psi_axes(pxx, pxy, pyy, k0):
    """Widths and curvatures of the principal axes of Im(Psi)."""
    M = np.array([[pxx, pxy], [pxy, pyy]], complex)
    lam, V = np.linalg.eigh(M.imag)
    W = np.sqrt(2.0 / np.clip(lam, 1e-12, None))
    curv = np.array([V[:, i] @ M.real @ V[:, i] for i in range(2)]) / k0
    return W, curv


def front_from_launch(beam: Beam, k0, launch, behind, n_ext):
    """Vacuum beam from ``behind`` m behind the antenna up to (not including)
    Scotty's first point, built from the launch width / curvature (circular).
    Returns the extended beam, or None if Scotty's first point is not on the
    launch line (output of a different launch?)."""
    L0, W, c = launch
    P0 = beam.points[0]
    t = np.real(beam.K[0]) / np.linalg.norm(np.real(beam.K[0]))
    D = float((P0 - L0) @ t)
    off = float(np.linalg.norm(P0 - L0 - D * t))
    if off > 2e-3 + 1e-3 * abs(D) or D < -1e-6:
        print(f"WARNING: Scotty's first point is {off:.3g} m off the launch line (or behind "
              "the antenna): the .h5 was not made with this launch -- using "
              "extend_first_from = 'entry' instead.")
        return None
    q0inv = 1.0 / (k0 * c + 2j / W ** 2)
    s0 = -float(behind)                             # start, measured from the antenna
    if s0 >= D - 1e-9:
        print(f"NOTE: extend_first = {behind:.4g} m would start the beam at or beyond the "
              f"plasma entry ({D:.4g} m from the antenna): the beam starts at the plasma entry.")
    sv = np.linspace(s0, D, n_ext + 1)[:-1] if D - s0 > 1e-9 else np.array([])
    if sv.size and s0 < 0 < D:
        sv = np.union1d(sv, [0.0])                  # a cross-section exactly at the antenna
    Wb, cb = _psi_axes(beam.psi_xx[0], beam.psi_xy[0], beam.psi_yy[0], k0)
    psi_vac = 1.0 / (q0inv + D / k0)
    W_vac = np.sqrt(2 / psi_vac.imag)
    if np.max(np.abs(Wb / W_vac - 1.0)) > 0.15:
        print(f"WARNING: the launch beam (W = {W * 1e3:.4g} mm, curvature {c:.4g} 1/m at the "
              f"antenna) is {W_vac * 1e3:.4g} mm wide at the plasma entry, but Scotty's beam "
              f"there is {Wb[0] * 1e3:.4g} x {Wb[1] * 1e3:.4g} mm: the launch width / curvature "
              "do not match this .h5 (set launch_beam_width / launch_beam_curvature in "
              "section 1 to the values Scotty was run with) -- using extend_first_from = "
              "'entry' instead.")
        return None
    where = (f"{behind:.4g} m behind the antenna" if behind > 0 else "at the antenna"
             if behind == 0 else f"{-behind:.4g} m ahead of the antenna (towards the plasma)")
    print(f"Front of the beam from the LAUNCH definition: antenna -> plasma entry {D:.4g} m, "
          f"starting {where}; width at the antenna {W * 1e3:.4g} mm "
          f"(curvature {c:.4g} 1/m).")
    print(f"  at the plasma entry: vacuum side W = {np.sqrt(2 / psi_vac.imag) * 1e3:.4g} mm, "
          f"1/R = {psi_vac.real / k0:.4g} 1/m (circular) | Scotty (plasma side, after "
          f"Psi_BC): W = {Wb[0] * 1e3:.4g} x {Wb[1] * 1e3:.4g} mm, "
          f"1/R = {cb[0]:.4g} / {cb[1]:.4g} 1/m")
    if sv.size == 0:
        return beam
    n = sv.size
    xh = beam.xhat[0] - (beam.xhat[0] @ t) * t
    xh /= np.linalg.norm(xh)
    yh = np.cross(t, xh)
    if yh @ beam.yhat[0] < 0:
        yh = -yh
    psi = 1.0 / (q0inv + sv / k0)
    front = dict(points=L0 + sv[:, None] * t, xhat=np.tile(xh, (n, 1)), yhat=np.tile(yh, (n, 1)),
                 psi_xx=psi, psi_xy=np.zeros(n, complex), psi_yy=psi,
                 distance=beam.distance[0] - (D - sv),
                 E=np.tile(beam.E[0], (n, 1)) if beam.E is not None else None,
                 K=np.tile(k0 * t, (n, 1)))
    cat = lambda f: (None if getattr(beam, f) is None or front[f] is None else
                     np.concatenate([front[f], getattr(beam, f)], axis=0))
    return Beam(**{f: cat(f) for f in front})


def extend_scotty_beam(beam: Beam, d, k0, cfg) -> Beam:
    """Front (launch definition or plasma-side Psi, see extend_first_from) and
    back (extend_last) vacuum parts of the Scotty beam."""
    if beam.psi_xx is None or beam.K is None or not k0:
        return extend_beam(beam, k0, cfg.extend_first, cfg.extend_last, cfg.extend_points)
    done = False
    if getattr(cfg, "extend_first_from", "entry") == "launch":
        launch = launch_definition(d, cfg)
        if launch is None:
            print("WARNING: launch width / curvature / position unknown (not in the .h5 and "
                  "not in section 1): using extend_first_from = 'entry'.")
        else:
            new = front_from_launch(beam, k0, launch, cfg.extend_first, cfg.extend_points)
            if new is not None:
                beam, done = new, True
    first = 0.0 if done else cfg.extend_first
    if first > 0 or cfg.extend_last > 0:
        beam = extend_beam(beam, k0, first, cfg.extend_last, cfg.extend_points)
    return beam


def report_beam_ends(beam: Beam, k0, d, cfg):
    """Print where the beam starts / ends and its 1/e semi-axes there; check
    that the start is the antenna with the launch width (circular)."""
    if beam.psi_xx is None or not k0:
        return
    print("Beam ends (1/e field semi-axes):")
    for name, i in (("start", 0), ("end  ", -1)):
        P = beam.points[i]
        W, _ = _psi_axes(beam.psi_xx[i], beam.psi_xy[i], beam.psi_yy[i], k0)
        print(f"  {name}: X,Y,Z = {P[0]:.5f}, {P[1]:.5f}, {P[2]:.5f} m  (R = {np.hypot(P[0], P[1]):.5f} m)"
              f"   W = {W[0] * 1e3:.3f} x {W[1] * 1e3:.3f} mm")
    launch = launch_definition(d, cfg)
    if launch is None:
        return
    L0, Wl, cl = launch
    t0 = np.real(beam.K[0]) / np.linalg.norm(np.real(beam.K[0])) if beam.K is not None else None
    s_start = float((beam.points[0] - L0) @ t0) if t0 is not None else 0.0
    off = (float(np.linalg.norm(beam.points[0] - L0 - s_start * t0)) if t0 is not None
           else np.inf)
    if s_start > 5e-3 and off < 5e-3:
        W0, _ = _psi_axes(beam.psi_xx[0], beam.psi_xy[0], beam.psi_yy[0], k0)
        psi_e = 1.0 / (1.0 / (k0 * cl + 2j / Wl ** 2) + s_start / k0)
        print(f"  the beam starts {s_start:.4g} m AHEAD of the antenna (towards the plasma): "
              f"W = {W0[0] * 1e3:.3f} x {W0[1] * 1e3:.3f} mm there (launched beam: "
              f"{np.sqrt(2 / psi_e.imag) * 1e3:.3f} mm, 1/R = {psi_e.real / k0:.4g} 1/m)")
        return
    dist = np.linalg.norm(beam.points - L0, axis=1)
    i = int(np.argmin(dist))
    W, _ = _psi_axes(beam.psi_xx[i], beam.psi_xy[i], beam.psi_yy[i], k0)
    print(f"  at the antenna (beam point {dist[i] * 1e3:.2f} mm from it): "
          f"W = {W[0] * 1e3:.3f} x {W[1] * 1e3:.3f} mm  (launch width {Wl * 1e3:.3f} mm)")
    if dist[i] > 5e-3:
        print("  WARNING: the beam does not pass through the antenna point.")
    elif abs(W[0] - Wl) > 0.01 * Wl or abs(W[1] - Wl) > 0.01 * Wl:
        print("  WARNING: the beam at the antenna is NOT the circular launch beam -- "
              "please send me this output.")


# ##########################################################################
# 5. CROSS-SECTION GEOMETRY
# ##########################################################################

def beam_cross_sections(beam: Beam, fallback_radius):
    """Return (A, B), each (N,3): conjugate semi-diameters of the elliptical
    1/e cross-section at every ray point, i.e.
        P(theta) = C + A cos(theta) + B sin(theta).
    A, B are principal semi-axes from diagonalising Im(Psi) in the local
    (xhat, yhat) frame, then re-phased by make_sections_continuous()."""
    n = beam.n
    if beam.psi_xx is None:
        A = fallback_radius * beam.xhat
        B = fallback_radius * beam.yhat
    else:
        M = np.empty((n, 2, 2))
        M[:, 0, 0] = beam.psi_xx.imag
        M[:, 0, 1] = M[:, 1, 0] = beam.psi_xy.imag
        M[:, 1, 1] = beam.psi_yy.imag
        lam, V = np.linalg.eigh(M)                    # batched, ascending
        if np.any(lam <= 0):
            print(f"WARNING: Im(Psi) not positive definite at {np.sum(np.any(lam <= 0, 1))} "
                  "points (cut-off / numerical noise?) -- width clipped there.")
        W = np.sqrt(2.0 / np.clip(lam, 1e-6, None))   # semi-axes (m)
        # 2D eigenvectors -> 3D: V[:, 0, k]*xhat + V[:, 1, k]*yhat
        ax0 = V[:, 0, 0, None] * beam.xhat + V[:, 1, 0, None] * beam.yhat
        ax1 = V[:, 0, 1, None] * beam.xhat + V[:, 1, 1, None] * beam.yhat
        A, B = W[:, 0, None] * ax0, W[:, 1, None] * ax1
    return make_sections_continuous(beam.points, A, B)


def make_sections_continuous(C, A, B):
    """Remove the arbitrary sign/rotation of each section's (A, B) pair.

    1. orientation: theta must turn the same way about the ray everywhere
       (gives outward STL normals). The first section is oriented so that
       A x B points along the ray; every following one so that its A x B
       agrees with the previous section's. (Comparing with neighbours rather
       than with a finite-difference tangent stays reliable where the ray
       bends sharply, e.g. at a reflection.)
    2. phase: rotate (A, B) -> (A cos a + B sin a, -A sin a + B cos a), the
       same ellipse with theta shifted by a, choosing a to best match the
       previous section (a 2x2 orthogonal Procrustes problem)."""
    A, B = A.copy(), B.copy()
    if np.dot(np.cross(A[0], B[0]), C[1] - C[0]) < 0:
        B[0] *= -1.0
    for i in range(1, len(C)):
        if np.dot(np.cross(A[i], B[i]), np.cross(A[i - 1], B[i - 1])) < 0:
            B[i] *= -1.0
        L = np.stack([A[i], B[i]], axis=1)            # 3x2
        Lp = np.stack([A[i - 1], B[i - 1]], axis=1)
        M = L.T @ Lp
        a = np.arctan2(M[1, 0] - M[0, 1], M[0, 0] + M[1, 1])
        c, s = np.cos(a), np.sin(a)
        A[i], B[i] = c * L[:, 0] + s * L[:, 1], -s * L[:, 0] + c * L[:, 1]
    return A, B


def ellipse_rings(C, A, B, theta):
    """(N, len(theta), 3) points on every cross-section."""
    return (C[:, None, :] + np.cos(theta)[None, :, None] * A[:, None, :]
            + np.sin(theta)[None, :, None] * B[:, None, :])


def section_widths(A, B):
    """Principal semi-axes (min, max) of each section (from A, B)."""
    G = np.stack([np.einsum("ij,ij->i", A, A), np.einsum("ij,ij->i", A, B),
                  np.einsum("ij,ij->i", B, B)], axis=1)
    tr, det = G[:, 0] + G[:, 2], G[:, 0] * G[:, 2] - G[:, 1] ** 2
    disc = np.sqrt(np.clip(tr ** 2 / 4 - det, 0, None))
    return np.sqrt(np.clip(tr / 2 - disc, 0, None)), np.sqrt(tr / 2 + disc)


def ring_selection(n, stride):
    idx = np.arange(0, n, max(stride, 1))
    return idx if idx[-1] == n - 1 else np.append(idx, n - 1)


# ##########################################################################
# 6. STL EXPORT
# ##########################################################################

def build_tube_mesh(rings, caps=True):
    """Triangulate rings (M, T, 3) into a tube, optionally capped by fans.
    Returns vertices (V,3) and faces (F,3); normals point outwards."""
    M, T, _ = rings.shape
    verts = rings.reshape(-1, 3)
    j = np.arange(T)
    jn = (j + 1) % T
    m = np.arange(M - 1)[:, None]
    a, b = m * T + j, m * T + jn
    c, d = a + T, b + T
    faces = [np.stack([a, b, d], -1).reshape(-1, 3),
             np.stack([a, d, c], -1).reshape(-1, 3)]
    if caps:
        c0, c1 = M * T, M * T + 1
        verts = np.vstack([verts, rings[0].mean(0), rings[-1].mean(0)])
        last = (M - 1) * T
        faces.append(np.stack([np.full(T, c0), jn, j], -1))
        faces.append(np.stack([np.full(T, c1), last + j, last + jn], -1))
    return verts, np.concatenate(faces).astype(np.int64)


def write_stl(path, verts, faces, ascii=False, name="scotty_beam"):
    """Binary (default) or ASCII STL, numpy only."""
    v = verts[faces]                                         # (F,3,3)
    nrm = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = nrm / np.where(ln == 0, 1, ln)
    if ascii:
        with open(path, "w") as f:
            f.write(f"solid {name}\n")
            for n_, tri in zip(nrm, v):
                f.write(f" facet normal {n_[0]:.6e} {n_[1]:.6e} {n_[2]:.6e}\n  outer loop\n")
                for p in tri:
                    f.write(f"   vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n")
                f.write("  endloop\n endfacet\n")
            f.write(f"endsolid {name}\n")
        return
    rec = np.zeros(len(faces), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = nrm, v
    with open(path, "wb") as f:
        f.write(b"Scotty beam envelope surface".ljust(80, b" "))
        f.write(np.uint32(len(faces)).tobytes())
        rec.tofile(f)


# ##########################################################################
# 7. NURBS EXPORT (IGES)
# ##########################################################################

# Exact unit circle as a rational quadratic NURBS with 9 control points
# (4 quarter arcs). Mapping (p, q) -> C + p A + q B turns it into the
# section ellipse (NURBS are invariant under affine maps).
_W45 = np.sqrt(0.5)
CIRCLE_PQ = np.array([[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0],
                      [-1, -1], [0, -1], [1, -1], [1, 0]], float)
CIRCLE_W = np.array([1, _W45, 1, _W45, 1, _W45, 1, _W45, 1])
CIRCLE_KNOTS = np.array([0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1])

# Planar elliptical end cap as ONE non-degenerate rational biquadratic
# patch whose 4 boundary edges are exactly the 4 quarter arcs above
# (indices into CIRCLE_PQ), so GiD can glue it to the tube patches.
CAP_IDX = np.array([[0, 7, 6],      # CAP_IDX[i, j]: u index i, v index j
                    [1, -1, 5],     # -1 = centre point
                    [2, 3, 4]])
CAP_W = np.outer([1, _W45, 1], [1, _W45, 1])


def section_control_points(C, A, B):
    """(N, 9, 3) control points of each section ellipse."""
    return (C[:, None, :] + CIRCLE_PQ[None, :, 0, None] * A[:, None, :]
            + CIRCLE_PQ[None, :, 1, None] * B[:, None, :])


def cap_control_points(C, A, B):
    """(3, 3, 3) control net of the planar cap of one section."""
    sec = section_control_points(C[None], A[None], B[None])[0]
    net = np.empty((3, 3, 3))
    for i in range(3):
        for j in range(3):
            k = CAP_IDX[i, j]
            net[i, j] = C if k < 0 else sec[k]
    return net


def bspline_basis(knots, p, x, n_ctrl):
    """Matrix N[k, i] = N_{i,p}(x_k) (Cox-de Boor), shape (len(x), n_ctrl).
    The right end of the knot range belongs to the last span."""
    x = np.atleast_1d(np.asarray(x, float))
    U = np.asarray(knots, float)
    N = np.zeros((len(x), len(U) - 1))
    for i in range(len(U) - 1):
        if U[i] < U[i + 1]:
            N[:, i] = (x >= U[i]) & (x < U[i + 1])
    last = np.nonzero(U[:-1] < U[1:])[0][-1]
    N[x >= U[-1], last] = 1.0
    for k in range(1, p + 1):
        Nn = np.zeros((len(x), len(U) - 1 - k))
        for i in range(len(U) - 1 - k):
            d1, d2 = U[i + k] - U[i], U[i + k + 1] - U[i + 1]
            if d1 > 0:
                Nn[:, i] += (x - U[i]) / d1 * N[:, i]
            if d2 > 0:
                Nn[:, i] += (U[i + k + 1] - x) / d2 * N[:, i + 1]
        N = Nn
    return N[:, :n_ctrl]


def interpolate_columns(Q, v):
    """Cubic (or lower) B-spline interpolation along the beam.

    Q : (S, ...) data at parameters v (S,), increasing in [0, 1].
    Returns (knots, degree, P) with P the control points (same shape as Q).
    The weights are constant along v, so interpolating the Cartesian
    control points directly keeps every section exact."""
    S = len(v)
    p = min(3, S - 1)
    inner = [np.mean(v[j:j + p]) for j in range(1, S - p)]
    knots = np.concatenate([np.zeros(p + 1), inner, np.ones(p + 1)])
    N = bspline_basis(knots, p, v, S)
    P = np.linalg.solve(N, Q.reshape(S, -1)).reshape(Q.shape)
    return knots, p, P


def nurbs_section_error(P, knots, p, v_rings, C, A, B, n_u=48):
    """Max radial distance (m) between the skinned surface and each exact
    cross-section, measured in that section's plane.

    For ring r, the surface iso-curve at v_r is a rational curve with
    control points  sum_j N_j(v_r) P[j]. Each curve point X is written as
    C + alpha A + beta B (+ out-of-plane part, which is just a slide along
    the surface and is ignored); its distance to the ellipse along the same
    radial line is |X_inplane - C| * |1 - 1/rho|, rho = sqrt(alpha^2+beta^2)."""
    Nv = bspline_basis(knots, p, v_rings, P.shape[0])            # (R, S)
    Pv = np.einsum("rs,sik->rik", Nv, P)                          # (R, 9, 3)
    u = np.linspace(0, 1, n_u, endpoint=False)
    Nu = bspline_basis(CIRCLE_KNOTS, 2, u, 9) * CIRCLE_W          # (U, 9)
    X = np.einsum("ui,rik->ruk", Nu, Pv) / Nu.sum(1)[None, :, None]
    L = np.stack([A, B], axis=2)                                  # (R,3,2)
    coef = np.einsum("rab,rub->rua", np.linalg.pinv(L), X - C[:, None, :])  # (R,U,2)
    inplane = np.einsum("rkb,rub->ruk", L, coef)
    rho = np.linalg.norm(coef, axis=2)
    err = np.linalg.norm(inplane, axis=2) * np.abs(1 - 1 / np.maximum(rho, 1e-12))
    return err.max(axis=1)


def fit_beam_nurbs(C, A, B, tol, n_init=30, n_max=600):
    """Choose cross-sections adaptively and skin a NURBS surface through them.

    The along-beam parameter of every Scotty point is its normalised arc
    length. Starting from ``n_init`` sections evenly spaced in arc length,
    the worst-fitted Scotty section in every gap is added until all sections
    are reproduced within ``tol`` (m) or ``n_max`` sections are used.
    Returns (knots_v, degree_v, P (S,9,3), selected indices, errors)."""
    n = len(C)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))])
    v_all = s / s[-1]
    sec = section_control_points(C, A, B)

    sel = np.unique(np.searchsorted(v_all, np.linspace(0, 1, min(n_init, n))).clip(0, n - 1))
    sel = np.unique(np.concatenate([[0], sel, [n - 1]]))
    while True:
        knots, p, P = interpolate_columns(sec[sel], v_all[sel])
        err = nurbs_section_error(P, knots, p, v_all, C, A, B)
        if err.max() <= tol or len(sel) >= min(n_max, n):
            break
        add = []
        for a, b in zip(sel[:-1], sel[1:]):
            if b - a > 1:
                k = a + 1 + np.argmax(err[a + 1:b])
                if err[k] > tol:
                    add.append(k)
        if not add:
            break
        sel = np.unique(np.concatenate([sel, add]))
    return knots, p, P, sel, err


def nurbs_beam_patches(C, A, B, tol, n_patches=4, caps=True, n_max=600, unit="m"):
    """Build the list of NURBS patches (dicts) describing the beam envelope.

    n_patches = 4 : the tube is split along the 4 quarter arcs (no closed /
                    seam surfaces -- safest for GiD meshing);
    n_patches = 1 : a single surface closed in u."""
    knots_v, p_v, P, sel, err = fit_beam_nurbs(C, A, B, tol, n_max=n_max)
    wmin, wmax = section_widths(A, B)
    k = np.argmax(err)
    print(f"NURBS: {len(sel)} cross-sections skinned (of {len(C)}); max deviation "
          f"{err[k]:.3g} {unit} at point {k} (local semi-axes {wmin[k]:.3g}-{wmax[k]:.3g} {unit}).")
    if err[k] > tol:
        print(f"WARNING: tolerance {tol:g} {unit} not reached -- raise --nurbs-max-sections "
              "or relax --nurbs-tol.")

    patches = []
    # P[j, i, :]: section j (along beam, v), circle control point i (u)
    if n_patches == 1:
        patches.append(dict(ctrl=P.transpose(1, 0, 2), w=np.repeat(CIRCLE_W[:, None], len(P), 1),
                            ku=CIRCLE_KNOTS, pu=2, kv=knots_v, pv=p_v, closed_u=True,
                            label="BEAMTUBE"))
    else:
        for q in range(4):
            cols = slice(2 * q, 2 * q + 3)
            patches.append(dict(ctrl=P[:, cols].transpose(1, 0, 2),
                                w=np.repeat(CIRCLE_W[cols, None], len(P), 1),
                                ku=np.array([0, 0, 0, 1, 1, 1.]), pu=2, kv=knots_v, pv=p_v,
                                closed_u=False, label=f"BEAMQ{q + 1}"))
    if caps:
        for end, i in (("CAPFIRST", 0), ("CAPLAST", -1)):
            patches.append(dict(ctrl=cap_control_points(C[i], A[i], B[i]), w=CAP_W,
                                ku=np.array([0, 0, 0, 1, 1, 1.]), pu=2,
                                kv=np.array([0, 0, 0, 1, 1, 1.]), pv=2,
                                closed_u=False, label=end))
    return patches


class IgesWriter:
    """Minimal IGES 5.3 writer: rational B-spline surfaces (entity 128) and
    curves (entity 126). Fixed 80-column format: S, G, D, P, T sections."""

    def __init__(self, unit="m"):
        self.unit = unit
        self.entities = []   # (type, form, label, [params], status)
        self.max_coord = 0.0

    @staticmethod
    def _r(x):
        return f"{float(x):.12E}"          # always contains a decimal point

    def _add(self, etype, form, label, prm, status="00000000"):
        """Append an entity; return its directory-entry (DE) number."""
        self.entities.append((etype, form, label, prm, status))
        return 2 * len(self.entities) - 1

    def add_surface(self, ctrl, w, ku, pu, kv, pv, closed_u=False, label="SURF",
                    status="00000000"):
        """ctrl (nu, nv, 3), w (nu, nv); u varies fastest in IGES order."""
        nu, nv = w.shape
        prm = [128, nu - 1, nv - 1, pu, pv, int(closed_u), 0,
               0 if np.ptp(w) > 0 else 1, 0, 0]
        prm += [self._r(x) for x in ku] + [self._r(x) for x in kv]
        prm += [self._r(w[i, j]) for j in range(nv) for i in range(nu)]
        prm += [self._r(c) for j in range(nv) for i in range(nu) for c in ctrl[i, j]]
        prm += [self._r(ku[pu]), self._r(ku[-pu - 1]), self._r(kv[pv]), self._r(kv[-pv - 1])]
        de = self._add(128, 0, label, prm, status)
        self.max_coord = max(self.max_coord, float(np.abs(ctrl).max()))
        return de

    def add_polyline_curve(self, pts, label="RAY"):
        """Exact polyline as a degree-1 B-spline (entity 126)."""
        n = len(pts)
        t = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        t /= t[-1]
        knots = np.concatenate([[0], t, [1]])
        prm = [126, n - 1, 1, 0, 0, 1, 0]
        prm += [self._r(x) for x in knots] + [self._r(1.0)] * n
        prm += [self._r(c) for p in pts for c in p]
        prm += [self._r(0), self._r(1), self._r(0), self._r(0), self._r(0)]
        self._add(126, 0, label, prm)
        self.max_coord = max(self.max_coord, float(np.abs(pts).max()))

    def add_planar_face(self, poly, label="FACE"):
        """Convex planar polygon (k,3) as a trimmed surface (entity 144): a
        bilinear plane patch (128) over the polygon's bounding rectangle,
        trimmed by the polygon given both in model space (102 of 110 lines)
        and in the patch's (u, v) parameter space; tied by entity 142."""
        poly = np.asarray(poly, float)
        c = poly.mean(0)
        nrm = np.cross(poly[1] - poly[0], poly[2] - poly[0])
        for i in range(2, len(poly) - 1):              # robust normal
            cand = np.cross(poly[i] - poly[0], poly[i + 1] - poly[0])
            if np.linalg.norm(cand) > np.linalg.norm(nrm):
                nrm = cand
        nrm /= np.linalg.norm(nrm)
        e1 = poly[1] - poly[0]
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(nrm, e1)
        a, b = (poly - c) @ e1, (poly - c) @ e2
        pad = 0.02 * max(np.ptp(a), np.ptp(b))
        a0, a1, b0, b1 = a.min() - pad, a.max() + pad, b.min() - pad, b.max() + pad
        corner = lambda x, y: c + x * e1 + y * e2
        ctrl = np.array([[corner(a0, b0), corner(a0, b1)],
                         [corner(a1, b0), corner(a1, b1)]])        # [u][v]
        k01 = np.array([0, 0, 1, 1.])
        srf = self.add_surface(ctrl, np.ones((2, 2)), k01, 1, k01, 1, label=label,
                               status="00010000")
        uv = np.stack([(a - a0) / (a1 - a0), (b - b0) / (b1 - b0), np.zeros(len(a))], 1)
        n = len(poly)
        model, param = [], []
        for i in range(n):
            j = (i + 1) % n
            model.append(self._add(110, 0, label, [110] + [self._r(x) for x in
                                                         (*poly[i], *poly[j])], "00010000"))
        for i in range(n):
            j = (i + 1) % n
            param.append(self._add(110, 0, label, [110] + [self._r(x) for x in
                                                         (*uv[i], *uv[j])], "00010500"))
        cm = self._add(102, 0, label, [102, n] + model, "00010000")
        cp = self._add(102, 0, label, [102, n] + param, "00010500")
        cos = self._add(142, 0, label, [142, 0, srf, cp, cm, 3], "00010000")
        self.max_coord = max(self.max_coord, float(np.abs(poly).max()))
        return self._add(144, 0, label, [144, srf, 1, 0, cos])

    @staticmethod
    def _hollerith(s):
        return f"{len(s)}H{s}"

    def write(self, path, description="Scotty beam envelope"):
        unit_flag = {"m": (6, "M"), "mm": (2, "MM"), "cm": (10, "CM")}[self.unit]
        now = datetime.datetime.now().strftime("%Y%m%d.%H%M%S")
        H = self._hollerith
        g = [H(","), H(";"), H("scotty_beam_to_gid"), H(os.path.basename(path)),
             H("scotty_beam_to_gid.py"), H("1.0"), "32", "38", "6", "308", "15",
             H("GiD"), "1.0", str(unit_flag[0]), H(unit_flag[1]), "1", "1.0",
             H(now), "1.0E-9", self._r(max(self.max_coord, 1.0)), H("Scotty"),
             H("-"), "11", "0", H(now)]
        g_lines = self._wrap(g, 72, terminator=";")
        s_lines = [description[:72]]

        d_lines, p_lines = [], []
        for k, (etype, form, label, prm, status) in enumerate(self.entities):
            de_seq = 2 * k + 1
            lines = self._wrap([str(x) for x in prm], 64, terminator=";")
            p_start = len(p_lines) + 1
            p_lines += [f"{ln:<64}{de_seq:>8}" for ln in lines]
            d_lines.append(f"{etype:>8}{p_start:>8}{0:>8}{0:>8}{0:>8}{0:>8}{0:>8}{0:>8}{status:>8}")
            d_lines.append(f"{etype:>8}{0:>8}{0:>8}{len(lines):>8}{form:>8}{'':>8}{'':>8}"
                           f"{label[:8]:>8}{k + 1:>8}")

        with open(path, "w", newline="\n") as f:
            for tag, lines in (("S", s_lines), ("G", g_lines), ("D", d_lines), ("P", p_lines)):
                for i, ln in enumerate(lines, 1):
                    f.write(f"{ln:<72}{tag}{i:>7}\n")
            t = f"S{len(s_lines):>7}G{len(g_lines):>7}D{len(d_lines):>7}P{len(p_lines):>7}"
            f.write(f"{t:<72}T{1:>7}\n")

    @staticmethod
    def _wrap(tokens, width, terminator):
        """Join tokens with ',' (last one followed by the terminator) into
        lines of at most ``width`` characters, never splitting a token."""
        lines, cur = [], ""
        for i, tok in enumerate(tokens):
            piece = tok + ("," if i < len(tokens) - 1 else terminator)
            if len(cur) + len(piece) > width:
                lines.append(cur)
                cur = ""
            cur += piece
        lines.append(cur)
        return lines


# ##########################################################################
# 8. TEXT EXPORTS
# ##########################################################################

def _fmt_c(z, tol=1e-9):
    """Real if the imaginary part is negligible, else Python 'a+bj'."""
    z = complex(z)
    if abs(z.imag) < tol * max(abs(z.real), 1.0):
        return f"{z.real:.8g}"
    return f"{z.real:.8g}{'+' if z.imag >= 0 else '-'}{abs(z.imag):.8g}j"


def write_central_ray_file(path, pts):
    np.savetxt(path, pts, fmt="%.8g", delimiter=",", header="x,y,z", comments="")


def write_vector_rows(path, arr, complex_values):
    with open(path, "w") as f:
        for row in arr:
            f.write(",".join(_fmt_c(c) if complex_values else f"{np.real(c):.8g}"
                             for c in row) + "\n")


def write_efield_endpoints(path, beam: Beam, scotty_first, scotty_last, scale):
    """E (complex) and k (real, Cartesian) at the first and last cross-section.
    In vacuum both are unchanged by the extension, so they are valid at the
    extended ends too; both positions are reported."""
    p = lambda v: ",".join(f"{c * scale:.8g}" for c in v)
    with open(path, "w") as f:
        for name, i, orig in (("First", 0, scotty_first), ("Last", -1, scotty_last)):
            f.write(f"# {name} surface  (position X,Y,Z = {p(beam.points[i])}; "
                    f"Scotty end point = {p(orig)})\n")
            f.write(",".join(_fmt_c(c) for c in beam.E[i]) + "\n")
            f.write(",".join(f"{np.real(c):.8g}" for c in beam.K[i]) + "\n")


# ##########################################################################
# 9. INTERACTIVE PLOT
# ##########################################################################

def _draw_beam_surface(ax, rings_fn, C, A, B, n_around, colour="C1"):
    """Beam surface through ALL cross-sections (no subsampling)."""
    th = np.linspace(0.0, 2.0 * np.pi, max(int(n_around), 8) + 1)   # closed rings
    R = rings_fn(C, A, B, th)                                          # (N, M+1, 3)
    ax.plot_surface(R[..., 0], R[..., 1], R[..., 2], rcount=R.shape[0], ccount=R.shape[1],
                    color=colour, alpha=0.55, linewidth=0, antialiased=False, shade=True)
    for i in (0, -1):                                                  # end sections
        ax.plot(*R[i].T, color="k", lw=0.8)
    return R


def _equal_axes(ax, pts):
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    c, r = 0.5 * (lo + hi), 0.5 * np.max(hi - lo)
    ax.set_xlim3d(c[0] - r, c[0] + r)
    ax.set_ylim3d(c[1] - r, c[1] + r)
    ax.set_zlim3d(c[2] - r, c[2] + r)
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass


def beam_plot(beam: Beam, A, B, n_around, envelope_faces=None):
    """Static 3D view: beam surface through every cross-section, central ray
    and (if built) the FEM envelope. The window can still be rotated/zoomed."""
    import matplotlib.pyplot as plt
    C = beam.points
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    R = _draw_beam_surface(ax, ellipse_rings, C, A, B, n_around)
    ax.plot(*C.T, color="k", lw=1.0)
    pts = [R.reshape(-1, 3)]
    if envelope_faces:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        colours = {"first": (0.2, 0.7, 0.2, 0.25), "last": (0.8, 0.2, 0.2, 0.25),
                   "lateral": (1.0, 0.6, 0.0, 0.06)}
        ax.add_collection3d(Poly3DCollection(
            [P for P, _ in envelope_faces], edgecolors="0.35", linewidths=0.4,
            facecolors=[colours.get(k, (1.0, 0.6, 0.0, 0.06)) for _, k in envelope_faces]))
        pts.append(np.vstack([P for P, _ in envelope_faces]))
    _equal_axes(ax, np.vstack(pts))
    ax.set_xlabel("X (tokamak R) [m]")
    ax.set_ylabel("Y (toroidal) [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title(f"Beam surface: {beam.n} cross-sections, ray length "
                 f"{beam.distance[-1] - beam.distance[0]:.3g} m")
    ax.view_init(elev=20, azim=-60)
    plt.tight_layout()
    plt.show()


# ##########################################################################
# 10. FEM ENVELOPE: CLOSED POLYHEDRON (PLANAR FACES) AROUND THE BEAM
# ##########################################################################
#
# The envelope is a chain of convex "cells" along the ray. Cell k lies
# between two joint planes (the planes of the Scotty cross-sections at two
# chosen ray points) and is closed laterally by ``envelope_sides`` planes:
#
#   * cell 0: the lateral planes start as a prism around the beam (normals
#     spread evenly around the first cross-section, aligned with its major
#     axis) and are then tilted to touch the beam;
#   * cell k > 0: lateral plane j is HINGED on the edge that face j of cell
#     k-1 left on the shared joint plane, and rotated about that edge until it
#     touches the beam. Neighbouring cells therefore share their joint polygon
#     exactly and every face is planar by construction.
#
# "Touches" means: every point of the beam surface in the cell is at a
# distance >= separation from the plane (support function of beam (+) ball).
# The beam is sampled with polygons CIRCUMSCRIBED to its elliptical sections,
# so the sampled hull contains the exact ellipses.
#
# Extra support planes (26 directions) are added where the hinged planes
# alone would leave a cell open, e.g. around a reflection.
#
# The first and last faces are the planes of the first and last beam
# cross-sections, i.e. they contain the beam's end caps.
#
# Joints between cells: first only at the bends of the central ray
# (Douglas-Peucker), merging cells until the chain is valid; then long cells
# are split (envelope_max_cell_length) and every added joint that spoils the
# chain is dropped again. "Valid" = cells closed, joint polygons matching,
# non-neighbouring cells disjoint, and every sampled beam point inside the
# envelope at >= separation from every side face (checked explicitly).

def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _plane_frame(n):
    """Two unit vectors spanning the plane with normal n."""
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    u = _unit(np.cross(n, a))
    return u, np.cross(n, u)


def _clip_polygon(poly, labels, n, h, label, eps):
    """Sutherland-Hodgman: keep the part of the planar polygon ``poly`` (list
    of 3D points; edge i = poly[i] -> poly[i+1] carries labels[i]) with
    n.x <= h. The new edge lying on the cutting plane gets ``label``."""
    m = len(poly)
    dist = [float(np.dot(n, p) - h) for p in poly]
    out, lab = [], []
    for i in range(m):
        p, q = poly[i], poly[(i + 1) % m]
        dp, dq = dist[i], dist[(i + 1) % m]
        if dp <= eps:
            out.append(p)
            lab.append(labels[i])
            if dq > eps:                                   # leaving
                t = dp / (dp - dq)
                out.append(p + t * (q - p))
                lab.append(label)
        elif dq <= eps:                                    # entering
            t = dp / (dp - dq)
            out.append(p + t * (q - p))
            lab.append(labels[i])
    keep_p, keep_l = [], []
    for p, l in zip(out, lab):                             # merge duplicates
        if keep_p and np.linalg.norm(p - keep_p[-1]) <= eps:
            keep_l[-1] = l
            continue
        keep_p.append(p)
        keep_l.append(l)
    while len(keep_p) > 1 and np.linalg.norm(keep_p[0] - keep_p[-1]) <= eps:
        keep_p.pop()
        keep_l.pop()
    if len(keep_p) < 3:
        return [], []
    return keep_p, keep_l


def _convex_face(planes, i, size, eps):
    """Face of the convex cell {n.x <= h for all planes} lying on plane i.
    Returns (vertices (k,3) ordered counter-clockwise seen from outside,
    labels: index of the neighbouring plane along each edge)."""
    n, h = planes[i]
    c = n * h
    u, v = _plane_frame(n)
    sq = [c + size * (su * u + sv * v) for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    poly, lab = sq, [-1] * 4                       # CCW about +n (outward)
    for j, (nj, hj) in enumerate(planes):
        if j != i:
            poly, lab = _clip_polygon(poly, lab, nj, hj, j, eps)
            if not poly:
                return np.zeros((0, 3)), []
    return np.array(poly), lab


def _poly_area(P, n):
    if len(P) < 3:
        return 0.0
    return 0.5 * float(np.dot(np.sum(np.cross(P, np.roll(P, -1, 0)), 0), n))


def _hinge_plane(q, e, r, P, d):
    """Plane through the line (q, e), outward normal as close as possible to
    r, touching the point set P with clearance d (n.(p-q) <= -d for all p).
    Returns (n, h) or None if impossible (line too close to the beam)."""
    n0 = _unit(r - np.dot(r, e) * e)
    m = np.cross(e, n0)
    D = P - q
    a, b = D @ n0, D @ m
    R = np.hypot(a, b)
    if np.any(R < d * (1 - 1e-9)):
        return None
    # each point allows the arc of normal angles  |phi - c| <= w
    c = np.arctan2(b, a) + np.pi
    w = np.pi - np.arccos(np.clip(-d / R, -1, 1))
    k = np.argmin(w)
    cr = (c - c[k] + np.pi) % (2 * np.pi) - np.pi
    lo, hi = np.max(cr - w), np.min(cr + w)
    if lo > hi:
        return None
    phi = max((lo + c[k], hi + c[k]), key=np.cos)    # tight side nearest r
    n = np.cos(phi) * n0 + np.sin(phi) * m
    return n, float(np.dot(n, q))


def _seg_simplify(C, tol, idx_lo, idx_hi, out):
    """Douglas-Peucker on the central ray (indices kept in ``out``)."""
    if idx_hi - idx_lo < 2:
        return
    a, b = C[idx_lo], C[idx_hi]
    ab = b - a
    L = np.linalg.norm(ab)
    seg = C[idx_lo + 1:idx_hi] - a
    if L > 0:
        dist = np.linalg.norm(np.cross(seg, ab / L), axis=1)
    else:
        dist = np.linalg.norm(seg, axis=1)
    k = int(np.argmax(dist))
    if dist[k] > tol:
        mid = idx_lo + 1 + k
        out.add(mid)
        _seg_simplify(C, tol, idx_lo, mid, out)
        _seg_simplify(C, tol, mid, idx_hi, out)


def envelope_candidate_joints(C, s, ray_tol, max_len, min_len):
    """Ray indices where a cell may start/end: Douglas-Peucker points of the
    central ray, long chords subdivided to <= max_len, points closer than
    min_len (along the ray) to the previous one dropped."""
    n = len(C)
    keep = {0, n - 1}
    _seg_simplify(C, ray_tol, 0, n - 1, keep)
    idx = sorted(keep)
    full = [idx[0]]
    for a, b in zip(idx[:-1], idx[1:]):
        nsub = int(np.ceil((s[b] - s[a]) / max_len)) if max_len > 0 else 1
        for t in range(1, nsub):
            full.append(int(np.searchsorted(s, s[a] + t * (s[b] - s[a]) / nsub)))
        full.append(b)
    out = [full[0]]
    for i in full[1:-1]:
        if s[i] - s[out[-1]] >= min_len and s[-1] - s[i] >= min_len:
            out.append(i)
    out.append(full[-1])
    return sorted(set(out))


_KDOP_DIRS = [np.array(v, float) / np.linalg.norm(v)
              for v in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).reshape(3, -1).T
              if np.any(v)]


def _cell_vertices_from_planes(planes, size, eps):
    """Vertices of the (possibly bounded by the big box) convex cell."""
    pts = [_convex_face(planes, i, size, eps)[0] for i in range(len(planes))]
    pts = [p for p in pts if len(p)]
    return np.vstack(pts) if pts else np.zeros((0, 3))


class EnvelopeError(RuntimeError):
    pass


def _build_cells(joints, C, A, B, nrm, s, samples, ring_of, d, n_sides, reach, size, eps):
    """Build the chain of cells for the given joint indices.
    Returns (cells, problems). Each cell: dict(planes=[(n,h)...] with
    planes[0]=entry, planes[1]=exit, faces=[...], lat=[ids], span=(i0,i1))."""
    cells, problems = [], []
    flat = samples.reshape(-1, 3)
    prev_exit = None
    for k in range(len(joints) - 1):
        i0, i1 = joints[k], joints[k + 1]
        n_in, n_out = nrm[i0], nrm[i1]
        h_in, h_out = float(np.dot(-n_in, C[i0])), float(np.dot(n_out, C[i1]))
        # beam points that this cell has to keep at distance >= d
        # (only rings of this cell and of its two neighbours, within reach
        # along the ray; any other part of the beam that comes close is
        # caught by the global checks in build_envelope)
        lo = max(np.searchsorted(s, s[i0] - reach[i0]), joints[max(k - 1, 0)])
        hi = min(np.searchsorted(s, s[i1] + reach[i1], side="right") - 1,
                 joints[min(k + 2, len(joints) - 1)])
        sel = (ring_of >= lo) & (ring_of <= hi)
        P, Pr = flat[sel], ring_of[sel]
        keep = (P @ -n_in <= h_in + d) & (P @ n_out <= h_out + d)
        P, Pr = P[keep], Pr[keep]
        if len(P) == 0:
            problems.append(("empty", k, k))
            P = C[i0:i1 + 1]
        lat = []
        if k == 0:
            # prism around the first section, aligned with its major axis
            ax = A[i0] if np.linalg.norm(A[i0]) >= np.linalg.norm(B[i0]) else B[i0]
            u = _unit(ax - np.dot(ax, n_in) * n_in)
            v = np.cross(n_in, u)
            ang = 2 * np.pi * np.arange(n_sides) / n_sides
            refs = [np.cos(t) * u + np.sin(t) * v for t in ang]
            # entry polygon: tight around the beam near the start; a face whose
            # hinge would then be infeasible falls back to the full support
            Ploc = P[np.abs(P @ n_in - np.dot(n_in, C[i0])) <= d + 1e-12]
            Ploc = Ploc if len(Ploc) else P
            hinges = []
            for r in refs:
                e = _unit(np.cross(n_in, r))
                for Ps in (Ploc, P):
                    hr = float(np.max(Ps @ r)) + d
                    q = np.linalg.lstsq(np.array([n_in, r]),
                                        np.array([np.dot(n_in, C[i0]), hr]), rcond=None)[0]
                    if _hinge_plane(q, e, r, P, d) is not None:
                        break
                hinges.append((q, e, r))
        else:
            hinges = []
            pc = prev_exit
            for j in pc["exit_lat"]:
                r = pc["planes"][j][0]
                e = np.cross(r, n_in)
                if np.linalg.norm(e) < 1e-9:
                    problems.append(("hinge", k, k))
                    continue
                e = _unit(e)
                M = np.array([n_in, r])
                q = np.linalg.lstsq(M, np.array([np.dot(n_in, C[i0]), pc["planes"][j][1]]),
                                    rcond=None)[0]
                hinges.append((q, e, r))
        planes = [(-n_in, h_in), (n_out, h_out)]
        # Points behind the entry plane (within d) are already enclosed by the
        # previous cell; for them only the distance to the FACE (a half-plane
        # starting at the hinge) matters, not the distance to the whole plane.
        behind = P @ n_in < np.dot(n_in, C[i0])
        P_front, P_back = P[~behind], P[behind]
        R_front = Pr[~behind]
        if len(P_front) == 0:
            P_front, R_front = P, Pr
        ok = True
        for q, e, r in hinges:
            pl = _hinge_plane(q, e, r, P_front, d)
            if pl is not None and len(P_back):
                nn = pl[0]
                f = np.cross(nn, e)
                f = f if np.dot(f, n_in) >= 0 else -f
                D = P_back - q
                over = D @ f > 0
                if np.any(D[over] @ nn > -d * (1 - 1e-9)):
                    pl = _hinge_plane(q, e, r, P, d)
            if pl is None:
                ok = False
                break
            planes.append(pl)
        if not ok:
            # Is it the cell's own stretch of beam, or another part of the
            # beam passing through this cell (sharp bend: the two legs of the
            # beam overlap)? In the latter case merge up to that part.
            own = (R_front >= i0) & (R_front <= i1)
            if np.any(~own) and all(_hinge_plane(q, e, r, P_front[own], d) is not None
                                    for q, e, r in hinges):
                far = R_front[~own]
                ring = far.max() if far.max() > i1 else far.min()
                kk = int(np.clip(np.searchsorted(joints, ring, side="right") - 1,
                                 0, len(joints) - 2))
                problems.append(("overlap", min(k, kk), max(k, kk)))
            else:
                problems.append(("hinge", k, k))
            cells.append(None)
            return cells, problems
        # extra support planes (26 directions in the cell frame) touching the
        # beam with clearance d. They close cells around sharp bends, where
        # the hinged planes alone leave the cell open. A plane is used only
        # if it does not cut the entry polygon (shared with the previous
        # cell) and actually cuts the cell built so far.
        chord = C[i1] - C[i0]
        fu = _unit(chord) if np.linalg.norm(chord) > eps else n_in
        fv = _unit(np.cross(fu, A[i0])) if np.linalg.norm(np.cross(fu, A[i0])) > 1e-9 \
            else _plane_frame(fu)[0]
        fw = np.cross(fu, fv)
        entry_poly = _convex_face(planes, 0, size, eps)[0]
        corner = _cell_vertices_from_planes(planes, size, eps)
        for g in _KDOP_DIRS:
            r = g[0] * fu + g[1] * fv + g[2] * fw
            hr = float(np.max(P @ r)) + d
            if len(entry_poly) and np.any(entry_poly @ r > hr + 10 * eps):
                continue
            if len(corner) and np.all(corner @ r <= hr + 10 * eps):
                continue
            planes.append((r, hr))
            corner = _cell_vertices_from_planes(planes, size, eps)
        faces = [_convex_face(planes, i, size, eps) for i in range(len(planes))]
        # lateral planes that still have an edge on the exit plane
        ex_poly, ex_lab = faces[1]
        exit_lat = []
        for t, l in enumerate(ex_lab):
            if l >= 2 and np.linalg.norm(ex_poly[(t + 1) % len(ex_poly)] - ex_poly[t]) > 10 * eps:
                exit_lat.append(l)
        cell = dict(planes=planes, faces=faces, span=(i0, i1), exit_lat=exit_lat,
                    n_lat=len(planes) - 2)
        # validity of this cell
        if len(faces[0][0]) < 3 or len(faces[1][0]) < 3 or len(exit_lat) < 3:
            problems.append(("degenerate", k, k))
        if any(l < 0 for f in faces for l in f[1]):
            problems.append(("unbounded", k, k))
        if prev_exit is not None:
            a_prev = _poly_area(prev_exit["faces"][1][0], prev_exit["planes"][1][0])
            a_here = _poly_area(faces[0][0], planes[0][0])
            if abs(a_prev - a_here) > 1e-6 * max(a_prev, 1e-30) + 1e-12:
                if 1 in faces[0][1]:          # exit plane cuts the entry polygon
                    problems.append(("crossing", k, k))
                else:
                    problems.append(("joint", k - 1, k))
        cells.append(cell)
        prev_exit = cell
    return cells, problems


def _cell_vertices(cell):
    return np.vstack([f[0] for f in cell["faces"] if len(f[0])])


def _cell_edges(cell):
    dirs = []
    for P, _ in cell["faces"]:
        if len(P):
            E = np.roll(P, -1, 0) - P
            L = np.linalg.norm(E, axis=1)
            dirs.append(E[L > 0] / L[L > 0, None])
    return np.vstack(dirs)


def _cells_overlap(c1, c2, tol):
    """Separating-axis test for two convex cells (True = interiors overlap)."""
    V1, V2 = _cell_vertices(c1), _cell_vertices(c2)
    E1, E2 = _cell_edges(c1), _cell_edges(c2)
    X = np.cross(E1[:, None, :], E2[None, :, :]).reshape(-1, 3)
    L = np.linalg.norm(X, axis=1)
    axes = np.vstack([[n for n, _ in c1["planes"]], [n for n, _ in c2["planes"]],
                      X[L > 1e-9] / L[L > 1e-9, None]])
    p1, p2 = V1 @ axes.T, V2 @ axes.T
    sep = (p1.max(0) <= p2.min(0) + tol) | (p2.max(0) <= p1.min(0) + tol)
    return not np.any(sep)


def _point_triangle_distance(P, a, b, c):
    """Distance from points P (n,3) to triangle (a,b,c) (Ericson, vectorised)."""
    ab, ac, ap = b - a, c - a, P - a
    d1, d2 = ap @ ab, ap @ ac
    bp = P - b
    d3, d4 = bp @ ab, bp @ ac
    cp = P - c
    d5, d6 = cp @ ab, cp @ ac
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    res = np.empty((len(P), 3))
    den = va + vb + vc
    den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    v, w = vb / den, vc / den
    res[:] = a + v[:, None] * ab + w[:, None] * ac           # interior
    m = (d1 <= 0) & (d2 <= 0); res[m] = a
    m = (d3 >= 0) & (d4 <= d3); res[m] = b
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    t = d1 / np.where(d1 - d3 == 0, 1, d1 - d3); res[m] = (a + t[:, None] * ab)[m]
    m = (d6 >= 0) & (d5 <= d6); res[m] = c
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    t = d2 / np.where(d2 - d6 == 0, 1, d2 - d6); res[m] = (a + t[:, None] * ac)[m]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6))
    res[m] = (b + t[:, None] * (c - b))[m]
    return np.linalg.norm(P - res, axis=1)


def envelope_boundary(cells):
    """Boundary faces of the chain: lateral faces of every cell, entry face of
    the first and exit face of the last. Returns list of (vertices, kind,
    cell index) with kind in {'lateral', 'first', 'last'}."""
    out = []
    for k, c in enumerate(cells):
        for i, (P, _) in enumerate(c["faces"]):
            if len(P) < 3:
                continue
            if i == 0 and k > 0 or i == 1 and k < len(cells) - 1:
                continue
            kind = "first" if i == 0 else "last" if i == 1 else "lateral"
            out.append((P, kind, k))
    return out


def _check_points(cells, samples, ring_of, joints, d, check_stride):
    """Return (worst clearance ratio, list of (point ring, cell)) failures."""
    flat = samples[::check_stride].reshape(-1, 3)
    rings = ring_of.reshape(samples.shape[:2])[::check_stride].reshape(-1)
    inside = np.zeros(len(flat), bool)
    owner = np.full(len(flat), -1)
    for k, c in enumerate(cells):
        ok = np.ones(len(flat), bool)
        for n, h in c["planes"]:
            ok &= flat @ n <= h + 1e-9 * (1 + abs(h))
        owner[ok & ~inside] = k
        inside |= ok
    fails = []
    bad = np.nonzero(~inside)[0]
    for p in bad[:1]:
        fails.append((int(rings[p]), int(np.searchsorted(joints, rings[p], side="right") - 1)))
    dmin = np.full(len(flat), np.inf)
    dcell = np.full(len(flat), -1)
    for P, kind, k in envelope_boundary(cells):
        if kind != "lateral":
            continue
        for t in range(1, len(P) - 1):
            dist = _point_triangle_distance(flat, P[0], P[t], P[t + 1])
            m = dist < dmin
            dmin[m], dcell[m] = dist[m], k
    ratio = dmin / d
    worst = int(np.argmin(ratio))
    if ratio[worst] < 1 - 1e-6:
        fails.append((int(rings[worst]), int(dcell[worst])))
    return float(ratio.min()), fails, int(np.sum(~inside))


def build_envelope(C, A, B, separation, n_sides=8, ray_tol=None, max_len=None,
                   n_theta=48, check_stride=1, max_iter=200, verbose=True):
    """Closed polyhedron (planar faces) around the beam with clearance
    ``separation`` from the beam's 1/e surface, in the units of C, A, B.

    Returns (faces, cells, info): faces = list of (vertices (k,3), kind)."""
    d = float(separation)
    if d <= 0:
        raise EnvelopeError("envelope separation must be > 0")
    if n_sides < 3:
        raise EnvelopeError("envelope_sides must be >= 3")
    n = len(C)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))])
    wmin, wmax = section_widths(A, B)
    nrm = _unit(np.cross(A, B))                    # section planes (forward)
    kfac = 1.0 / np.cos(np.pi / n_theta)           # circumscribed polygons
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    samples = ellipse_rings(C, A * kfac, B * kfac, theta)
    ring_of = np.repeat(np.arange(n), n_theta)
    reach = 2.0 * (wmax * kfac + d)                # along-ray window per point
    span = np.ptp(samples.reshape(-1, 3), axis=0).max() + 4 * d
    size = 20.0 * span
    eps = 1e-9 * span
    ray_tol = 0.5 * d if ray_tol is None else ray_tol
    max_len = (np.median(wmax) + d) if max_len is None else max_len
    min_len = 0.25 * (np.median(wmax) + d)
    args = (C, A, B, nrm, s, samples, ring_of, d, n_sides, reach, size, eps)

    def validate(joints):
        """(cells, problem or None, worst clearance ratio)."""
        cells, problems = _build_cells(joints, *args)
        if problems:
            return cells, problems[0], None
        for i in range(len(cells)):
            for j in range(i + 2, len(cells)):
                if _cells_overlap(cells[i], cells[j], eps * 10):
                    return cells, ("overlap", i, j), None
        worst, fails, _ = _check_points(cells, samples, ring_of, np.array(joints), d,
                                        check_stride)
        if fails:
            r, c = fails[0]
            a = int(np.clip(np.searchsorted(joints, r, side="right") - 1, 0, len(cells) - 1))
            return cells, ("clearance", min(a, c), max(a, c)), worst
        return cells, None, worst

    # Phase 1: joints only at the bends of the ray (Douglas-Peucker points);
    # merge cells until the chain is valid.
    joints = envelope_candidate_joints(C, s, ray_tol, 0, min_len)
    if verbose:
        print(f"Envelope: separation {d:g}, {n_sides} lateral faces, "
              f"{len(joints) - 1} cell(s) between the bends of the ray")
    for it in range(max_iter):
        cells, prob, worst = validate(joints)
        if prob is None:
            break
        kind, a, b = prob
        if len(joints) <= 2:
            raise EnvelopeError(f"could not build a valid envelope ({kind}); try a larger "
                                "separation or more lateral faces")
        # merge: remove the joints strictly inside cells a..b. For a problem
        # in a single cell, a hinge that cannot be placed means the cell
        # starts too close to a sharp bend -> merge with the previous cell;
        # otherwise merge with the next one.
        if b > a:
            rm = joints[a + 1:b + 1]
        elif kind == "hinge" and a > 0:
            rm = [joints[a]]
        elif kind == "crossing" and a + 1 < len(joints) - 1:
            rm = [joints[a + 1]]
        else:
            # remove the end of cell a whose plane is most tilted with
            # respect to the cell's chord (never the first/last joint)
            ends = [x for x in (a, a + 1) if 0 < x < len(joints) - 1]
            if not ends:
                raise EnvelopeError(f"could not build a valid envelope ({kind}); try a "
                                    "larger separation or more lateral faces")
            chord = _unit(C[joints[a + 1]] - C[joints[a]])
            rm = [joints[min(ends, key=lambda x: np.dot(nrm[joints[x]], chord))]]
        if verbose > 1:
            print(f"  {kind} problem in cell(s) {a}-{b}: merging (removing joint(s) {rm})")
        joints = [j for j in joints if j not in rm]
    else:
        raise EnvelopeError("envelope: too many merge iterations")
    base, base_cells, base_worst = list(joints), cells, worst

    # Phase 2: split long cells (<= max_len) so the envelope follows the beam
    # width more closely; drop the added joints wherever they spoil the chain.
    if max_len > 0:
        added = set(envelope_candidate_joints(C, s, ray_tol, max_len, min_len)) - set(base)
        added = {j for j in added if not any(abs(s[j] - s[b]) < min_len for b in base)}
        joints = sorted(set(base) | added)
        for it in range(max_iter):
            cells, prob, worst = validate(joints)
            if prob is None:
                break
            kind, a, b = prob
            ends = {joints[x] for x in (a, a + 1, b, b + 1) if x < len(joints)}
            rm = sorted(j for j in ends if j in added)
            if not rm:
                lo, hi = joints[max(a - 1, 0)], joints[min(b + 2, len(joints) - 1)]
                rm = [j for j in joints if lo <= j <= hi and j in added]
            if verbose > 1:
                print(f"  refine: {kind} problem in cell(s) {a}-{b}: dropping joint(s) {rm}")
            if not rm:
                joints, cells, worst = base, base_cells, base_worst
                break
            added -= set(rm)
            joints = [j for j in joints if j not in rm]
        else:
            joints, cells, worst = base, base_cells, base_worst

    faces = [(P, kind) for P, kind, _ in envelope_boundary(cells)]
    info = dict(cells=len(cells), faces=len(faces), clearance_ratio=worst, joints=joints)
    if verbose:
        vol = envelope_volume(faces)
        edge = min(np.linalg.norm(np.roll(P, -1, 0) - P, axis=1).min() for P, _ in faces)
        print(f"Envelope: {len(cells)} cell(s), {len(faces)} planar faces, volume {vol:.4g}, "
              f"min clearance {worst * d:.4g} (requested {d:g}), shortest edge {edge:.3g}")
    return faces, cells, info


def envelope_volume(faces):
    v = 0.0
    for P, _ in faces:
        for t in range(1, len(P) - 1):
            v += np.dot(P[0], np.cross(P[t], P[t + 1])) / 6.0
    return v


def envelope_triangles(faces):
    """Fan-triangulate the (convex, outward-ordered) faces -> verts, tris."""
    verts, tris = [], []
    for P, _ in faces:
        b = len(verts)
        verts.extend(P)
        tris.extend([b, b + t, b + t + 1] for t in range(1, len(P) - 1))
    return np.array(verts), np.array(tris, dtype=np.int64)


def write_envelope_iges(path, faces, unit):
    """One trimmed planar surface per face. Labels: ENVFIRST / ENVLAST for
    the end faces (planes of the first / last beam cross-section, i.e. the
    beam's end caps) and ENVSIDE for the others."""
    iges = IgesWriter(unit=unit)
    names = {"first": "ENVFIRST", "last": "ENVLAST", "lateral": "ENVSIDE"}
    for P, kind in faces:
        iges.add_planar_face(P, label=names[kind])
    iges.write(path, description="Polyhedral FEM envelope around the Scotty beam")


# ##########################################################################
# 11. ERMES <-> SCOTTY, O-X CONVERSION DIAGNOSTICS, ERMES FEM COST
# ##########################################################################

def _unit3(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def scotty_launch_from_ermes(point_xyz, k_xyz, w0, waist_distance, freq):
    """ERMES launch (Cartesian point, direction k0, waist w0 and its distance
    ahead of the launch point) -> Scotty launch:
        (R, zeta, Z), poloidal and toroidal TORBEAM angles [deg],
        launch_beam_width [m], launch_beam_curvature [1/m].

    Scotty builds the launch wavevector as
        K_R = -K cos(tor) cos(pol), K_zeta/R = -K sin(tor) cos(pol),
        K_Z = -K sin(pol),
    with (R, zeta, Z) the LOCAL cylindrical frame at the launch point, so k0
    is first rotated into that frame. Gaussian beam (z = distance past the
    waist, zR = pi w0^2 / lambda):  Psi = k0 / (z - i zR)  ->
        W = w0 sqrt(1 + (z/zR)^2),  curvature = Re(Psi)/k0 = z / (z^2 + zR^2).
    """
    x, y, z = map(float, point_xyz)
    k = _unit3(k_xyz)
    R, zeta = np.hypot(x, y), np.arctan2(y, x)
    c, s = np.cos(zeta), np.sin(zeta)
    k_R, k_phi, k_Z = k[0] * c + k[1] * s, -k[0] * s + k[1] * c, k[2]
    pol = -np.degrees(np.arcsin(np.clip(k_Z, -1.0, 1.0)))
    tor = np.degrees(np.arctan2(-k_phi, -k_R))
    lam = SPEED_OF_LIGHT / freq
    zR = np.pi * w0 ** 2 / lam
    zz = -float(waist_distance)
    W = w0 * np.sqrt(1.0 + (zz / zR) ** 2)
    curv = zz / (zz ** 2 + zR ** 2) + 0.0
    return (R, zeta, z), pol, tor, W, curv, zR


def torbeam_angles_to_k(launch_position, pol_deg, tor_deg):
    """Inverse check: Scotty launch -> Cartesian unit k."""
    R, zeta, _ = launch_position
    p, t = np.radians(pol_deg), np.radians(tor_deg)
    v = cyl_to_cart_vector(np.atleast_1d(-np.cos(t) * np.cos(p)),
                           np.atleast_1d(-np.sin(t) * np.cos(p)),
                           np.atleast_1d(-np.sin(p)), np.atleast_1d(zeta))
    return v[0]


def apply_ermes_input(cfg):
    """Replace the Scotty launch parameters in cfg by the ERMES-equivalent ones."""
    pos, pol, tor, W, curv, zR = scotty_launch_from_ermes(
        cfg.ermes_launch_point, cfg.ermes_k0, cfg.ermes_beam_width,
        cfg.ermes_waist_distance, cfg.frequency)
    cfg.launch_position = tuple(float(v) for v in pos)
    cfg.poloidal_angle = round(float(pol), 6)
    cfg.toroidal_angle = round(float(tor), 6)
    cfg.launch_beam_width = float(W)
    cfg.launch_beam_curvature = float(curv)
    k_back = torbeam_angles_to_k(cfg.launch_position, cfg.poloidal_angle, cfg.toroidal_angle)
    err = np.linalg.norm(k_back - _unit3(cfg.ermes_k0))
    print("ERMES -> Scotty launch")
    print(f"  ERMES point X,Y,Z [m]    : {', '.join(f'{v:.6g}' for v in cfg.ermes_launch_point)}")
    print(f"  ERMES k0 (unit)          : {', '.join(f'{v:.6f}' for v in _unit3(cfg.ermes_k0))}")
    print(f"  ERMES w0 [m], waist ahead: {cfg.ermes_beam_width:.6g}, {cfg.ermes_waist_distance:.6g}"
          f"   (Rayleigh length {zR:.4g} m)")
    print(f"  Scotty R, zeta, Z        : {pos[0]:.7f} m, {pos[1]:.7f} rad "
          f"({np.degrees(pos[1]):.4f} deg), {pos[2]:.7f} m")
    print(f"  Scotty pol, tor (TORBEAM): {cfg.poloidal_angle:.4f} deg, {cfg.toroidal_angle:.4f} deg")
    print(f"  Scotty width, curvature  : {W:.6g} m, {curv:.6g} 1/m")
    print(f"  check k from angles      : {', '.join(f'{v:.6f}' for v in k_back)}  (|error| {err:.1e})")


def psi_waist_parameters(psi_xx, psi_xy, psi_yy, k0):
    """For the two principal directions of the complex 2x2 Psi at one point:
    local width W, waist w0 and distance from the point to the waist along the
    beam (> 0: waist ahead). From 1/Psi = (z - i zR)/k0."""
    M = np.array([[psi_xx, psi_xy], [psi_xy, psi_yy]], dtype=complex)
    out = []
    for lam in np.linalg.eigvals(M):
        inv = 1.0 / lam
        z, zR = k0 * inv.real, -k0 * inv.imag
        W = np.sqrt(2.0 / lam.imag) if lam.imag > 0 else np.inf
        w0 = np.sqrt(2.0 * zR / k0) if zR > 0 else np.nan
        out.append((W, w0, -z))
    return out


def ermes_polarisation(E, k, epol_ref=None):
    """Complex E (Cartesian, Scotty convention) -> ERMES form:
    (Epol, |Epar|, phase_par, |Eper|, phase_per), with Epar along Epol and
    Eper along k x Epol, normalised and with phase_par = 0. If epol_ref is
    None, Epol is the major axis of the polarisation ellipse."""
    k = _unit3(np.real(k))
    E = np.asarray(E, dtype=complex)
    E = E - (E @ k) * k
    if epol_ref is not None:
        u = np.asarray(epol_ref, dtype=float)
        u = _unit3(u - (u @ k) * k)
    else:
        u = _unit3(np.real(E * np.exp(-0.5j * np.angle(np.sum(E * E)))))
    v = np.cross(k, u)
    a, b = E @ u, E @ v
    nrm = np.sqrt(abs(a) ** 2 + abs(b) ** 2)
    g = np.exp(-1j * np.angle(a)) / nrm
    a, b = a * g, b * g
    return u, abs(a), 0.0, abs(b), float(np.angle(b))


def ermes_field_vector(k, epol, epar, eper):
    """ERMES (k0, Epol, (|Epar|, phase), (|Eper|, phase)) -> complex E."""
    k, u = _unit3(k), _unit3(epol)
    return (epar[0] * np.exp(1j * epar[1]) * u
            + eper[0] * np.exp(1j * eper[1]) * np.cross(k, u))


def polarisation_overlaps(E_scotty, E_ermes):
    """|<E_ermes, E_scotty>| (same time convention) and the same with E_ermes
    conjugated (opposite time convention: e^{+jwt} vs e^{-iwt})."""
    a = np.asarray(E_scotty, complex) / np.linalg.norm(E_scotty)
    b = np.asarray(E_ermes, complex) / np.linalg.norm(E_ermes)
    return abs(np.vdot(b, a)), abs(np.sum(b * a))


def write_ermes_launch(path, cfg, beam, k0, omega):
    """Scotty beam at its first point (after any vacuum extension) in ERMES
    form, plus comparison with the ERMES polarisation in cfg."""
    E, K = beam.E[0], np.real(beam.K[0])
    k = _unit3(K)
    if hasattr(cfg, "propagation_length"):          # beam_to_gid.py (vacuum beam)
        head = ["# ERMES launch data of the vacuum Gaussian beam (beam_to_gid.py),",
                "# at the beam start written to GiD."]
    else:
        head = ["# ERMES launch data from Scotty. E and k are Scotty's at the plasma entry;",
                "# in vacuum they are the same at the beam start written to GiD."]
    lines = head + [
             f"Freq    = {omega / (2 * np.pi):.6g}",
             f"Point   = {','.join(f'{v:.9g}' for v in beam.points[0])}   "
             "(first point of the beam written to GiD)",
             f"k0      = {','.join(f'{v:.6f}' for v in k)}"]
    use_ermes = bool(getattr(cfg, "ermes_input", False))
    ref = cfg.ermes_epol if (use_ermes and cfg.ermes_epol is not None) else None
    for label, r in (("major axis of the ellipse", None),
                     ("your ERMES Epol", ref)):
        if label != "major axis of the ellipse" and r is None:
            continue
        u, ap, pp, ar, pr = ermes_polarisation(E, k, r)
        lines += [f"# Epol = {label}",
                  f"Epol    = {','.join(f'{v:.6f}' for v in u)}",
                  f"|Epar|  = {ap:.6f},  Phase = {pp:.6g}",
                  f"|Eper|  = {ar:.6f},  Phase = {pr:.6g}"]
    if beam.psi_xx is not None:
        for j, (W, w0, ahead) in enumerate(psi_waist_parameters(
                beam.psi_xx[0], beam.psi_xy[0], beam.psi_yy[0], k0)):
            lines.append(f"# principal direction {j + 1}: width here {W:.6g} m, "
                         f"waist w0 {w0:.6g} m at {ahead + 0.0:.6g} m ahead")
    if (use_ermes and cfg.ermes_epol is not None and cfg.ermes_epar is not None
            and cfg.ermes_eper is not None):
        k_e = cfg.ermes_k0 if cfg.ermes_k0 is not None else k
        Ee = ermes_field_vector(k_e, cfg.ermes_epol, cfg.ermes_epar, cfg.ermes_eper)
        same, conj = polarisation_overlaps(E, Ee)
        lines += ["# Overlap |<E_ERMES, E_Scotty>| (1 = same polarisation; power in the",
                  "# Scotty mode = overlap^2):",
                  f"#   same time convention     : {same:.4f}  (power {same ** 2:.4f})",
                  f"#   opposite time convention : {conj:.4f}  (power {conj ** 2:.4f})",
                  "# If the second is the larger one, ERMES and Scotty use opposite",
                  "# time conventions: use the complex conjugate of the phases above."]
        ang = np.degrees(np.arccos(np.clip(abs(k @ _unit3(k_e)), -1, 1)))
        if ang > 0.05:
            lines.append(f"# NOTE: Scotty k and ERMES k0 differ by {ang:.3f} deg.")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(l for l in lines if not l.startswith("# ERMES launch") and
                    not l.startswith("# in vacuum")))
    print(f"Wrote ERMES launch data: {path}")


# ---------------------------- plasma input readers ---------------------------

def _read_torbeam_ascii_topfile(path):
    """TORBEAM ASCII topfile -> dict(R, Z, Br, Bt, Bz, pol_flux), grids [Z][R].
    Sections are found by their label lines (X/R-coordinates, Z-coordinates,
    B_R, B_t, B_Z, psi), as Scotty reads them; the numbers of a section run
    until the next text line."""
    sections, label, nums = [], None, []
    with open(path) as f:
        for line in f:
            parts = line.replace(",", " ").split()
            if not parts:
                continue
            try:
                vals = [float(t) for t in parts]
            except ValueError:
                if label is not None:
                    sections.append((label, nums))
                label, nums = line.strip().lower(), []
                continue
            nums.extend(vals)
    if label is not None:
        sections.append((label, nums))

    def find(*keys, exclude=("number of", "grid points", "psi_sep", "outside radius")):
        tests = (lambda lab, k: lab == k, lambda lab, k: lab.startswith(k),
                 lambda lab, k: k in lab)
        for test in tests:                          # exact > starts with > contains
            for k in keys:                          # keys in order of priority
                for lab, v in sections:
                    if test(lab, k) and not any(e in lab for e in exclude):
                        return np.asarray(v, float)
        raise KeyError(f"{path}: no section labelled {keys} "
                       f"(labels: {[lab for lab, _ in sections]})")
    R = find("x-coord", "r-coord", "radial grid", "radial coord")
    Z = find("z-coord", "vertical grid", "vertical coord")
    out = dict(R=R, Z=Z)
    for key, names in (("Br", ("b_r",)), ("Bt", ("b_t", "b_phi")), ("Bz", ("b_z",)),
                       ("pol_flux", ("psi", "pol_flux", "poloidal flux"))):
        a = find(*names)
        if a.size != R.size * Z.size:
            raise ValueError(f"{path}: {key} has {a.size} values, expected "
                             f"{R.size} x {Z.size}")
        out[key] = a.reshape(Z.size, R.size)        # R fastest, as Scotty reads it
    return out


def find_topfile(folder, suffix):
    """topfile<suffix>.json or the TORBEAM ASCII topfile<suffix>, or None."""
    for name in (f"topfile{suffix}.json", f"topfile{suffix}", f"topfile{suffix}.txt",
                 f"topfile{suffix}.dat"):
        fp = os.path.join(folder, name)
        if os.path.isfile(fp):
            return fp
    return None


def load_topfile(path):
    """topfile<suffix>.json or TORBEAM ASCII topfile: returns f(R, Z) ->
    (B_R, B_t, B_Z, psi_n, dpsi/dR, dpsi/dZ). Grids stored [Z][R] as Scotty
    expects (also accepts [R][Z])."""
    import json
    from scipy.interpolate import RectBivariateSpline
    if path.lower().endswith(".json"):
        with open(path) as f:
            data = json.load(f)
    else:
        data = _read_torbeam_ascii_topfile(path)
    low = {k.lower(): k for k in data}

    def get(*names):
        for n in names:
            if n.lower() in low:
                return np.asarray(data[low[n.lower()]], dtype=float)
        raise KeyError(f"{path}: none of {names} found (keys: {list(data)})")

    R, Z = get("R").ravel(), get("Z").ravel()

    def grid(*names):
        a = get(*names)
        if a.ndim == 1 and a.size == R.size * Z.size:
            a = a.reshape(Z.size, R.size)          # flattened, R fastest (eqdsk_to_scotty, OMFIT)
        if a.shape == (Z.size, R.size):
            a = a.T
        elif a.shape != (R.size, Z.size):
            raise ValueError(f"{path}: {names[0]} has shape {a.shape}, "
                             f"expected ({Z.size}, {R.size})")
        if R[0] > R[-1]:
            a = a[::-1]
        if Z[0] > Z[-1]:
            a = a[:, ::-1]
        return a

    Rs, Zs = np.sort(R), np.sort(Z)
    sp = {key: RectBivariateSpline(Rs, Zs, grid(*names))
          for key, names in (("BR", ("Br", "B_R")), ("Bt", ("Bt", "B_t", "Bphi")),
                             ("BZ", ("Bz", "B_Z")),
                             ("psi", ("pol_flux", "poloidal_flux", "psi")))}

    def field(Rp, Zp):
        return (sp["BR"].ev(Rp, Zp), sp["Bt"].ev(Rp, Zp), sp["BZ"].ev(Rp, Zp),
                sp["psi"].ev(Rp, Zp), sp["psi"].ev(Rp, Zp, dx=1),
                sp["psi"].ev(Rp, Zp, dy=1))
    return field


def load_ne_profile(path, coordinate="rho", units=1e19, psi_zero=None):
    """ne<suffix>.dat -> f(psi_n) -> (ne [m^-3], dne/dpsi_n). Lines that do not
    hold two numbers (e.g. the point count) are skipped."""
    rows = []
    with open(path) as f:
        for line in f:
            try:
                vals = [float(t) for t in line.replace(",", " ").split()]
            except ValueError:
                continue
            if len(vals) >= 2:
                rows.append(vals[:2])
    x, ne = np.array(rows).T
    psi = x ** 2 if coordinate == "rho" else x
    o = np.argsort(psi)
    psi, ne = psi[o], ne[o] * units
    dne = np.gradient(ne, psi)

    def prof(p):
        p = np.asarray(p, dtype=float)
        n = np.interp(p, psi, ne, right=0.0)
        dn = np.interp(p, psi, dne, right=0.0)
        if psi_zero is not None:
            n, dn = np.where(p >= psi_zero, 0.0, n), np.where(p >= psi_zero, 0.0, dn)
        return np.clip(n, 0.0, None), dn
    return prof


# ---------------------------- O-X conversion ---------------------------------

def _gauss_avg(a, mu, sigma):
    """< exp(-a x^2) > for x ~ Normal(mu, sigma^2)."""
    if not np.isfinite(a):          # L_n = inf (no density gradient): window closed
        return 0.0
    q = 1.0 + 2.0 * a * sigma ** 2
    return np.exp(-a * mu ** 2 / q) / np.sqrt(q)


def ox_conversion_diagnostics(d, K, cfg, w0_eff, path, note=None):
    """N_par, N_y, L_n, Mjolhus T along the Scotty ray (before extension);
    write a table to ``path`` and print the values at the cut-off (or at the
    point of highest density if the ray never reaches X = 1)."""
    suffix = cfg.input_filename_suffix
    top = find_topfile(cfg.magnetic_data_path, suffix)
    nef = os.path.join(cfg.ne_data_path, f"ne{suffix}.dat")
    if top is None:
        print(f"WARNING: O-X diagnostics skipped: no topfile{suffix}(.json) in "
              f"{os.path.abspath(cfg.magnetic_data_path)}.")
        return None
    if not os.path.isfile(nef):
        print(f"WARNING: O-X diagnostics skipped: {os.path.abspath(nef)} not found.")
        return None
    print(f"O-X diagnostics: B and psi from {top}, ne from {nef}")
    if d["omega"] is None or K is None:
        print("WARNING: O-X diagnostics skipped: no frequency or wavevector in the .h5.")
        return None
    field = load_topfile(top)
    prof = load_ne_profile(nef, cfg.ox_ne_radial_coordinate, cfg.ox_ne_units,
                           cfg.poloidal_flux_zero_density)

    omega = d["omega"]
    k0 = omega / SPEED_OF_LIGHT
    P = d["beam"]
    R, zeta, Z = np.hypot(P[:, 0], P[:, 1]), np.arctan2(P[:, 1], P[:, 0]), P[:, 2]
    BR, Bt, BZ, psi, dpR, dpZ = field(R, Z)
    B = cyl_to_cart_vector(BR, Bt, BZ, zeta)
    G = cyl_to_cart_vector(dpR, np.zeros_like(dpR), dpZ, zeta)
    Bm, Gm = np.linalg.norm(B, axis=1), np.linalg.norm(G, axis=1)
    bh, nh = B / Bm[:, None], G / Gm[:, None]
    yh = np.cross(bh, nh)
    yh /= np.linalg.norm(yh, axis=1)[:, None]
    N = np.real(K) / k0
    Npar, Ny, Nn = (np.sum(N * v, axis=1) for v in (bh, yh, nh))
    ne, dne = prof(psi)
    nc = EPS0 * E_MASS * omega ** 2 / E_CHARGE ** 2
    X = ne / nc
    Y = E_CHARGE * Bm / (E_MASS * omega)
    Nopt = np.sqrt(Y / (1.0 + Y))
    with np.errstate(divide="ignore", invalid="ignore"):
        Ln = np.where(np.abs(dne) > 0, ne / (np.abs(dne) * Gm), np.inf)
    a = np.pi * k0 * Ln * np.sqrt(Y / 2.0) * 2.0 * (1.0 + Y)
    b = np.pi * k0 * Ln * np.sqrt(Y / 2.0)
    with np.errstate(invalid="ignore", over="ignore"):
        T = np.exp(-a * (np.abs(Npar) - Nopt) ** 2 - b * Ny ** 2)

    cols = dict(distance=d["distance"], R=R, Z=Z, psi_n=psi, ne=ne, X=X, Y=Y, B=Bm,
                N_par=Npar, N_y=Ny, N_normal=Nn, N_par_opt=Nopt, L_n=Ln, T_mjolhus=T)
    with open(path, "w") as f:
        f.write("# " + " ".join(f"{c:>14s}" for c in cols) + "\n")
        for row in zip(*cols.values()):
            f.write("  " + " ".join(f"{v:14.7g}" for v in row) + "\n")

    # evaluation point: first X = 1 crossing, else maximum X
    i = int(np.argmax(X))
    up = np.nonzero((X[:-1] < 1.0) & (X[1:] >= 1.0))[0]
    if up.size:
        j = up[0]
        t = (1.0 - X[j]) / (X[j + 1] - X[j])
        at = lambda v: (1 - t) * v[j] + t * v[j + 1]
        where = f"X = 1 crossing, distance {at(d['distance']):.4g} m"
    else:
        at = lambda v: v[i]
        why = "ray turns before the cut-off" if X[i] < 1.0 else "no X = 1 crossing found along the ray"
        where = f"maximum X = {X[i]:.3f} ({why}), distance {d['distance'][i]:.4g} m"
    y_, np_, ny_, no_, ln_ = at(Y), at(Npar), at(Ny), at(Nopt), at(Ln)
    a_, b_ = at(a), at(b)
    dN = abs(np_) - no_
    T_pw = float(np.exp(-a_ * dN ** 2 - b_ * ny_ ** 2))
    sigma = 1.0 / (k0 * w0_eff) if w0_eff and np.isfinite(w0_eff) else 0.0
    T_beam = _gauss_avg(a_, dN, sigma) * _gauss_avg(b_, ny_, sigma)
    T_best = _gauss_avg(a_, 0.0, sigma) * _gauss_avg(b_, 0.0, sigma)

    print(f"O-X conversion (Mjolhus), evaluated at {where}:")
    print(f"  Y = {y_:.4f}   N_par,opt = {no_:.4f}   L_n = {ln_:.4g} m")
    print(f"  N_par = {np_:+.4f}  (|N_par| - N_par,opt = {dN:+.4f})   N_y = {ny_:+.4f}")
    print(f"  window 1/e half-widths: dN_par = {1 / np.sqrt(a_):.4f}, dN_y = {1 / np.sqrt(b_):.4f}")
    print(f"  plane-wave T                          = {T_pw:.4f}")
    if sigma:
        print(f"  beam-averaged T (w0 = {w0_eff * 1e3:.4g} mm, sigma_N = {sigma:.4f}) = {T_beam:.4f}")
        print(f"  best beam-averaged T with this w0 (aim exactly at the window) = {T_best:.4f}")
    if X.max() < 0.5:
        print("  WARNING: the ray never gets near the O cut-off (X_max < 0.5): check "
              "mode_flag (+1 = O on MAST-U) and the launch angles.")
    print("  NOTE: " + (note or "beam tracing is not valid at the cut-off; these are "
                         "estimates. N_par/N_y are Scotty's local values there."))
    print(f"Wrote O-X diagnostics along the ray: {path}")
    return dict(T_plane=T_pw, T_beam=T_beam, T_best=T_best, N_par=np_, N_y=ny_,
                N_par_opt=no_, Y=y_, L_n=ln_)


# ---------------------------- ERMES FEM cost ---------------------------------

def _human_bytes(b):
    for u in ("B", "kB", "MB", "GB", "TB", "PB"):
        if b < 1000 or u == "PB":
            return f"{b:.3g} {u}"
        b /= 1000.0


def _human_time(s):
    for u, f in (("s", 60), ("min", 60), ("h", 24), ("days", 365)):
        if s < f:
            return f"{s:.3g} {u}"
        s /= f
    return f"{s:.3g} years"


# Tetrahedra per h^3 for a good-quality mesh of mean edge h (Delaunay of a
# jittered close-packed point set: 8.6-10.6; regular tetrahedra: 8.5)
TETS_PER_H3 = (7.0, 9.0, 11.0)


def fem_cost_parameters(cfg, prefix="fem_", **override):
    """Collect the FEM-cost settings from a cfg namespace (names prefix+...)."""
    g = lambda n, dflt: getattr(cfg, prefix + n, dflt)
    p = dict(elements_per_wavelength=g("elements_per_wavelength", 10.0),
             element_order=g("element_order", 1), dofs_per_node=g("dofs_per_node", 3),
             vertices_per_tet=g("vertices_per_tet", 0.164), solver=g("solver", "direct"),
             fill_constant=g("fill_constant", 15.0), flop_constant=g("flop_constant", 60.0),
             cores=g("cores", 128), gflops_per_core=g("gflops_per_core", 10.0),
             iterations=g("iterations", 1000), refractive_index=g("refractive_index", 1.0))
    p.update(override)
    return p


def fem_cost_lines(V, area, k0, p, n_ray_max=1.0, title=None):
    """Mesh and solver cost of an ERMES model of volume V [m^3] and surface
    area [m^2], meshed with tetrahedra of mean edge
    h = lambda0 / (elements_per_wavelength * |N|).  Returns (lines, info).
    Counts: T ~ 9 V/h^3 (7..11), Euler relations F ~ 2T, E ~ V + T.
    Direct solver: nested-dissection scalings fill ~ Ndof^(4/3), flops ~ Ndof^2."""
    lam0 = 2 * np.pi / k0
    Nidx = p["refractive_index"] if p["refractive_index"] > 0 else max(1.0, n_ray_max)
    epw = float(p["elements_per_wavelength"])
    h = lam0 / (epw * Nidx)
    lo, T, hi = (c * V / h ** 3 for c in TETS_PER_H3)
    Vv = p["vertices_per_tet"] * T
    E = Vv + T
    F = 2 * T
    if p["element_order"] == 1:
        nodes, pairs = Vv, E
    else:                                   # 10-node tetrahedra: vertex + edge nodes
        nodes, pairs = Vv + E, 3 * E + 6 * F + 3 * T
    dpn = p["dofs_per_node"]
    Ndof = dpn * nodes
    nnz = dpn ** 2 * (nodes + 2 * pairs)
    tri = area / (np.sqrt(3) / 4 * h ** 2) if area else 0.0
    mat_mem = nnz * 20.0 + Ndof * 16.0       # complex AIJ, 32-bit indices
    rate = p["cores"] * p["gflops_per_core"] * 1e9
    bar = "=" * 74
    L = [bar, (title or " ESTIMATE OF THE ERMES MESH FOR THE FEM ENVELOPE").center(74), bar,
         f"  frequency / wavelength    : {SPEED_OF_LIGHT / lam0 / 1e9:.4g} GHz / {lam0 * 1e3:.4g} mm",
         f"  domain volume / surface   : {V:.4g} m^3 / {area:.4g} m^2"
         f"  (= {V / (lam0 / Nidx) ** 3:.3e} local lambda^3)",
         f"  mean tetrahedron edge h   : lambda/({epw:g} x |N| = {Nidx:.3g}) = {h * 1e3:.4g} mm",
         f"  TETRAHEDRA                : ~ {T:.3e}   (range {lo:.2e} .. {hi:.2e})",
         f"  boundary triangles        : ~ {tri:.3e}",
         f"  mesh vertices / edges     : ~ {Vv:.3e} / {E:.3e}",
         f"  nodes (order {p['element_order']})           : ~ {nodes:.3e}",
         f"  UNKNOWNS (complex)        : ~ {Ndof:.3e}   ({dpn} per node)",
         f"  matrix non-zeros          : ~ {nnz:.3e}  ({nnz / Ndof:.0f} per row)",
         f"  matrix memory             : ~ {_human_bytes(mat_mem)}"]
    if p["solver"] == "direct":
        fill = p["fill_constant"] * Ndof ** (4.0 / 3.0)
        flops = p["flop_constant"] * Ndof ** 2
        L += [f"  direct: factor entries    : ~ {fill:.3e}",
              f"  direct: memory            : ~ {_human_bytes(1.2 * 16 * fill + mat_mem)}",
              f"  direct: factorisation     : ~ {flops:.2e} flop -> "
              f"{_human_time(flops / rate)} on {p['cores']} cores "
              f"x {p['gflops_per_core']:g} GFLOP/s"]
    else:
        restart = 30
        mem = mat_mem * 3 + (restart + 5) * 16.0 * Ndof
        f_it = 3 * 8 * nnz + 8 * 2 * restart * Ndof
        L += [f"  iterative: memory (matrix + preconditioner ~2x + GMRES({restart}))"
              f" : ~ {_human_bytes(mem)}",
              f"  iterative: {p['iterations']} iterations : "
              f"{_human_time(p['iterations'] * f_it / rate)} at the flop rate "
              "(memory-bandwidth bound in practice: expect slower)"]
    other = ", ".join(f"lambda/{e:g}: {T * (e / epw) ** 3:.2e}" for e in (6, 8, 12, 15)
                      if e != epw)
    L += [f"  tetrahedra, other sizes   : {other}",
          "  scaling with (elements per wavelength x |N|): unknowns ^3, "
          "direct memory ^4, direct time ^6",
          "  order of magnitude only: calibrate the fill / flop constants with one "
          "real ERMES run"]
    info = dict(h=h, tets=T, nodes=nodes, ndof=Ndof, nnz=nnz, n_index=Nidx)
    return L, info


def faces_area(faces):
    """Total area of planar polygon faces [(P, kind), ...]."""
    return sum(0.5 * np.linalg.norm(np.sum(np.cross(P[1:-1] - P[0], P[2:] - P[0]), axis=0))
               for P, _ in faces)


def fem_cost_estimate(faces, k0, n_ray_max, cfg, path, plasma=True):
    """ERMES cost if the polyhedral envelope ``faces`` (metres) is the FEM
    domain; settings from the fem_* parameters of ``cfg``. Prints and writes."""
    p = fem_cost_parameters(cfg)
    L, info = fem_cost_lines(abs(envelope_volume(faces)), faces_area(faces), k0, p, n_ray_max)
    if plasma and info["n_index"] <= 1.0 and getattr(cfg, "ox_diagnostics", False):
        L.insert(-1, "  NOTE: for O-X conversion the slow X wave near the UHR has |N| >> 1: "
                     "set fem_refractive_index for that region")
    L.append("=" * 74)
    print("\n" + "\n".join(L))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"Wrote FEM estimate: {path}")
    return info


# ##########################################################################
# 12. COMMAND LINE / MAIN
# ##########################################################################

def _parameter_comments():
    """Read the comment written next to each parameter in section 1 (used
    as its --help text) and the section titles (used as help groups)."""
    info, group, current = {}, "Parameters", None
    try:
        with open(os.path.abspath(__file__)) as f:
            lines = f.read().split("\n")
    except OSError:
        return info
    start = next((i for i, l in enumerate(lines) if l.startswith("_names_before_parameters")), 0)
    for line in lines[start + 1:]:
        if line.startswith("# (end of user parameters"):
            break
        s = line.strip()
        if s.startswith("# ----"):
            group, current = s.strip("#- ").strip(), None
            continue
        head = line.split("#", 1)[0]
        comment = line.split("#", 1)[1].strip() if "#" in line else ""
        if "=" in head and not line.startswith(" "):
            current = head.split("=", 1)[0].strip()
            info[current] = [group, comment]
        elif current and s.startswith("#") and not head.strip():
            info[current][1] += "; " + comment
        else:
            current = None
    return info


class _HelpFormatter(argparse.RawDescriptionHelpFormatter,
                     argparse.ArgumentDefaultsHelpFormatter):
    pass


def build_parser():
    """Command line built automatically from the section-1 parameters."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=_HelpFormatter)
    comments = _parameter_comments()
    groups = {}
    defaults = default_parameters()
    for name in USER_PARAMETERS:
        default = getattr(defaults, name)
        grp, text = comments.get(name, ("Parameters", ""))
        g = groups.get(grp) or groups.setdefault(grp, ap.add_argument_group(grp))
        flags = ["--" + name.replace("_", "-")] + PARAMETER_ALIASES.get(name, [])
        kw = dict(dest=name, default=default, help=text or " ")
        if isinstance(default, bool):
            kw["action"] = argparse.BooleanOptionalAction
        elif name == "launch_position":
            kw.update(nargs=3, type=float, metavar=("R", "ZETA", "Z"))
        elif name in PARAMETER_NARGS:
            kw.update(nargs=PARAMETER_NARGS[name], type=float, metavar="V")
        else:
            kw["type"] = PARAMETER_TYPES.get(name, type(default))
            if name in PARAMETER_CHOICES:
                kw["choices"] = PARAMETER_CHOICES[name]
            else:
                kw["metavar"] = "VALUE"
        g.add_argument(*flags, **kw)

    g = ap.add_argument_group("Tools")
    g.add_argument("--show-params", action="store_true",
                   help="Print every parameter value in use and stop")
    g.add_argument("--list-datasets", action="store_true",
                   help="List the datasets inside the .h5 file and stop")
    # Kept for command lines written for the old scotty2gid.py
    ap.add_argument("h5_positional", nargs="?", help=argparse.SUPPRESS)
    ap.add_argument("--no-stl", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--run-only", action="store_true", help=argparse.SUPPRESS)
    return ap


def parse_args(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    cfg = build_parser().parse_args(argv)
    if cfg.h5_positional:
        cfg.scotty_h5_file = cfg.h5_positional
    # A file given on the command line means "export this file" ...
    file_given = cfg.h5_positional or any(a.startswith(("--h5", "--scotty-h5-file")) for a in argv)
    # ... unless --run-scotty is also given explicitly
    if file_given and not any(a == "--run-scotty" for a in argv):
        cfg.run_scotty = False
    if cfg.no_stl:
        cfg.output_format = {"stl": "none", "both": "nurbs"}.get(cfg.output_format,
                                                                 cfg.output_format)
    if cfg.run_only:
        cfg.export_geometry = False
    for name in USER_PARAMETERS:
        if isinstance(getattr(cfg, name), list):
            setattr(cfg, name, tuple(getattr(cfg, name)))
    return cfg


def show_parameters(cfg):
    """Print every parameter with the value in use (marks changed ones)."""
    defaults = default_parameters()
    for name in USER_PARAMETERS:
        val = getattr(cfg, name)
        mark = "   <- changed" if val != getattr(defaults, name) else ""
        print(f"  {name:<28} = {val!r}{mark}")


def main(cfg=None):
    """Run everything. ``cfg`` defaults to the command line (which itself
    defaults to the values in section 1)."""
    cfg = cfg or parse_args()
    if cfg.ermes_input:
        apply_ermes_input(cfg)
    if cfg.show_params:
        show_parameters(cfg)
        return

    # ---- stage 1: get the Scotty output (run it, or use an existing file) ----
    if cfg.run_scotty:
        h5 = run_scotty_simulation(cfg)
    else:
        h5 = cfg.scotty_h5_file or newest_scotty_output(cfg.scotty_output_path)
        print(f"Using existing Scotty output: {h5}")
    if not cfg.export_geometry:
        return
    if not os.path.isfile(h5):
        sys.exit(f"File not found: {h5}")

    if cfg.list_datasets:
        import h5py
        with h5py.File(h5, "r") as f:
            print("\n".join(list_all_datasets(f, with_shapes=True)))
        return

    # ---- stage 2: load, extend, build cross-sections ----
    d = load_scotty_beam(h5)
    rotate_scotty_output(d, scotty_zeta_rotation(cfg))
    E, K = cartesian_fields(d)
    if E is None and (cfg.write_efield or cfg.write_field_arrays):
        print("WARNING: e_hat / K / q_R / q_zeta missing -- no E-field files.")
    beam = Beam(points=d["beam"], xhat=d["xhat"], yhat=d["yhat"], psi_xx=d["psi_xx"],
                psi_xy=d["psi_xy"], psi_yy=d["psi_yy"], distance=d["distance"], E=E, K=K)
    scotty_first, scotty_last = beam.points[0].copy(), beam.points[-1].copy()
    k0 = d["omega"] / SPEED_OF_LIGHT if d["omega"] else None

    os.makedirs(cfg.out_dir, exist_ok=True)
    stem = cfg.file_stem or os.path.splitext(os.path.basename(h5))[0].replace(".", "_")
    out = lambda name: os.path.join(cfg.out_dir, name)

    # launched waist for the beam's angular spectrum: from the LAUNCH definition
    # (Scotty's Psi at its first point is on the plasma side of the Psi_BC jump)
    w0_eff = None
    launch = launch_definition(d, cfg) if k0 else None
    if launch is not None:
        _, Wl, cl = launch
        psi_l = k0 * cl + 2j / Wl ** 2
        w0s = [w0 for _, w0, _ in psi_waist_parameters(psi_l, 0.0, psi_l, k0)]
        if all(np.isfinite(w0s)):
            w0_eff = float(np.sqrt(w0s[0] * w0s[1]))
    if w0_eff is None and k0 and beam.psi_xx is not None:
        w0s = [w0 for _, w0, _ in psi_waist_parameters(
            beam.psi_xx[0], beam.psi_xy[0], beam.psi_yy[0], k0)]
        if all(np.isfinite(w0s)):
            w0_eff = float(np.sqrt(w0s[0] * w0s[1]))
    n_ray_max = (float(np.max(np.linalg.norm(np.real(K), axis=1)) / k0)
                 if (K is not None and k0) else 1.0)
    if cfg.ox_diagnostics:
        try:
            ox_conversion_diagnostics(d, K, cfg, w0_eff, out(f"{stem}_OX_along_ray.dat"))
        except Exception as err:          # never lose the geometry output
            print(f"WARNING: O-X diagnostics failed: {err}")

    beam = extend_scotty_beam(beam, d, k0, cfg)
    report_beam_ends(beam, k0, d, cfg)

    if beam.E is not None and k0:
        write_ermes_launch(out(f"{stem}_ermes_launch.txt"), cfg, beam, k0, d["omega"])

    ray_len = np.sum(np.linalg.norm(np.diff(beam.points, axis=0), axis=1))
    A, B = beam_cross_sections(beam, fallback_radius=0.02 * ray_len or 0.01)
    wmin, wmax = section_widths(A, B)
    print(f"Beam: {beam.n} cross-sections, ray length {ray_len:.4g} m, "
          f"semi-axes {wmin.min():.3g} .. {wmax.max():.3g} m")
    if wmin.min() < 0.05 * np.median(wmin):
        i = np.argmin(wmin)
        print(f"NOTE: the beam is pinched to {wmin[i]:.3g} m at point {i} "
              f"(distance {beam.distance[i]:.4g} m) -- check Scotty there.")
    if wmax.max() > ray_len:
        print("WARNING: some beam widths exceed the ray length (Im(Psi) ~ 0, e.g. near a "
              "cut-off); the surface will look unphysical there.")

    env_faces = None
    if cfg.envelope:
        try:
            env_faces, _, _ = build_envelope(
                beam.points, A, B, cfg.envelope_separation, n_sides=cfg.envelope_sides,
                ray_tol=cfg.envelope_ray_tolerance or None,
                max_len=None if cfg.envelope_max_cell_length == 0
                else max(cfg.envelope_max_cell_length, 0.0))
        except EnvelopeError as err:
            print(f"WARNING: no FEM envelope written: {err}")
    if env_faces and cfg.fem_cost_estimate and k0:
        fem_cost_estimate(env_faces, k0, n_ray_max, cfg, out(f"{stem}_ermes_fem_estimate.txt"))

    theta = np.linspace(0, 2 * np.pi, cfg.stl_points_per_section, endpoint=False)
    if cfg.show_plot:
        beam_plot(beam, A, B, cfg.plot_points_around, env_faces)

    # ---- stage 3: write outputs ----
    sc = UNIT_SCALE[cfg.units]
    C_s, A_s, B_s = beam.points * sc, A * sc, B * sc

    if cfg.output_format in ("stl", "both"):
        idx = ring_selection(beam.n, cfg.stl_stride)
        verts, faces = build_tube_mesh(ellipse_rings(C_s[idx], A_s[idx], B_s[idx], theta),
                                       cfg.caps)
        path = out(f"{stem}_beam_surface.stl")
        write_stl(path, verts, faces, ascii=cfg.stl_ascii)
        print(f"Wrote STL ({len(idx)} sections, {len(faces)} triangles, {cfg.units}): {path}")

    if cfg.output_format in ("nurbs", "both"):
        patches = nurbs_beam_patches(C_s, A_s, B_s, tol=cfg.nurbs_tolerance * sc,
                                     n_patches=cfg.nurbs_patches, caps=cfg.caps,
                                     n_max=cfg.nurbs_max_sections, unit=cfg.units)
        iges = IgesWriter(unit=cfg.units)
        for pt in patches:
            iges.add_surface(pt["ctrl"], pt["w"], pt["ku"], pt["pu"], pt["kv"], pt["pv"],
                             closed_u=pt["closed_u"], label=pt["label"])
        if cfg.nurbs_include_ray:
            iges.add_polyline_curve(C_s)
        path = out(f"{stem}_beam_surface.igs")
        iges.write(path)
        print(f"Wrote NURBS/IGES ({len(patches)} surfaces, {cfg.units}): {path}")

    if env_faces:
        faces_s = [(P * sc, kind) for P, kind in env_faces]
        if cfg.output_format in ("stl", "both"):
            path = out(f"{stem}_envelope.stl")
            v, t = envelope_triangles(faces_s)
            write_stl(path, v, t, ascii=cfg.stl_ascii, name="fem_envelope")
            print(f"Wrote envelope STL ({len(t)} triangles, {cfg.units}): {path}")
        if cfg.output_format in ("nurbs", "both"):
            path = out(f"{stem}_envelope.igs")
            write_envelope_iges(path, faces_s, cfg.units)
            print(f"Wrote envelope IGES ({len(faces_s)} planar faces, {cfg.units}): {path}")

    if cfg.write_central_ray:
        path = out(f"{stem}_central_ray.dat")
        write_central_ray_file(path, C_s)
        print(f"Wrote central ray ({beam.n} points): {path}")

    if beam.E is not None and cfg.write_efield:
        write_efield_endpoints(out(cfg.efield_file), beam, scotty_first, scotty_last, sc)
        print(f"Wrote E-field/k end-point summary: {out(cfg.efield_file)}")

    if beam.E is not None and cfg.write_field_arrays:
        write_vector_rows(out(f"{stem}_Efield_along_ray.dat"), beam.E, True)
        write_vector_rows(out(f"{stem}_Kvector_along_ray.dat"), beam.K, False)
        print(f"Wrote per-point E-field and k ({beam.n} rows each).")


if __name__ == "__main__":
    main()
