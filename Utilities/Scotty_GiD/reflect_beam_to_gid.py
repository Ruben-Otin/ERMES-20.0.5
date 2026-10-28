#!/usr/bin/env python3
"""
reflect_beam_to_gid.py
======================

Same job as scotty_to_gid.py / beam_to_gid.py -- turn a Gaussian beam into
geometry GiD can import, plus a polyhedral FEM volume around it for ERMES --
but the beam is REFLECTED by perfectly conducting (PEC) surfaces it meets:

  1. BEAM     The incident beam is either
                * "scotty": the Scotty beam (Scotty is run, or an existing
                  scotty_output*.h5 is read) with its vacuum extensions,
                  exactly as scotty_to_gid.py builds it, or
                * "vacuum": the free-space Gaussian beam of beam_to_gid.py.
              Everything that defines the beam is taken from section 1 of
              THOSE scripts. Edit it there, or give it on this command line
              (any option this script does not know is passed on to them):
                  python reflect_beam_to_gid.py --h5 FILE.h5 --extend-last 0.6
                  python reflect_beam_to_gid.py --beam-source vacuum --length 1.5
  2. REFLECT  The central ray is followed from the launch point. Where it
              first crosses one of the surfaces of section 1, the beam is cut
              there, reflected, and continued as a Gaussian beam in vacuum for
              the REST of its length: the total path length is the same as
              without surfaces (Scotty beam + extensions, or
              propagation_length). Repeated up to max_reflections times.
  3. EXPORT   Per leg (piece of beam between reflections): the beam surface as
              STL and/or NURBS (IGES). One closed polyhedron with planar faces
              around the whole folded beam (FEM volume for ERMES), whose faces
              lying on the mirrors are labelled ENVPEC<k>. Central ray, E and k
              files, a reflection report and the reflecting surfaces as STL.

This script imports scotty_to_gid.py (always) and beam_to_gid.py (for the
vacuum beam): keep it in the same folder as them.

Reflection of a Gaussian beam on a curved PEC surface
-----------------------------------------------------
Conventions as in Scotty / scotty_to_gid.py: field ~ exp(i k s + i w.Psi.w/2),
w the transverse displacement in the (x_hat, y_hat) frame, 1/e amplitude
contour w^T Im(Psi) w = 2, Re(Psi) = k/R (> 0 diverging).

At the hit point X the surface has unit normal n (taken on the side the beam
comes from) and, in an orthonormal tangent frame (e1, e2), height
h(u) = 1/2 u^T K u along n (K > 0: concave as seen by the beam). Writing the
incident phase on the surface and asking the reflected phase to match it to
second order in u (the PEC condition E_tan,inc + E_tan,ref = 0 on the surface):

    ray            t_r = M t_i,              M = I - 2 n n^T
    basis          x_r = M x_i,  y_r = -M y_i   (keeps (x, y, t) right-handed)
    beam tensor    Psi_r = D [ Psi_i - 2 k cos(theta) J^-T K J^-1 ] D
                   J_ab = x_a . e_b  (x_1 = x_hat, x_2 = y_hat),  D = diag(1,-1)
    polarisation   E_r = -M E_i      (tangential E flips, normal E stays)
    wavevector     K_r = M K_i

Im(Psi) (the beam width at the mirror) is unchanged; Re(Psi) receives the
mirror's focusing power. For a sphere of radius R this gives the classic
tangential / sagittal focal lengths R cos(theta)/2 and R / (2 cos(theta)).
Flat mirror: Psi_r = D Psi_i D, the beam just continues "unfolded".
After the reflection the beam is propagated in vacuum with
Psi^-1(s) = Psi^-1(0) + (s/k) I (as the vacuum extensions of scotty_to_gid.py).
NOTE: if the hit lies inside the part of the beam traced by Scotty (i.e. in
the plasma), the reflected beam is still propagated in VACUUM (a warning is
printed); re-run Scotty from the reflection point if that matters.

Reflecting surfaces (section 1)
-------------------------------
``surfaces`` is a list of dicts (lengths in m):
    dict(type="plane",    point=(x,y,z), normal=(nx,ny,nz))
            optional: radius=r (disk)  or  size=(a,b), u_axis=(..) (rectangle
            of half-sizes a along u_axis and b across it)
    dict(type="sphere",   center=(x,y,z), radius=R)
            optional: axis=(..), max_angle=deg (spherical cap around axis)
    dict(type="cylinder", point=(x,y,z), axis=(ax,ay,az), radius=R)
            optional: length=L (from point along axis; default infinite)
    dict(type="mesh",     file="mirror.stl")   any file GiD exports (Files >
            Export), read by gid_surface_readers.py (same folder):
              IGES .igs  - NATIVE, no extra package: the trimmed NURBS surfaces
                           exactly as GiD wrote them, unit from the file
              STEP .stp  - gmsh (pip install gmsh)
              DXF, VTK (.vtk), GiD mesh (.msh), Exodus II (.exo), OBJ/OFF/PLY,
              STL - native;  .vtu / other meshio formats - meshio;
              Rhino .3dm - rhino3dm;  HDF5 dump, Surfaces Json - best effort
              ACIS .sat  - not readable (export IGES from GiD instead)
            Volume meshes are reduced to their skin. Optional keys:
              layers=["Inner_Tube", ...]  only these GiD layers / mesh blocks
              units="mm"     unit of MESH files (default surface_file_units);
                             IGES / STEP carry their own unit
              mesh_size=0.01 CAD surfaces: triangle size [m] (surface_mesh_size)
              reader="gmsh"  IGES through gmsh instead of the native reader
            (type="file" is the same; check what is read with
             python gid_surface_readers.py FILE -o check.stl)
    dict(type="revolved", contour="wall.dat" | "eqdsk:path/file.eqdsk")
            surface of revolution about Z of an (R, Z) contour (2 columns,
            or the limiter of a G-EQDSK file); optional phi_range=(deg0, deg1),
            n_phi=360
Every type also takes  name="..."  and  bounds=((xmin,ymin,zmin),(xmax,ymax,zmax))
(hits outside the box are ignored). Quick ones from the command line:
    --plane PX PY PZ NX NY NZ   --sphere CX CY CZ R   --cylinder PX PY PZ AX AY AZ R
    --surface-files mirror.stl ../MAST-U-Geo.stp Inner.igs
For a mesh the normal and the curvature K at the hit are obtained from a
quadric fitted to the mesh vertices around the hit (radius ~ beam footprint,
see curvature_fit_radius), so faceting does not scatter the reflected beam.
If too few vertices are available the facet normal is used and K = 0.

ERMES input, O-X conversion, ERMES cost
---------------------------------------
The ERMES launch input (ermes_input, ermes_k0, ...) and the O-X diagnostics
(ox_diagnostics, plasma files) are settings of the BEAM script's section 1
(scotty_to_gid.py or beam_to_gid.py) and work here too:
  * ermes_input = True: the incident beam is defined by the ERMES data,
  * <stem>_ermes_launch.txt: the beam start in ERMES form (Epol, |Epar|,
    |Eper|, phases), compared with the ERMES polarisation,
  * <stem>_OX_along_ray.dat: N_par, N_y, Mjolhus T along the whole folded ray
    (Scotty source: Scotty's refracted K on the traced part; vacuum source or
    after a reflection: straight lines, no refraction).
The ERMES estimate at the end gives tetrahedra, unknowns, matrix size, memory
and solve time for the FEM envelope (ermes_* parameters of section 1).

FEM envelope with reflections
-----------------------------
scotty_to_gid.py's polyhedral envelope (a chain of convex cells) is built for
every leg separately. Each leg is then cut by the tangent plane of the mirror
at its reflection point(s) (moved behind the mirror by envelope_pec_offset),
and the exact boundary of the UNION of all legs is computed (planar faces,
possibly non-convex / with holes, conforming edges). The legs overlap near
each mirror, so the union is one closed volume:
    ENVFIRST  : plane of the first beam cross-section (beam start / port)
    ENVLAST   : plane of the last beam cross-section
    ENVPEC<k> : on mirror k (assign the PEC condition)
    ENVSIDE   : everything else
CURVED mirrors (envelope_pec_conform = True, default): the ENVPEC face is then
moved onto the REAL surface. Every side face that touches it is extended (or
trimmed) in its own plane until it meets the surface, so side faces stay planar
and the volume stays closed; the ENVPEC face becomes the piece of the mirror
they enclose (STL: triangles on the surface; IGES: a bicubic B-spline patch
fitted to the surface, trimmed by those edges). This needs the mirror to cover
the whole ENVPEC face as seen along its normal: if the envelope is wider than
the mirror (e.g. a thin centre column with a large envelope_separation), the
face is left on the tangent plane with a warning -- reduce envelope_separation.
With envelope_pec_conform = False the ENVPEC face is the tangent plane; the
report gives the sagitta of the mirror over the footprint.

Parameters
----------
Section 1 below, or --name value on the command line ("-" for "_").
    python reflect_beam_to_gid.py --show-params     # values in use (+ source)
    python reflect_beam_to_gid.py --help            # options of THIS script
    python reflect_beam_to_gid.py --source-help     # options of the beam script

Usage examples
--------------
    # existing Scotty output, flat mirror at R = 0.3 m facing outwards
    python reflect_beam_to_gid.py --h5 ./Scotty_results/scotty_output_freq50.0_pol10.0_tor0.0000.h5 \\
        --extend-last 1.0 --plane 0.3 0 0 1 0 0

    # vacuum beam, 2 m, against a spherical mirror, no plot window
    python reflect_beam_to_gid.py --beam-source vacuum --length 2 \\
        --sphere 0.0 0 -0.6 1.5 --no-show-plot

    # reflect on the MAST-U CAD (needs gmsh), mesh it at 1.5 cm
    python reflect_beam_to_gid.py --surface-files ./MAST-U-Geo.stp --surface-mesh-size 0.015

    # reflect on geometry exported from GiD as IGES (no extra package needed)
    python reflect_beam_to_gid.py --surface-files ./Inner.igs --surface-mesh-size 0.01

Requires Python >= 3.9, numpy (+ h5py and scotty for the Scotty beam,
+ matplotlib for the plot, + scipy for clean trimmed CAD edges; gmsh /
meshio / rhino3dm / h5py only for the formats listed above that need them).
"""

from __future__ import annotations

import argparse
import copy
import importlib
import math
import os
import re
import struct
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

###############################################################################
# 1. USER PARAMETERS
###############################################################################
# Change the values here.
# Any of them can ALSO be changed for a single run from the command line:
#     --name value          (the same name, with "-" instead of "_")
#     --name / --no-name    for True/False parameters
# The parameters that define the BEAM are in section 1 of scotty_to_gid.py
# (beam_source = "scotty") or beam_to_gid.py (beam_source = "vacuum"); they can
# be given on this command line too (they are passed on to those scripts).
###############################################################################

_names_before_parameters = set(globals())

# ---------------------------- What to do -------------------------------------

beam_source = "vacuum"         # "scotty": Scotty beam + vacuum extensions, built as in
                               #           scotty_to_gid.py (settings: ITS section 1)
                               # "vacuum": free-space Gaussian beam of beam_to_gid.py
                               #           (settings: ITS section 1)
                               
show_plot = True               # True / False: 3D plot (beam surface, surfaces,
                               # reflections, FEM envelope); waits until closed

# ---------------------------- Reflecting (PEC) surfaces ----------------------

# List of basic surfaces - multiple [dict(..),dict(...)]
surfaces = [dict(type="cylinder", point=(0, 0, -2), axis=(0, 0, 1), radius=1, length=4, name="Column")]     
#surfaces = [dict(type="plane", point=(0.2, 0.0, 0.0), normal=(1, 0, 0), radius=1, name="flat")]
#surfaces = [dict(type="sphere", center=(0.0, 0.0, 0.0), radius=1.5, name="focusing_mirror", axis=(0, 0, 1), max_angle=180)]                            

# Surface from files - formats: [ .igs, .stp, .dxf, .vtk, .msh, .exo, .stl, ... ]
#surfaces = [dict(type="mesh", file="./Column.igs", name="mirror")]
#surfaces = [dict(type="mesh", file="./Column_Coils.iges", name="Inner", layers=["Column","Coils"])]
#surfaces = [dict(type="mesh", file="./Coils.igs", name="Coils", layers=["Coils"], mesh_size=0.01)]
#surfaces = [dict(type="mesh", file="./Wall_mesh.msh", name="wall", units="mm")]
#surfaces = [dict(type="mesh", file="./Wall.vtk", name="wall")]
#surfaces = [dict(type="revolved", contour="eqdsk:./Plasma/plasma_equilibrium.eqdsk", name="Vessel")]

# No surfaces - Uncomment this if none of the above is selected
#surfaces = []

surface_files = []             # extra files used as surfaces: any GiD export
                               # parameter for the commnad line, surfaces above 
                               # for setting surfaces in code 
                               # format (.igs .stp .dxf .vtk .msh .exo .stl ...)
                             
surface_file_units = "m"       # "m" | "cm" | "mm": unit of MESH / contour files
                               # (IGES and STEP are converted from their own unit)
                               
surface_mesh_size = 0.02       # CAD (IGES/STEP/3dm): target triangle size [m] (finer
                               # = more accurate hits/curvature, slower)
                               
surface_cad_reader = "auto"    # "auto" / "native": IGES read natively (no package)
                               # "gmsh": IGES through gmsh (STEP always uses gmsh)
                               
max_reflections = 10           # max number of reflections; at the next surface after
                               # the last one the beam is stopped
                               
skip_backward_extension = True # True : look for surfaces only from the launch point
                               #        (Scotty extend_first part never reflects)
                               # False: from the very first point of the beam
                               
ignore_hits_before = 0.0       # ignore surfaces closer than this to the launch point,
                               # along the beam [m] (e.g. the antenna itself)
                               
curvature_fit_radius = 0.0     # mesh surfaces: radius of the quadric fit around the
                               # hit [m]; 0 = automatic (beam footprint)

# ---------------------------- Reflected beam ---------------------------------

leg_point_spacing = 0.0        # spacing of the cross-sections of the reflected legs
                               # [m]; 0 = same as the incident beam (20..2000 per leg)

# ---------------------------- Geometry output (GiD) --------------------------

output_format = "both"          # "stl"   : triangle mesh   (GiD: Import > STL)
                                # "nurbs" : NURBS / planar faces (GiD: Import > IGES)
                                # "both"  | "none"
                               
units     = "m"                 # "m" | "cm" | "mm"  unit of all written geometry
caps      = True                # True : close both ends of each beam leg | False: open
out_dir   = "./GiD_input/"      # folder for ALL files written (created if needed)
file_stem = None                # start of the output file names
                                # None: from the .h5 name / vacuum beam + "_refl"
                               
envelope_factor = 1.0           # beam surface radius = factor x W (1/e amplitude);
                                # the FEM envelope keeps its separation from it
                               
stl_points_per_section = 64     # STL: points around each cross-section
stl_stride = 1                  # STL: use every N-th cross-section (1 = all)
stl_ascii = False               # STL: True = ASCII file | False = binary file
nurbs_tolerance = 1e-4          # NURBS: max deviation from any cross-section [m]
nurbs_max_sections = 600        # NURBS: max cross-sections the surface goes through
nurbs_patches = 4               # NURBS: 4 = tube in 4 quarter surfaces | 1 = one surface
nurbs_include_ray = False       # NURBS: True = add the central ray as a curve

write_surfaces = True           # True / False: <stem>_surfaces.stl, the reflecting
                                # surfaces near the beam (to check them in GiD)

# ---------------------------- FEM envelope (volume for ERMES) ----------------

envelope = True                # True : closed polyhedron around the folded beam
                               # (<stem>_envelope.igs / .stl, following output_format)
                               
envelope_separation = 0.050    # minimum distance from the beam surface to the side
                               # faces of the envelope [m]
                               
envelope_sides = 10            # side faces around the beam in each cell (>= 3)
envelope_max_cell_length = 0.0 # longest cell along the beam [m]; 0 = automatic,
                               # -1 = no subdivision (fewest faces)
                               
envelope_ray_tolerance = 0.0   # max deviation of the ray from a straight line in a
                               # cell [m]; 0 = separation / 2
                               
envelope_pec_offset = 0.0      # move the cut at each mirror BEHIND the mirror's
                               # tangent plane by this [m] (0 = cut on the plane)
                               
envelope_pec_conform = True    # True : the ENVPEC faces follow the REAL mirror surface
                               #        (side faces extended in their own planes down
                               #        to it; curved B-spline face in the IGES)
                               # False: ENVPEC faces on the tangent plane (planar)
                               
envelope_pec_mesh_size = 0.0   # point spacing along the curved ENVPEC edges / inside
                               # the face [m]; 0 = automatic (perimeter / 150)
                               
ermes_elements_per_wavelength = 10.0  # ERMES estimate printed at the end: mean
                                      # tetrahedron edge h = lambda / this value
                               
ermes_element_order = 1        # 1 | 2: nodal Lagrange elements (3 unknowns per node)
ermes_refractive_index = 1.0   # |N| in lambda_local = lambda0/|N| (> 1 if the domain
                               # holds short-wavelength plasma regions, e.g. slow X)
                               
ermes_solver = "direct"        # "direct" (MUMPS-like) | "iterative" (Krylov)
ermes_fill_constant = 15.0     # direct: factor entries ~ C * Ndof^(4/3)  (calibrate)
ermes_flop_constant = 60.0     # direct: factorisation flops ~ C * Ndof^2 (calibrate)
ermes_cores = 128              # cores used by the solver
ermes_gflops_per_core = 10.0   # sustained GFLOP/s per core
ermes_iterations = 1000        # iterative: Krylov iterations assumed for the time

# ---------------------------- Other output files -----------------------------

write_central_ray = True       # True / False: <stem>_central_ray.dat (all legs)
write_efield = True            # True / False: E-field and k at both beam ends

efield_file = "E_final.dat"    # name of that end-point file

write_field_arrays = True      # True / False: E and k at every ray point (all legs)
write_ermes_launch = True      # True / False: <stem>_ermes_launch.txt (beam start in
                               # ERMES form, polarisation compared with ermes_epol...)
                               
write_reflections = True       # True / False: <stem>_reflections.txt (hit points,
                               # normals, curvature, Psi / E before and after, ...)
                               
plot_points_around = 36        # plot: points around each cross-section (the surface
                               # always uses ALL cross-sections)

###############################################################################
# (end of user parameters -- nothing below needs to be edited)
###############################################################################

USER_PARAMETERS = [n for n in globals()
                   if n not in _names_before_parameters and not n.startswith("_")]
_PARAMETER_DEFAULTS = {n: globals()[n] for n in USER_PARAMETERS}

PARAMETER_CHOICES = {
    "beam_source": ["scotty", "vacuum"],
    "output_format": ["stl", "nurbs", "both", "none"],
    "units": ["m", "cm", "mm"],
    "surface_file_units": ["m", "cm", "mm"],
    "surface_cad_reader": ["auto", "native", "gmsh"],
    "nurbs_patches": [1, 4],
    "ermes_solver": ["direct", "iterative"],
    "ermes_element_order": [1, 2],
}
PARAMETER_ALIASES = {
    "output_format": ["--format"],
    "show_plot": ["--plot"],
    "nurbs_tolerance": ["--nurbs-tol"],
}
PARAMETER_TYPES = {"file_stem": str}
NO_COMMAND_LINE = {"surfaces"}          # edited in section 1 only

SPEED_OF_LIGHT = 299_792_458.0
UNIT_SCALE = {"m": 1.0, "cm": 100.0, "mm": 1000.0}
MIN_LEG_LENGTH = 1e-4                    # m: hits closer to a leg start are ignored
_HERE = os.path.dirname(os.path.abspath(__file__))


def _import_first(names, what):
    """Import the first module of ``names`` found (this script's folder first)."""
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    for n in names:
        try:
            return importlib.import_module(n)
        except ModuleNotFoundError as err:
            if err.name != n:
                raise
    sys.exit(f"Cannot import {what} ({' / '.join(n + '.py' for n in names)}); "
             f"put it in the same folder as this script ({_HERE}).")


s2g = _import_first(("scotty_to_gid", "e_scotty_to_gid", "scotty_beam_to_gid"),
                    "scotty_to_gid.py")
Beam = s2g.Beam


def default_parameters():
    """Section-1 values as a namespace (for use from another script)."""
    ns = argparse.Namespace(**{k: (list(v) if isinstance(v, list) else v)
                               for k, v in _PARAMETER_DEFAULTS.items()})
    ns.show_params = ns.source_help = False
    ns.plane = ns.sphere = ns.cylinder = None
    return ns


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _plane_frame(n):
    return s2g._plane_frame(np.asarray(n, float))


# ##########################################################################
# 2. REFLECTING SURFACES
# ##########################################################################

class Surface:
    """Base class. A surface finds the first crossing of a polyline and gives
    the local geometry (normal, tangent frame, curvature) at a point."""
    kind = "surface"

    def __init__(self, name, bounds=None):
        self.name = name
        self.bounds = None if bounds is None else np.asarray(bounds, float).reshape(2, 3)

    def _in_bounds(self, X):
        if self.bounds is None:
            return np.ones(len(X), bool)
        lo, hi = self.bounds
        return np.all((X >= lo - 1e-12) & (X <= hi + 1e-12), axis=1)

    # analytic surfaces: (m, k) candidate parameters along each segment
    def _roots(self, P0, D):
        raise NotImplementedError

    def _extent_ok(self, X):
        return np.ones(len(X), bool)

    def first_hit(self, P0, P1):
        """Earliest crossing of the segments P0[i] -> P1[i], taken in order.
        Returns (i, t, X, data) with X = P0[i] + t (P1[i]-P0[i]), or None."""
        D = P1 - P0
        T = self._roots(P0, D)
        m, k = T.shape
        ok = np.isfinite(T) & (T >= -1e-12) & (T <= 1 + 1e-12)
        if not ok.any():
            return None
        X = P0[:, None, :] + np.where(ok, T, 0)[..., None] * D[:, None, :]
        flat = X.reshape(-1, 3)
        okf = ok.reshape(-1).copy()
        okf[okf] &= self._extent_ok(flat[okf]) & self._in_bounds(flat[okf])
        if not okf.any():
            return None
        key = np.where(okf.reshape(m, k), np.arange(m)[:, None] + np.clip(T, 0, 1), np.inf)
        i, j = np.unravel_index(np.argmin(key), key.shape)
        return int(i), float(np.clip(T[i, j], 0, 1)), X[i, j].copy(), None

    def rough_normal(self, X, data):
        return self.local_geometry(X, data, 0.0)[0]

    def local_geometry(self, X, data, fit_radius):
        """(n, e1, e2, K, note): unit normal (any side), orthonormal tangent
        frame, 2x2 curvature tensor such that the surface is
        X + u1 e1 + u2 e2 + 1/2 u^T K u n."""
        raise NotImplementedError

    def display_mesh(self, box):
        """(V, T) triangles to draw / export inside box (2, 3)."""
        return np.zeros((0, 3)), np.zeros((0, 3), int)


class PlaneSurface(Surface):
    kind = "plane"

    def __init__(self, name, point, normal, radius=None, size=None, u_axis=None, bounds=None):
        super().__init__(name, bounds)
        self.p = np.asarray(point, float)
        self.n = _unit(normal)
        self.radius = radius
        self.size = None if size is None else np.asarray(size, float).reshape(2)
        u = None
        if u_axis is not None:
            u = np.asarray(u_axis, float) - np.dot(u_axis, self.n) * self.n
            u = _unit(u) if np.linalg.norm(u) > 1e-9 else None
        self.u = u if u is not None else _plane_frame(self.n)[0]
        self.v = np.cross(self.n, self.u)

    def _roots(self, P0, D):
        den = D @ self.n
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (np.dot(self.p, self.n) - P0 @ self.n) / den
        t[np.abs(den) < 1e-300] = np.nan
        return t[:, None]

    def _extent_ok(self, X):
        r = X - self.p
        ok = np.ones(len(X), bool)
        if self.radius is not None:
            ok &= np.linalg.norm(r - np.outer(r @ self.n, self.n), axis=1) <= self.radius
        if self.size is not None:
            ok &= (np.abs(r @ self.u) <= self.size[0]) & (np.abs(r @ self.v) <= self.size[1])
        return ok

    def local_geometry(self, X, data, fit_radius):
        return self.n.copy(), self.u.copy(), self.v.copy(), np.zeros((2, 2)), "exact (plane)"

    def display_mesh(self, box):
        c = box.mean(0)
        c = c - np.dot(c - self.p, self.n) * self.n
        half = 0.5 * np.linalg.norm(box[1] - box[0])
        if self.radius is not None:
            c, a = self.p, np.full(2, self.radius)
            ang = np.linspace(0, 2 * np.pi, 49)[:-1]
            ring = c + self.radius * (np.cos(ang)[:, None] * self.u + np.sin(ang)[:, None] * self.v)
            V = np.vstack([c, ring])
            T = np.array([[0, 1 + i, 1 + (i + 1) % 48] for i in range(48)])
            return V, T
        if self.size is not None:
            c, (a, b) = self.p, self.size
        else:
            a = b = half
        V = np.array([c + sa * a * self.u + sb * b * self.v
                      for sa, sb in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
        return V, np.array([[0, 1, 2], [0, 2, 3]])


class SphereSurface(Surface):
    kind = "sphere"

    def __init__(self, name, center, radius, axis=None, max_angle=None, bounds=None):
        super().__init__(name, bounds)
        self.c = np.asarray(center, float)
        self.R = float(radius)
        self.axis = None if axis is None else _unit(axis)
        self.max_angle = max_angle

    def _roots(self, P0, D):
        W = P0 - self.c
        a = np.einsum("ij,ij->i", D, D)
        b = 2 * np.einsum("ij,ij->i", D, W)
        cc = np.einsum("ij,ij->i", W, W) - self.R ** 2
        disc = b * b - 4 * a * cc
        sq = np.sqrt(np.where(disc >= 0, disc, np.nan))
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.stack([(-b - sq) / (2 * a), (-b + sq) / (2 * a)], 1)

    def _extent_ok(self, X):
        if self.axis is None or self.max_angle is None:
            return np.ones(len(X), bool)
        cosang = ((X - self.c) / self.R) @ self.axis
        return cosang >= np.cos(np.deg2rad(self.max_angle)) - 1e-12

    def local_geometry(self, X, data, fit_radius):
        n = _unit(X - self.c)
        e1, e2 = _plane_frame(n)
        return n, e1, e2, -np.eye(2) / self.R, "exact (sphere)"

    def display_mesh(self, box):
        th = np.linspace(0, np.pi, 37)
        ph = np.linspace(0, 2 * np.pi, 73)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        if self.axis is not None:
            a = self.axis
            u, v = _plane_frame(a)
        else:
            a, u, v = np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
        dirs = (np.cos(TH)[..., None] * a + np.sin(TH)[..., None] *
                (np.cos(PH)[..., None] * u + np.sin(PH)[..., None] * v))
        V = (self.c + self.R * dirs).reshape(-1, 3)
        return _grid_triangles(V, len(th), len(ph), box, self._extent_ok)


class CylinderSurface(Surface):
    kind = "cylinder"

    def __init__(self, name, point, axis, radius, length=None, bounds=None):
        super().__init__(name, bounds)
        self.p = np.asarray(point, float)
        self.e = _unit(axis)
        self.R = float(radius)
        self.length = length

    def _roots(self, P0, D):
        W = P0 - self.p
        Dp = D - np.outer(D @ self.e, self.e)
        Wp = W - np.outer(W @ self.e, self.e)
        a = np.einsum("ij,ij->i", Dp, Dp)
        b = 2 * np.einsum("ij,ij->i", Dp, Wp)
        cc = np.einsum("ij,ij->i", Wp, Wp) - self.R ** 2
        disc = b * b - 4 * a * cc
        sq = np.sqrt(np.where((disc >= 0) & (a > 1e-300), disc, np.nan))
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.stack([(-b - sq) / (2 * a), (-b + sq) / (2 * a)], 1)

    def _extent_ok(self, X):
        if self.length is None:
            return np.ones(len(X), bool)
        z = (X - self.p) @ self.e
        return (z >= -1e-12) & (z <= self.length + 1e-12)

    def local_geometry(self, X, data, fit_radius):
        r = X - self.p
        r = r - np.dot(r, self.e) * self.e
        n = _unit(r)
        e1 = self.e.copy()
        e2 = np.cross(n, e1)
        return n, e1, e2, np.diag([0.0, -1.0 / self.R]), "exact (cylinder)"

    def display_mesh(self, box):
        if self.length is not None:
            z0, z1 = 0.0, self.length
        else:
            corners = np.array([[box[i, 0], box[j, 1], box[k, 2]]
                                for i in (0, 1) for j in (0, 1) for k in (0, 1)])
            z = (corners - self.p) @ self.e
            z0, z1 = z.min(), z.max()
        zz = np.linspace(z0, z1, 25)
        ph = np.linspace(0, 2 * np.pi, 73)
        u, v = _plane_frame(self.e)
        Z, PH = np.meshgrid(zz, ph, indexing="ij")
        V = (self.p + Z[..., None] * self.e + self.R *
             (np.cos(PH)[..., None] * u + np.sin(PH)[..., None] * v)).reshape(-1, 3)
        return _grid_triangles(V, len(zz), len(ph), box, None)


def _grid_triangles(V, n0, n1, box, extent_ok):
    idx = np.arange(n0 * n1).reshape(n0, n1)
    a, b = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel()
    c, d = idx[1:, :-1].ravel(), idx[1:, 1:].ravel()
    T = np.vstack([np.stack([a, b, d], 1), np.stack([a, d, c], 1)])
    pad = 0.05 * np.linalg.norm(box[1] - box[0])
    inside = np.all((V >= box[0] - pad) & (V <= box[1] + pad), axis=1)
    if extent_ok is not None:
        inside &= extent_ok(V)
    T = T[np.all(inside[T], axis=1)]
    return V, T


class MeshSurface(Surface):
    """Triangulated surface (STL, meshed CAD, revolved contour)."""
    kind = "mesh"

    def __init__(self, name, V, T, bounds=None, source=""):
        super().__init__(name, bounds)
        tri = V[T]
        e1, e2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
        cr = np.cross(e1, e2)
        area2 = np.linalg.norm(cr, axis=1)
        good = area2 > 1e-14 * max(np.ptp(V, axis=0).max(), 1e-9) ** 2
        self.V, self.T = V, T[good]
        tri = V[self.T]
        self.v0 = tri[:, 0]
        self.e1, self.e2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
        self.fn = _unit(np.cross(self.e1, self.e2))
        self.tmin, self.tmax = tri.min(1), tri.max(1)
        self.cen = tri.mean(1)
        self.edge = float(np.median(np.linalg.norm(self.e1, axis=1))) if len(tri) else 0.0
        self.source = source

    def first_hit(self, P0, P1, chunk=32):
        m = len(P0)
        tol = 1e-9 * max(1.0, float(np.abs(P0).max()))
        for c0 in range(0, m, chunk):
            c1 = min(c0 + chunk, m)
            lo = np.minimum(P0[c0:c1], P1[c0:c1]).min(0) - tol
            hi = np.maximum(P0[c0:c1], P1[c0:c1]).max(0) + tol
            cand = np.nonzero(np.all(self.tmax >= lo, 1) & np.all(self.tmin <= hi, 1))[0]
            if not len(cand):
                continue
            for i in range(c0, c1):
                slo = np.minimum(P0[i], P1[i]) - tol
                shi = np.maximum(P0[i], P1[i]) + tol
                sub = cand[np.all(self.tmax[cand] >= slo, 1) & np.all(self.tmin[cand] <= shi, 1)]
                if not len(sub):
                    continue
                res = self._segment_hits(P0[i], P1[i] - P0[i], sub)
                if res is not None:
                    t, j = res
                    return i, t, P0[i] + t * (P1[i] - P0[i]), int(j)
        return None

    def _segment_hits(self, p0, d, sub):
        """Moller-Trumbore against the triangles ``sub``: earliest (t, tri)."""
        e1, e2, v0 = self.e1[sub], self.e2[sub], self.v0[sub]
        pvec = np.cross(d, e2)
        det = np.einsum("ij,ij->i", e1, pvec)
        scale = np.linalg.norm(d) * np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1)
        ok = np.abs(det) > 1e-12 * scale
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        tvec = p0 - v0
        u = np.einsum("ij,ij->i", tvec, pvec) * inv
        qvec = np.cross(tvec, e1)
        v = (qvec @ d) * inv
        t = np.einsum("ij,ij->i", e2, qvec) * inv
        eps = 1e-10
        ok &= (u >= -eps) & (v >= -eps) & (u + v <= 1 + eps) & (t >= 0) & (t <= 1)
        if not ok.any():
            return None
        idx = np.nonzero(ok)[0]
        X = p0 + t[idx, None] * d
        keep = self._in_bounds(X)
        if not keep.any():
            return None
        idx = idx[keep]
        k = idx[np.argmin(t[idx])]
        return float(t[k]), sub[k]

    def rough_normal(self, X, data):
        return self.fn[data].copy()

    def local_geometry(self, X, data, fit_radius):
        nf = self.fn[data]
        r0 = max(fit_radius, 2.5 * self.edge, 1e-6)
        for attempt in range(3):
            r = r0 * 2 ** attempt
            res = self._quadric_fit(X, data, nf, r)
            if res is not None:
                return res
        e1, e2 = _plane_frame(nf)
        return (nf.copy(), e1, e2, np.zeros((2, 2)),
                f"FACET normal, K = 0 (not enough mesh vertices within {r:.3g} m to fit "
                "the curvature; refine the mesh or set curvature_fit_radius)")

    def _quadric_fit(self, X, j, nf, r):
        # triangles near X, connected to the hit triangle (same sheet)
        near = np.nonzero(np.linalg.norm(self.cen - X, axis=1) <= r + 2 * self.edge)[0]
        near_set = set(near.tolist()) | {j}
        by_vertex = defaultdict(list)
        for t in near_set:
            for v in self.T[t]:
                by_vertex[v].append(t)
        comp, stack = {j}, [j]
        while stack:
            t = stack.pop()
            for v in self.T[t]:
                for t2 in by_vertex[v]:
                    if t2 not in comp:
                        comp.add(t2)
                        stack.append(t2)
        vid = np.unique(self.T[list(comp)].ravel())
        P = self.V[vid]
        dist = np.linalg.norm(P - X, axis=1)
        P, dist = P[dist <= r], dist[dist <= r]
        if len(P) < 8:
            return None
        e1, e2 = _plane_frame(nf)
        d = P - X
        u1, u2, h = d @ e1, d @ e2, d @ nf
        G = np.stack([np.ones_like(u1), u1, u2, 0.5 * u1 ** 2, u1 * u2, 0.5 * u2 ** 2], 1)
        w = np.sqrt(np.clip(1 - (dist / r) ** 2, 0.05, 1))
        Gs = G * w[:, None] / np.array([1, r, r, r * r, r * r, r * r])  # scaled columns
        sv = np.linalg.svd(Gs, compute_uv=False)
        if sv[-1] < 1e-6 * sv[0]:
            return None
        coef = np.linalg.lstsq(G * w[:, None], h * w, rcond=None)[0]
        _, g1, g2, a, b, c = coef
        # surface z = f(u): normal, first and second fundamental forms at u = 0
        g = np.array([g1, g2])
        q = math.sqrt(1 + g @ g)
        n = _unit(nf - g1 * e1 - g2 * e2)
        Tm = np.stack([e1 + g1 * nf, e2 + g2 * nf], 1)            # 3x2 tangents
        L = np.array([[a, b], [b, c]]) / q
        Q, Rm = np.linalg.qr(Tm)
        Ri = np.linalg.inv(Rm)
        K = Ri.T @ L @ Ri
        K = 0.5 * (K + K.T)
        t1, t2 = Q[:, 0], Q[:, 1]
        if np.dot(np.cross(t1, t2), n) < 0:
            t2 = -t2
            K[0, 1] = K[1, 0] = -K[0, 1]
        tilt = np.degrees(np.arccos(np.clip(np.dot(n, nf), -1, 1)))
        return n, t1, t2, K, (f"quadric fit to {len(P)} mesh vertices within {r:.3g} m "
                              f"(normal {tilt:.2f} deg from the facet normal)")

    def display_mesh(self, box):
        pad = 0.05 * np.linalg.norm(box[1] - box[0])
        keep = np.all((self.cen >= box[0] - pad) & (self.cen <= box[1] + pad), axis=1)
        return self.V, self.T[keep]


# ---------------------------- loading surfaces ------------------------------

def merge_vertices(tris):
    """(nt, 3, 3) triangle soup -> (V, T) with coincident vertices merged."""
    P = tris.reshape(-1, 3).astype(float)
    if not len(P):
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    tol = 1e-7 * max(float(np.ptp(P, axis=0).max()), 1e-9)
    q = np.round(P / tol).astype(np.int64)
    _, first, inv = np.unique(q, axis=0, return_index=True, return_inverse=True)
    return P[first], inv.reshape(-1, 3)


def read_stl(path):
    """ASCII or binary STL -> (nt, 3, 3) triangles (file units)."""
    with open(path, "rb") as f:
        data = f.read()
    if len(data) >= 84:
        n = struct.unpack("<I", data[80:84])[0]
        if 84 + 50 * n == len(data):
            rec = np.frombuffer(data, dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)),
                                                     ("a", "<u2")]), count=n, offset=84)
            return rec["v"].astype(float)
    txt = data.decode("latin-1")
    num = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
    v = re.findall(r"vertex\s+" + num + r"\s+" + num + r"\s+" + num, txt)
    if not v or len(v) % 3:
        sys.exit(f"Could not read STL file {path}")
    return np.array(v, float).reshape(-1, 3, 3)


def _readers():
    return _import_first(("gid_surface_readers", "e_gid_surface_readers"),
                         "gid_surface_readers.py")


def read_surface_file(path, cfg, region, units=None, mesh_size=None, reader=None, layers=None):
    """Any GiD export format -> (nt, 3, 3) triangles in m (gid_surface_readers.py)."""
    gsr = _readers()
    try:
        return gsr.read_surface_file(
            path, mesh_size=cfg.surface_mesh_size if mesh_size is None else mesh_size,
            region=region, units=units or cfg.surface_file_units,
            cad_reader=reader or getattr(cfg, "surface_cad_reader", "auto"), layers=layers)
    except gsr.SurfaceReadError as err:
        sys.exit(f"Surface file {path}: {err}")


def read_contour(spec, scale):
    """(R, Z) contour [m] from a 2-column file or 'eqdsk:<file>' (limiter)."""
    if str(spec).lower().startswith("eqdsk:"):
        eqmod = _import_first(("eqdsk_to_scotty", "e_eqdsk_to_scotty"), "eqdsk_to_scotty.py")
        eq = eqmod.read_geqdsk(spec.split(":", 1)[1])
        RZ = np.stack([np.asarray(eq.rlim, float), np.asarray(eq.zlim, float)], 1)
        if len(RZ) < 3:
            sys.exit(f"{spec}: the G-EQDSK file has no limiter / wall contour.")
        return RZ
    RZ = np.loadtxt(spec, comments=("#", "!", "%"), usecols=(0, 1)) / scale
    return RZ


def revolved_mesh(RZ, phi_range=(0.0, 360.0), n_phi=360):
    """Surface of revolution about Z of the (R, Z) contour -> (V, T)."""
    RZ = np.asarray(RZ, float)
    if np.linalg.norm(RZ[0] - RZ[-1]) > 1e-9:
        RZ = np.vstack([RZ, RZ[:1]])                       # close the contour
    full = abs((phi_range[1] - phi_range[0]) % 360) < 1e-9 and phi_range[1] != phi_range[0]
    nph = int(n_phi) if full else int(n_phi) + 1
    ph = np.deg2rad(phi_range[0] + (phi_range[1] - phi_range[0]) * np.arange(nph) / int(n_phi))
    nc = len(RZ)
    V = np.stack([RZ[None, :, 0] * np.cos(ph)[:, None], RZ[None, :, 0] * np.sin(ph)[:, None],
                  np.broadcast_to(RZ[None, :, 1], (nph, nc))], -1).reshape(-1, 3)
    idx = np.arange(nph * nc).reshape(nph, nc)
    nxt = np.roll(idx, -1, axis=0) if full else idx[1:]
    cur = idx if full else idx[:-1]
    a, b = cur[:, :-1].ravel(), cur[:, 1:].ravel()
    c, d = nxt[:, :-1].ravel(), nxt[:, 1:].ravel()
    T = np.vstack([np.stack([a, b, d], 1), np.stack([a, d, c], 1)])
    return V, T


def _region_filter(tris, region):
    if region is None or not len(tris):
        return tris
    lo, hi = tris.min(1), tris.max(1)
    keep = np.all(hi >= region[0], 1) & np.all(lo <= region[1], 1)
    return tris[keep]


def make_surface(spec, cfg, region):
    spec = dict(spec)
    typ = str(spec.pop("type", "")).lower()
    name = spec.pop("name", None)
    bounds = spec.pop("bounds", None)
    try:
        if typ == "plane":
            return PlaneSurface(name or "plane", bounds=bounds, **spec)
        if typ == "sphere":
            return SphereSurface(name or "sphere", bounds=bounds, **spec)
        if typ == "cylinder":
            return CylinderSurface(name or "cylinder", bounds=bounds, **spec)
    except TypeError as err:
        sys.exit(f"Surface {typ!r}: {err}")
    if not typ and "file" in spec:
        typ = "mesh"
    if typ in ("mesh", "file", "cad"):
        path = spec.pop("file", None)
        if not path or not os.path.isfile(path):
            sys.exit(f"Mesh surface: file not found: {path!r}")
        tris = read_surface_file(path, cfg, region, units=spec.pop("units", None),
                                 mesh_size=spec.pop("mesh_size", None),
                                 reader=spec.pop("reader", None), layers=spec.pop("layers", None))
        if spec:
            print(f"WARNING: surface '{name}': unknown keys ignored: {sorted(spec)}")
        n_all = len(tris)
        tris = _region_filter(tris, region)
        V, T = merge_vertices(tris)
        print(f"Surface '{name or os.path.basename(path)}': {path}, {n_all} triangles, "
              f"{len(T)} within reach of the beam")
        if not len(T):
            print(f"WARNING: surface '{name or path}' is ignored: none of it is within reach "
                  "of the beam (check its units / layers / position)")
            return None
        return MeshSurface(name or os.path.basename(path), V, T, bounds=bounds, source=path)
    if typ == "revolved":
        RZ = read_contour(spec.pop("contour"), UNIT_SCALE[spec.pop("units", cfg.surface_file_units)])
        V, T = revolved_mesh(RZ, spec.pop("phi_range", (0.0, 360.0)), spec.pop("n_phi", 360))
        tris = _region_filter(V[T], region)
        V, T = merge_vertices(tris)
        print(f"Surface '{name or 'revolved'}': contour of {len(RZ)} points revolved, "
              f"{len(T)} triangles within reach of the beam")
        return MeshSurface(name or "revolved", V, T, bounds=bounds, source="revolved contour")
    sys.exit(f"Unknown surface type {typ!r} (plane, sphere, cylinder, mesh/file, revolved)")


def surface_specs(cfg):
    specs = [dict(s) for s in (cfg.surfaces or [])]
    for f in cfg.surface_files or []:
        specs.append(dict(type="mesh", file=f))
    for p in cfg.plane or []:
        specs.append(dict(type="plane", point=p[:3], normal=p[3:]))
    for p in cfg.sphere or []:
        specs.append(dict(type="sphere", center=p[:3], radius=p[3]))
    for p in cfg.cylinder or []:
        specs.append(dict(type="cylinder", point=p[:3], axis=p[3:6], radius=p[6]))
    for k, s in enumerate(specs):
        s.setdefault("name", f"{s.get('type', 'surface')}{k + 1}")
    return specs


# ##########################################################################
# 3. THE INCIDENT BEAM (Scotty or vacuum)
# ##########################################################################

def source_module(cfg):
    if cfg.beam_source == "scotty":
        return s2g
    return _import_first(("beam_to_gid", "e_beam_to_gid", "vacuum_beam_to_gid"),
                         "beam_to_gid.py")


def incident_beam(cfg, scfg):
    """Build the incident beam with the source script's own functions.
    Returns dict(beam, k0, launch_index, scotty_range, ends, stem, straight)."""
    if cfg.beam_source == "scotty":
        if scfg.run_scotty:
            h5 = s2g.run_scotty_simulation(scfg)
        else:
            h5 = scfg.scotty_h5_file or s2g.newest_scotty_output(scfg.scotty_output_path)
            print(f"Using existing Scotty output: {h5}")
        if not os.path.isfile(h5):
            sys.exit(f"File not found: {h5}")
        d = s2g.load_scotty_beam(h5)
        if hasattr(s2g, "rotate_scotty_output"):
            s2g.rotate_scotty_output(d, s2g.scotty_zeta_rotation(scfg))
        E, K = s2g.cartesian_fields(d)
        if E is None and (cfg.write_efield or cfg.write_field_arrays):
            print("WARNING: e_hat / K / q_R / q_zeta missing -- no E-field files.")
        beam = Beam(points=d["beam"], xhat=d["xhat"], yhat=d["yhat"], psi_xx=d["psi_xx"],
                    psi_xy=d["psi_xy"], psi_yy=d["psi_yy"], distance=d["distance"], E=E, K=K)
        first, last, n_sc = beam.points[0].copy(), beam.points[-1].copy(), beam.n
        k0 = d["omega"] / SPEED_OF_LIGHT if d["omega"] else None
        beam = s2g.extend_scotty_beam(beam, d, k0, scfg)
        if hasattr(s2g, "report_beam_ends"):
            s2g.report_beam_ends(beam, k0, d, scfg)
        i0 = int(np.argmin(np.linalg.norm(beam.points - first, axis=1)))
        stem = (cfg.file_stem or scfg.file_stem or
                os.path.splitext(os.path.basename(h5))[0].replace(".", "_") + "_refl")
        return dict(beam=beam, k0=k0, launch_index=i0, scotty_range=(i0, i0 + n_sc - 1),
                    ends=(first, last), stem=stem, straight=False)

    vb = source_module(cfg)
    if scfg.propagation_length <= 0:
        sys.exit("propagation_length must be > 0.")
    gb = vb.gaussian_beam_from_config(scfg)
    vb.describe_beam(gb, scfg.propagation_length)
    b = vb.sample_beam(gb, scfg.propagation_length, scfg.n_points)
    beam = Beam(points=b.points, xhat=b.xhat, yhat=b.yhat, psi_xx=b.psi_xx, psi_xy=b.psi_xy,
                psi_yy=b.psi_yy, distance=b.distance, E=b.E, K=b.K)
    stem = cfg.file_stem or (f"vacuum_beam_{scfg.frequency / 1e9:g}GHz_"
                             f"L{scfg.propagation_length:g}m").replace(".", "_") + "_refl"
    return dict(beam=beam, k0=gb.k, launch_index=0, scotty_range=None,
                ends=(beam.points[0].copy(), beam.points[-1].copy()), stem=stem, straight=True)


# ##########################################################################
# 4. REFLECTIONS
# ##########################################################################

@dataclass
class RayState:
    """Everything needed to continue the beam from one point."""
    point: np.ndarray
    t: np.ndarray                       # unit direction of the ray
    xhat: np.ndarray
    yhat: np.ndarray
    psi: Optional[np.ndarray]           # 2x2 complex in (xhat, yhat), or None
    E: Optional[np.ndarray]
    K: Optional[np.ndarray]
    distance: float


@dataclass
class Reflection:
    number: int
    surface: str
    point: np.ndarray
    normal: np.ndarray                  # unit, towards the incoming beam
    t_in: np.ndarray
    t_out: np.ndarray
    incidence_deg: float
    e1: np.ndarray
    e2: np.ndarray
    curvature: np.ndarray               # 2x2 in (e1, e2), w.r.t. ``normal``
    geometry_note: str
    distance: float                     # along the beam (beam's distance coordinate)
    path_length: float                  # geometric path length from the start
    psi_in: Optional[np.ndarray] = None
    psi_out: Optional[np.ndarray] = None
    E_in: Optional[np.ndarray] = None
    E_out: Optional[np.ndarray] = None
    footprint: Optional[np.ndarray] = None   # semi-axes of the 1/e spot on the surface
    sagitta: float = 0.0
    focal_lengths: Optional[np.ndarray] = None
    warnings: list = field(default_factory=list)
    surface_obj: object = None               # the Surface that was hit
    pec_conformed: bool = False              # ENVPEC face follows the real surface


def psi_matrix(beam, i):
    if beam.psi_xx is None:
        return None
    return np.array([[beam.psi_xx[i], beam.psi_xy[i]], [beam.psi_xy[i], beam.psi_yy[i]]],
                    dtype=complex)


def propagate_psi(psi, s, k0):
    """Vacuum law Psi^-1(s) = Psi^-1(0) + (s/k0) I; s array -> (len(s), 2, 2)."""
    s = np.atleast_1d(np.asarray(s, float))
    if psi is None or k0 is None:
        return None
    Minv = np.linalg.inv(psi)[None] + (s / k0)[:, None, None] * np.eye(2)[None]
    out = np.linalg.inv(Minv)
    return 0.5 * (out + out.transpose(0, 2, 1))


def orthonormal_basis(xhat, yhat, t):
    """Gram-Schmidt (x, y) against t, keeping the sign of y."""
    x = _unit(xhat - np.dot(xhat, t) * t)
    y = yhat - np.dot(yhat, t) * t - np.dot(yhat, x) * x
    return x, _unit(y)


def state_at(beam, j, X, t, k0):
    """State at X on segment j (vacuum propagation from point j)."""
    ds = float(np.linalg.norm(X - beam.points[j]))
    psi = psi_matrix(beam, j)
    if psi is not None:
        psi = propagate_psi(psi, ds, k0)[0] if k0 else psi
    x, y = orthonormal_basis(beam.xhat[j], beam.yhat[j], t)
    E = None if beam.E is None else beam.E[j].copy()
    K = None if beam.K is None else beam.K[j].copy()
    return RayState(point=np.asarray(X, float), t=t, xhat=x, yhat=y, psi=psi, E=E, K=K,
                    distance=float(beam.distance[j] + ds))


def beam_from_states(points, xhat, yhat, psi, distance, E, K):
    n = len(points)
    if psi is None:
        pxx = pxy = pyy = None
    else:
        pxx, pxy, pyy = psi[:, 0, 0], psi[:, 0, 1], psi[:, 1, 1]
    return Beam(points=np.asarray(points, float), xhat=np.asarray(xhat, float),
                yhat=np.asarray(yhat, float), psi_xx=pxx, psi_xy=pxy, psi_yy=pyy,
                distance=np.asarray(distance, float),
                E=None if E is None else np.asarray(E), K=None if K is None else np.asarray(K))


def vacuum_leg(st: RayState, length, n_points, k0):
    """Straight Gaussian beam in vacuum from ``st`` over ``length``."""
    s = np.linspace(0.0, length, max(int(n_points), 2))
    n = len(s)
    psi = propagate_psi(st.psi, s, k0) if st.psi is not None else None
    K = None
    if st.K is not None:
        K = np.tile(st.K, (n, 1))
    elif k0:
        K = np.tile(k0 * st.t, (n, 1))
    return beam_from_states(st.point + s[:, None] * st.t, np.tile(st.xhat, (n, 1)),
                            np.tile(st.yhat, (n, 1)), psi, st.distance + s,
                            None if st.E is None else np.tile(st.E, (n, 1)), K)


def _take(beam, sl):
    g = lambda a: None if a is None else a[sl]
    return Beam(points=beam.points[sl], xhat=beam.xhat[sl], yhat=beam.yhat[sl],
                psi_xx=g(beam.psi_xx), psi_xy=g(beam.psi_xy), psi_yy=g(beam.psi_yy),
                distance=beam.distance[sl], E=g(beam.E), K=g(beam.K))


def _append_state(beam, st: RayState):
    """Beam + one more point described by ``st``."""
    cat = lambda a, v: None if a is None or v is None else np.concatenate([a, [v]])
    psi = st.psi
    return Beam(points=np.vstack([beam.points, st.point]),
                xhat=np.vstack([beam.xhat, st.xhat]), yhat=np.vstack([beam.yhat, st.yhat]),
                psi_xx=cat(beam.psi_xx, None if psi is None else psi[0, 0]),
                psi_xy=cat(beam.psi_xy, None if psi is None else psi[0, 1]),
                psi_yy=cat(beam.psi_yy, None if psi is None else psi[1, 1]),
                distance=np.concatenate([beam.distance, [st.distance]]),
                E=cat(beam.E, st.E), K=cat(beam.K, st.K))


def truncate_leg(beam, j, st):
    """Points 0..j of ``beam`` followed by the hit point (state st)."""
    part = _take(beam, slice(0, j + 1))
    if np.linalg.norm(st.point - beam.points[j]) < 1e-12 and j > 0:
        part = _take(beam, slice(0, j))
    return _append_state(part, st)


def densify(beam, min_points, k0):
    """Subdivide segments (vacuum propagation) until the beam has min_points."""
    if beam.n >= min_points:
        return beam
    m = int(np.ceil((min_points - 1) / max(beam.n - 1, 1)))
    out = None
    for j in range(beam.n - 1):
        seg = beam.points[j + 1] - beam.points[j]
        L = np.linalg.norm(seg)
        t = seg / L if L > 0 else _unit(np.cross(beam.xhat[j], beam.yhat[j]))
        for f in np.arange(m) / m:
            st = state_at(beam, j, beam.points[j] + f * seg, t, k0)
            st.xhat, st.yhat = beam.xhat[j], beam.yhat[j]
            out = beam_from_states([st.point], [st.xhat], [st.yhat],
                                   None if st.psi is None else st.psi[None], [st.distance],
                                   None if st.E is None else [st.E],
                                   None if st.K is None else [st.K]) if out is None \
                else _append_state(out, st)
    return _append_state(out, RayState(point=beam.points[-1], t=t, xhat=beam.xhat[-1],
                                       yhat=beam.yhat[-1], psi=psi_matrix(beam, -1),
                                       E=None if beam.E is None else beam.E[-1],
                                       K=None if beam.K is None else beam.K[-1],
                                       distance=float(beam.distance[-1])))


def arclength(points):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])


def first_surface_hit(beam, surfaces, j0, X0, straight):
    """Earliest crossing of the beam's central ray, from X0 (on segment j0)
    on, with any surface. Returns (surface, segment j, X, data) or None."""
    if not surfaces:
        return None
    if straight:
        P0, P1 = X0[None], beam.points[-1][None]
    else:
        P0 = np.vstack([X0, beam.points[j0 + 1:-1]])
        P1 = beam.points[j0 + 1:]
    best = None
    for srf in surfaces:
        h = srf.first_hit(P0, P1)
        if h is None:
            continue
        i, t, X, data = h
        key = (i + t) if not straight else t
        if best is None or key < best[0]:
            best = (key, srf, i, X, data)
    if best is None:
        return None
    _, srf, i, X, data = best
    if straight:
        s = arclength(beam.points)
        sx = np.linalg.norm(X - beam.points[0])
        j = int(np.clip(np.searchsorted(s, sx, side="right") - 1, 0, beam.n - 2))
    else:
        j = j0 + i
    return srf, j, X, data


def reflect_state(st: RayState, geom, k0):
    """Apply the PEC reflection (see the top of the file). Returns
    (new state, oriented normal, curvature, cos(theta), psi correction)."""
    n, e1, e2, Kc = geom
    n, Kc = np.asarray(n, float), np.asarray(Kc, float)
    if np.dot(n, st.t) > 0:                       # normal towards the incoming beam
        n, Kc = -n, -Kc
        e2 = -e2                                  # keep (e1, e2, n) right-handed
        Kc[0, 1] = Kc[1, 0] = -Kc[0, 1]
    cos_in = float(-np.dot(n, st.t))
    M = np.eye(3) - 2 * np.outer(n, n)
    t_r = _unit(M @ st.t)
    x, y = st.xhat, st.yhat
    psi_out, dpsi = None, None
    D = np.diag([1.0, -1.0])
    if st.psi is not None:
        J = np.array([[x @ e1, x @ e2], [y @ e1, y @ e2]])
        if abs(np.linalg.det(J)) < 1e-6 or k0 is None:
            dpsi = np.zeros((2, 2))
        else:
            Ji = np.linalg.inv(J)
            dpsi = -2.0 * k0 * cos_in * (Ji.T @ Kc @ Ji)
        psi_out = D @ (st.psi + dpsi) @ D
        psi_out = 0.5 * (psi_out + psi_out.T)
    new = RayState(point=st.point.copy(), t=t_r, xhat=M @ x, yhat=-(M @ y), psi=psi_out,
                   E=None if st.E is None else -(M @ st.E),
                   K=None if st.K is None else M @ st.K, distance=st.distance)
    return new, n, e1, e2, Kc, cos_in, dpsi


def beam_widths(psi):
    lam = np.linalg.eigvalsh(np.imag(psi))
    return np.sqrt(2.0 / np.clip(lam, 1e-12, None))


def trace_reflections(src, surfaces, cfg):
    """Follow the beam, reflect it on the surfaces. Returns (legs, reflections)."""
    beam, k0 = src["beam"], src["k0"]
    s_all = arclength(beam.points)
    total = s_all[-1]
    seg_len = np.diff(s_all)
    spacing = cfg.leg_point_spacing or float(np.median(seg_len[seg_len > 0]))
    j0 = src["launch_index"] if cfg.skip_backward_extension else 0
    s_start = s_all[j0] + max(cfg.ignore_hits_before, 0.0)
    legs, refl = [], []
    if not surfaces or s_start >= total:
        return [beam], refl
    j = int(np.clip(np.searchsorted(s_all, s_start, side="right") - 1, 0, beam.n - 2))
    f = (s_start - s_all[j]) / max(seg_len[j], 1e-300)
    X_start = beam.points[j] + f * (beam.points[j + 1] - beam.points[j])
    current, straight, used = beam, src["straight"], 0.0

    for r in range(cfg.max_reflections + 1):
        hit = first_surface_hit(current, surfaces, j, X_start, straight)
        if hit is None:
            legs.append(current)
            break
        srf, jh, X, data = hit
        seg = current.points[jh + 1] - current.points[jh]
        t_in = _unit(seg)
        st = state_at(current, jh, X, t_in, k0)
        leg = densify(truncate_leg(current, jh, st), 8, k0)
        legs.append(leg)
        used += arclength(leg.points)[-1]
        if r == cfg.max_reflections:
            print(f"WARNING: max_reflections = {cfg.max_reflections} reached: the beam is "
                  f"STOPPED on surface '{srf.name}' at {np.round(X, 4)} (path {used:.4g} m "
                  f"of {total:.4g} m).")
            break
        warn = []
        if src["scotty_range"] is not None and r == 0:
            lo, hi = src["scotty_range"]
            if lo <= jh < hi:
                warn.append("hit inside the part of the beam traced by Scotty (plasma?): "
                            "the reflected beam is propagated in VACUUM")
        # beam width at the hit -> fit radius for mesh curvature
        fit_r = cfg.curvature_fit_radius
        if not fit_r:
            wmax = beam_widths(st.psi).max() * cfg.envelope_factor if st.psi is not None else 0.02
            cos0 = abs(np.dot(srf.rough_normal(X, data), t_in))
            fit_r = wmax / max(cos0, 0.2)
        geom = srf.local_geometry(X, data, fit_r)
        note = geom[4]
        new, n, e1, e2, Kc, cos_in, dpsi = reflect_state(st, geom[:4], k0)
        if cos_in < 0.02:
            warn.append(f"grazing incidence (theta = {np.degrees(np.arccos(cos_in)):.2f} deg): "
                        "the Gaussian-beam reflection is inaccurate")
        if isinstance(srf, MeshSurface) and np.dot(new.t, srf.fn[data] *
                                                   np.sign(-np.dot(srf.fn[data], t_in))) < 0:
            warn.append("the reflected ray points into the mesh facet that was hit "
                        "(fitted normal differs a lot from the facet): check the mesh")
        rec = Reflection(number=r + 1, surface=srf.name, point=X.copy(), normal=n, t_in=t_in,
                         t_out=new.t, incidence_deg=float(np.degrees(np.arccos(np.clip(cos_in, -1, 1)))),
                         e1=e1, e2=e2, curvature=Kc, geometry_note=note, distance=st.distance,
                         path_length=used, psi_in=st.psi, psi_out=new.psi, E_in=st.E,
                         E_out=new.E, warnings=warn)
        rec.surface_obj = srf
        if st.psi is not None:
            J = np.array([[st.xhat @ e1, st.xhat @ e2], [st.yhat @ e1, st.yhat @ e2]])
            G = J.T @ np.imag(st.psi) * cfg.envelope_factor ** -2 @ J
            fp = np.sqrt(2.0 / np.clip(np.linalg.eigvalsh(G), 1e-12, None))
            rec.footprint = fp
            rec.sagitta = 0.5 * float(np.abs(np.linalg.eigvalsh(Kc)).max()) * fp.max() ** 2
            mu = np.linalg.eigvalsh(np.real(dpsi)) if dpsi is not None else np.zeros(2)
            with np.errstate(divide="ignore"):
                rec.focal_lengths = np.where(np.abs(mu) > 1e-12, -k0 / mu, np.inf) if k0 else None
        refl.append(rec)
        print(f"Reflection {r + 1}: surface '{srf.name}' at ({X[0]:.5g}, {X[1]:.5g}, {X[2]:.5g}) m, "
              f"incidence {rec.incidence_deg:.2f} deg, path {used:.4g} m; {note}")
        for w in warn:
            print(f"  WARNING: {w}")
        remaining = total - used
        if remaining <= MIN_LEG_LENGTH:
            break
        n_pts = int(np.clip(np.ceil(remaining / spacing) + 1, 20, 2000))
        current = vacuum_leg(new, remaining, n_pts, k0)
        straight = True
        j = 0
        X_start = new.point + min(MIN_LEG_LENGTH, 0.5 * remaining) * new.t
    return legs, refl


def side_contacts(legs, leg_AB, surfaces, reflections, d, n_theta=16, every=3):
    """Places where the beam (1/e x factor) or its envelope (+ separation)
    crosses a surface although the central ray does not hit it there (the
    part of the beam that would be reflected is NOT modelled). Returns a list
    of (leg, distance along the leg [m], surface name, 'beam' | 'envelope')."""
    found = []
    th = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    for k, (leg, (A, B)) in enumerate(zip(legs, leg_AB)):
        s = arclength(leg.points)
        near_refl = []
        for r in reflections:
            if r.number - 1 in (k - 1, k):
                fp = (r.footprint.max() if r.footprint is not None else 0.05) + d
                near_refl.append((r.point, fp / max(np.cos(np.deg2rad(r.incidence_deg)), 0.2)))
        last = {}
        for i in range(0, leg.n, every):
            C = leg.points[i]
            if any(np.linalg.norm(C - p) < rad for p, rad in near_refl):
                continue
            ring = (np.cos(th)[:, None] * A[i] + np.sin(th)[:, None] * B[i])
            rad = np.linalg.norm(ring, axis=1, keepdims=True)
            for what, R_ in (("beam", ring), ("envelope", ring * (1 + d / np.maximum(rad, 1e-12)))):
                P0 = np.repeat(C[None], n_theta, 0)
                for srf in surfaces:
                    if srf.first_hit(P0, C + R_) is not None:
                        key = (srf.name, what)
                        if key not in last or s[i] - last[key][-1] > 0.05:
                            found.append([k, s[i], s[i], srf.name, what])
                            last[key] = found[-1]
                            last[key].append(s[i])
                        else:
                            last[key][2] = last[key][-1] = s[i]
                        break
    # keep 'envelope' entries only where the beam itself is clear
    out = []
    for k, s0, s1, name, what, _ in found:
        if what == "envelope" and any(f[0] == k and f[3] == name and f[4] == "beam" and
                                      f[1] <= s1 and f[2] >= s0 for f in found):
            continue
        out.append((k, s0, s1, name, what))
    return out


def join_legs(legs):
    """One Beam with all legs (hit points not duplicated) + start index of
    every leg in it."""
    parts, starts, n = [], [], 0
    for k, leg in enumerate(legs):
        sl = slice(0, None) if k == 0 else slice(1, None)
        parts.append(_take(leg, sl))
        starts.append(max(n - 1, 0) if k else 0)
        n += parts[-1].n
    cat = lambda f: (None if any(getattr(p, f) is None for p in parts)
                     else np.concatenate([getattr(p, f) for p in parts]))
    return Beam(**{f: cat(f) for f in ("points", "xhat", "yhat", "psi_xx", "psi_xy",
                                       "psi_yy", "distance", "E", "K")}), starts


# ##########################################################################
# 5. FEM ENVELOPE OF THE FOLDED BEAM: per-leg chains, mirror cuts, union
# ##########################################################################

def _cell_faces(planes, size, eps):
    return [s2g._convex_face(planes, i, size, eps) for i in range(len(planes))]


def clip_chain_at_plane(cells, n, h, label, from_end, X, r_local, size, eps, force_cap=False):
    """Cut the cells of one leg by the half-space n.x <= h, starting from the
    leg end that touches the mirror and continuing as long as the joint
    polygon towards the rest of the chain still reaches behind the plane."""
    m = len(cells)
    behind = lambda P: len(P) > 0 and bool(np.any(P @ n > h + 10 * eps))
    # cells that reach behind the plane close to the hit point
    near = []
    for k in range(m):
        V = np.vstack([P for P, _ in cells[k]["faces"] if len(P)])
        bad = V[V @ n > h + 10 * eps]
        if len(bad) and np.min(np.linalg.norm(bad - X, axis=1)) < r_local:
            near.append(k)
    if from_end:
        lo = min(near + [m - 1])
        while lo > 0 and behind(cells[lo]["faces"][0][0]):   # entry joint still behind
            lo -= 1
        rng = list(range(lo, m))
    else:
        hi = max(near + [0])
        while hi < m - 1 and behind(cells[hi]["faces"][1][0]):
            hi += 1
        rng = list(range(0, hi + 1))
    ptol = 1e-7 * size / 20.0
    span_ = size / 20.0
    # leg-end cap almost on the mirror plane (near-normal incidence): use the
    # mirror plane itself as the cap (avoids micron-thin wedge faces)
    ke, ie = (m - 1, 1) if from_end else (0, 0)
    ce = cells[ke]
    nc, hc = ce["planes"][ie]
    Pc = ce["faces"][ie][0]
    if force_cap or (np.dot(nc, n) > 0.999 and len(Pc) >= 3
                     and np.max(np.abs(Pc @ n - h)) < 1e-4 * span_):
        # force_cap: the leg runs on to the cut plane (no notch between the legs
        # behind the mirror's tangent plane, needed when the cut is moved behind)
        ce["planes"] = list(ce["planes"])
        ce["planes"][ie] = (np.asarray(n, float), float(h))
        ce["faces"] = _cell_faces(ce["planes"], size, eps)
    for k in sorted(set(rng)):
        c = cells[k]
        span_ = size / 20.0
        same = [i for i, (ni, hi_) in enumerate(c["planes"])
                if np.dot(ni, n) > 0.99 and len(c["faces"][i][0]) >= 3 and
                np.max(np.abs(c["faces"][i][0] @ n - h)) < 1e-6 * span_]
        if same:                 # the cut coincides with a face (normal incidence)
            c.setdefault("pec", {})[same[0]] = label
            continue
        c["planes"] = list(c["planes"]) + [(np.asarray(n, float), float(h))]
        c.setdefault("pec", {})[len(c["planes"]) - 1] = label
        c["faces"] = _cell_faces(c["planes"], size, eps)
    return sorted(set(rng))


def chain_boundary(cells, leg, first_kind, last_kind):
    out = []
    last = len(cells) - 1
    for k, c in enumerate(cells):
        for i, (P, _) in enumerate(c["faces"]):
            if len(P) < 3:
                continue
            if (i == 0 and k > 0) or (i == 1 and k < last):
                continue
            kind = c.get("pec", {}).get(i) or (first_kind if i == 0 else
                                                 last_kind if i == 1 else "lateral")
            out.append(dict(P=np.asarray(P, float), kind=kind, plane=c["planes"][i], leg=leg))
    return out


class _Pool:
    """Merge points closer than tol (grid hashing)."""

    def __init__(self, tol):
        self.tol = tol
        self.pts = []
        self.grid = defaultdict(list)

    def add(self, p):
        key = tuple(np.floor(p / self.tol).astype(np.int64))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in self.grid.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if np.linalg.norm(self.pts[i] - p) <= self.tol:
                            return i
        self.pts.append(np.asarray(p, float))
        self.grid[key].append(len(self.pts) - 1)
        return len(self.pts) - 1

    def array(self):
        return np.array(self.pts)


def _inside_cell(c, F, cell_leg, face_leg, ptol, ntol=1e-7):
    """Is point c (on face F) inside the convex cell? Coplanar planes: same
    orientation -> inside only if the face belongs to a later leg (so one of
    two coincident faces survives); opposite orientation -> inside."""
    nF, hF = F["plane"]
    for n, h in cell_leg["planes"]:
        if np.linalg.norm(n - nF) < ntol and abs(h - hF) < ptol:
            if face_leg < cell_leg["leg"]:
                return False
            continue
        if np.linalg.norm(n + nF) < ntol and abs(h + hF) < ptol:
            continue
        if np.dot(n, c) > h + ptol:
            return False
    return True


def _face_minus_cells(F, cells, eps, ptol):
    """Convex fragments of face F lying outside every cell in ``cells``
    (arrangement of the lines bounding F n cell, then classification)."""
    P = F["P"]
    nF = F["plane"][0]
    lines = []
    hit_cells = []
    for c in cells:
        poly, lab = [p for p in P], [-1] * len(P)
        for idx, (n, h) in enumerate(c["planes"]):
            if np.linalg.norm(np.cross(n, nF)) < 1e-9:      # plane parallel to F
                if abs(np.dot(n, P[0]) - h) <= ptol:          # coplanar: see _inside_cell
                    continue
                if np.dot(n, P[0]) > h:                       # F outside this cell
                    poly = []
                    break
                continue
            poly, lab = s2g._clip_polygon(poly, lab, n, h, idx, eps)
            if not poly:
                break
        if not poly:
            continue
        hit_cells.append(c)
        for l in set(lab):
            if l >= 0:
                lines.append(c["planes"][l])
    if not hit_cells:
        return None
    frags = [[p for p in P]]
    for n, h in lines:
        new = []
        for g in frags:
            d = np.array([np.dot(n, p) - h for p in g])
            if d.max() <= eps or d.min() >= -eps:
                new.append(g)
                continue
            a, _ = s2g._clip_polygon(g, [0] * len(g), n, h, 0, eps)
            b, _ = s2g._clip_polygon(g, [0] * len(g), -n, -h, 0, eps)
            new.extend(x for x in (a, b) if len(x) >= 3)
        frags = new
    keep = [g for g in frags
            if not any(_inside_cell(np.mean(g, 0), F, c, F["leg"], ptol) for c in hit_cells)]
    return keep


def _loops_from_fragments(frags, pool, nF):
    """Boundary loops (vertex ids) of the union of convex fragments (CCW
    about nF); interior edges shared by two fragments cancel."""
    edges = defaultdict(int)
    for g in frags:
        ids = [pool.add(p) for p in g]
        ids = [v for i, v in enumerate(ids) if v != ids[i - 1]]
        for a, b in zip(ids, ids[1:] + ids[:1]):
            if a != b:
                edges[(a, b)] += 1
    rem = []
    for (a, b), k in edges.items():
        k2 = edges.get((b, a), 0)
        rem.extend([(a, b)] * max(k - k2, 0))
    succ = defaultdict(list)
    for a, b in rem:
        succ[a].append(b)
    loops = []
    while True:
        start = next((a for a, l in succ.items() if l), None)
        if start is None:
            break
        loop, cur = [start], start
        for _ in range(100000):
            b = succ[cur].pop(0)
            if b == start:
                break
            loop.append(b)
            cur = b
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def _drop_collinear(loop, X, tol):
    changed = True
    while changed and len(loop) > 3:
        changed = False
        for i in range(len(loop)):
            p, v, q = X[loop[i - 1]], X[loop[i]], X[loop[(i + 1) % len(loop)]]
            a, b = v - p, q - v
            la, lb = np.linalg.norm(a), np.linalg.norm(b)
            if la < tol or lb < tol or np.linalg.norm(np.cross(a, b)) <= 1e-9 * la * lb \
                    and np.dot(a, b) > 0:
                del loop[i]
                changed = True
                break
    return loop


def _signed_area(loop, X, n):
    P = X[loop]
    return 0.5 * float(np.dot(np.sum(np.cross(P, np.roll(P, -1, 0)), 0), n))


def _point_in_loop_2d(q, poly):
    x, y = q
    inside = False
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xc = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if xc > x:
                inside = not inside
    return inside


def union_boundary(leg_cells, leg_faces, span):
    """Exact boundary of the union of all legs' cells. Returns a list of
    dict(outer=(k,3), holes=[(m,3)...], kind, normal)."""
    eps = 1e-9 * span
    ptol = 1e-7 * span
    all_cells = []
    for leg, cells in enumerate(leg_cells):
        for c in cells:
            V = np.vstack([P for P, _ in c["faces"] if len(P)])
            all_cells.append(dict(planes=c["planes"], leg=leg, lo=V.min(0), hi=V.max(0)))
    pool = _Pool(1e-6 * span)
    raw = []                                          # (loops, kind, normal)
    n_cut = 0
    for F in leg_faces:
        lo, hi = F["P"].min(0) - ptol, F["P"].max(0) + ptol
        cand = [c for c in all_cells if c["leg"] != F["leg"]
                and np.all(c["hi"] >= lo) and np.all(c["lo"] <= hi)]
        frags = _face_minus_cells(F, cand, eps, ptol) if cand else None
        if frags is None:
            frags = [list(F["P"])]
        else:
            n_cut += 1
        if not frags:
            continue
        loops = _loops_from_fragments(frags, pool, F["plane"][0])
        if loops:
            raw.append((loops, F["kind"], F["plane"][0], F["leg"]))
    X = pool.array()
    # simplify, then make every edge conforming (insert T-junction vertices)
    for loops, _, _, _ in raw:
        for i, l in enumerate(loops):
            loops[i] = _drop_collinear(list(l), X, 1e-6 * span)
    used = sorted({v for loops, _, _, _ in raw for l in loops for v in l})
    U = X[used]
    tol = 1e-6 * span
    for loops, _, _, _ in raw:
        for i, l in enumerate(loops):
            out = []
            for a, b in zip(l, l[1:] + l[:1]):
                out.append(a)
                pa, pb = X[a], X[b]
                d = pb - pa
                L2 = float(d @ d)
                if L2 <= 0:
                    continue
                t = (U - pa) @ d / L2
                dist = np.linalg.norm(U - pa - t[:, None] * d, axis=1)
                sel = np.nonzero((t > 1e-9) & (t < 1 - 1e-9) & (dist < tol))[0]
                sel = [used[k] for k in sel[np.argsort(t[sel])] if used[k] not in (a, b)]
                out.extend(sel)
            loops[i] = out
    faces = []
    for loops, kind, nrm, leg in raw:
        outers = [l for l in loops if _signed_area(l, X, nrm) > 0]
        holes = [l for l in loops if _signed_area(l, X, nrm) <= 0]
        u, v = _plane_frame(nrm)
        to2 = lambda l: np.stack([X[l] @ u, X[l] @ v], 1)
        groups = {i: [] for i in range(len(outers))}
        for hl in holes:
            owner = next((i for i, o in enumerate(outers)
                          if _point_in_loop_2d(to2(hl)[0], to2(o))), None)
            if owner is not None:
                groups[owner].append(hl)
        for i, o in enumerate(outers):
            faces.append(dict(outer=X[o], holes=[X[h] for h in groups[i]], kind=kind, leg=leg,
                              normal=nrm))
    return faces, n_cut


# ---------------------------- triangulation of planar faces ----------------

def _seg_cross(p, q, a, b, eps):
    def orient(o, s, t):
        return (s[0] - o[0]) * (t[1] - o[1]) - (s[1] - o[1]) * (t[0] - o[0])
    d1, d2 = orient(a, b, p), orient(a, b, q)
    d3, d4 = orient(p, q, a), orient(p, q, b)
    return ((d1 > eps and d2 < -eps) or (d1 < -eps and d2 > eps)) and \
           ((d3 > eps and d4 < -eps) or (d3 < -eps and d4 > eps))


def triangulate_face(outer, holes, nrm):
    """Ear-clipping triangulation of a planar polygon with holes.
    Returns triangles as (T, 3, 3) points (CCW about nrm)."""
    u, v = _plane_frame(nrm)
    pts3 = [outer] + list(holes)
    allp = np.vstack(pts3)
    P = np.stack([allp @ u, allp @ v], 1)
    scale = max(np.ptp(P, axis=0).max(), 1e-12)
    eps = 1e-12 * scale * scale
    poly = list(range(len(outer)))
    off = len(outer)
    hole_idx = []
    for h in holes:
        hole_idx.append(list(range(off, off + len(h))))
        off += len(h)
    hole_idx.sort(key=lambda idx: -P[idx, 0].max())
    for k, idx in enumerate(hole_idx):
        # a hole touching the boundary at a vertex: splice it in there
        pinch = None
        for jl, hv in enumerate(idx):
            dd = np.linalg.norm(P[poly] - P[hv], axis=1)
            if dd.min() <= 1e-9 * scale:
                pinch = (jl, int(np.argmin(dd)))
                break
        if pinch is not None:
            jl, kpos = pinch
            poly = poly[:kpos + 1] + idx[jl + 1:] + idx[:jl] + [idx[jl]] + poly[kpos + 1:]
            continue
        edges = [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]
        for other in hole_idx[k + 1:]:
            edges += [(other[i], other[(i + 1) % len(other)]) for i in range(len(other))]
        edges += [(idx[i], idx[(i + 1) % len(idx)]) for i in range(len(idx))]
        nb = {}
        for lst in [poly, idx]:
            for i, w in enumerate(lst):
                nb.setdefault(w, []).extend([lst[i - 1], lst[(i + 1) % len(lst)]])

        def collinear(o, w, z):
            d1, d2 = P[w] - P[o], P[z] - P[o]
            return abs(d1[0] * d2[1] - d1[1] * d2[0]) <= 1e-9 * np.linalg.norm(d1) * np.linalg.norm(d2)

        pairs = sorted(((np.linalg.norm(P[w] - P[h]), jl, w) for jl, h in enumerate(idx)
                        for w in set(poly)), key=lambda x: x[0])
        best = None
        for strict in (True, False):
            for _, jl, w in pairs:
                h = idx[jl]
                if any(_seg_cross(P[h], P[w], P[a], P[b], eps)
                       for a, b in edges if h not in (a, b) and w not in (a, b)):
                    continue
                if strict and (any(collinear(h, w, z) for z in nb[h]) or
                               any(collinear(w, h, z) for z in nb[w]) or
                               any(np.linalg.norm(P[z] - P[h]) > 1e-12 and
                                   np.linalg.norm(P[z] - P[w]) > 1e-12 and collinear(h, w, z) and
                                   0 < np.dot(P[z] - P[h], P[w] - P[h]) < np.dot(P[w] - P[h], P[w] - P[h])
                                   for z in poly + idx)):
                    continue
                best, hm_l, hm = w, jl, h
                break
            if best is not None:
                break
        if best is None:
            hm_l = int(np.argmax(P[idx, 0]))
            hm = idx[hm_l]
            best = min(set(poly), key=lambda w: np.linalg.norm(P[w] - P[hm]))
        kpos = poly.index(best)
        seq = idx[hm_l:] + idx[:hm_l] + [hm]
        poly = poly[:kpos + 1] + seq + [best] + poly[kpos + 1:]
    tris = []
    V = poly[:]
    cross = lambda a, b, c: ((P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) -
                             (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0]))

    def in_tri(w, a, b, c):
        for z in (a, b, c):
            if np.linalg.norm(P[w] - P[z]) <= 1e-9 * scale:
                return False
        return cross(a, b, w) >= -eps and cross(b, c, w) >= -eps and cross(c, a, w) >= -eps

    guard = 0
    while len(V) > 3 and guard < 20 * len(poly) + 100:
        guard += 1
        found = False
        for i in range(len(V)):
            a, b, c = V[i - 1], V[i], V[(i + 1) % len(V)]
            if cross(a, b, c) <= eps:
                continue
            if any(in_tri(w, a, b, c) for w in V if w not in (a, b, c)):
                continue
            tris.append((a, b, c))
            del V[i]
            found = True
            break
        if not found:                    # allow other vertices ON the ear's edges
            for i in range(len(V)):
                a, b, c = V[i - 1], V[i], V[(i + 1) % len(V)]
                if cross(a, b, c) <= eps:
                    continue
                if any(cross(a, b, w) > eps and cross(b, c, w) > eps and cross(c, a, w) > eps
                       for w in V if w not in (a, b, c)):
                    continue
                tris.append((a, b, c))
                del V[i]
                found = True
                break
        if not found:                    # degenerate remainder: clip the flattest corner
            i = int(np.argmax([cross(V[i - 1], V[i], V[(i + 1) % len(V)]) for i in range(len(V))]))
            tris.append((V[i - 1], V[i], V[(i + 1) % len(V)]))
            del V[i]
    if len(V) == 3:
        tris.append(tuple(V))
    return allp[np.array(tris)] if tris else np.zeros((0, 3, 3))


def faces_to_triangles(faces):
    tris = [triangulate_face(f["outer"], f["holes"], f["normal"]) for f in faces]
    tris = [t for t in tris if len(t)]
    return np.concatenate(tris) if tris else np.zeros((0, 3, 3))


def mesh_check(tris, span):
    """(volume, number of edges not shared by exactly 2 opposite triangles)."""
    V, T = merge_vertices(tris)
    vol = float(np.sum(np.einsum("ij,ij->i", V[T[:, 0]], np.cross(V[T[:, 1]], V[T[:, 2]])))) / 6
    cnt = defaultdict(int)
    for a, b, c in T:
        for e in ((a, b), (b, c), (c, a)):
            if e[0] != e[1]:
                cnt[e] += 1
    bad = sum(1 for (a, b), k in cnt.items() if k != cnt.get((b, a), 0) or k != 1)
    return vol, bad


def repair_small_gaps(faces, pec_planes, span, max_size_frac=1e-3):
    """Close tiny holes left by floating-point slivers of the union: chain
    the unmatched (boundary) edges into loops and fill each loop that is
    smaller than max_size_frac x span with triangles. Returns (#filled, #big)."""
    tris = [triangulate_face(f["outer"], f["holes"], f["normal"]) for f in faces]
    tris = np.concatenate([t for t in tris if len(t)])
    V, T = merge_vertices(tris)
    cnt = defaultdict(int)
    for a, b, c in T:
        for e in ((a, b), (b, c), (c, a)):
            if e[0] != e[1]:
                cnt[e] += 1
    succ = defaultdict(list)          # reversed boundary edges, oriented for filling
    for (a, b), k in cnt.items():
        for _ in range(k - cnt.get((b, a), 0)):
            succ[b].append(a)
    filled = big = 0
    tol = 1e-6 * span
    while True:
        start = next((a for a, l in succ.items() if l), None)
        if start is None:
            break
        loop, cur = [start], start
        while succ[cur]:
            nxt = succ[cur].pop()
            if nxt == start:
                break
            loop.append(nxt)
            cur = nxt
        P = V[loop]
        if len(P) < 3:
            continue
        if np.ptp(P, axis=0).max() > max_size_frac * span:
            big += 1
            continue
        kind = "lateral"
        for label, n, h in pec_planes:
            if np.all(np.abs(P @ n - h) < 10 * tol):
                kind = label
        for t in range(1, len(P) - 1):
            tri = np.array([P[0], P[t], P[t + 1]])
            nrm = np.cross(tri[1] - tri[0], tri[2] - tri[0])
            if np.linalg.norm(nrm) <= 1e-300:
                continue
            faces.append(dict(outer=tri, holes=[], kind=kind, leg=-1, normal=_unit(nrm)))
        filled += 1
    return filled, big


# ---------------------------- ENVPEC faces on the real mirror surface -------

class _ConformError(Exception):
    pass


def _surface_hit_near(srf, p, d, back, front):
    """Crossing of the line p + t d (|d| = 1) with the surface, t in [-back,
    front], closest to p. Returns t or None."""
    best = None
    for sgn, L in ((1.0, front), (-1.0, back)):
        if L <= 0:
            continue
        h = srf.first_hit(p[None], (p + sgn * L * d)[None])
        if h is not None:
            t = sgn * float(np.linalg.norm(h[2] - p))
            if best is None or abs(t) < abs(best):
                best = t
    return best


def _simple_polygon_2d(P):
    n = len(P)
    if n < 3:
        return False
    scale = max(np.ptp(P, axis=0).max(), 1e-12)
    eps = 1e-12 * scale * scale
    seg = [(P[i], P[(i + 1) % n]) for i in range(n)]
    lo = np.minimum(P, np.roll(P, -1, 0))
    hi = np.maximum(P, np.roll(P, -1, 0))
    for i in range(n):
        cand = np.nonzero(np.all(hi >= lo[i], 1) & np.all(lo <= hi[i], 1))[0]
        for j in cand:
            if j <= i + 1 or (i == 0 and j == n - 1):
                continue
            if _seg_cross(seg[i][0], seg[i][1], seg[j][0], seg[j][1], eps):
                return False
    return True


def _incircle(a, b, c, d):
    m = np.array([[a[0] - d[0], a[1] - d[1], (a[0] - d[0]) ** 2 + (a[1] - d[1]) ** 2],
                  [b[0] - d[0], b[1] - d[1], (b[0] - d[0]) ** 2 + (b[1] - d[1]) ** 2],
                  [c[0] - d[0], c[1] - d[1], (c[0] - d[0]) ** 2 + (c[1] - d[1]) ** 2]])
    return np.linalg.det(m)


def _orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _cdt_by_flips(loops2, G2):
    """Constrained triangulation of a polygon with holes (loops2) plus
    interior points G2, robust to degenerate (hair-like) boundaries: ear
    clipping, then Lawson flips and point insertion; boundary edges are never
    flipped. Returns (P2, triangles) with P2 = boundary points + kept G2."""
    loops2 = [np.asarray(l, float) for l in loops2]
    area = lambda l: 0.5 * np.sum(l[:, 0] * np.roll(l[:, 1], -1) - np.roll(l[:, 0], -1) * l[:, 1])
    k0 = int(np.argmax([abs(area(l)) for l in loops2]))
    outer = loops2[k0] if area(loops2[k0]) > 0 else loops2[k0][::-1]
    holes = [l if area(l) < 0 else l[::-1] for k, l in enumerate(loops2) if k != k0]
    B = np.vstack(loops2)
    lookup = {tuple(np.round(p, 14)): i for i, p in enumerate(B)}
    T3 = triangulate_face(np.column_stack([outer, np.zeros(len(outer))]),
                          [np.column_stack([h, np.zeros(len(h))]) for h in holes],
                          np.array([0.0, 0.0, 1.0]))
    tris = []
    for t in T3:
        ids = [lookup.get(tuple(np.round(p[:2], 14))) for p in t]
        if None in ids:
            raise _ConformError("could not triangulate the curved ENVPEC face")
        tris.append(ids)
    P = [p for p in B]
    cons = set()
    off = 0
    for l in loops2:
        m = len(l)
        for i in range(m):
            cons.add(frozenset((off + i, off + (i + 1) % m)))
        off += m
    tris = [t if _orient(P[t[0]], P[t[1]], P[t[2]]) >= 0 else [t[0], t[2], t[1]] for t in tris]

    def edge_map():
        em = {}
        for ti, t in enumerate(tris):
            if t is None:
                continue
            for k in range(3):
                em[(t[k], t[(k + 1) % 3])] = ti
        return em

    def legalize(stack, em):
        guard = 0
        while stack and guard < 200000:
            guard += 1
            a, b = stack.pop()
            if frozenset((a, b)) in cons:
                continue
            t1, t2 = em.get((a, b)), em.get((b, a))
            if t1 is None or t2 is None or tris[t1] is None or tris[t2] is None:
                continue
            c = next(v for v in tris[t1] if v not in (a, b))
            d = next(v for v in tris[t2] if v not in (a, b))
            # t1 = (a, b, c) ccw, t2 = (b, a, d) ccw; flip if d inside circle(a, b, c)
            if _incircle(P[a], P[b], P[c], P[d]) <= 1e-30:
                continue
            if _orient(P[c], P[d], P[b]) <= 0 or _orient(P[d], P[c], P[a]) <= 0:
                continue                             # flip would fold
            for e in ((a, b), (b, c), (c, a), (b, a), (a, d), (d, b)):
                em.pop(e, None)
            tris[t1], tris[t2] = [c, d, b], [d, c, a]
            for ti in (t1, t2):
                t = tris[ti]
                for k in range(3):
                    em[(t[k], t[(k + 1) % 3])] = ti
            stack += [(b, c), (c, a), (a, d), (d, b)]

    em = edge_map()
    legalize([e for e in list(em)], em)
    Ta = None
    for g in np.asarray(G2, float).reshape(-1, 2):
        live = [ti for ti, t in enumerate(tris) if t is not None]
        A = np.array([P[tris[ti][0]] for ti in live])
        Bv = np.array([P[tris[ti][1]] for ti in live])
        C = np.array([P[tris[ti][2]] for ti in live])
        w0 = (Bv[:, 0] - A[:, 0]) * (g[1] - A[:, 1]) - (Bv[:, 1] - A[:, 1]) * (g[0] - A[:, 0])
        w1 = (C[:, 0] - Bv[:, 0]) * (g[1] - Bv[:, 1]) - (C[:, 1] - Bv[:, 1]) * (g[0] - Bv[:, 0])
        w2 = (A[:, 0] - C[:, 0]) * (g[1] - C[:, 1]) - (A[:, 1] - C[:, 1]) * (g[0] - C[:, 0])
        tot = w0 + w1 + w2
        inside = np.nonzero((w0 > 1e-6 * tot) & (w1 > 1e-6 * tot) & (w2 > 1e-6 * tot))[0]
        if not len(inside):
            continue                                 # on an edge / outside: skip it
        ti = live[inside[0]]
        a, b, c = tris[ti]
        P.append(g)
        n = len(P) - 1
        for e in ((a, b), (b, c), (c, a)):
            em.pop(e, None)
        tris[ti] = [a, b, n]
        tris.append([b, c, n])
        tris.append([c, a, n])
        for tj in (ti, len(tris) - 2, len(tris) - 1):
            t = tris[tj]
            for k in range(3):
                em[(t[k], t[(k + 1) % 3])] = tj
        legalize([(a, b), (b, c), (c, a)], em)
    tri = np.array([t for t in tris if t is not None], int)
    return np.array(P), tri


def _conforming_delaunay(loops2, loops3, G2, G3, max_iter=30):
    """Delaunay triangulation of the polygon (loops2, with 3D twins loops3)
    plus interior points G2 / G3 that contains every boundary segment:
    encroached segments are split at their midpoint (Ruppert), encroaching
    interior points dropped. Returns (loops2, loops3, P2, P3, triangles,
    splits) -- splits = [(a3, b3, mid3)] to apply to the neighbouring faces."""
    from scipy.spatial import Delaunay, cKDTree
    L2 = [list(map(np.asarray, l)) for l in loops2]
    L3 = [list(map(np.asarray, l)) for l in loops3]
    G2, G3 = np.asarray(G2, float).reshape(-1, 2), np.asarray(G3, float).reshape(-1, 3)
    scale = max(np.ptp(np.vstack(loops2), axis=0).max(), 1e-12)
    splits = []
    for _ in range(max_iter):
        B2 = np.vstack([np.array(l) for l in L2])
        nb = len(B2)
        tree = cKDTree(np.vstack([B2, G2]) if len(G2) else B2)
        enc, drop, off = [], set(), 0
        for li, l in enumerate(L2):
            m = len(l)
            for i in range(m):
                a, b = l[i], l[(i + 1) % m]
                r = 0.5 * np.linalg.norm(b - a)
                for j in tree.query_ball_point(0.5 * (a + b), r * (1 - 1e-9)):
                    if j >= nb:
                        drop.add(j - nb)
                    elif j not in (off + i, off + (i + 1) % m) and r > 1e-9 * scale:
                        enc.append((li, i))
                        break
            off += m
        if drop:
            keep = np.array([k not in drop for k in range(len(G2))])
            G2, G3 = G2[keep], G3[keep]
        if not enc:
            break
        for li, i in sorted(set(enc), reverse=True):
            l2, l3 = L2[li], L3[li]
            j = (i + 1) % len(l2)
            mid3 = 0.5 * (l3[i] + l3[j])
            splits.append((l3[i], l3[j], mid3))
            l2.insert(i + 1, 0.5 * (l2[i] + l2[j]))
            l3.insert(i + 1, mid3)
    loops2 = [np.array(l) for l in L2]
    loops3 = [np.array(l) for l in L3]
    B2, B3 = np.vstack(loops2), np.vstack(loops3)
    P2 = np.vstack([B2, G2]) if len(G2) else B2
    P3 = np.vstack([B3, G3]) if len(G3) else B3
    tri = Delaunay(P2).simplices
    tri = tri[_readers()._points_in_loops(P2[tri].mean(1), loops2)]
    edges = Counter(frozenset((t[k], t[(k + 1) % 3])) for t in tri for k in range(3))
    off, ok = 0, True
    for l in loops2:
        m = len(l)
        for i in range(m):
            if edges.get(frozenset((off + i, off + (i + 1) % m)), 0) != 1:
                ok = False
        off += m
    if not ok:
        # degenerate boundary (e.g. hair-thin parts): constrained flips instead
        P2c, tri = _cdt_by_flips(loops2, G2)
        nb = len(B2)
        keep = {tuple(np.round(g, 14)): k for k, g in enumerate(G2)}
        G3c = np.array([G3[keep[tuple(np.round(g, 14))]] for g in P2c[nb:]]).reshape(-1, 3)
        P2 = P2c
        P3 = np.vstack([B3, G3c]) if len(G3c) else B3
    return loops2, loops3, P2, P3, tri, splits


def _apply_splits(faces, splits, tol):
    """Insert the split points into every face loop edge they belong to."""
    for a, b, m in splits:
        for f in faces:
            if "curved" in f:
                continue
            ls = [f["outer"]] + list(f["holes"])
            changed = False
            for k, L in enumerate(ls):
                da = np.linalg.norm(L - a, axis=1)
                db = np.linalg.norm(L - b, axis=1)
                ia, ib = np.nonzero(da < tol)[0], np.nonzero(db < tol)[0]
                for i in ia:
                    n = len(L)
                    j = (i + 1) % n
                    if j in ib or (i - 1) % n in ib:
                        pos = i + 1 if j in ib else i
                        ls[k] = np.insert(L, pos, m, axis=0)
                        changed = True
                        break
            if changed:
                f["outer"], f["holes"] = ls[0], ls[1:]


def _curved_patch(srf, o, e1, e2, nf, loops3, loops2, h, lmax, outward):
    """Triangles of the surface over the region bounded by loops2 (2D coords
    in the frame o, e1, e2) + a height grid for the IGES B-spline patch."""
    try:
        from scipy.spatial import Delaunay, cKDTree
    except ImportError:
        raise _ConformError("needs scipy (pip install scipy)")
    gsr = _readers()
    B2 = np.vstack(loops2)
    B3 = np.vstack(loops3)
    lo, hi = B2.min(0), B2.max(0)

    def heights(ab):
        H = np.full(len(ab), np.nan)
        for k, (a, b) in enumerate(ab):
            p = o + a * e1 + b * e2
            t = _surface_hit_near(srf, p, nf, lmax, lmax)
            if t is not None:
                H[k] = t
        return H

    gx = np.arange(lo[0] + 0.5 * h, hi[0], h)
    gy = np.arange(lo[1] + 0.5 * h, hi[1], h)
    G = np.stack(np.meshgrid(gx, gy, indexing="ij"), -1).reshape(-1, 2)
    if len(G):
        G = G[gsr._points_in_loops(G, loops2)]
    if len(G):
        dist, _ = cKDTree(B2).query(G)
        G = G[dist > 0.5 * h]
    Hin = heights(G) if len(G) else np.zeros(0)
    if np.isnan(Hin).any():
        raise _ConformError("the surface has holes / ends under the ENVPEC face")
    G3 = (o + G[:, :1] * e1 + G[:, 1:] * e2 + Hin[:, None] * nf) if len(G) else np.zeros((0, 3))
    loops2, loops3, P2, P3, tri, splits = _conforming_delaunay(loops2, loops3, G, G3)
    B2 = np.vstack(loops2)
    B3 = np.vstack(loops3)
    T3 = P3[tri]
    nrm = np.cross(T3[:, 1] - T3[:, 0], T3[:, 2] - T3[:, 0])
    flip = nrm @ outward < 0
    T3[flip] = T3[flip][:, ::-1]
    # height grid over the bounding box (for the IGES patch)
    pad = 0.02 * max(np.ptp(B2, axis=0).max(), 1e-9)
    n = [int(np.clip(math.ceil((hi[i] - lo[i] + 2 * pad) / h), 4, 80)) + 1 for i in range(2)]
    ga = np.linspace(lo[0] - pad, hi[0] + pad, n[0])
    gb = np.linspace(lo[1] - pad, hi[1] + pad, n[1])
    GG = np.stack(np.meshgrid(ga, gb, indexing="ij"), -1).reshape(-1, 2)
    HG = heights(GG)
    ok = ~np.isnan(HG)
    known2 = np.vstack([GG[ok], B2])
    knownh = np.concatenate([HG[ok], (B3 - o) @ nf])
    if (~ok).any():                                 # outside the mirror: nearest height
        _, j = cKDTree(known2).query(GG[~ok])
        HG[~ok] = knownh[j]
    curved = dict(origin=o, e1=e1, e2=e2, n=nf, ga=ga, gb=gb, H=HG.reshape(n[0], n[1]),
                  loops2=loops2)
    return T3, curved, loops3, splits


def _local_surface(srf, o, frame, lo, hi):
    """Mesh surfaces: a copy with only the triangles in the box lo..hi of the
    local frame (rows of ``frame``) -- many ray casts are made against it."""
    if not isinstance(srf, MeshSurface):
        return srf
    pad = 2 * srf.edge
    L = (srf.V - o) @ np.asarray(frame).T
    inside = np.all((L >= np.asarray(lo) - pad) & (L <= np.asarray(hi) + pad), axis=1)
    T = srf.T[np.any(inside[srf.T], axis=1)]
    loc = MeshSurface(srf.name, srf.V, T, bounds=srf.bounds, source=srf.source)
    return loc


def _pec_faces(faces, rec):
    return [i for i, f in enumerate(faces) if f["kind"] == f"pec{rec.number}"]


def _pec_heights(faces, rec, srf, span, off):
    """Heights of the mirror (along its normal, from the hit point) over the
    ENVPEC face(s) of ``rec``: (gmin, gmax, size). Raises _ConformError."""
    idx = _pec_faces(faces, rec)
    if not idx:
        raise _ConformError("no ENVPEC face")
    nf = _unit(rec.normal)
    e1, e2 = _plane_frame(nf)
    o = rec.point
    loops = [L for i in idx for L in [faces[i]["outer"]] + list(faces[i]["holes"])]
    P = np.vstack(loops)
    ab = np.stack([(P - o) @ e1, (P - o) @ e2], 1)
    size = float(np.ptp(ab, axis=0).max())
    lo, hi = ab.min(0), ab.max(0)
    st = size / 20.0
    G = np.stack(np.meshgrid(np.arange(lo[0], hi[0] + st, st), np.arange(lo[1], hi[1] + st, st),
                             indexing="ij"), -1).reshape(-1, 2)
    loops2 = [np.stack([(L - o) @ e1, (L - o) @ e2], 1) for L in loops]
    G = G[_readers()._points_in_loops(G, loops2)] if len(G) else G
    W = 2 * size + 2 * abs(off) + 1e-6 * span
    srf = _local_surface(srf, o, [e1, e2, nf], np.r_[lo - size, -W], np.r_[hi + size, W])
    g = []
    for a, b in np.vstack([ab, G]):
        t = _surface_hit_near(srf, o + a * e1 + b * e2, nf, W, W)
        if t is None:
            raise _ConformError(
                f"the surface '{srf.name}' does not cover the whole ENVPEC face ({size:.3g} m "
                "wide): the envelope is wider than the part of the surface facing the beam "
                "-- reduce envelope_separation")
        g.append(t)
    return float(min(g)), float(max(g)), size


def _chain_loops(edges):
    succ = defaultdict(list)
    for a, b in edges:
        succ[a].append(b)
    loops = []
    while True:
        st = next((a for a, l in succ.items() if l), None)
        if st is None:
            return loops
        L, cur = [st], st
        for _ in range(1000000):
            if not succ[cur]:
                raise _ConformError("an edge of the curved ENVPEC face is not closed")
            b = succ[cur].pop(0)
            if b == st:
                break
            L.append(b)
            cur = b
        if len(L) >= 3:
            loops.append(L)


_TRIM_DELTA = 2e-5          # cut this far (x span) behind the surface
_TRIM_SNAP = 1e-4            # vertices with phi below this (x span) are put on the cut
_TRIM_MERGE = 0.0           # cut points closer than this (x span) to a vertex use it


def _trim_by_surface(faces, rec, srf, cfg, span, off, hinfo):
    """Remove the part of the envelope lying behind the mirror surface: every
    face near the mirror is clipped (in its own plane) by the surface, and the
    ENVPEC face is replaced by the piece of the surface those cuts enclose."""
    gsr = _readers()
    try:
        from scipy.spatial import Delaunay, cKDTree
    except ImportError:
        raise _ConformError("needs scipy (pip install scipy)")
    gmin, gmax, size = hinfo
    nf = _unit(rec.normal)
    e1, e2 = _plane_frame(nf)
    o = rec.point
    label = f"pec{rec.number}"
    pec_idx = _pec_faces(faces, rec)
    per = sum(np.linalg.norm(np.roll(faces[i]["outer"], -1, 0) - faces[i]["outer"], axis=1).sum()
              for i in pec_idx)
    hs = cfg.envelope_pec_mesh_size if cfg.envelope_pec_mesh_size > 0 else per / 150.0
    hs = max(hs, 1e-6 * span)
    mar = max(0.05 * size, 2 * hs)
    Pp = np.vstack([faces[i]["outer"] for i in pec_idx])
    abp = np.stack([(Pp - o) @ e1, (Pp - o) @ e2], 1)
    box_lo = np.r_[abp.min(0) - 0.5 * size, -off - mar]
    box_hi = np.r_[abp.max(0) + 0.5 * size, gmax + mar]
    W = (box_hi[2] - box_lo[2]) + size
    srf = _local_surface(srf, o, [e1, e2, nf], box_lo - np.r_[0, 0, W], box_hi + np.r_[0, 0, W])
    BIG = 1e30
    pool = _Pool(1e-8 * span)
    delta = _TRIM_DELTA * span
    snap = _TRIM_SNAP * span
    merge_d = _TRIM_MERGE * span
    X = pool.pts

    def local(x):
        d = np.asarray(x) - o
        return np.array([d @ e1, d @ e2, d @ nf])

    def phi_point(x):
        a, b, z = local(x)
        if np.any(np.r_[a, b, z] < box_lo) or np.any(np.r_[a, b, z] > box_hi):
            return BIG
        t = _surface_hit_near(srf, o + a * e1 + b * e2, nf, W, W)
        # (+ delta: the cut lies a hair behind the surface, so that no vertex of
        # the envelope is ever exactly on it)
        return BIG if t is None else z - t + delta

    phic = {}

    def phi(i):
        if i not in phic:
            val = phi_point(X[i])
            # a vertex (almost) on the surface is taken as ON the cut: the cut then
            # passes through it instead of producing micro-edges next to it
            phic[i] = -1e-300 if 0 < val < snap else val
        return phic[i]

    floops = [[[pool.add(p) for p in L] for L in [f["outer"]] + list(f["holes"])] for f in faces]
    orig = set(range(len(X)))
    subc, crossc, crossing = {}, {}, set()

    def edge_points(a, b):
        key = (min(a, b), max(a, b))
        if key not in subc:
            pa, pb = X[key[0]], X[key[1]]
            n = max(1, int(math.ceil(np.linalg.norm(pb - pa) / hs)))
            subc[key] = [key[0]] + [pool.add(pa + (pb - pa) * k / n) for k in range(1, n)] + [key[1]]
        s = subc[key]
        return s if key[0] == a else s[::-1]

    def cross(i, j):                                # phi(i) > 0 >= phi(j)
        key = frozenset((i, j))
        if key not in crossc:
            pi_, pj = X[i], X[j]
            lo_, hi_ = 0.0, 1.0
            fi, fj = phi(i), phi(j)
            for _ in range(40):
                m = 0.5 * (lo_ + hi_)
                if hi_ - lo_ < 1e-9:
                    break
                fm = phi_point(pi_ + (pj - pi_) * m)
                if fm > 0:
                    lo_ = m
                else:
                    hi_ = m
            p = pi_ + (pj - pi_) * hi_
            # a cut point next to an existing vertex IS that vertex (no micro-edges)
            Lij = np.linalg.norm(pj - pi_)
            if hi_ * Lij < merge_d:
                k = i
            elif (1 - hi_) * Lij < merge_d:
                k = j
            else:
                k = pool.add(p)
            crossc[key] = k
            crossing.add(k)
            phic.setdefault(k, 0.0)
        return crossc[key]

    def in_box(P3):
        L = (P3 - o) @ np.stack([e1, e2, nf], 1)
        return np.all(L.max(0) >= box_lo) and np.all(L.min(0) <= box_hi)

    new_faces, iso = [], []
    for fi, f in enumerate(faces):
        allp = np.vstack([f["outer"]] + list(f["holes"]))
        if not in_box(allp):
            new_faces.append(f)
            continue
        nrm = _unit(f["normal"])
        u, v = _plane_frame(nrm)
        uv_sign = 1.0 if np.cross(u, v) @ nrm > 0 else -1.0
        sub_loops = []
        for L in floops[fi]:
            s = []
            for a, b in zip(L, L[1:] + L[:1]):
                s += edge_points(a, b)[:-1]
            sub_loops.append(s)
        bids = [i for s in sub_loops for i in s]
        to2 = lambda ids: np.array([[X[i] @ u, X[i] @ v] for i in ids])
        B2 = to2(bids)
        loops2 = [to2(s) for s in sub_loops]
        lo2, hi2 = B2.min(0), B2.max(0)
        gx = np.arange(lo2[0] + 0.5 * hs, hi2[0], hs)
        gy = np.arange(lo2[1] + 0.5 * hs, hi2[1], hs)
        G = np.stack(np.meshgrid(gx, gy, indexing="ij"), -1).reshape(-1, 2)
        if len(G):
            G = G[gsr._points_in_loops(G, loops2)]
        if len(G):
            dd, _ = cKDTree(B2).query(G)
            G = G[dd > 0.5 * hs]
        # keep only interior points inside the slab (elsewhere the face is not cut)
        w0 = (X[bids[0]] @ nrm)
        Gp = [G[k, 0] * u + G[k, 1] * v + w0 * nrm for k in range(len(G))]
        Gp = [p for p in Gp if phi_point(p) < BIG]
        iids = [pool.add(p) for p in Gp]
        ph = np.array([phi(i) for i in bids + iids])
        if ph.min() > 0:
            new_faces.append(f)
            continue
        if fi in pec_idx and ph.max() <= 0:
            continue                                # removed entirely
        if fi in pec_idx:
            raise _ConformError(f"part of the ENVPEC face stays in front of '{srf.name}'")
        ids = bids + iids
        P2 = to2(ids)
        tri = None
        try:
            tri = Delaunay(P2).simplices
            tri = tri[gsr._points_in_loops(P2[tri].mean(1), loops2)]
            ec = Counter(frozenset((t[k], t[(k + 1) % 3])) for t in tri for k in range(3))
            off_ = 0
            for s in sub_loops:
                m = len(s)
                if any(ec.get(frozenset((off_ + k, off_ + (k + 1) % m)), 0) != 1 for k in range(m)):
                    tri = None
                    break
                off_ += m
        except Exception:
            tri = None
        if tri is None:                             # boundary only (ear clipping)
            T3 = triangulate_face(np.array([X[i] for i in sub_loops[0]]),
                                  [np.array([X[i] for i in s]) for s in sub_loops[1:]], nrm)
            tris = [[pool.add(p) for p in t] for t in T3]   # already CCW about nrm
        else:
            tris = []
            for t in tri:
                a, b, c = (ids[k] for k in t)
                p0, p1, p2 = P2[t[0]], P2[t[1]], P2[t[2]]
                ccw = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])
                tris.append([a, b, c] if ccw * uv_sign > 0 else [a, c, b])
        kept, cnt = [], Counter()
        for t in tris:
            poly, isx = [], []
            for k in range(3):
                p, q = t[k], t[(k + 1) % 3]
                if phi(p) > 0:
                    poly.append(p)
                    isx.append(False)
                if (phi(p) > 0) != (phi(q) > 0):
                    poly.append(cross(p, q) if phi(p) > 0 else cross(q, p))
                    isx.append(True)
            cps = [k for k in range(len(poly)) if isx[k]]
            if len(cps) == 2:                       # the cut through this triangle
                k0, k1 = cps
                if (k0 - k1) % len(poly) == 1:
                    k0, k1 = k1, k0
                if poly[k0] != poly[k1]:
                    iso.append((poly[k0], poly[k1]))
            poly = [x for k, x in enumerate(poly) if x != poly[k - 1]]
            if len(poly) < 3:
                continue
            for k in range(1, len(poly) - 1):
                tr = (poly[0], poly[k], poly[k + 1])
                kept.append(tr)
                for e in ((tr[0], tr[1]), (tr[1], tr[2]), (tr[2], tr[0])):
                    cnt[e] += 1
        bnd = [e for e, k in cnt.items() if k > cnt.get((e[1], e[0]), 0)]
        if not bnd:
            continue
        loops = _chain_loops(bnd)
        Xa = np.array(X)
        keep_ids = orig | crossing
        loops = [[i for i in L if i in keep_ids] for L in loops]
        loops = [L for L in loops if len(L) >= 3]
        outers = [L for L in loops if _signed_area(L, Xa, nrm) > 0]
        holes = [L for L in loops if _signed_area(L, Xa, nrm) <= 0]
        to2b = lambda L: np.stack([Xa[L] @ u, Xa[L] @ v], 1)
        for oi, L in enumerate(outers):
            hl = [h for h in holes if len(outers) == 1 or _point_in_loop_2d(to2b(h)[0], to2b(L))]
            g = dict(f)
            g["outer"], g["holes"] = Xa[L], [Xa[h] for h in hl]
            new_faces.append(g)
    # the new ENVPEC face: the piece of the surface enclosed by the cuts
    Xa = np.array(X)
    ploops = _chain_loops([(b, a) for a, b in iso])
    if not ploops:
        raise _ConformError("no intersection of the envelope with the surface found")
    outward = faces[pec_idx[0]]["normal"]
    loops3 = [Xa[L] for L in ploops]
    outers = [k for k, L in enumerate(ploops) if _signed_area(L, Xa, outward) > 0]
    holes = [k for k in range(len(ploops)) if k not in outers]
    for k in outers:
        hl = holes if len(outers) == 1 else [
            h for h in holes if _point_in_loop_2d(
                np.array([(loops3[h][0] - o) @ e1, (loops3[h][0] - o) @ e2]),
                np.stack([(loops3[k] - o) @ e1, (loops3[k] - o) @ e2], 1))]
        L3 = [loops3[k]] + [loops3[h] for h in hl]
        L2 = [np.stack([(L - o) @ e1, (L - o) @ e2], 1) for L in L3]
        for l2 in L2:
            if not _simple_polygon_2d(l2):
                raise _ConformError("the curved ENVPEC face folds over itself")
        g = dict(faces[pec_idx[0]])
        g["tris"], g["curved"], L3, splits = _curved_patch(srf, o, e1, e2, nf, L3, L2, hs, W,
                                                           outward)
        g["outer"], g["holes"] = L3[0], L3[1:]
        g["surface"] = srf.name
        new_faces = [dict(f) for f in new_faces]
        _apply_splits(new_faces, splits, 1e-9 * span)
        new_faces.append(g)
    # must still be one closed volume
    parts = [f["tris"] if "tris" in f else triangulate_face(f["outer"], f["holes"], f["normal"])
             for f in new_faces]
    vol, bad = mesh_check(np.concatenate([p for p in parts if len(p)]), span)
    if bad:
        raise _ConformError(f"the trimmed envelope is not closed ({bad} open edges)")
    return new_faces


def _envelope_union(leg_cells0, leg_AB, reflections, offsets, span, d, extend=()):
    leg_cells = copy.deepcopy(leg_cells0)
    d_size, eps = 20.0 * span, 1e-9 * span
    for r, rec in enumerate(reflections):
        n = -rec.normal                                  # points behind the mirror
        h = float(np.dot(n, rec.point)) + offsets[r]
        wmax = max(np.linalg.norm(leg_AB[r][0][-1]), np.linalg.norm(leg_AB[r][1][-1]))
        r_local = 3.0 * (wmax + d) / max(np.cos(np.deg2rad(rec.incidence_deg)), 0.1)
        label = f"pec{rec.number}"
        fc = r in extend
        clip_chain_at_plane(leg_cells[r], n, h, label, True, rec.point, r_local, d_size, eps, fc)
        if r + 1 < len(leg_cells):
            clip_chain_at_plane(leg_cells[r + 1], n, h, label, False, rec.point, r_local,
                                d_size, eps, fc)
    leg_faces = []
    for k, cells in enumerate(leg_cells):
        leg_faces += chain_boundary(cells, k, "first" if k == 0 else "cap",
                                    "last" if k == len(leg_cells) - 1 else "cap")
    faces, n_cut = union_boundary(leg_cells, leg_faces, span)
    pec_planes = [(f"pec{rec.number}", rec.normal, float(np.dot(rec.normal, rec.point))
                   - offsets[r]) for r, rec in enumerate(reflections)]
    filled, big = repair_small_gaps(faces, pec_planes, span)
    if filled:
        print(f"Envelope: closed {filled} tiny sliver gap(s) of the union (< {1e-3 * span:.2g} m)")
    if big:
        print(f"WARNING: envelope: {big} gap(s) larger than {1e-3 * span:.2g} m left open")
    for f in faces:
        if f["kind"] == "cap":
            f["kind"] = "lateral"
    return faces, n_cut


def build_folded_envelope(legs, leg_AB, reflections, cfg):
    """Per-leg envelopes (scotty_to_gid.build_envelope), cut by the mirrors,
    then the union; ENVPEC faces optionally moved onto the real surfaces.
    Returns (faces, info) or raises EnvelopeError."""
    d = cfg.envelope_separation
    leg_cells = []
    for k, (leg, (A, B)) in enumerate(zip(legs, leg_AB)):
        print(f"Envelope of leg {k} ({leg.n} cross-sections):")
        _, cells, _ = s2g.build_envelope(
            leg.points, A, B, d, n_sides=cfg.envelope_sides,
            ray_tol=cfg.envelope_ray_tolerance or None,
            max_len=None if cfg.envelope_max_cell_length == 0
            else max(cfg.envelope_max_cell_length, 0.0))
        leg_cells.append(cells)
    allv = np.vstack([P for cells in leg_cells for c in cells for P, _ in c["faces"] if len(P)])
    span = float(np.ptp(allv, axis=0).max()) + 4 * d
    offsets = [cfg.envelope_pec_offset] * len(reflections)
    faces, n_cut = _envelope_union(leg_cells, leg_AB, reflections, offsets, span, d)
    if cfg.envelope_pec_conform and reflections:
        faces = conform_pec_faces(faces, leg_cells, leg_AB, reflections, cfg, span, offsets)
    return faces, dict(span=span, cut_faces=n_cut, leg_volumes=_leg_volumes(leg_cells))


def _leg_volumes(leg_cells):
    """Volume of each leg's own envelope (before cutting at the mirrors)."""
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return None
    out = []
    for cells in leg_cells:
        v = 0.0
        for c in cells:
            P = np.vstack([F for F, _ in c["faces"] if len(F)])
            try:
                v += ConvexHull(P).volume
            except Exception:
                pass
        out.append(v)
    return out


def ermes_estimate(vol, lam, cfg, leg_vols=None, area=0.0):
    """Lines of text: tetrahedra / unknowns / memory / time ERMES will need for
    the volume (shared estimate of scotty_to_gid.py, ermes_* settings)."""
    p = s2g.fem_cost_parameters(cfg, prefix="ermes_")
    L, _ = s2g.fem_cost_lines(vol, area, 2 * np.pi / lam, p)
    if leg_vols:
        tot = sum(leg_vols) or 1.0
        parts = ", ".join(f"leg {k}: {100 * v / tot:.0f}%" for k, v in enumerate(leg_vols))
        L.insert(-3, f"  share of the volume       : {parts}")
    L.append("=" * 74)
    return L


def conform_pec_faces(faces, leg_cells, leg_AB, reflections, cfg, span, offsets):
    """ENVPEC faces follow the real reflecting surfaces (envelope_pec_conform).
    The legs are first cut by a plane BEHIND the mirror (deep enough for the
    whole face), then everything behind the real surface is trimmed away."""
    user = list(offsets)
    todo = []
    for r, rec in enumerate(reflections):
        srf = rec.surface_obj
        if srf is None:
            continue
        if srf.kind == "plane" and offsets[r] == 0:
            rec.pec_conformed = True                # already exact
            continue
        todo.append(r)
    if not todo:
        return faces

    def warn(r, err):
        print(f"WARNING: ENVPEC{reflections[r].number} left PLANAR (tangent plane of "
              f"'{reflections[r].surface_obj.name}'): {err}")

    faces, _ = _envelope_union(leg_cells, leg_AB, reflections, offsets, span,
                               cfg.envelope_separation, set(todo))
    hinfo = {}
    for _ in range(4):                              # find a deep enough cut
        changed = False
        for r in list(todo):
            rec = reflections[r]
            try:
                gmin, gmax, size = _pec_heights(faces, rec, rec.surface_obj, span, offsets[r])
            except _ConformError as err:
                warn(r, err)
                todo.remove(r)
                offsets[r] = user[r]
                changed = True                      # rebuild without the extended caps
                continue
            old = hinfo.get(r, (np.inf, -np.inf, 0.0))
            hinfo[r] = (min(old[0], gmin), max(old[1], gmax), max(old[2], size))
            need = -gmin + max(0.03 * size, 1e-6 * span)
            if need > offsets[r] + 1e-9 * span:
                offsets[r] = 1.1 * need
                changed = True
        if not changed:
            break
        faces, _ = _envelope_union(leg_cells, leg_AB, reflections, offsets, span,
                                   cfg.envelope_separation, set(todo))
    # trim; if a mirror fails, retry it with the cut plane further behind
    # (a different cut changes the near-degenerate intersections), then give up
    depth = [1.0, 1.7, 2.8]
    base = list(offsets)
    level = {r: 0 for r in todo}
    out = faces
    for _ in range(3 * len(depth)):
        failed = []
        out = faces
        for r in todo:
            try:
                out = _trim_by_surface(out, reflections[r], reflections[r].surface_obj, cfg,
                                       span, offsets[r], hinfo[r])
            except _ConformError as err:
                failed.append((r, err))
        if not failed:
            break
        for r, err in failed:
            level[r] += 1
            if level[r] >= len(depth):
                warn(r, err)
                todo.remove(r)
                offsets[r] = user[r]
            else:
                offsets[r] = base[r] * depth[level[r]]
        faces, _ = _envelope_union(leg_cells, leg_AB, reflections, offsets, span,
                                   cfg.envelope_separation, set(todo))
        for r in list(todo):
            try:
                g0, g1, sz = _pec_heights(faces, reflections[r], reflections[r].surface_obj,
                                          span, offsets[r])
                old = hinfo[r]
                hinfo[r] = (min(old[0], g0), max(old[1], g1), max(old[2], sz))
            except _ConformError as err:
                warn(r, err)
                todo.remove(r)
                offsets[r] = user[r]
                faces, _ = _envelope_union(leg_cells, leg_AB, reflections, offsets, span,
                                           cfg.envelope_separation, set(todo))
        out = faces
    for r in todo:
        rec = reflections[r]
        rec.pec_conformed = True
        idx = _pec_faces(out, rec)
        nt = sum(len(out[i]["tris"]) for i in idx)
        dz = max(np.ptp(np.concatenate([(L - rec.point) @ rec.normal
                                        for L in [out[i]["outer"]] + list(out[i]["holes"])]))
                 for i in idx)
        print(f"Envelope: ENVPEC{rec.number} follows surface '{rec.surface_obj.name}' "
              f"({nt} triangles, height variation {dz:.3g} m)")
    return out


# ##########################################################################
# 6. OUTPUT
# ##########################################################################

class IgesWriterHoles(s2g.IgesWriter):
    """scotty_to_gid's IGES writer + planar faces that may be non-convex and
    have holes (entity 144 with inner boundaries)."""

    def add_planar_face_holes(self, outer, holes, nrm, label="FACE"):
        outer = np.asarray(outer, float)
        nrm = _unit(nrm)
        c = outer.mean(0)
        L = np.linalg.norm(np.roll(outer, -1, 0) - outer, axis=1)
        k = int(np.argmax(L))
        e1 = _unit(outer[(k + 1) % len(outer)] - outer[k])
        e1 = _unit(e1 - np.dot(e1, nrm) * nrm)
        e2 = np.cross(nrm, e1)
        a, b = (outer - c) @ e1, (outer - c) @ e2
        pad = 0.02 * max(np.ptp(a), np.ptp(b))
        a0, a1, b0, b1 = a.min() - pad, a.max() + pad, b.min() - pad, b.max() + pad
        corner = lambda x, y: c + x * e1 + y * e2
        ctrl = np.array([[corner(a0, b0), corner(a0, b1)], [corner(a1, b0), corner(a1, b1)]])
        k01 = np.array([0, 0, 1, 1.])
        srf = self.add_surface(ctrl, np.ones((2, 2)), k01, 1, k01, 1, label=label,
                               status="00010000")

        def loop(poly):
            poly = np.asarray(poly, float)
            uv = np.stack([((poly - c) @ e1 - a0) / (a1 - a0),
                           ((poly - c) @ e2 - b0) / (b1 - b0), np.zeros(len(poly))], 1)
            n = len(poly)
            model = [self._add(110, 0, label, [110] + [self._r(x) for x in
                                                       (*poly[i], *poly[(i + 1) % n])], "00010000")
                     for i in range(n)]
            param = [self._add(110, 0, label, [110] + [self._r(x) for x in
                                                       (*uv[i], *uv[(i + 1) % n])], "00010500")
                     for i in range(n)]
            cm = self._add(102, 0, label, [102, n] + model, "00010000")
            cp = self._add(102, 0, label, [102, n] + param, "00010500")
            self.max_coord = max(self.max_coord, float(np.abs(poly).max()))
            return self._add(142, 0, label, [142, 0, srf, cp, cm, 3], "00010000")

        cos_out = loop(outer)
        cos_in = [loop(h) for h in holes]
        return self._add(144, 0, label, [144, srf, 1, len(cos_in), cos_out] + cos_in)

    def add_curved_face(self, f, sc, label="FACE"):
        """ENVPEC face on the real mirror: bicubic B-spline height patch over
        the tangent plane (fitted to the surface), trimmed by the face loops.
        Returns the max deviation of the patch at the face boundary [m]."""
        cv = f["curved"]
        gsr = _readers()
        ga, gb, H = cv["ga"], cv["gb"], cv["H"]
        o, e1, e2, nn = cv["origin"], cv["e1"], cv["e2"], cv["n"]
        p = 3

        def knots(m):
            inner = np.linspace(0.0, 1.0, m - p + 1)
            return np.concatenate([[0.0] * p, inner, [1.0] * p])
        nu = int(np.clip(len(ga) // 2 + 2, p + 1, 40))
        nv = int(np.clip(len(gb) // 2 + 2, p + 1, 40))
        ku, kv = knots(nu), knots(nv)
        xi = np.array([ku[i + 1:i + p + 1].mean() for i in range(nu)])
        eta = np.array([kv[j + 1:j + p + 1].mean() for j in range(nv)])
        da, db = ga[-1] - ga[0], gb[-1] - gb[0]
        Nu = gsr.bspline_basis(ku, p, (ga - ga[0]) / da)
        Nv = gsr.bspline_basis(kv, p, (gb - gb[0]) / db)
        C = np.linalg.pinv(Nu) @ H @ np.linalg.pinv(Nv).T
        A = ga[0] + xi * da
        B = gb[0] + eta * db
        ctrl = (o[None, None] + A[:, None, None] * e1 + B[None, :, None] * e2
                + C[..., None] * nn) * sc
        srf = self.add_surface(ctrl, np.ones((nu, nv)), ku, p, kv, p, label=label,
                               status="00010000")
        loops3 = [f["outer"]] + list(f["holes"])
        dev = 0.0

        def loop(poly, l2):
            nonlocal dev
            poly = np.asarray(poly, float) * sc
            uv = np.stack([(l2[:, 0] - ga[0]) / da, (l2[:, 1] - gb[0]) / db,
                           np.zeros(len(l2))], 1)
            hfit = np.einsum("mi,ij,mj->m", gsr.bspline_basis(ku, p, uv[:, 0]), C,
                             gsr.bspline_basis(kv, p, uv[:, 1]))
            dev = max(dev, float(np.abs(hfit - (poly / sc - o) @ nn).max()))
            n = len(poly)
            model = [self._add(110, 0, label, [110] + [self._r(x) for x in
                                                       (*poly[i], *poly[(i + 1) % n])], "00010000")
                     for i in range(n)]
            param = [self._add(110, 0, label, [110] + [self._r(x) for x in
                                                       (*uv[i], *uv[(i + 1) % n])], "00010500")
                     for i in range(n)]
            cm = self._add(102, 0, label, [102, n] + model, "00010000")
            cp = self._add(102, 0, label, [102, n] + param, "00010500")
            self.max_coord = max(self.max_coord, float(np.abs(poly).max()))
            return self._add(142, 0, label, [142, 0, srf, cp, cm, 3], "00010000")

        cos = [loop(L3, L2) for L3, L2 in zip(loops3, cv["loops2"])]
        self._add(144, 0, label, [144, srf, 1, len(cos) - 1, cos[0]] + cos[1:])
        return dev


FACE_LABELS = {"first": "ENVFIRST", "last": "ENVLAST", "lateral": "ENVSIDE"}


def face_label(kind):
    if kind.startswith("pec"):
        return "ENVPEC" + kind[3:]
    return FACE_LABELS.get(kind, "ENVSIDE")


def write_efield_ends(path, beam, ends, scale, source):
    p = lambda v: ",".join(f"{c * scale:.8g}" for c in v)
    with open(path, "w") as f:
        for name, i, orig in (("First", 0, ends[0]), ("Last", -1, ends[1])):
            extra = f"; {source} end point = {p(orig)}" if source == "Scotty" else ""
            f.write(f"# {name} surface  (position X,Y,Z = {p(beam.points[i])}{extra})\n")
            f.write(",".join(s2g._fmt_c(c) for c in beam.E[i]) + "\n")
            f.write(",".join(f"{np.real(c):.8g}" for c in beam.K[i]) + "\n")


def write_reflection_report(path, reflections, legs, cfg, unit, scale):
    v = lambda a: "(" + ", ".join(f"{x * scale:.6g}" for x in a) + ")"
    u = lambda a: "(" + ", ".join(f"{x:.6f}" for x in a) + ")"
    c = lambda a: "(" + ", ".join(s2g._fmt_c(x) for x in a) + ")"
    with open(path, "w") as f:
        f.write(f"# Reflections of the beam on PEC surfaces (lengths in {unit})\n")
        f.write(f"# legs: {len(legs)}; reflections: {len(reflections)}\n")
        for k, leg in enumerate(legs):
            L = arclength(leg.points)[-1]
            f.write(f"# leg {k}: {leg.n} points, length {L * scale:.6g}, from "
                    f"{v(leg.points[0])} to {v(leg.points[-1])}\n")
        for r in reflections:
            f.write(f"\n[reflection {r.number}]  surface '{r.surface}'  "
                    f"(envelope faces ENVPEC{r.number})\n")
            f.write(f"point               = {v(r.point)}\n")
            f.write(f"normal (to beam)    = {u(r.normal)}\n")
            f.write(f"direction in / out  = {u(r.t_in)} / {u(r.t_out)}\n")
            f.write(f"incidence angle     = {r.incidence_deg:.4f} deg (from the normal)\n")
            f.write(f"path length to here = {r.path_length * scale:.6g}\n")
            kk = np.linalg.eigvalsh(r.curvature)
            rad = ", ".join("inf" if abs(x) < 1e-12 else f"{1 / x * scale:.6g}" for x in kk)
            f.write(f"principal curvature radii = {rad}  (> 0 concave seen from the beam)\n")
            f.write(f"surface geometry    = {r.geometry_note}\n")
            if r.footprint is not None:
                f.write(f"beam footprint semi-axes on the surface = "
                        f"{', '.join(f'{x * scale:.6g}' for x in r.footprint)}\n")
                f.write(f"mirror sagitta over the footprint       = {r.sagitta * scale:.4g}"
                        + ("  (ENVPEC face follows the real surface)\n" if r.pec_conformed
                           else "  (deviation of the ENVPEC plane from the real surface)\n"))
            if r.focal_lengths is not None:
                fl = ", ".join("inf" if not np.isfinite(x) else f"{x * scale:.6g}"
                               for x in r.focal_lengths)
                f.write(f"focal lengths of the mirror for this beam = {fl}  (> 0 focusing)\n")
            if r.psi_in is not None:
                f.write(f"Psi in  (xx, xy, yy) [1/m] = {c([r.psi_in[0, 0], r.psi_in[0, 1], r.psi_in[1, 1]])}\n")
                f.write(f"Psi out (xx, xy, yy) [1/m] = {c([r.psi_out[0, 0], r.psi_out[0, 1], r.psi_out[1, 1]])}\n")
                w = beam_widths(r.psi_in) * cfg.envelope_factor
                f.write(f"beam semi-axes at the mirror = {', '.join(f'{x * scale:.6g}' for x in w)}\n")
            if r.E_in is not None:
                f.write(f"E in  = {c(r.E_in)}\nE out = {c(r.E_out)}\n")
            for w in r.warnings:
                f.write(f"WARNING: {w}\n")


def write_beam_leg(cfg, k, leg, A, B, theta, sc, out, stem):
    C_s, A_s, B_s = leg.points * sc, A * sc, B * sc
    if cfg.output_format in ("stl", "both"):
        idx = s2g.ring_selection(leg.n, cfg.stl_stride)
        verts, faces = s2g.build_tube_mesh(s2g.ellipse_rings(C_s[idx], A_s[idx], B_s[idx], theta),
                                           cfg.caps)
        path = out(f"{stem}_leg{k}_beam_surface.stl")
        s2g.write_stl(path, verts, faces, ascii=cfg.stl_ascii, name=f"beam_leg{k}")
        print(f"Wrote leg {k} STL ({len(idx)} sections, {len(faces)} triangles): {path}")
    if cfg.output_format in ("nurbs", "both"):
        patches = s2g.nurbs_beam_patches(C_s, A_s, B_s, tol=cfg.nurbs_tolerance * sc,
                                         n_patches=cfg.nurbs_patches, caps=cfg.caps,
                                         n_max=cfg.nurbs_max_sections, unit=cfg.units)
        iges = s2g.IgesWriter(unit=cfg.units)
        for pt in patches:
            iges.add_surface(pt["ctrl"], pt["w"], pt["ku"], pt["pu"], pt["kv"], pt["pv"],
                             closed_u=pt["closed_u"], label=pt["label"])
        if cfg.nurbs_include_ray:
            iges.add_polyline_curve(C_s)
        path = out(f"{stem}_leg{k}_beam_surface.igs")
        iges.write(path, description=f"Beam envelope, leg {k}")
        print(f"Wrote leg {k} NURBS/IGES ({len(patches)} surfaces): {path}")


# ##########################################################################
# 7. INTERACTIVE PLOT
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


def beam_plot(beam, A, B, n_around, legs, starts, reflections, surf_meshes,
              env_tris=None, env_faces=None):
    """Static 3D view: beam surface through every cross-section (one colour per
    leg), central rays, reflection points, reflecting surfaces and the FEM
    envelope. The window can still be rotated/zoomed."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    loop = lambda r: np.vstack([r, r[:1]])
    C = beam.points
    fig = plt.figure(figsize=(10, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    bounds = list(starts) + [beam.n]
    pts = []
    for k in range(len(bounds) - 1):
        a, b = bounds[k], min(bounds[k + 1] + 1, beam.n)      # legs share their joint
        if b - a < 2:
            continue
        R = _draw_beam_surface(ax, s2g.ellipse_rings, C[a:b], A[a:b], B[a:b], n_around,
                               colour=f"C{(k % 9) + (k > 0)}" if k else "C1")
        pts.append(R.reshape(-1, 3))
    for k, leg in enumerate(legs):
        ax.plot(*leg.points.T, color="k", lw=0.9)
    span = np.ptp(C, axis=0).max()
    for r in reflections:
        ax.scatter(*r.point, color="k", s=25)
        ax.quiver(*r.point, *(0.12 * span * r.normal), color="k", lw=1)
        ax.text(*(r.point + 0.03 * span * r.normal), f"R{r.number}", fontsize=8)
    for V, T in surf_meshes:
        if len(T):
            T = T[np.random.default_rng(0).permutation(len(T))[:4000]] if len(T) > 4000 else T
            ax.add_collection3d(Poly3DCollection(
                V[T], facecolors=(0.5, 0.5, 0.6, 0.18), edgecolors=(0.4, 0.4, 0.5, 0.3),
                linewidths=0.2))
    if env_tris is not None and len(env_tris):
        colours = []
        for f in env_faces:
            kind = f["kind"]
            col = ((0.2, 0.7, 0.2, 0.25) if kind == "first" else (0.8, 0.2, 0.2, 0.25)
                   if kind == "last" else (0.1, 0.3, 0.9, 0.3) if kind.startswith("pec")
                   else (1.0, 0.6, 0.0, 0.06))
            colours += [col] * f["ntri"]
        ax.add_collection3d(Poly3DCollection(env_tris, facecolors=colours, edgecolors="none"))
        for f in env_faces:
            for e in [f["outer"]] + list(f["holes"]):
                ax.plot(*loop(e).T, color="0.35", lw=0.4)
        pts.append(np.asarray(env_tris).reshape(-1, 3))
    _equal_axes(ax, np.vstack(pts) if pts else C)
    ax.set_xlabel("X (tokamak R) [m]")
    ax.set_ylabel("Y (toroidal) [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title(f"Beam surface: {beam.n} cross-sections, {len(legs)} leg(s), "
                 f"{len(reflections)} reflection(s)")
    ax.view_init(elev=20, azim=-60)
    plt.tight_layout()
    plt.show()


# ##########################################################################
# 8. COMMAND LINE / MAIN
# ##########################################################################

def _parameter_comments():
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
        elif current and s.startswith("#") and not head.strip() and line.startswith(" "):
            info[current][1] += "; " + comment
        else:
            current = None
    return info


class _HelpFormatter(argparse.RawDescriptionHelpFormatter,
                     argparse.ArgumentDefaultsHelpFormatter):
    pass


def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=_HelpFormatter, allow_abbrev=False,
        epilog="Options not listed here are passed to scotty_to_gid.py / beam_to_gid.py "
               "(see --source-help).")
    comments = _parameter_comments()
    groups = {}
    defaults = default_parameters()
    for name in USER_PARAMETERS:
        if name in NO_COMMAND_LINE:
            continue
        default = getattr(defaults, name)
        grp, text = comments.get(name, ("Parameters", ""))
        g = groups.get(grp) or groups.setdefault(grp, ap.add_argument_group(grp))
        flags = ["--" + name.replace("_", "-")] + PARAMETER_ALIASES.get(name, [])
        kw = dict(dest=name, default=default, help=(text or " ").replace("%", "%%"))
        if isinstance(default, bool):
            kw["action"] = argparse.BooleanOptionalAction
        elif isinstance(default, list):
            kw.update(nargs="*", type=str, metavar="FILE")
        else:
            kw["type"] = PARAMETER_TYPES.get(name, type(default))
            if name in PARAMETER_CHOICES:
                kw["choices"] = PARAMETER_CHOICES[name]
            else:
                kw["metavar"] = "VALUE"
        g.add_argument(*flags, **kw)
    g = ap.add_argument_group("Quick surfaces (added to section 1 'surfaces')")
    g.add_argument("--plane", nargs=6, type=float, action="append",
                   metavar=("PX", "PY", "PZ", "NX", "NY", "NZ"), help="infinite plane")
    g.add_argument("--sphere", nargs=4, type=float, action="append",
                   metavar=("CX", "CY", "CZ", "R"), help="sphere")
    g.add_argument("--cylinder", nargs=7, type=float, action="append",
                   metavar=("PX", "PY", "PZ", "AX", "AY", "AZ", "R"), help="infinite cylinder")
    g = ap.add_argument_group("Tools")
    g.add_argument("--show-params", action="store_true",
                   help="Print every parameter in use (this script + beam script) and stop")
    g.add_argument("--source-help", action="store_true",
                   help="Show the options of the beam script (scotty_to_gid / beam_to_gid)")
    return ap


def parse_args(argv=None):
    """Returns (cfg of this script, cfg of the beam script)."""
    argv = sys.argv[1:] if argv is None else argv
    # argparse takes negative numbers in scientific notation (-1e-05) for
    # options: write them as plain decimals (also for the beam script)
    sci = re.compile(r"^-(\d+\.?\d*|\.\d+)[eE][-+]?\d+$")
    argv = [np.format_float_positional(float(a), trim="-") if sci.match(a) else a
            for a in argv]
    cfg, rest = build_parser().parse_known_args(argv)
    cfg.surfaces = list(_PARAMETER_DEFAULTS["surfaces"])
    mod = source_module(cfg)
    if cfg.source_help:
        mod.build_parser().print_help()
        sys.exit(0)
    scfg = mod.parse_args(rest)
    return cfg, scfg


def show_parameters(cfg, scfg, mod):
    defaults = default_parameters()
    print(f"reflect_beam_to_gid.py:")
    for name in USER_PARAMETERS:
        val = getattr(cfg, name)
        mark = "   <- changed" if val != getattr(defaults, name) else ""
        print(f"  {name:<28} = {val!r}{mark}")
    for extra in ("plane", "sphere", "cylinder"):
        if getattr(cfg, extra):
            print(f"  --{extra:<26} {getattr(cfg, extra)}")
    print(f"\n{os.path.basename(mod.__file__)} (beam settings):")
    mod.show_parameters(scfg)


def main(cfg=None, scfg=None):
    if cfg is None:
        cfg, scfg = parse_args()
    mod = source_module(cfg)
    if scfg is None:
        scfg = mod.default_parameters()
    if getattr(scfg, "ermes_input", False) and hasattr(mod, "apply_ermes_input"):
        mod.apply_ermes_input(scfg)
    if cfg.show_params:
        show_parameters(cfg, scfg, mod)
        return

    # ---- the incident beam ----
    src = incident_beam(cfg, scfg)
    beam = src["beam"]
    s_all = arclength(beam.points)
    print(f"Incident beam: {beam.n} cross-sections, path length {s_all[-1]:.4g} m")

    # ---- surfaces (meshes only kept within reach of the beam) ----
    specs = surface_specs(cfg)
    j0 = src["launch_index"] if cfg.skip_backward_extension else 0
    reach = s_all[-1] - s_all[j0] + 0.1
    region = np.array([beam.points[j0] - reach, beam.points[j0] + reach])
    surfs = [m for m in (make_surface(s, cfg, region) for s in specs) if m is not None]
    if not surfs:
        print("NOTE: no reflecting surfaces defined (section 1 'surfaces' / 'surface_files', "
              "or --plane / --sphere / --cylinder / --surface-files): no reflections.")

    # ---- reflections ----
    legs, reflections = trace_reflections(src, surfs, cfg)
    full, starts = join_legs(legs)
    print(f"Beam after reflections: {len(legs)} leg(s), {len(reflections)} reflection(s), "
          f"total path {arclength(full.points)[-1]:.4g} m")

    leg_AB = []
    for leg in legs:
        L = arclength(leg.points)[-1]
        A, B = s2g.beam_cross_sections(leg, fallback_radius=0.02 * L or 0.01)
        leg_AB.append((A * cfg.envelope_factor, B * cfg.envelope_factor))
    A_full = np.vstack([leg_AB[0][0]] + [ab[0][1:] for ab in leg_AB[1:]])
    B_full = np.vstack([leg_AB[0][1]] + [ab[1][1:] for ab in leg_AB[1:]])
    wmin, wmax = s2g.section_widths(A_full, B_full)
    print(f"Beam semi-axes {wmin.min():.3g} .. {wmax.max():.3g} m "
          f"(envelope_factor {cfg.envelope_factor:g})")

    if surfs:
        contacts = side_contacts(legs, leg_AB, surfs, reflections,
                                 cfg.envelope_separation if cfg.envelope else 0.0)
        for k, s0, s1, name, what in contacts:
            if what == "beam":
                print(f"WARNING: leg {k}, {s0:.3g}-{s1:.3g} m along it: the BEAM (1/e x "
                      f"envelope_factor) crosses surface '{name}' although the central ray "
                      "does not hit it there -- that part of the beam is not reflected here.")
            else:
                print(f"NOTE: leg {k}, {s0:.3g}-{s1:.3g} m along it: the FEM envelope "
                      f"(separation {cfg.envelope_separation:g} m) crosses surface '{name}' "
                      "(the beam itself is clear) -- trim the volume in GiD if needed.")

    # ---- FEM envelope ----
    env_faces, env_tris, env_info, env_vol = None, None, None, None
    if cfg.envelope:
        try:
            env_faces, info = build_folded_envelope(legs, leg_AB, reflections, cfg)
            env_info = info
            parts = [f["tris"] if "tris" in f else
                     triangulate_face(f["outer"], f["holes"], f["normal"]) for f in env_faces]
            for f, t in zip(env_faces, parts):
                f["ntri"] = len(t)
            env_tris = np.concatenate([t for t in parts if len(t)])
            vol, bad = mesh_check(env_tris, info["span"])
            env_vol = vol
            n_pec = sum(f["kind"].startswith("pec") for f in env_faces)
            n_cur = sum("curved" in f for f in env_faces)
            print(f"FEM envelope: {len(env_faces)} faces ({n_pec} on mirrors, {n_cur} of them "
                  f"curved), volume "
                  f"{vol:.4g} m^3, {len(legs)} leg(s) joined"
                  + ("" if bad == 0 else f"  WARNING: {bad} unmatched edges (not watertight)"))
        except s2g.EnvelopeError as err:
            print(f"WARNING: no FEM envelope written: {err}")
            env_faces = env_tris = None

    # ---- plot ----
    theta = np.linspace(0, 2 * np.pi, cfg.stl_points_per_section, endpoint=False)
    box = np.array([full.points.min(0), full.points.max(0)])
    pad = 0.15 * np.ptp(full.points, axis=0).max() + 0.05
    box = box + np.array([[-pad] * 3, [pad] * 3])
    surf_meshes = [s.display_mesh(box) for s in surfs]
    if cfg.show_plot:
        beam_plot(full, A_full, B_full, cfg.plot_points_around, legs, starts,
                         reflections, surf_meshes, env_tris, env_faces)

    # ---- outputs ----
    os.makedirs(cfg.out_dir, exist_ok=True)
    stem = src["stem"]
    out = lambda name: os.path.join(cfg.out_dir, name)
    sc = UNIT_SCALE[cfg.units]
    for k, (leg, (A, B)) in enumerate(zip(legs, leg_AB)):
        write_beam_leg(cfg, k, leg, A, B, theta, sc, out, stem)

    if env_faces:
        if cfg.output_format in ("stl", "both"):
            path = out(f"{stem}_envelope.stl")
            V, T = merge_vertices(env_tris * sc)
            s2g.write_stl(path, V, T, ascii=cfg.stl_ascii, name="fem_envelope")
            print(f"Wrote envelope STL ({len(T)} triangles, {cfg.units}): {path}")
        if cfg.output_format in ("nurbs", "both"):
            path = out(f"{stem}_envelope.igs")
            iges = IgesWriterHoles(unit=cfg.units)
            for f in env_faces:
                if "curved" in f:
                    dev = iges.add_curved_face(f, sc, label=face_label(f["kind"]))
                    print(f"  {face_label(f['kind'])}: B-spline patch on '{f['surface']}', "
                          f"max deviation at its edges {dev:.2g} m")
                    continue
                iges.add_planar_face_holes(f["outer"] * sc, [h * sc for h in f["holes"]],
                                           f["normal"], label=face_label(f["kind"]))
            iges.write(path, description="Polyhedral FEM envelope around the reflected beam")
            print(f"Wrote envelope IGES ({len(env_faces)} faces, {cfg.units}): {path}")

    if cfg.write_surfaces and surfs:
        tris = [V[T] for V, T in surf_meshes if len(T)]
        if tris:
            V, T = merge_vertices(np.concatenate(tris) * sc)
            path = out(f"{stem}_surfaces.stl")
            s2g.write_stl(path, V, T, ascii=cfg.stl_ascii, name="pec_surfaces")
            print(f"Wrote reflecting surfaces near the beam ({len(T)} triangles): {path}")

    if cfg.write_central_ray:
        path = out(f"{stem}_central_ray.dat")
        s2g.write_central_ray_file(path, full.points * sc)
        print(f"Wrote central ray ({full.n} points, all legs): {path}")
    if full.E is not None and full.K is not None:
        if cfg.write_efield:
            write_efield_ends(out(cfg.efield_file), full, src["ends"], sc,
                              "Scotty" if cfg.beam_source == "scotty" else "vacuum")
            print(f"Wrote E-field/k end-point summary: {out(cfg.efield_file)}")
        if cfg.write_field_arrays:
            s2g.write_vector_rows(out(f"{stem}_Efield_along_ray.dat"), full.E, True)
            s2g.write_vector_rows(out(f"{stem}_Kvector_along_ray.dat"), full.K, False)
            print(f"Wrote per-point E-field and k ({full.n} rows each).")
    k0 = src.get("k0")
    if k0 and full.E is not None and full.K is not None:
        omega = k0 * SPEED_OF_LIGHT
        if cfg.write_ermes_launch and hasattr(s2g, "write_ermes_launch"):
            s2g.write_ermes_launch(out(f"{stem}_ermes_launch.txt"), scfg, full, k0, omega)
        if getattr(scfg, "ox_diagnostics", False) and hasattr(s2g, "ox_conversion_diagnostics"):
            # waist of the launched beam: at the beam START (vacuum part built from
            # the launch definition; Scotty's first point is on the plasma side of
            # the Psi_BC jump, so its Psi does not give the launched waist)
            i0 = 0
            w0_eff = None
            if full.psi_xx is not None:
                w0s = [w for _, w, _ in s2g.psi_waist_parameters(
                    full.psi_xx[i0], full.psi_xy[i0], full.psi_yy[i0], k0)]
                if all(np.isfinite(w0s)):
                    w0_eff = float(np.sqrt(w0s[0] * w0s[1]))
            note = ("Scotty's refracted K on the traced part; straight vacuum lines on the "
                    "extensions and after reflections." if cfg.beam_source == "scotty" else
                    "STRAIGHT lines, no refraction (vacuum k0 in the plasma); run Scotty "
                    "for the refracted N_par / N_y.")
            try:
                s2g.ox_conversion_diagnostics(
                    dict(beam=full.points, distance=arclength(full.points), omega=omega),
                    full.K, scfg, w0_eff, out(f"{stem}_OX_along_ray.dat"), note=note)
            except Exception as err:
                print(f"WARNING: O-X diagnostics failed: {err}")
    if cfg.write_reflections and reflections:
        path = out(f"{stem}_reflections.txt")
        write_reflection_report(path, reflections, legs, cfg, cfg.units, sc)
        print(f"Wrote reflection report: {path}")

    # ---- size of the ERMES problem ----
    if env_vol and k0:
        T = np.asarray(env_tris, float)
        area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]),
                                    axis=1).sum()
        lines = ermes_estimate(env_vol, 2 * np.pi / k0, cfg,
                               env_info.get("leg_volumes") if env_info else None, area)
        print("\n" + "\n".join(lines))
        try:
            with open(out(f"{stem}_ermes_estimate.txt"), "w") as f:
                f.write("\n".join(lines) + "\n")
        except OSError:
            pass
    elif cfg.envelope:
        print("\nNo ERMES estimate: " + ("no FEM envelope" if not env_vol
                                         else "frequency unknown (no omega in the source)"))
    return dict(legs=legs, reflections=reflections, envelope=env_faces, beam=full)


if __name__ == "__main__":
    main()
