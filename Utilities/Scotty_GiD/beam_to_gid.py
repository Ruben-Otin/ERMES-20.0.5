#!/usr/bin/env python3
"""
beam_to_gid.py
=====================

Propagate a Gaussian beam through VACUUM (free space) from a user-given
starting point and turn it into geometry that GiD can import -- the same
outputs as scotty_to_gid.py, but without Scotty and without a plasma:

  * STL   : triangulated "beam envelope" (tube of elliptical cross-sections,
            capped at both ends).  GiD: Files > Import > STL.
  * NURBS : the same envelope as exact rational B-spline surfaces (IGES
            entity 128) in an .igs file.  GiD: Files > Import > IGES.
  * text  : central ray, E-field / k at both ends (same format as
            scotty_to_gid.py), E and k along the ray, and a summary of the
            beam parameters (widths, curvatures, Psi, Gouy phase) at both ends.

What you give (section 1)
-------------------------
  * frequency
  * start point (Cartesian X, Y, Z, or cylindrical R, zeta, Z)
  * direction of propagation (a vector, or TORBEAM-like poloidal/toroidal
    angles)
  * beam profile, either
       beam_spec = "waist"  : waist radius w0 and the signed distance from the
                              start point to the waist, along the beam, or
       beam_spec = "launch" : width and curvature AT the start point (the
                              same quantities as Scotty's launch_beam_width /
                              launch_beam_curvature)
    each given for the two principal axes, so astigmatic / elliptical beams
    (with a rotation of the ellipse) are allowed
  * polarisation (Jones vector in the transverse basis)
  * propagation_length: the beam is computed from the start point (s = 0) to
    s = propagation_length, measured along the propagation direction.

Physics (paraxial Gaussian beam in free space)
----------------------------------------------
Same conventions as Scotty: field ~ exp(i k s + i w.Psi.w / 2) with w the
transverse displacement, so that the 1/e amplitude contour is
w^T Im(Psi) w = 2.  Along each principal axis i

    Psi_i(s) = k / ((s - d_i) - i zR_i),      zR_i = k w0_i^2 / 2

    W_i(s)   = w0_i sqrt(1 + ((s - d_i)/zR_i)^2)        (1/e amplitude radius)
    1/R_i(s) = (s - d_i) / ((s - d_i)^2 + zR_i^2)       (> 0: diverging)

which is exactly the law Psi^-1(s) = Psi^-1(0) + (s/k) I used by
scotty_to_gid.py for its vacuum extensions.  In vacuum the ray is a straight
line, k = k0 t_hat is constant and the (paraxial) polarisation is constant.

Transverse basis
----------------
x_hat = component of ``transverse_x_hint`` perpendicular to the direction,
y_hat = t_hat x x_hat, so (x_hat, y_hat, t_hat) is right-handed.  The beam's
principal axes are x_hat, y_hat rotated by ``ellipse_rotation`` about t_hat,
and the polarisation Jones vector is given in (x_hat, y_hat).

Parameters
----------
All parameters are plain assignments in section 1 ("USER PARAMETERS"), each
with its possible values in a comment. Edit them there, or change any of
them for one run from the command line with  --name value  (same name, "-"
instead of "_").
    python vacuum_beam_to_gid.py --show-params     # values in use
    python vacuum_beam_to_gid.py --help            # all flags with defaults

Usage examples
--------------
    # parameters of section 1
    python vacuum_beam_to_gid.py

    # 1 m of propagation, start at (2.2, 0, -0.1), going in -X, 3 cm waist
    # 0.6 m in front of the start point, no plot window
    python vacuum_beam_to_gid.py --start-point 2.2 0 -0.1 --direction -1 0 0 \\
        --propagation-length 1.0 --waist-width 0.03 --waist-distance 0.6 \\
        --no-show-plot

    # Scotty-style launch: width and curvature at the start point
    python vacuum_beam_to_gid.py --beam-spec launch --launch-width 0.04 \\
        --launch-curvature -0.9

    # elliptical beam, rotated 30 deg, circular polarisation, NURBS only
    python vacuum_beam_to_gid.py --waist-width 0.02 0.04 --ellipse-rotation 30 \\
        --polarization rcp --output-format nurbs

Requires Python >= 3.9 and numpy (+ matplotlib for the plot).
"""

from __future__ import annotations

import argparse
import datetime
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
# e.g.  python vacuum_beam_to_gid.py --frequency 32e9 --start-point 1.0 0 1.4 \
#              --propagation-length 0.8 --output-format both --no-show-plot
# "python vacuum_beam_to_gid.py --show-params" prints the values in use.
###############################################################################

_names_before_parameters = set(globals())

# ---------------------------- What to do -------------------------------------

show_plot = True               # True / False: interactive 3D plot of the beam
                               # (the script waits until the window is closed)

# ---------------------------- Wave -------------------------------------------

frequency = 50e9               # frequency [Hz]  (lambda = c / frequency)

# ---------------------------- Start point and direction ----------------------

start_point = (2.2, 0.0, -0.1) # starting point of the beam [m]
start_point_cylindrical = True # True : start_point = (R [m], zeta [rad], Z [m])
                               #        (as Scotty's launch_position)
                               # False: start_point = (X, Y, Z) [m]
                               # X = R cos(zeta), Y = R sin(zeta)
direction_from_angles = True   # True : direction from poloidal/toroidal_angle
                               # False: direction from the vector 'direction'
direction = (-1.0, 0.0, 0.0)   # propagation direction (X, Y, Z), any length
poloidal_angle = 20.0          # [deg] TORBEAM-like: 0,0 = towards -R (inwards);
                               # poloidal > 0 tilts the beam towards -Z
toroidal_angle = 20.0           # [deg] toroidal > 0 tilts the beam towards -zeta
                               # t_hat = -(cos p cos t) R_hat - (cos p sin t) zeta_hat
                               #         - (sin p) Z_hat
propagation_length = 2.0       # length of beam to compute [m], measured from
                               # start_point along the propagation direction
n_points = 401                 # number of cross-sections along the beam

# ---------------------------- Beam profile -----------------------------------
# Two values = the two principal axes (1, 2); one value on the command line
# (e.g. --waist-width 0.03) means a circular beam.

beam_spec = "waist"            # "waist" : waist_width + waist_distance
                               # "launch": launch_width + launch_curvature
waist_width = (0.03, 0.03)     # w0: 1/e AMPLITUDE radius at the waist [m]
                               # (1/e^2 intensity radius, as usual for Gaussian beams)
waist_distance = (0.5, 0.5)    # signed distance start point -> waist along the
                               # beam [m]; > 0: waist in front of the start point,
                               # < 0: behind it (beam already diverging)
launch_width = (0.04, 0.04)    # W at the start point, 1/e amplitude radius [m]
launch_curvature = (-1.0, -1.0) # 1/R at the start point [1/m]; < 0 converging
                               # (focusing), > 0 diverging, 0 = waist at start
                               # (same meaning as Scotty's launch_beam_curvature)
ellipse_rotation = 0.0         # [deg] rotation of principal axis 1 from x_hat
                               # towards y_hat (about the propagation direction)
transverse_x_hint = (0.0, 0.0, 1.0)  # defines x_hat: its part perpendicular to the
                               # direction (if parallel, another axis is used)

# ---------------------------- Polarisation -----------------------------------

polarization = "x"             # Jones vector in (x_hat, y_hat):
                               # "x" | "y" | "+45" | "-45" | "rcp" | "lcp"
                               # "linear:<deg>"  linear, <deg> from x_hat to y_hat
                               # "<ex>,<ey>"     complex, e.g. "1,1j" or "0.6,0.8"
                               # (normalised automatically)
time_convention = "-iwt"       # "-iwt": fields ~ exp(-i w t) (physics, Scotty)
                               # "+jwt": fields ~ exp(+j w t) (engineering)
                               # only used to give rcp/lcp the right handedness
                               # (IEEE: rcp rotates clockwise seen from the source)

# ---------------------------- Geometry output (GiD) --------------------------

envelope_factor = 1.0          # envelope radius = factor x W (1/e amplitude radius)
                               # 1.0: 86.5 % of the power inside, 1.5: 98.9 %,
                               # 2.0: 99.97 %
output_format = "both"         # "stl"   : triangle mesh   (GiD: Import > STL)
                               # "nurbs" : NURBS surfaces  (GiD: Import > IGES)
                               # "both"  | "none"
units = "m"                    # "m" | "cm" | "mm"  unit of all written geometry
caps = True                    # True : close both ends of the beam | False: open tube
out_dir = "./GiD_input/"       # folder for ALL files written;
                               # created automatically if it does not exist
file_stem = None               # start of the output file names
                               # None: vacuum_beam_<freq>GHz_L<length>m
stl_points_per_section = 64    # STL: points around each cross-section
stl_stride = 1                 # STL: use every N-th cross-section (1 = all)
stl_ascii = False              # STL: True = ASCII file | False = binary file
nurbs_tolerance = 1e-4         # NURBS: max deviation from any cross-section [m]
nurbs_max_sections = 600       # NURBS: max cross-sections the surface goes through
nurbs_patches = 4              # NURBS: 4 = tube in 4 quarter surfaces (best for GiD)
                               #        1 = one closed surface
nurbs_include_ray = False      # NURBS: True = add the central ray as a curve

# ---------------------------- Other output files -----------------------------

write_central_ray = True       # True / False: <stem>_central_ray.dat (x,y,z per line)
write_efield = True            # True / False: E-field and k at both beam ends
efield_file = "Efield_output_data.dat"  # name of that end-point file
write_field_arrays = True      # True / False: E-field and k at every ray point
                               # (<stem>_Efield_along_ray.dat, <stem>_Kvector_along_ray.dat)
write_beam_parameters = True   # True / False: <stem>_beam_parameters.txt (widths,
                               # curvatures, Psi, Gouy phase, waists at both ends)
plot_every = 20                # plot: spacing (ray points) of the extra ellipses

###############################################################################
# (end of user parameters -- nothing below needs to be edited)
###############################################################################

USER_PARAMETERS = [n for n in globals()
                   if n not in _names_before_parameters and not n.startswith("_")]
_PARAMETER_DEFAULTS = {n: globals()[n] for n in USER_PARAMETERS}   # snapshot

# Values that are checked when given on the command line
PARAMETER_CHOICES = {
    "beam_spec": ["waist", "launch"],
    "time_convention": ["-iwt", "+jwt"],
    "output_format": ["stl", "nurbs", "both", "none"],
    "units": ["m", "cm", "mm"],
    "nurbs_patches": [1, 4],
}
# Extra / shorter command-line names
PARAMETER_ALIASES = {
    "frequency": ["--freq"],
    "poloidal_angle": ["--pol"],
    "toroidal_angle": ["--tor"],
    "propagation_length": ["--length"],
    "polarization": ["--polarisation"],
    "output_format": ["--format"],
    "show_plot": ["--plot"],
    "write_efield": ["--efield"],
    "write_field_arrays": ["--field-arrays"],
    "write_central_ray": ["--central-ray", "--txt"],
    "nurbs_tolerance": ["--nurbs-tol"],
}
PARAMETER_TYPES = {"file_stem": str, "polarization": str}
PARAMETER_METAVARS = {
    "start_point": ("P1", "P2", "P3"),
    "direction": ("DX", "DY", "DZ"),
    "transverse_x_hint": ("HX", "HY", "HZ"),
}
# Two-axis parameters: one value on the command line = both axes
TWO_AXIS_PARAMETERS = ("waist_width", "waist_distance", "launch_width", "launch_curvature")

SPEED_OF_LIGHT = 299_792_458.0  # m/s
UNIT_SCALE = {"m": 1.0, "cm": 100.0, "mm": 1000.0}


def default_parameters():
    """The values of section 1 as a namespace. From another script:
        import vacuum_beam_to_gid as v
        cfg = v.default_parameters(); cfg.propagation_length = 0.5; v.main(cfg)"""
    return argparse.Namespace(**dict(_PARAMETER_DEFAULTS), show_params=False)


# ##########################################################################
# 2. GAUSSIAN BEAM IN VACUUM
# ##########################################################################

@dataclass
class Beam:
    """Everything known at each point of the ray."""
    points: np.ndarray                 # (N,3) central ray, Cartesian
    xhat: np.ndarray                   # (N,3) local transverse basis
    yhat: np.ndarray                   # (N,3)
    psi_xx: Optional[np.ndarray]       # (N,) complex, in the (xhat, yhat) frame
    psi_xy: Optional[np.ndarray]
    psi_yy: Optional[np.ndarray]
    distance: np.ndarray               # (N,) distance along the ray from the start
    E: Optional[np.ndarray] = None     # (N,3) complex polarisation, Cartesian
    K: Optional[np.ndarray] = None     # (N,3) wavevector, Cartesian [1/m]

    @property
    def n(self):
        return self.points.shape[0]


@dataclass
class GaussianBeamVacuum:
    """Simple-astigmatic Gaussian beam on a straight axis.

    Principal axis i (i = 0, 1) has unit vector u[i], waist at signed
    distance d[i] from the start point and Rayleigh range zR[i]."""
    k: float                           # wavenumber [1/m]
    start: np.ndarray                  # (3,) start point
    t: np.ndarray                      # (3,) unit propagation direction
    xhat: np.ndarray                   # (3,) transverse basis
    yhat: np.ndarray
    rot: float                         # principal-axis rotation [rad]
    d: np.ndarray                      # (2,) waist positions [m]
    zR: np.ndarray                     # (2,) Rayleigh ranges [m]
    jones: np.ndarray                  # (2,) complex, in (xhat, yhat)

    @property
    def u(self):
        c, s = np.cos(self.rot), np.sin(self.rot)
        return np.array([c * self.xhat + s * self.yhat, -s * self.xhat + c * self.yhat])

    @property
    def w0(self):
        return np.sqrt(2.0 * self.zR / self.k)

    def psi_principal(self, s):
        """(len(s), 2) complex Psi along the principal axes."""
        s = np.atleast_1d(np.asarray(s, float))[:, None]
        return self.k / ((s - self.d) - 1j * self.zR)

    def widths(self, s):
        """(len(s), 2) 1/e amplitude radii W_i(s)."""
        s = np.atleast_1d(np.asarray(s, float))[:, None]
        return self.w0 * np.sqrt(1.0 + ((s - self.d) / self.zR) ** 2)

    def curvatures(self, s):
        """(len(s), 2) wavefront curvature 1/R_i [1/m] (> 0 diverging)."""
        s = np.atleast_1d(np.asarray(s, float))[:, None]
        return (s - self.d) / ((s - self.d) ** 2 + self.zR ** 2)

    def gouy(self, s):
        """Gouy phase [rad] relative to the waist(s), simple astigmatism."""
        s = np.atleast_1d(np.asarray(s, float))[:, None]
        return 0.5 * np.sum(np.arctan((s - self.d) / self.zR), axis=1)

    def psi_xy_frame(self, s):
        """Psi_xx, Psi_xy, Psi_yy in the (xhat, yhat) frame."""
        P = self.psi_principal(s)
        c, sn = np.cos(self.rot), np.sin(self.rot)
        return (P[:, 0] * c * c + P[:, 1] * sn * sn,
                (P[:, 0] - P[:, 1]) * c * sn,
                P[:, 0] * sn * sn + P[:, 1] * c * c)

    @property
    def E_cart(self):
        return self.jones[0] * self.xhat + self.jones[1] * self.yhat


def unit_vector(v, what="vector"):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n == 0 or not np.isfinite(n):
        sys.exit(f"The {what} {tuple(v)} has zero length.")
    return v / n


def start_point_cartesian(cfg):
    p = np.asarray(cfg.start_point, float)
    if cfg.start_point_cylindrical:
        R, zeta, Z = p
        return np.array([R * np.cos(zeta), R * np.sin(zeta), Z])
    return p


def propagation_direction(cfg, start):
    """Unit propagation direction (Cartesian)."""
    if not cfg.direction_from_angles:
        return unit_vector(cfg.direction, "direction")
    p, t = np.deg2rad([cfg.poloidal_angle, cfg.toroidal_angle])
    zeta = np.arctan2(start[1], start[0]) if np.hypot(start[0], start[1]) > 0 else 0.0
    k_R, k_phi, k_Z = -np.cos(p) * np.cos(t), -np.cos(p) * np.sin(t), -np.sin(p)
    c, s = np.cos(zeta), np.sin(zeta)
    return unit_vector([k_R * c - k_phi * s, k_R * s + k_phi * c, k_Z])


def transverse_basis(t, hint):
    """x_hat = part of ``hint`` perpendicular to t; y_hat = t x x_hat."""
    candidates = [np.asarray(hint, float), np.array([0, 0, 1.0]),
                  np.array([0, 1.0, 0]), np.array([1.0, 0, 0])]
    for i, h in enumerate(candidates):
        if np.linalg.norm(h) == 0:
            continue
        h = h / np.linalg.norm(h)
        x = h - np.dot(h, t) * t
        if np.linalg.norm(x) > 1e-3:
            if i > 0:
                print(f"NOTE: transverse_x_hint is (almost) parallel to the direction; "
                      f"using {tuple(candidates[i])} instead to define x_hat.")
            x = x / np.linalg.norm(x)
            return x, np.cross(t, x)
    raise RuntimeError("could not build a transverse basis")


def jones_vector(spec, time_convention="-iwt"):
    """Normalised complex Jones vector (ex, ey) from the polarization string."""
    s = str(spec).strip().lower().replace(" ", "")
    r2 = 1 / np.sqrt(2)
    # IEEE: rcp rotates clockwise seen from the source, i.e. from x_hat towards
    # y_hat (right-hand rule about t_hat). With exp(-iwt): x + i y.
    sgn = 1.0 if time_convention == "-iwt" else -1.0
    presets = {"x": (1, 0), "y": (0, 1), "+45": (r2, r2), "-45": (r2, -r2),
               "rcp": (r2, sgn * 1j * r2), "lcp": (r2, -sgn * 1j * r2)}
    if s in presets:
        v = np.array(presets[s], complex)
    elif s.startswith("linear:"):
        a = np.deg2rad(float(s.split(":", 1)[1]))
        v = np.array([np.cos(a), np.sin(a)], complex)
    elif "," in s:
        try:
            v = np.array([complex(x) for x in s.split(",")], complex)
        except ValueError:
            sys.exit(f"Cannot read polarization '{spec}' (use e.g. '1,1j').")
        if v.size != 2:
            sys.exit(f"polarization '{spec}' must have exactly two components.")
    else:
        sys.exit(f"Unknown polarization '{spec}'. Use x, y, +45, -45, rcp, lcp, "
                 "linear:<deg> or '<ex>,<ey>'.")
    n = np.linalg.norm(v)
    if n == 0:
        sys.exit("polarization vector has zero length.")
    return v / n


def two(values, name):
    v = np.asarray(values, float).reshape(-1)
    if v.size == 1:
        v = np.repeat(v, 2)
    if v.size != 2:
        sys.exit(f"{name} needs one or two values, got {v.size}.")
    return v


def gaussian_beam_from_config(cfg) -> GaussianBeamVacuum:
    k = 2 * np.pi * cfg.frequency / SPEED_OF_LIGHT
    start = start_point_cartesian(cfg)
    t = propagation_direction(cfg, start)
    xhat, yhat = transverse_basis(t, cfg.transverse_x_hint)
    clean = lambda v: np.where(np.abs(v) < 1e-14, 0.0, v) + 0.0     # no "-0"
    t, xhat, yhat = clean(t), clean(xhat), clean(yhat)

    if cfg.beam_spec == "waist":
        w0 = two(cfg.waist_width, "waist_width")
        if np.any(w0 <= 0):
            sys.exit("waist_width must be > 0.")
        d = two(cfg.waist_distance, "waist_distance")
        zR = k * w0 ** 2 / 2
    else:  # "launch": Psi(0) = k * curvature + 2i / W^2, as in Scotty
        W = two(cfg.launch_width, "launch_width")
        if np.any(W <= 0):
            sys.exit("launch_width must be > 0.")
        c = two(cfg.launch_curvature, "launch_curvature")
        inv = 1.0 / (k * c + 2j / W ** 2)          # = (-d - i zR) / k
        d, zR = -k * inv.real, -k * inv.imag

    return GaussianBeamVacuum(k=k, start=start, t=t, xhat=xhat, yhat=yhat,
                              rot=np.deg2rad(cfg.ellipse_rotation), d=d, zR=zR,
                              jones=jones_vector(cfg.polarization, cfg.time_convention))


def sample_beam(gb: GaussianBeamVacuum, length, n_points) -> Beam:
    """Cross-sections from s = 0 to s = length along the beam."""
    s = np.linspace(0.0, length, max(int(n_points), 2))
    n = len(s)
    pxx, pxy, pyy = gb.psi_xy_frame(s)
    return Beam(points=gb.start + s[:, None] * gb.t,
                xhat=np.tile(gb.xhat, (n, 1)), yhat=np.tile(gb.yhat, (n, 1)),
                psi_xx=pxx, psi_xy=pxy, psi_yy=pyy, distance=s,
                E=np.tile(gb.E_cart, (n, 1)), K=np.tile(gb.k * gb.t, (n, 1)))


def envelope_sections(gb: GaussianBeamVacuum, beam: Beam, factor=1.0):
    """(A, B): conjugate semi-diameters of each envelope cross-section,
    P(theta) = C + A cos(theta) + B sin(theta), with A, B along the principal
    axes (constant directions in vacuum, so no re-phasing is needed)."""
    W = factor * gb.widths(beam.distance)
    u = gb.u
    return W[:, 0, None] * u[0], W[:, 1, None] * u[1]


def describe_beam(gb: GaussianBeamVacuum, length):
    lam = 2 * np.pi / gb.k
    print(f"Wavelength {lam * 1e3:.4g} mm, k = {gb.k:.6g} 1/m")
    print(f"Start point (X,Y,Z) = ({gb.start[0]:.6g}, {gb.start[1]:.6g}, {gb.start[2]:.6g}) m")
    print(f"Direction t_hat     = ({gb.t[0]:.6f}, {gb.t[1]:.6f}, {gb.t[2]:.6f})")
    print(f"x_hat               = ({gb.xhat[0]:.6f}, {gb.xhat[1]:.6f}, {gb.xhat[2]:.6f})")
    print(f"y_hat               = ({gb.yhat[0]:.6f}, {gb.yhat[1]:.6f}, {gb.yhat[2]:.6f})")
    end = gb.start + length * gb.t
    print(f"End point   (X,Y,Z) = ({end[0]:.6g}, {end[1]:.6g}, {end[2]:.6g}) m")
    for i in range(2):
        theta = 2.0 / (gb.k * gb.w0[i])
        print(f"Axis {i + 1}: w0 = {gb.w0[i]:.4g} m, waist at s = {gb.d[i]:.4g} m, "
              f"zR = {gb.zR[i]:.4g} m, far-field half-angle = {np.rad2deg(theta):.3g} deg")
        if theta > 0.3:
            print(f"WARNING: axis {i + 1}: divergence {np.rad2deg(theta):.3g} deg -- the "
                  "paraxial Gaussian-beam model is not accurate for such a tight waist.")
    jx, jy = gb.jones
    print(f"Polarisation (x_hat, y_hat) = ({_fmt_c(jx)}, {_fmt_c(jy)})")


# ##########################################################################
# 5. CROSS-SECTION GEOMETRY
# ##########################################################################

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


def write_stl(path, verts, faces, ascii=False, name="vacuum_beam"):
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
        f.write(b"Gaussian beam envelope surface (vacuum)".ljust(80, b" "))
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

    The along-beam parameter of every point is its normalised arc
    length. Starting from ``n_init`` sections evenly spaced in arc length,
    the worst-fitted section in every gap is added until all sections
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
        self.entities = []   # (type, form, label, [params])
        self.max_coord = 0.0

    @staticmethod
    def _r(x):
        return f"{float(x):.12E}"          # always contains a decimal point

    def add_surface(self, ctrl, w, ku, pu, kv, pv, closed_u=False, label="SURF"):
        """ctrl (nu, nv, 3), w (nu, nv); u varies fastest in IGES order."""
        nu, nv = w.shape
        prm = [128, nu - 1, nv - 1, pu, pv, int(closed_u), 0,
               0 if np.ptp(w) > 0 else 1, 0, 0]
        prm += [self._r(x) for x in ku] + [self._r(x) for x in kv]
        prm += [self._r(w[i, j]) for j in range(nv) for i in range(nu)]
        prm += [self._r(c) for j in range(nv) for i in range(nu) for c in ctrl[i, j]]
        prm += [self._r(ku[pu]), self._r(ku[-pu - 1]), self._r(kv[pv]), self._r(kv[-pv - 1])]
        self.entities.append((128, 0, label, prm))
        self.max_coord = max(self.max_coord, float(np.abs(ctrl).max()))

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
        self.entities.append((126, 0, label, prm))
        self.max_coord = max(self.max_coord, float(np.abs(pts).max()))

    @staticmethod
    def _hollerith(s):
        return f"{len(s)}H{s}"

    def write(self, path, description="Gaussian beam envelope"):
        unit_flag = {"m": (6, "M"), "mm": (2, "MM"), "cm": (10, "CM")}[self.unit]
        now = datetime.datetime.now().strftime("%Y%m%d.%H%M%S")
        H = self._hollerith
        g = [H(","), H(";"), H("vacuum_beam_to_gid"), H(os.path.basename(path)),
             H("vacuum_beam_to_gid.py"), H("1.0"), "32", "38", "6", "308", "15",
             H("GiD"), "1.0", str(unit_flag[0]), H(unit_flag[1]), "1", "1.0",
             H(now), "1.0E-9", self._r(max(self.max_coord, 1.0)), H("user"),
             H("-"), "11", "0", H(now)]
        g_lines = self._wrap(g, 72, terminator=";")
        s_lines = [description[:72]]

        d_lines, p_lines = [], []
        for k, (etype, form, label, prm) in enumerate(self.entities):
            de_seq = 2 * k + 1
            lines = self._wrap([str(x) for x in prm], 64, terminator=";")
            p_start = len(p_lines) + 1
            p_lines += [f"{ln:<64}{de_seq:>8}" for ln in lines]
            d_lines.append(f"{etype:>8}{p_start:>8}{0:>8}{0:>8}{0:>8}{0:>8}{0:>8}{0:>8}{'00000000':>8}")
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



def write_efield_endpoints(path, beam: Beam, scale):
    """E (complex) and k (real, Cartesian) at the first and last cross-section,
    same layout as scotty_to_gid.py (comment, E line, k line per surface)."""
    p = lambda v: ",".join(f"{c * scale:.8g}" for c in v)
    with open(path, "w") as f:
        for name, i in (("First", 0), ("Last", -1)):
            f.write(f"# {name} surface  (position X,Y,Z = {p(beam.points[i])}; "
                    f"vacuum Gaussian beam, s = {beam.distance[i] * scale:.8g})\n")
            f.write(",".join(_fmt_c(c) for c in beam.E[i]) + "\n")
            f.write(",".join(f"{np.real(c):.8g}" for c in beam.K[i]) + "\n")


def write_beam_parameters(path, gb: GaussianBeamVacuum, beam: Beam, scale, unit,
                          factor, cfg):
    """Human-readable summary of the beam at both ends (lengths in ``unit``;
    Psi, k and curvature are converted consistently: 1/unit)."""
    v3 = lambda v: f"({v[0]:.8g}, {v[1]:.8g}, {v[2]:.8g})"
    L = lambda x: x * scale
    lines = [
        "# Gaussian beam in vacuum -- written by vacuum_beam_to_gid.py "
        f"({datetime.datetime.now():%Y-%m-%d %H:%M})",
        "# Convention: field ~ exp(i k s + i w.Psi.w/2); 1/e amplitude contour w^T Im(Psi) w = 2",
        f"# Lengths in {unit}; Psi, k, 1/R in 1/{unit}",
        f"frequency_Hz          = {cfg.frequency:.10g}",
        f"wavelength            = {L(2 * np.pi / gb.k):.8g}",
        f"k0                    = {gb.k / scale:.8g}",
        f"propagation_length    = {L(beam.distance[-1]):.8g}",
        f"start_point           = {v3(L(gb.start))}",
        f"end_point             = {v3(L(beam.points[-1]))}",
        f"t_hat                 = {v3(gb.t)}",
        f"x_hat                 = {v3(gb.xhat)}",
        f"y_hat                 = {v3(gb.yhat)}",
        f"principal_axis_1      = {v3(gb.u[0])}",
        f"principal_axis_2      = {v3(gb.u[1])}",
        f"jones_(x_hat,y_hat)   = {_fmt_c(gb.jones[0])}, {_fmt_c(gb.jones[1])}",
        f"E_hat_cartesian       = {', '.join(_fmt_c(c) for c in gb.E_cart)}",
        f"envelope_factor       = {factor:g}   (geometry radius = factor * W)",
    ]
    for i in range(2):
        lines += [
            f"w0_axis{i + 1}              = {L(gb.w0[i]):.8g}",
            f"waist_s_axis{i + 1}         = {L(gb.d[i]):.8g}   (distance from start point)",
            f"waist_point_axis{i + 1}     = {v3(L(gb.start + gb.d[i] * gb.t))}",
            f"rayleigh_range_axis{i + 1}  = {L(gb.zR[i]):.8g}",
        ]
    for name, idx in (("FIRST", 0), ("LAST", -1)):
        s = beam.distance[idx]
        W, c = gb.widths(s)[0], gb.curvatures(s)[0]
        P = gb.psi_principal(s)[0]
        lines += [
            f"",
            f"[{name} surface]  s = {L(s):.8g}",
            f"position              = {v3(L(beam.points[idx]))}",
            f"W_axis1, W_axis2      = {L(W[0]):.8g}, {L(W[1]):.8g}   (1/e amplitude radius)",
            f"envelope_radius_1,2   = {L(factor * W[0]):.8g}, {L(factor * W[1]):.8g}",
            f"curvature_1, _2 (1/R) = {c[0] / scale:.8g}, {c[1] / scale:.8g}   (> 0 diverging)",
            f"Psi_axis1, Psi_axis2  = {_fmt_c(P[0] / scale)}, {_fmt_c(P[1] / scale)}",
            f"Psi_xx, Psi_xy, Psi_yy = {_fmt_c(beam.psi_xx[idx] / scale)}, "
            f"{_fmt_c(beam.psi_xy[idx] / scale)}, {_fmt_c(beam.psi_yy[idx] / scale)}",
            f"gouy_phase_rad        = {gb.gouy(s)[0]:.8g}",
        ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ##########################################################################
# 9. INTERACTIVE PLOT
# ##########################################################################

def interactive_plot(beam: Beam, A, B, theta, every):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import CheckButtons, Slider

    C = beam.points
    loop = lambda r: np.vstack([r, r[:1]])
    ring = lambda i: loop(ellipse_rings(C[i:i + 1], A[i:i + 1], B[i:i + 1], theta)[0])

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    plt.subplots_adjust(bottom=0.22)
    ax.plot(*C.T, color="k", lw=1.2, label="Central ray")
    pts = ring(0)
    ell, = ax.plot(*pts.T, color="C1", lw=2, label="Beam cross-section")
    dot, = ax.plot(*C[:1].T, "o", color="C1")
    ctx = [ax.plot(*ring(i).T, color="C0", lw=0.6, alpha=0.35, visible=False)[0]
           for i in range(0, beam.n, max(every, 1))]

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.legend(loc="upper left")
    ax.view_init(elev=20, azim=-60)
    lim = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    c, r = lim.mean(1), 0.5 * np.ptp(lim, axis=1).max()
    ax.set_xlim3d(c[0] - r, c[0] + r)
    ax.set_ylim3d(c[1] - r, c[1] + r)
    ax.set_zlim3d(c[2] - r, c[2] + r)

    slider = Slider(plt.axes([0.2, 0.08, 0.6, 0.03]), "Ray index", 0, beam.n - 1,
                    valinit=0, valstep=1)

    def update(_):
        i = int(slider.val)
        p = ring(i)
        ell.set_data(p[:, 0], p[:, 1])
        ell.set_3d_properties(p[:, 2])
        dot.set_data([C[i, 0]], [C[i, 1]])
        dot.set_3d_properties([C[i, 2]])
        ax.set_title(f"index {i}/{beam.n - 1} | s = {beam.distance[i]:.4g} m")
        fig.canvas.draw_idle()

    check = CheckButtons(plt.axes([0.02, 0.02, 0.22, 0.06]), ["Show all cross-sections"], [False])

    def toggle(_):
        vis = check.get_status()[0]
        for ln in ctx:
            ln.set_visible(vis)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    check.on_clicked(toggle)
    update(0)
    plt.show()



# ##########################################################################
# 10. COMMAND LINE / MAIN
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
        kw = dict(dest=name, default=default, help=(text or " ").replace("%", "%%"))
        if isinstance(default, bool):
            kw["action"] = argparse.BooleanOptionalAction
        elif name in TWO_AXIS_PARAMETERS:
            kw.update(nargs="+", type=float, metavar="VALUE")
        elif isinstance(default, tuple):
            kw.update(nargs=len(default), type=float,
                      metavar=PARAMETER_METAVARS.get(name, "VALUE"))
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
    return ap


def parse_args(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    cfg = build_parser().parse_args(argv)
    for name in USER_PARAMETERS:
        val = getattr(cfg, name)
        if isinstance(val, list):
            if name in TWO_AXIS_PARAMETERS and len(val) == 1:
                val = val * 2
            setattr(cfg, name, tuple(val))
    # giving a direction vector on the command line means "use it"
    if any(a == "--direction" for a in argv) and not any(
            a in ("--direction-from-angles", "--no-direction-from-angles") for a in argv):
        cfg.direction_from_angles = False
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
    if cfg.show_params:
        show_parameters(cfg)
        return
    if cfg.propagation_length <= 0:
        sys.exit("propagation_length must be > 0.")

    # ---- the beam ----
    gb = gaussian_beam_from_config(cfg)
    describe_beam(gb, cfg.propagation_length)
    beam = sample_beam(gb, cfg.propagation_length, cfg.n_points)
    A, B = envelope_sections(gb, beam, cfg.envelope_factor)
    wmin, wmax = section_widths(A, B)
    print(f"Beam: {beam.n} cross-sections over {cfg.propagation_length:.4g} m, "
          f"envelope semi-axes {wmin.min():.3g} .. {wmax.max():.3g} m "
          f"(envelope_factor {cfg.envelope_factor:g})")

    theta = np.linspace(0, 2 * np.pi, cfg.stl_points_per_section, endpoint=False)
    if cfg.show_plot:
        interactive_plot(beam, A, B, theta, cfg.plot_every)

    # ---- outputs ----
    os.makedirs(cfg.out_dir, exist_ok=True)
    stem = cfg.file_stem or (f"vacuum_beam_{cfg.frequency / 1e9:g}GHz_"
                             f"L{cfg.propagation_length:g}m").replace(".", "_")
    out = lambda name: os.path.join(cfg.out_dir, name)
    sc = UNIT_SCALE[cfg.units]
    C_s, A_s, B_s = beam.points * sc, A * sc, B * sc

    if cfg.output_format in ("stl", "both"):
        idx = ring_selection(beam.n, cfg.stl_stride)
        verts, faces = build_tube_mesh(ellipse_rings(C_s[idx], A_s[idx], B_s[idx], theta),
                                       cfg.caps)
        path = out(f"{stem}_beam_surface.stl")
        write_stl(path, verts, faces, ascii=cfg.stl_ascii, name="vacuum_beam")
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
        iges.write(path, description="Gaussian beam envelope in vacuum")
        print(f"Wrote NURBS/IGES ({len(patches)} surfaces, {cfg.units}): {path}")

    if cfg.write_central_ray:
        path = out(f"{stem}_central_ray.dat")
        write_central_ray_file(path, C_s)
        print(f"Wrote central ray ({beam.n} points): {path}")

    if cfg.write_efield:
        write_efield_endpoints(out(cfg.efield_file), beam, sc)
        print(f"Wrote E-field/k end-point summary: {out(cfg.efield_file)}")

    if cfg.write_field_arrays:
        write_vector_rows(out(f"{stem}_Efield_along_ray.dat"), beam.E, True)
        write_vector_rows(out(f"{stem}_Kvector_along_ray.dat"), beam.K, False)
        print(f"Wrote per-point E-field and k ({beam.n} rows each).")

    if cfg.write_beam_parameters:
        path = out(f"{stem}_beam_parameters.txt")
        write_beam_parameters(path, gb, beam, sc, cfg.units, cfg.envelope_factor, cfg)
        print(f"Wrote beam parameters: {path}")


if __name__ == "__main__":
    main()
