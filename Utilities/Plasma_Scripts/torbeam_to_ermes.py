#!/usr/bin/env python3
"""
torbeam_to_ermes.py
===================
Electron density and magnetic field files for an ERMES 3D model, evaluated
at every node of the GiD mesh (``<problem>-1.dat``), from a Scotty/TORBEAM
topfile (B, psi) or a G-EQDSK file, and a density profile ne(rho_pol).
Optionally adds a synthetic turbulent density field that emulates tokamak
edge turbulence, and field-aligned blobs / holes.

This is also the engine of eqdsk_to_ermes_turb.py (which sets
equilibrium = "eqdsk" and the file names used by problem_generator.py).

Equilibrium (section 1, ``equilibrium``)
----------------------------------------
"topfile" : Scotty/TORBEAM topfile (JSON or ASCII) with B and psi_n.
"eqdsk"   : EFIT G-EQDSK file, read with eqdsk_to_scotty.py (same folder):
            B_R = -(1/R) dpsi/dZ, B_Z = (1/R) dpsi/dR, B_t = F(psi)/R (vacuum F
            outside the plasma), psi_n = (psi - psi_axis)/(psi_LCFS - psi_axis).
            Nodes with psi_n < 1 outside the LCFS polygon (private-flux region
            below an X-point) are treated as outside the plasma.
The density is NOT in an EQDSK: it always comes from ne_file.

Output (sorted by node ID; same format as gen2Dfullwavefile.py):
    ne file  : "<node_id> <ne [m^-3]>"
    mag file : "<node_id> <Bx> <By> <Bz>"   [T], Cartesian, mesh frame
    (optional) dn file : "<node_id> <delta n / n>"   (turbulence only)

Geometry (section 1, ``geometry``)
----------------------------------
"slab"     : as gen2Dfullwavefile.py: x = R, y = Z, the toroidal direction is
             -z (straight). ne and B do not depend on z; the turbulence does.
"toroidal" : full 3D torus. The tokamak axis is the Cartesian axis
             ``tokamak_axis`` (default "z", with x, y horizontal):
             R = sqrt(x^2 + y^2), phi = atan2(y, x), and
             B = B_R R_hat + B_t phi_hat + B_Z Z_hat (phi_hat counter-clockwise
             seen from +axis; B_t as signed in the topfile).

Turbulence model (``turbulence = True``)
---------------------------------------
Edge/core drift-wave turbulence is, to a good approximation, a random density
field that is
  * field-aligned  : very long correlation along B (L_par ~ metres, k_par ~ 0),
  * short across B : radial and binormal correlation lengths of ~1-5 cm,
  * stronger at the edge: delta n / n ~ 0.1-1 % in the core, ~5-10 % at the
    pedestal and 10-50 % (intermittent, positively skewed) in the SOL.
It is built as a sum of N random-phase Fourier modes
    F(x) = sqrt(2/N) sum_j cos(k_r,j x_r + k_b,j x_b + k_par,j x_par + phase_j)
in local field-aligned coordinates at each node:
    x_r   : radial distance [m] (flux-surface label rho_pol mapped to metres
            at the outboard midplane), so the eddies follow the flux surfaces,
    x_b   : binormal distance [m], r * (theta - p * phi), with
            p = d theta / d phi the local field-line pitch (constant along a
            field line, so the field is aligned with B),
    x_par : distance along B [m] (~ R * phi).
The wavenumbers are drawn from a Gaussian spectrum, which gives the
correlation exp(-(dr/L_r)^2 - (db/L_b)^2 - (dpar/L_par)^2); an optional
binormal peak k_b0 adds drift-wave-like structure. F has zero mean and unit
variance; the density is
    "gaussian"  : ne = ne0 (1 + A(rho) F)                 (clipped at 0)
    "lognormal" : ne = ne0 exp(s F - s^2 / 2), s^2 = ln(1 + A^2)
                  (positive, same mean and rms, positively skewed like the
                   intermittent SOL blobs)
with A(rho) the rms delta n / n profile ``turb_amplitude_profile``.
B is not perturbed (delta B / B ~ 1e-4 is negligible for the wave).

Blobs and holes (``blobs = True``)
----------------------------------
Intermittent edge transport is carried by filaments: coherent, field-aligned
structures of excess density ("blobs"), born near the separatrix and moving out
through the SOL, a few cm across and metres long along B; and density deficits
("holes") left just inside the separatrix. They are added on top of the
turbulence as discrete objects in the same field-aligned coordinates:
    blob : ne += dn_k exp(-(dr/w_r)^2 - (db/w_b)^2 - (dpar/L_par)^2)
    hole : ne *= 1 - h_k exp(-(...)^2)
with random positions, log-normal sizes and amplitudes. The blob density is
dn_k = blob_amplitude * n_ref * exp(-(r_k - r_sep)/blob_decay_length), with n_ref
the profile density at blob_reference_rho, so blobs keep carrying density far
into the SOL, where they dominate (delta n / n >> 1, as in the far SOL). Their
number follows the packing fraction (fraction of the radial x binormal area of
their region covered by blobs; SOL measurements: ~5-30 %).
The field-aligned coordinates use a local pitch, which is accurate for a
domain that spans a limited poloidal / toroidal range (an ERMES model).

Usage
-----
Edit section 1 and run ``python torbeam_to_ermes.py``, or override any
parameter on the command line:  --name value  (``--help``, ``--show-params``),
e.g.  python torbeam_to_ermes.py --equilibrium eqdsk --turbulence
Typical chain: eqdsk_to_scotty.py writes Scotty_input/topfile_sample, which is
the default ``topfile`` here.
Requires numpy, scipy (pandas speeds up writing; matplotlib for the plot).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import numpy as np

# ##########################################################################
# 1. USER PARAMETERS
# ##########################################################################

# ---------------------------- Input files ------------------------------------

mesh_file = "./MASTU-2D-t20_p20.gid/MASTU-2D-t20_p20-1.dat"  # GiD mesh nodes: No[ID] = p(x,y,z);
mesh_units = "m"               # "m" | "cm" | "mm": units of the mesh coordinates
equilibrium = "topfile"        # "topfile": B and psi_n from ``topfile``
                               # "eqdsk"  : B and psi_n from the G-EQDSK ``eqdsk_file``
topfile = "./Scotty_input/topfile_sample"  # Scotty JSON (2D or flattened) or
                               # TORBEAM ASCII topfile (eqdsk_to_scotty.py)
eqdsk_file = "./EQDSK_Samples/mast-u-sample.eqdsk"  # G-EQDSK, only read if equilibrium = "eqdsk"
                               # (eqdsk_to_ermes_turb.py passes its own EQDSK_FILE)
eqdsk_psi_in_weber = False     # True if psi in the EQDSK is in Wb (not Wb/rad)
eqdsk_flip_bp = False          # True to flip B_R and B_Z (other COCOS convention)
eqdsk_use_lcfs = True          # True: psi_n < 1 outside the LCFS polygon (private
                               # flux region) is treated as outside the plasma
ne_file = "./EQDSK_Samples/rho_ne.txt"  # ne profile, either
                               #  * two columns "x  ne" (an optional first line
                               #    with the point count is skipped; Scotty format)
                               #  * one column "ne", with x read from ne_x_file
ne_x_file = "./EQDSK_Samples/rho_flux.txt"  # radial coordinate for a ONE-column ne_file ("" = none)
                               # (as PROFILE_X_FILE in eqdsk_to_scotty.py)
ne_radial_coordinate = "rho"   # x in the ne files: "rho" = sqrt(psi_n) | "psi" = psi_n
ne_units = 1.0                 # factor from the ne file values to m^-3 (rho_ne.txt is in m^-3)
                               # (Scotty ne<sfx>.dat: 1e19; raw profiles in m^-3: 1)
ne_outside = "constant"        # ne beyond the last point of the profile:
                               # "constant": its last value (a flat SOL tail continues)
                               # "zero"    : 0
sol_ne_separatrix = 0.0        # SOL model: ne at the separatrix [m^-3]; 0 = off (ne
                               # as in ne_file, which may be 0 outside rho = 1)
sol_decay_length = 0.01        # SOL model: e-folding length of ne outside the
                               # separatrix [m] (outboard-midplane distance)
sol_blend_width = 0.02         # SOL model: rho_pol width inside the separatrix over
                               # which the profile blends into ne_sep (smooth join)

# ---------------------------- Geometry ---------------------------------------

geometry = "toroidal"          # "slab"     : x = R, y = Z, toroidal = -z (straight),
                               #              as gen2Dfullwavefile.py
                               # "toroidal" : full 3D torus, axis = tokamak_axis
tokamak_axis = "z"             # "toroidal": Cartesian axis along the tokamak axis
                               # "x" | "y" | "z" (the other two, in cyclic order,
                               # are the horizontal plane)
flip_bt = False                # True: reverse the toroidal field (B_t -> -B_t)
spline_order = 3               # 1..5: spline order for psi, B (topfile) and ne(psi)

# ---------------------------- Output files -----------------------------------

out_dir = "./Plasma"           # folder for the output files
ne_out = "ne_3D.dat"           # electron density file   (node, ne [m^-3])
mag_out = "Be_3D.dat"          # magnetic field file     (node, Bx, By, Bz [T])
write_turbulence_file = True   # True: also <dn_out> (node, delta n / n), turbulence only
dn_out = "dn_3D.dat"           # turbulence file
plot = True                    # True: <out_dir>/plasma_check.png (ne, delta n / n, |B|
                               # on a cut through the mesh and the profiles)

# ---------------------------- Turbulence -------------------------------------

turbulence = False             # True : add the synthetic turbulent density field
turb_amplitude_profile = [(0.0, 0.003), (0.7, 0.005), (0.9, 0.02), (0.97, 0.05),
                          (1.0, 0.10), (1.2, 0.30)]
                               # rms delta n / n vs rho_pol: (rho, A) pairs,
                               # linearly interpolated, constant outside
turb_Lr = 0.015                # radial correlation length [m]   (edge: 1-3 cm)
turb_Lb = 0.030                # binormal correlation length [m] (edge: 2-5 cm)
turb_Lpar = 5.0                # parallel correlation length [m] (field-aligned:
                               # metres; 0 = exactly k_par = 0)
turb_kb_peak = 0.0             # binormal wavenumber of the spectral peak [1/m];
                               # 0 = Gaussian spectrum centred at k = 0
                               # (e.g. 2*pi/0.05 for 5 cm drift-wave structures)
turb_modes = 256               # number of random Fourier modes (more = smoother
                               # statistics, slower: time ~ nodes x modes)
turb_statistics = "lognormal"  # "lognormal" (positive, skewed, SOL-like) |
                               # "gaussian" (symmetric, clipped at ne >= 0)
turb_seed = 12345              # random seed (same seed = same turbulence)

# ---------------------------- Blobs and holes (filaments) --------------------

blobs = False                  # True : add field-aligned blobs (and holes, below)
blob_packing_fraction = 0.15   # fraction of the radial x binormal area of the blob
                               # region covered by blobs (SOL: ~0.05-0.3)
blob_rho_range = (0.98, 1.40)  # rho_pol range of the blob centres
blob_width_r = 0.012           # radial 1/e half-width [m]   (SOL blobs: 0.5-3 cm)
blob_width_b = 0.020           # binormal 1/e half-width [m] (SOL blobs: 1-4 cm)
blob_length_par = 10.0         # parallel 1/e half-length [m] (filaments: metres)
blob_size_spread = 0.3         # log-normal spread of the sizes (0 = all equal)
blob_amplitude = 1.0           # peak density excess / n_ref at the separatrix (SOL: ~0.5-2)
blob_amplitude_spread = 0.5    # log-normal spread of the amplitudes
blob_decay_length = 0.04       # e-folding of the blob amplitude with radial distance
                               # outside the separatrix [m] (slower than the profile:
                               # blobs carry density into the far SOL)
blob_reference_rho = 1.0       # n_ref = profile density at this rho (blob_ne_reference = 0);
                               # 1.0 = separatrix density (profiles ending at 0 at
                               # rho = 1: use e.g. 0.98, or blob_ne_reference)
blob_ne_reference = 0.0        # n_ref [m^-3]; 0 = from the profile at blob_reference_rho
hole_fraction = 1.0            # hole packing relative to the blob packing (0 = no holes)
hole_rho_range = (0.90, 1.00)  # rho_pol range of the hole centres
hole_amplitude = 0.3           # peak relative density deficit of the holes (0..1)
blob_seed = 2024               # random seed for the blobs and holes

# ---------------------------- Performance ------------------------------------

chunk_nodes = 200_000          # nodes evaluated at once (memory ~ chunk x modes x 8 B)

# ##########################################################################
# 2. INTERNAL SETTINGS
# ##########################################################################

USER_PARAMETERS = [n for n in list(globals()) if not n.startswith("_") and n not in (
    "annotations", "argparse", "json", "os", "re", "sys", "time", "np")]
PARAMETER_CHOICES = {"mesh_units": ["m", "cm", "mm"], "ne_radial_coordinate": ["rho", "psi"],
                     "ne_outside": ["constant", "zero"],
                     "equilibrium": ["topfile", "eqdsk"],
                     "geometry": ["slab", "toroidal"], "tokamak_axis": ["x", "y", "z"],
                     "turb_statistics": ["lognormal", "gaussian"]}
PARAMETER_NARGS = {"blob_rho_range": 2, "hole_rho_range": 2}
UNIT_SCALE = {"m": 1.0, "cm": 0.01, "mm": 0.001}
NODE_RE = re.compile(rb"No\[(\d+)\]\s*=\s*p\(\s*([^,\s]+)\s*,\s*([^,\s]+)\s*,\s*([^)\s]+)\s*\)")


# ##########################################################################
# 3. INPUT READERS
# ##########################################################################

def read_mesh_nodes(path, block_bytes=64 << 20):
    """GiD .dat nodes 'No[ID] = p(x,y,z);' -> (ids int64, xyz (n,3) float64).
    Streams the file in blocks (meshes of tens of millions of nodes)."""
    ids, xyz = [], []
    t0 = time.time()
    with open(path, "rb") as f:
        rest = b""
        while True:
            block = f.read(block_bytes)
            if not block:
                break
            block = rest + block
            cut = block.rfind(b"\n")
            if cut < 0:
                rest = block
                continue
            rest, block = block[cut + 1:], block[:cut]
            m = NODE_RE.findall(block)
            if m:
                a = np.array(m, dtype=object)
                ids.append(a[:, 0].astype(np.int64))
                xyz.append(a[:, 1:].astype(np.float64))
        if rest:
            m = NODE_RE.findall(rest)
            if m:
                a = np.array(m, dtype=object)
                ids.append(a[:, 0].astype(np.int64))
                xyz.append(a[:, 1:].astype(np.float64))
    if not ids:
        raise RuntimeError(f"No nodes 'No[ID] = p(x,y,z);' found in {path}.")
    ids, xyz = np.concatenate(ids), np.concatenate(xyz)
    print(f"Read {ids.size:,} nodes from {path} ({time.time() - t0:.1f} s)")
    return ids, xyz


def _read_torbeam_ascii(path):
    """TORBEAM ASCII topfile -> dict(R, Z, Br, Bt, Bz, pol_flux) with 2D [Z][R]."""
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
    excl = ("number of", "grid points", "psi_sep", "outside radius")

    def find(*keys):
        for test in (lambda l, k: l == k, lambda l, k: l.startswith(k), lambda l, k: k in l):
            for k in keys:
                for lab, v in sections:
                    if test(lab, k) and not any(e in lab for e in excl):
                        return np.asarray(v, float)
        raise KeyError(f"{path}: no section {keys}")
    R, Z = find("x-coord", "r-coord"), find("z-coord")
    out = dict(R=R, Z=Z)
    for key, names in (("Br", ("b_r",)), ("Bt", ("b_t",)), ("Bz", ("b_z",)),
                       ("pol_flux", ("psi", "poloidal flux"))):
        out[key] = find(*names).reshape(Z.size, R.size)
    return out


def load_topfile(path, k):
    """Topfile -> dict of RectBivariateSpline (R, Z) for psi, BR, Bt, BZ + grid.
    The magnetic axis is taken at the minimum of psi_n on the grid."""
    from scipy.interpolate import RectBivariateSpline
    if path.lower().endswith(".json"):
        with open(path) as f:
            data = json.load(f)
    else:
        data = _read_torbeam_ascii(path)
    low = {kk.lower(): kk for kk in data}
    get = lambda n: np.asarray(data[low[n.lower()]], float)
    R, Z = get("R").ravel(), get("Z").ravel()

    def grid(name):
        a = get(name)
        if a.ndim == 1 and a.size == R.size * Z.size:
            a = a.reshape(Z.size, R.size)                 # flattened, R fastest
        if a.shape == (Z.size, R.size):
            a = a.T                                       # stored [Z][R] (Scotty)
        if a.shape != (R.size, Z.size):
            raise ValueError(f"{path}: {name} has shape {a.shape}")
        return a
    sp = {key: RectBivariateSpline(R, Z, grid(name), kx=k, ky=k, s=0)
          for key, name in (("psi", "pol_flux"), ("BR", "Br"), ("Bt", "Bt"), ("BZ", "Bz"))}
    P = grid("pol_flux")
    i, j = np.unravel_index(np.argmin(P), P.shape)
    sp.update(R=R, Z=Z, R_axis=R[i], Z_axis=Z[j])
    return sp


def _numeric_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            try:
                v = [float(t) for t in line.replace(",", " ").split()]
            except ValueError:
                continue
            if v:
                rows.append(v)
    return rows


def _import_eqdsk_module():
    """eqdsk_to_scotty.py from the folder of this script."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        import eqdsk_to_scotty
    except ModuleNotFoundError as err:
        if err.name != "eqdsk_to_scotty":
            raise
        sys.exit(f"equilibrium = 'eqdsk' needs eqdsk_to_scotty.py in {here}.")
    return eqdsk_to_scotty


def load_eqdsk(path, k, psi_in_weber=False, flip_bp=False):
    """G-EQDSK -> same dict as load_topfile (+ the LCFS polygon), using the
    field calculation of eqdsk_to_scotty.py."""
    from scipy.interpolate import RectBivariateSpline
    if not path or not os.path.isfile(path):
        sys.exit(f"equilibrium = 'eqdsk' but eqdsk_file = '{path}' does not exist.")
    E = _import_eqdsk_module()
    eq = E.read_geqdsk(path)
    BR, Bt, BZ, psi_n = E.compute_fields(eq, psi_in_weber=psi_in_weber, flip_bp=flip_bp,
                                         interp_order=k)
    R, Z = np.asarray(eq.R, float), np.asarray(eq.Z, float)
    sp = {key: RectBivariateSpline(R, Z, a, kx=k, ky=k, s=0)
          for key, a in (("psi", psi_n), ("BR", BR), ("Bt", Bt), ("BZ", BZ))}
    sp.update(R=R, Z=Z, R_axis=float(eq.rmaxis), Z_axis=float(eq.zmaxis), lcfs=None)
    if np.size(eq.rbbbs) >= 3:
        sp["lcfs"] = np.column_stack([eq.rbbbs, eq.zbbbs])
    print(f"EQDSK {path}: {eq.nw} x {eq.nh} grid, I_p = {eq.current / 1e6:.3f} MA, "
          f"B_t(R0 = {eq.rcentr:.3f} m) = {eq.bcentr:.3f} T, "
          f"LCFS polygon {'with ' + str(len(eq.rbbbs)) + ' points' if sp['lcfs'] is not None else 'missing'}")
    return sp


def load_ne_profile(path, coordinate, units, k, x_file="", outside="constant"):
    """ne file -> f(psi_n) [m^-3], clipped at >= 0; beyond the data: the last
    value (outside = "constant") or 0 ("zero").
    Two columns (x, ne) with an optional count line, or one column (ne) with
    x from x_file."""
    from scipy.interpolate import UnivariateSpline
    rows = _numeric_rows(path)
    if sum(len(r) >= 2 for r in rows) >= 0.5 * len(rows):
        x, ne = np.array([r[:2] for r in rows if len(r) >= 2]).T
    else:
        if not x_file:
            raise ValueError(f"{path} has one column: give its radial coordinate in ne_x_file")
        ne = np.array([r[0] for r in rows])
        x = np.array([r[0] for r in _numeric_rows(x_file)])
        if ne.size == x.size + 1 and float(ne[0]).is_integer() and ne[0] == x.size:
            ne = ne[1:]                                   # leading point count
        if ne.size != x.size:
            raise ValueError(f"{path} ({ne.size} values) and {x_file} ({x.size}) differ")
    psi = x ** 2 if coordinate == "rho" else x
    psi, idx = np.unique(psi, return_index=True)          # sorted, no duplicates
    ne = ne[idx]
    o = np.arange(psi.size)
    spl = UnivariateSpline(psi[o], ne[o], k=k, s=0, ext=3 if outside == "constant" else 1)
    return lambda p: np.clip(spl(p), 0.0, None) * units, (psi.min(), psi.max())


# ##########################################################################
# 4. GEOMETRY: MESH FRAME <-> (R, phi, Z)
# ##########################################################################

AXIS_ORDER = {"z": (0, 1, 2), "x": (1, 2, 0), "y": (2, 0, 1)}   # (h1, h2, vertical), cyclic


def to_cylindrical(xyz, cfg, R_ref):
    """Mesh coordinates -> R, Z, phi [rad] and the toroidal arc s_phi [m]."""
    if cfg.geometry == "slab":
        R, Z = xyz[:, 0], xyz[:, 1]
        s_phi = -xyz[:, 2]                    # toroidal direction = -z
        return R, Z, s_phi / R_ref, s_phi
    a, b, c = AXIS_ORDER[cfg.tokamak_axis]
    h1, h2, Z = xyz[:, a], xyz[:, b], xyz[:, c]
    R, phi = np.hypot(h1, h2), np.arctan2(h2, h1)
    return R, Z, phi, R * phi


def b_to_mesh(BR, Bt, BZ, phi, cfg):
    """(B_R, B_t, B_Z) -> Cartesian components in the mesh frame (n, 3)."""
    out = np.empty((BR.size, 3))
    if cfg.geometry == "slab":
        out[:, 0], out[:, 1], out[:, 2] = BR, BZ, -Bt      # as gen2Dfullwavefile.py
        return out
    a, b, c = AXIS_ORDER[cfg.tokamak_axis]
    cp, sp = np.cos(phi), np.sin(phi)
    out[:, a] = BR * cp - Bt * sp
    out[:, b] = BR * sp + Bt * cp
    out[:, c] = BZ
    return out


def midplane_radius_map(top):
    """rho_pol -> outboard-midplane distance from the magnetic axis [m]."""
    Rs = np.linspace(top["R_axis"], top["R"].max(), 2000)
    rho = np.sqrt(np.clip(top["psi"].ev(Rs, np.full_like(Rs, top["Z_axis"])), 0, None))
    keep = np.concatenate([[True], np.diff(rho) > 0])
    keep &= np.maximum.accumulate(rho) == rho          # monotonic part
    rho, r = rho[keep], Rs[keep] - top["R_axis"]

    def f(x):
        x = np.asarray(x, float)
        y = np.interp(x, rho, r)
        hi = x > rho[-1]                                # linear extrapolation outside
        if hi.any():
            slope = (r[-1] - r[-2]) / (rho[-1] - rho[-2])
            y[hi] = r[-1] + slope * (x[hi] - rho[-1])
        return y
    return f


# ##########################################################################
# 5. TURBULENCE
# ##########################################################################

class Turbulence:
    """Field-aligned random-phase Fourier field F (zero mean, unit variance)."""

    def __init__(self, cfg):
        rng = np.random.default_rng(cfg.turb_seed)
        n = int(cfg.turb_modes)
        # correlation exp(-(d/L)^2)  <->  k ~ Normal(0, sqrt(2)/L)
        kr = rng.normal(0.0, np.sqrt(2) / cfg.turb_Lr, n)
        kb = rng.normal(0.0, np.sqrt(2) / cfg.turb_Lb, n)
        if cfg.turb_kb_peak > 0:                       # drift-wave-like binormal peak
            kb += np.where(rng.random(n) < 0.5, -1.0, 1.0) * cfg.turb_kb_peak
        kp = (rng.normal(0.0, np.sqrt(2) / cfg.turb_Lpar, n) if cfg.turb_Lpar > 0
              else np.zeros(n))
        self.K = np.stack([kr, kb, kp], 1)             # (n, 3)
        self.phase = rng.uniform(0, 2 * np.pi, n)
        self.amp = np.sqrt(2.0 / n)
        prof = np.asarray(cfg.turb_amplitude_profile, float)
        self.rho_tab, self.A_tab = prof[:, 0], prof[:, 1]
        self.statistics = cfg.turb_statistics

    def field(self, X):
        """X (m, 3) = (x_r, x_b, x_par) -> F (m,)."""
        return np.cos(X @ self.K.T + self.phase) @ np.full(len(self.phase), self.amp)

    def amplitude(self, rho):
        return np.interp(rho, self.rho_tab, self.A_tab)

    def apply(self, ne0, F, rho):
        A = self.amplitude(rho)
        if self.statistics == "lognormal":
            s = np.sqrt(np.log1p(A ** 2))
            ne = ne0 * np.exp(s * F - 0.5 * s ** 2)
        else:
            ne = np.clip(ne0 * (1.0 + A * F), 0.0, None)
        with np.errstate(divide="ignore", invalid="ignore"):
            dn = np.where(ne0 > 0, ne / ne0 - 1.0, 0.0)
        return ne, dn


class Blobs:
    """Field-aligned Gaussian filaments: blobs (density excess) and holes
    (relative deficit), placed at random in the field-aligned coordinates."""

    def __init__(self, cfg, X_sample, r_of_rho, n_ref):
        rng = np.random.default_rng(cfg.blob_seed)
        lo, hi = X_sample.min(axis=0), X_sample.max(axis=0)
        wr, wb, lp = cfg.blob_width_r, cfg.blob_width_b, cfg.blob_length_par
        pad = np.array([3 * wr, 3 * wb, 3 * lp])
        lo_p, hi_p = lo - pad, hi + pad                     # also blobs centred just outside
        self.r_sep = float(r_of_rho(np.array([1.0]))[0])
        self.n_ref = n_ref

        def place(rho_range, number_scale):
            r_lo, r_hi = (float(r_of_rho(np.array([x]))[0]) for x in rho_range)
            r_lo, r_hi = max(r_lo, lo_p[0]), min(r_hi, hi_p[0])
            if r_hi <= r_lo:
                return None
            area = (r_hi - r_lo) * (hi_p[1] - lo_p[1])
            n_par = max(1.0, (hi_p[2] - lo_p[2]) / (2 * lp))
            n = int(round(number_scale * area / (np.pi * wr * wb) * n_par))
            if n < 1:
                return None
            c = np.column_stack([rng.uniform(r_lo, r_hi, n), rng.uniform(lo_p[1], hi_p[1], n),
                                 rng.uniform(lo_p[2], hi_p[2], n)])
            sp = cfg.blob_size_spread
            w = np.column_stack([wr * rng.lognormal(-sp ** 2 / 2, sp, n) if sp > 0 else np.full(n, wr),
                                 wb * rng.lognormal(-sp ** 2 / 2, sp, n) if sp > 0 else np.full(n, wb),
                                 np.full(n, lp)])
            return c, w

        b = place(cfg.blob_rho_range, cfg.blob_packing_fraction)
        self.blob_c, self.blob_w = b if b else (np.zeros((0, 3)), np.zeros((0, 3)))
        sa = cfg.blob_amplitude_spread
        nb = len(self.blob_c)
        amp = cfg.blob_amplitude * (rng.lognormal(-sa ** 2 / 2, sa, nb) if sa > 0 else np.ones(nb))
        decay = np.exp(-np.clip(self.blob_c[:, 0] - self.r_sep, 0, None) / cfg.blob_decay_length)
        self.blob_dn = amp * n_ref * decay
        h = (place(cfg.hole_rho_range, cfg.blob_packing_fraction * cfg.hole_fraction)
             if cfg.hole_fraction > 0 else None)
        self.hole_c, self.hole_w = h if h else (np.zeros((0, 3)), np.zeros((0, 3)))
        nh = len(self.hole_c)
        self.hole_h = np.clip(cfg.hole_amplitude * (rng.lognormal(-sa ** 2 / 2, sa, nh)
                                                    if sa > 0 else np.ones(nh)), 0, 0.95)

    @staticmethod
    def _shape(X, c, w):
        """Gaussian of every structure at the nodes X (only nodes within 4 widths)."""
        out = []
        for ck, wk in zip(c, w):
            d = np.abs(X - ck) / wk
            m = np.nonzero((d[:, 0] < 4) & (d[:, 1] < 4) & (d[:, 2] < 4))[0]
            if m.size:
                out.append((m, np.exp(-np.sum(d[m] ** 2, axis=1))))
            else:
                out.append(None)
        return out

    def apply(self, ne, X):
        ne = ne.copy()
        for g, dnk in zip(self._shape(X, self.blob_c, self.blob_w), self.blob_dn):
            if g is not None:
                ne[g[0]] += dnk * g[1]
        for g, hk in zip(self._shape(X, self.hole_c, self.hole_w), self.hole_h):
            if g is not None:
                ne[g[0]] *= 1.0 - hk * g[1]
        return np.clip(ne, 0.0, None)

    def summary(self):
        return (f"Blobs: {len(self.blob_c)} filaments (peak excess {self.blob_dn.min():.2e}.."
                f"{self.blob_dn.max():.2e} m^-3, n_ref = {self.n_ref:.2e} m^-3), "
                f"{len(self.hole_c)} holes (peak deficit {100 * self.hole_h.min():.0f}.."
                f"{100 * self.hole_h.max():.0f} %)" if len(self.hole_c) else
                f"Blobs: {len(self.blob_c)} filaments (peak excess {self.blob_dn.min():.2e}.."
                f"{self.blob_dn.max():.2e} m^-3, n_ref = {self.n_ref:.2e} m^-3), no holes")


def field_aligned_coords(R, Z, phi, s_phi, BR, Bt, BZ, rho, top, r_of_rho, th_c, ph_c,
                         R_ref, geometry):
    """(x_r, x_b, x_par) [m] at each node: radial (flux label in metres),
    binormal r (theta - p phi) with the local field-line pitch p = d theta/d phi
    (constant along a field line), and parallel ~ R phi."""
    dR, dZ = R - top["R_axis"], Z - top["Z_axis"]
    r_geo = np.hypot(dR, dZ)
    th_abs = np.arctan2(dZ, dR)
    th = np.angle(np.exp(1j * (th_abs - th_c)))              # relative to the domain, wrapped
    ph = (np.angle(np.exp(1j * (phi - ph_c))) if geometry == "toroidal" else s_phi / R_ref)
    B_th = -BR * np.sin(th_abs) + BZ * np.cos(th_abs)        # poloidal field along theta_hat
    R_tor = R if geometry == "toroidal" else R_ref          # slab: straight toroidal axis
    with np.errstate(divide="ignore", invalid="ignore"):
        pitch = np.where(np.abs(Bt) > 0, B_th * R_tor / (r_geo * Bt), 0.0)
    return np.stack([r_of_rho(rho), r_geo * (th - pitch * ph), R_tor * ph], 1)


# ##########################################################################
# 6. EVALUATION
# ##########################################################################

def load_equilibrium(cfg):
    """(top, (ne_of_psi, psi_range)) from the files in cfg."""
    k = int(cfg.spline_order)
    if cfg.equilibrium == "eqdsk":
        top = load_eqdsk(cfg.eqdsk_file, k, cfg.eqdsk_psi_in_weber, cfg.eqdsk_flip_bp)
    else:
        top = load_topfile(cfg.topfile, k)
        top["lcfs"] = None
    return top, load_ne_profile(cfg.ne_file, cfg.ne_radial_coordinate, cfg.ne_units, k,
                                cfg.ne_x_file, getattr(cfg, "ne_outside", "constant"))


def evaluate(ids, xyz, cfg, loaded=None):
    """ne, B (and delta n / n) at every node. ``loaded`` = load_equilibrium(cfg)
    to reuse an equilibrium already read."""
    top, (ne_of_psi, psi_range) = loaded if loaded is not None else load_equilibrium(cfg)
    lcfs = None
    if top.get("lcfs") is not None and cfg.eqdsk_use_lcfs:
        from matplotlib.path import Path as MplPath
        lcfs = MplPath(top["lcfs"])
        # tolerance: the boundary polygon is coarse; points within ~1.5 grid cells
        # of it count as inside (only the private-flux region is excluded)
        lcfs_tol = 1.5 * max(np.diff(top["R"]).max(), np.diff(top["Z"]).max())
    private = 0
    Rg, Zg = top["R"], top["Z"]
    print(f"Equilibrium grid: R {Rg.min():.3f}..{Rg.max():.3f} m, Z {Zg.min():.3f}..{Zg.max():.3f} m, "
          f"magnetic axis (R, Z) = ({top['R_axis']:.4f}, {top['Z_axis']:.4f}) m")
    R_ref = float(np.median(xyz[:, 0])) if cfg.geometry == "slab" else 1.0
    if not hasattr(cfg, "chunk_nodes") or cfg.chunk_nodes < 1:
        cfg.chunk_nodes = 200_000
    turb = Turbulence(cfg) if cfg.turbulence else None
    use_blobs = bool(getattr(cfg, "blobs", False))
    use_sol = cfg.sol_ne_separatrix > 0
    r_of_rho = midplane_radius_map(top) if (turb or use_sol or use_blobs) else None
    if use_sol:
        r_sep = float(r_of_rho(np.array([1.0]))[0])
        print(f"SOL model: ne_sep = {cfg.sol_ne_separatrix:.3e} m^-3, decay length "
              f"{cfg.sol_decay_length} m, blended over rho {1 - cfg.sol_blend_width:.3f}..1")

    n = ids.size
    ne = np.empty(n)
    B = np.empty((n, 3))
    dn = np.zeros(n) if (turb or use_blobs) else None
    rho_all = np.empty(n)
    outside = 0
    # reference angles of the domain (keep theta, phi continuous inside it)
    R0, Z0, phi0, _ = to_cylindrical(xyz[:min(n, 100000)], cfg, R_ref)
    th_c = np.angle(np.mean(np.exp(1j * np.arctan2(Z0 - top["Z_axis"], R0 - top["R_axis"]))))
    ph_c = np.angle(np.mean(np.exp(1j * phi0)))

    def node_fields(xyz_part):
        """R, Z, phi, s_phi, B components, psi_n and outside-grid mask."""
        R, Z, phi, s_phi = to_cylindrical(xyz_part, cfg, R_ref)
        out = (R < Rg.min()) | (R > Rg.max()) | (Z < Zg.min()) | (Z > Zg.max())
        Rc, Zc = np.clip(R, Rg.min(), Rg.max()), np.clip(Z, Zg.min(), Zg.max())
        psi = top["psi"].ev(Rc, Zc)
        n_priv = 0
        if lcfs is not None:                    # private-flux region: outside the plasma
            cand = np.nonzero(psi < 1.0)[0]
            if cand.size:
                pts = np.column_stack([R[cand], Z[cand]])
                outp = ~(lcfs.contains_points(pts, radius=lcfs_tol) |
                         lcfs.contains_points(pts, radius=-lcfs_tol))
                psi[cand[outp]] = np.maximum(psi[cand[outp]], 1.0)
                n_priv = int(outp.sum())
        BR, Bt, BZ = top["BR"].ev(Rc, Zc), top["Bt"].ev(Rc, Zc), top["BZ"].ev(Rc, Zc)
        if cfg.flip_bt:
            Bt = -Bt
        return R, Z, phi, s_phi, BR, Bt, BZ, psi, out, n_priv

    blobs = None
    if use_blobs:
        sub = xyz if n <= 200_000 else xyz[np.random.default_rng(0).choice(n, 200_000, replace=False)]
        R, Z, phi, s_phi, BR, Bt, BZ, psi, out, _ = node_fields(sub)
        rho_s = np.sqrt(np.clip(psi, 0, None))
        Xs = field_aligned_coords(R, Z, phi, s_phi, BR, Bt, BZ, rho_s, top, r_of_rho,
                                  th_c, ph_c, R_ref, cfg.geometry)
        n_ref = cfg.blob_ne_reference or float(ne_of_psi(np.array([cfg.blob_reference_rho ** 2]))[0])
        blobs = Blobs(cfg, Xs, r_of_rho, n_ref)
        print(blobs.summary())
    t0 = time.time()
    for a in range(0, n, cfg.chunk_nodes):
        b = min(a + cfg.chunk_nodes, n)
        R, Z, phi, s_phi, BR, Bt, BZ, psi, out, n_priv = node_fields(xyz[a:b])
        outside += int(out.sum())
        private += n_priv
        n0 = ne_of_psi(psi)
        if use_sol:
            rr = np.sqrt(np.clip(psi, 0, None))
            sol = cfg.sol_ne_separatrix * np.exp(-np.clip(r_of_rho(rr) - r_sep, 0, None)
                                                 / cfg.sol_decay_length)
            w = np.clip((rr - (1.0 - cfg.sol_blend_width)) / max(cfg.sol_blend_width, 1e-9), 0, 1)
            n0 = np.where(rr >= 1.0, sol, (1 - w) * n0 + w * cfg.sol_ne_separatrix)
        n0[out] = 0.0
        B[a:b] = b_to_mesh(BR, Bt, BZ, phi, cfg)
        rho = np.sqrt(np.clip(psi, 0, None))
        rho_all[a:b] = rho
        if turb or blobs:
            X = field_aligned_coords(R, Z, phi, s_phi, BR, Bt, BZ, rho, top, r_of_rho,
                                     th_c, ph_c, R_ref, cfg.geometry)
            nbg = n0.copy()                     # background profile
            if turb:
                n0, _ = turb.apply(n0, turb.field(X), rho)
            if blobs:
                n0 = blobs.apply(n0, X)
                n0[out] = 0.0
            # delta n / n relative to the background (floored at 1 % of n_ref with blobs,
            # which also put density where the background profile has none)
            den = np.maximum(nbg, 0.01 * blobs.n_ref if blobs else 0.0)
            ok = den > 0
            dn[a:b] = 0.0
            dn[a:b][ok] = (n0[ok] - nbg[ok]) / den[ok]
        ne[a:b] = n0
        if b == n or (a // cfg.chunk_nodes) % 20 == 19:
            print(f"  {b:,} / {n:,} nodes ({time.time() - t0:.1f} s)")
    if private:
        print(f"NOTE: {private:,} nodes have psi_n < 1 outside the LCFS (private-flux region): "
              "treated as outside the plasma (psi_n set to 1).")
    if outside:
        print(f"WARNING: {outside:,} nodes are outside the equilibrium grid: B there is the value "
              "at the nearest grid edge and ne = 0.")
    return ne, B, dn, rho_all, psi_range


# ##########################################################################
# 7. OUTPUT
# ##########################################################################

def write_columns(path, ids, cols, fmt="%.8e"):
    """'id v1 v2 ...' per line, fast for tens of millions of rows."""
    t0 = time.time()
    data = np.asarray(cols, float)
    data = data[:, None] if data.ndim == 1 else data
    try:
        import pandas as pd
        df = pd.DataFrame(data)
        df.insert(0, "id", ids)
        df.to_csv(path, sep=" ", header=False, index=False, float_format=fmt,
                  chunksize=1_000_000)
    except ImportError:
        with open(path, "w") as f:
            for a in range(0, len(ids), 1_000_000):
                np.savetxt(f, np.column_stack([ids[a:a + 1_000_000], data[a:a + 1_000_000]]),
                           fmt=["%d"] + [fmt] * data.shape[1])
    print(f"Wrote {path} ({len(ids):,} rows, {time.time() - t0:.1f} s)")


def check_plot(path, xyz, ne, B, dn, rho, cfg):
    """Cut through the mesh (ne, delta n / n, |B|) and the profiles."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    if cfg.geometry == "slab":
        c = xyz[:, 2]
        h1, h2, lab = xyz[:, 0], xyz[:, 1], ("x = R [m]", "y = Z [m]")
    else:
        a, b, cax = AXIS_ORDER[cfg.tokamak_axis]
        phi = np.arctan2(xyz[:, b], xyz[:, a])
        c = np.angle(np.exp(1j * (phi - np.angle(np.mean(np.exp(1j * phi))))))
        h1, h2, lab = np.hypot(xyz[:, a], xyz[:, b]), xyz[:, cax], ("R [m]", "Z [m]")
    mid = np.median(c)
    tol = 0.02 * (np.ptp(c) or 1.0)
    sel = np.nonzero(np.abs(c - mid) <= tol)[0]
    if sel.size > 300_000:
        sel = rng.choice(sel, 300_000, replace=False)
    panels = [("ne [m^-3]", ne), ("|B| [T]", np.linalg.norm(B, axis=1))]
    if dn is not None:
        panels.insert(1, ("delta n / n", dn))
    fig, axs = plt.subplots(1, len(panels) + 1, figsize=(5 * (len(panels) + 1), 4.6))
    for ax, (name, v) in zip(axs, panels):
        vmax = np.percentile(np.abs(v[sel]), 99.5) if name == "delta n / n" else None
        sc = ax.scatter(h1[sel], h2[sel], c=v[sel], s=1, cmap="RdBu_r" if vmax else "viridis",
                        vmin=-vmax if vmax else None, vmax=vmax)
        fig.colorbar(sc, ax=ax, shrink=0.85)
        ax.set_title(f"{name} (cut at the mid {'z' if cfg.geometry == 'slab' else 'phi'})")
        ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1]); ax.set_aspect("equal")
    ax = axs[-1]
    s2 = rng.choice(ne.size, min(ne.size, 200_000), replace=False)
    ax.plot(rho[s2], ne[s2], ",", alpha=0.4, label="ne at nodes")
    ax.set_xlabel("rho_pol"); ax.set_ylabel("ne [m^-3]")
    if dn is not None:
        ax2 = ax.twinx()
        bins = np.linspace(rho.min(), rho.max(), 40)
        idx = np.digitize(rho, bins)
        rms = [np.sqrt(np.mean(dn[idx == i] ** 2)) if np.any(idx == i) else np.nan
               for i in range(1, len(bins))]
        ax2.plot(0.5 * (bins[1:] + bins[:-1]), rms, "r-", label="rms delta n / n (nodes)")
        if cfg.turbulence:
            ax2.plot(bins, np.interp(bins, *np.asarray(cfg.turb_amplitude_profile).T), "k--",
                     label="target A(rho), turbulence")
        ax2.set_ylabel("rms delta n / n", color="r")
        ax2.legend(loc="upper right", fontsize=8)
    ax.set_title("profiles")
    plt.tight_layout()
    plt.savefig(path, dpi=110)
    plt.close(fig)
    print(f"Wrote {path}")


# ##########################################################################
# 8. COMMAND LINE / MAIN
# ##########################################################################

def default_parameters():
    g = globals()
    return argparse.Namespace(**{n: g[n] for n in USER_PARAMETERS})


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="ne and B files for ERMES from a GiD mesh, "
                                "a topfile and ne(rho), with optional turbulence.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    for name in USER_PARAMETERS:
        dflt = globals()[name]
        flag = "--" + name.replace("_", "-")
        if isinstance(dflt, bool):
            p.add_argument(flag, action=argparse.BooleanOptionalAction, default=dflt)
        elif name in PARAMETER_NARGS:
            p.add_argument(flag, type=float, nargs=PARAMETER_NARGS[name], default=dflt)
        elif name == "turb_amplitude_profile":
            p.add_argument(flag, type=float, nargs="+", default=dflt,
                           help="rho A rho A ... pairs")
        else:
            p.add_argument(flag, type=type(dflt), default=dflt,
                           choices=PARAMETER_CHOICES.get(name))
    p.add_argument("--show-params", action="store_true", help="print the parameters and exit")
    cfg = p.parse_args(argv)
    prof = cfg.turb_amplitude_profile
    if prof and not isinstance(prof[0], (tuple, list)):
        if len(prof) % 2:
            p.error("--turb-amplitude-profile needs rho A pairs")
        cfg.turb_amplitude_profile = [(prof[i], prof[i + 1]) for i in range(0, len(prof), 2)]
    return cfg


def main(cfg=None):
    cfg = cfg or parse_args()
    if getattr(cfg, "show_params", False):
        for n in USER_PARAMETERS:
            print(f"  {n:28s} = {getattr(cfg, n)!r}")
        return
    ids, xyz = read_mesh_nodes(cfg.mesh_file)
    xyz *= UNIT_SCALE[cfg.mesh_units]
    order = np.argsort(ids, kind="stable")
    if np.any(order != np.arange(ids.size)):
        ids, xyz = ids[order], xyz[order]
    if np.any(np.diff(ids) == 0):
        print("WARNING: repeated node IDs in the mesh file.")
    print(f"Geometry: {cfg.geometry}" + (f" (tokamak axis = {cfg.tokamak_axis})"
                                         if cfg.geometry == "toroidal" else " (x = R, y = Z, toroidal = -z)"))
    if getattr(cfg, "blobs", False):
        print(f"Blobs: packing {cfg.blob_packing_fraction}, rho {cfg.blob_rho_range}, widths "
              f"{cfg.blob_width_r} x {cfg.blob_width_b} m, L_par {cfg.blob_length_par} m, "
              f"amplitude {cfg.blob_amplitude} n_ref, holes {cfg.hole_fraction}/blob, seed {cfg.blob_seed}")
    if cfg.turbulence:
        print(f"Turbulence: {cfg.turb_statistics}, L_r = {cfg.turb_Lr} m, L_b = {cfg.turb_Lb} m, "
              f"L_par = {cfg.turb_Lpar} m, k_b0 = {cfg.turb_kb_peak} 1/m, {cfg.turb_modes} modes, "
              f"seed {cfg.turb_seed}")
    ne, B, dn, rho, psi_range = evaluate(ids, xyz, cfg)
    print(f"ne: max {ne.max():.4e} m^-3, {np.count_nonzero(ne):,} nodes with ne > 0; "
          f"rho_pol at nodes {rho.min():.3f}..{rho.max():.3f} "
          f"(profile covers psi_n {psi_range[0]:.3f}..{psi_range[1]:.3f})")
    if dn is not None:
        m = ne > 0
        print(f"delta n / n at nodes with ne > 0: rms {np.sqrt(np.mean(dn[m] ** 2)):.4f}, "
              f"min {dn[m].min():+.3f}, max {dn[m].max():+.3f}")
    os.makedirs(cfg.out_dir, exist_ok=True)
    out = lambda name: os.path.join(cfg.out_dir, name)
    write_columns(out(cfg.ne_out), ids, ne)
    write_columns(out(cfg.mag_out), ids, B)
    if dn is not None and cfg.write_turbulence_file:
        write_columns(out(cfg.dn_out), ids, dn)
    if cfg.plot:
        try:
            check_plot(out("plasma_check.png"), xyz, ne, B, dn, rho, cfg)
        except Exception as err:                 # never lose the data files
            print(f"WARNING: plot failed: {err}")


if __name__ == "__main__":
    main()
