#!/usr/bin/env python3
"""
gid_surface_readers.py
======================

Reads the geometry / mesh files GiD can export (Files > Export) and returns
the SURFACES as a triangle soup in metres, for reflect_beam_to_gid.py (PEC
reflecting surfaces). Can also be run on its own to convert a file to STL and
check in GiD what will be used as reflector:

    python gid_surface_readers.py Inner.igs -o check.stl --mesh-size 0.01
    python gid_surface_readers.py Inner.igs --layers Inner_Tube -o tube.stl
    python gid_surface_readers.py model.msh --info          # layers / blocks only

Formats (GiD export menu -> what is read)
------------------------------------------
  IGES            .igs .iges   NATIVE (numpy only): trimmed NURBS surfaces
                               (128 + 144/142 or 143/141), GiD layers (406),
                               units from the file. gmsh used instead if
                               reader="gmsh".
  STEP            .stp .step   gmsh (pip install gmsh); also .brep
  DXF             .dxf         NATIVE: 3DFACE, SOLID, polyface / polygon-mesh
                               POLYLINE, MESH; layers (group code 8)
  ACIS            .sat .sab    no free reader exists -> export IGES/STEP instead
  Rhinoceros      .3dm         rhino3dm (pip install rhino3dm): meshes, render
                               meshes of breps, else untrimmed NURBS faces
  VTK mesh        .vtk         NATIVE: legacy ASCII / BINARY, POLYDATA or
                               UNSTRUCTURED_GRID (volume cells -> their skin)
                  .vtu .vtp    meshio (pip install meshio)
  GiD mesh        .msh (GiD)   NATIVE: all element types, MESH blocks as layers
                  .msh (gmsh)  gmsh or meshio
  Meshio          any          meshio; .obj .off .ply also NATIVE
  Exodus II       .exo .e .g   NATIVE via scipy (netCDF3) or h5py (netCDF4);
                               element blocks as layers
  DumpHDF5        .h5 .hdf5    h5py, best effort (coordinates + connectivity
                               found by shape); prints the file tree otherwise
  Surfaces Json   .json        best effort (vertex + face arrays found by name)
  STL             .stl         NATIVE (ASCII / binary)
Volume meshes (tetra, hexa, prism, pyramid) are reduced to their boundary
faces. Not geometry: Nodes->Raster, Text data report, Calculation file.

Units: IGES (and STEP through gmsh) carry their unit and are converted to m
automatically; for all mesh formats give ``units`` ("m", "cm", "mm", ...).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import struct
import sys
from collections import Counter, defaultdict

import numpy as np

UNIT_TO_M = {"m": 1.0, "cm": 0.01, "mm": 1e-3, "um": 1e-6, "km": 1e3,
             "in": 0.0254, "ft": 0.3048}


class SurfaceReadError(Exception):
    """Raised with a message meant for the user."""


def _say(verbose, *a):
    if verbose:
        print(*a)


# ##########################################################################
# Element / cell helpers
# ##########################################################################

_FACES = {
    "tri": [(0, 1, 2)],
    "quad": [(0, 1, 2, 3)],
    "tet": [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)],
    "hex": [(0, 3, 2, 1), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)],
    "wedge": [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "pyr": [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
_NCORNER = {"tri": 3, "quad": 4, "tet": 4, "hex": 8, "wedge": 6, "pyr": 5}


def _split_faces(A):
    A = np.asarray(A, np.int64)
    if A.shape[1] == 3:
        return [A]
    return [A[:, [0, 1, 2]], A[:, [0, 2, 3]]]


def cells_to_triangles(X, blocks):
    """X (n, 3) points; blocks: list of (kind, conn) with kind in tri, quad,
    tet, hex, wedge, pyr (conn (m, k) 0-based, extra (high-order) nodes
    ignored) or "poly" (list of index lists). Volume cells give their boundary
    faces. Returns (nt, 3, 3)."""
    X = np.asarray(X, float)
    surf, vol, tris = [], defaultdict(list), []
    for kind, conn in blocks:
        if kind == "poly":
            for poly in conn:
                poly = list(poly)
                tris += [[poly[0], poly[i], poly[i + 1]] for i in range(1, len(poly) - 1)]
            continue
        if kind not in _NCORNER:
            continue
        conn = np.asarray(conn, np.int64)
        if conn.ndim != 2 or not len(conn) or conn.shape[1] < _NCORNER[kind]:
            continue
        conn = conn[:, :_NCORNER[kind]]
        for f in _FACES[kind]:
            if kind in ("tri", "quad"):
                surf.append(conn[:, f])
            else:
                vol[len(f)].append(conn[:, f])
    out = [np.asarray(tris, np.int64).reshape(-1, 3)]
    for A in surf:
        out += _split_faces(A)
    for k, arrs in vol.items():
        A = np.vstack(arrs)
        _, inv, cnt = np.unique(np.sort(A, 1), axis=0, return_inverse=True, return_counts=True)
        out += _split_faces(A[cnt[inv.ravel()] == 1])
    T = np.vstack(out)
    if not len(T):
        return np.zeros((0, 3, 3))
    if T.max() >= len(X) or T.min() < 0:
        raise SurfaceReadError("element connectivity refers to nodes that do not exist")
    return X[T]


def _layer_ok(name, layers):
    if not layers:
        return True
    names = {str(l).strip().lower() for l in layers}
    return str(name).strip().lower() in names


def write_stl(path, tris, name="surfaces"):
    """Binary STL of a (nt, 3, 3) array."""
    tris = np.asarray(tris, np.float32)
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(tris), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)),
                                              ("a", "<u2")]))
    rec["n"], rec["v"] = n, tris
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        f.write(rec.tobytes())


# ##########################################################################
# STL
# ##########################################################################

def read_stl(path):
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
        raise SurfaceReadError(f"could not read STL file {path}")
    return np.array(v, float).reshape(-1, 3, 3)


# ##########################################################################
# NURBS evaluation
# ##########################################################################

def bspline_basis(U, p, t):
    """B-spline basis functions (len(t), len(U)-p-1) of degree p at t."""
    U = np.asarray(U, float)
    t = np.clip(np.atleast_1d(np.asarray(t, float)), U[p], U[-p - 1])
    nk = len(U)
    k = np.clip(np.searchsorted(U, t, side="right") - 1, p, nk - p - 2)
    N = np.zeros((len(t), nk - 1))
    N[np.arange(len(t)), k] = 1.0
    tt = t[:, None]
    for d in range(1, p + 1):
        i = np.arange(nk - 1 - d)
        den1, den2 = U[i + d] - U[i], U[i + d + 1] - U[i + 1]
        a = np.where(den1 > 0, (tt - U[i]) / np.where(den1 > 0, den1, 1), 0.0)
        b = np.where(den2 > 0, (U[i + d + 1] - tt) / np.where(den2 > 0, den2, 1), 0.0)
        N = a * N[:, i] + b * N[:, i + 1]
    return N


class NurbsCurve:
    def __init__(self, U, p, W, P, t0, t1):
        self.U, self.p = np.asarray(U, float), int(p)
        self.W, self.P = np.asarray(W, float), np.asarray(P, float)
        self.t0, self.t1 = float(t0), float(t1)

    def __call__(self, t):
        N = bspline_basis(self.U, self.p, t) * self.W
        return (N @ self.P) / N.sum(1, keepdims=True)

    def spans(self):
        kn = self.U[(self.U > self.t0) & (self.U < self.t1)]
        return np.unique(np.concatenate([[self.t0], kn, [self.t1]]))


class NurbsSurface:
    def __init__(self, U, V, p, q, W, P, dom):
        self.U, self.V = np.asarray(U, float), np.asarray(V, float)
        self.p, self.q = int(p), int(q)
        self.W = np.asarray(W, float)                      # (nu, nv)
        self.WP = self.W[..., None] * np.asarray(P, float)  # (nu, nv, 3)
        self.P = np.asarray(P, float)
        self.dom = np.asarray(dom, float)                  # u0, u1, v0, v1

    def __call__(self, u, v):
        u, v = np.atleast_1d(u), np.atleast_1d(v)
        out = np.empty((len(u), 3))
        nu, nv = self.W.shape
        WPf = self.WP.reshape(nu, nv * 3)
        for c0 in range(0, len(u), 4000):
            c1 = min(c0 + 4000, len(u))
            Nu = bspline_basis(self.U, self.p, u[c0:c1])
            Nv = bspline_basis(self.V, self.q, v[c0:c1])
            num = np.einsum("mj,mjc->mc", Nv, (Nu @ WPf).reshape(-1, nv, 3))
            den = np.einsum("mj,mj->m", Nv, Nu @ self.W)
            out[c0:c1] = num / den[:, None]
        return out

    def grid(self, us, vs):
        Nu, Nv = bspline_basis(self.U, self.p, us), bspline_basis(self.V, self.q, vs)
        num = np.einsum("ai,bj,ijc->abc", Nu, Nv, self.WP)
        den = np.einsum("ai,bj,ij->ab", Nu, Nv, self.W)
        return num / den[..., None]

    def knots_in_domain(self, which):
        U, a, b = (self.U, *self.dom[:2]) if which == 0 else (self.V, *self.dom[2:])
        return np.unique(np.concatenate([[a], U[(U > a) & (U < b)], [b]]))

    def project(self, X, g=41):
        """(u, v) of the points of the surface closest to X (m, 3)."""
        u0, u1, v0, v1 = self.dom
        us, vs = np.linspace(u0, u1, g), np.linspace(v0, v1, g)
        G = self.grid(us, vs).reshape(-1, 3)
        UV = np.stack(np.meshgrid(us, vs, indexing="ij"), -1).reshape(-1, 2)
        idx = np.empty(len(X), int)
        for c0 in range(0, len(X), 500):
            d = ((X[c0:c0 + 500, None, :] - G[None]) ** 2).sum(-1)
            idx[c0:c0 + 500] = d.argmin(1)
        uv = UV[idx].copy()
        hu, hv = 1e-6 * max(u1 - u0, 1e-12), 1e-6 * max(v1 - v0, 1e-12)
        for _ in range(15):
            S = self(uv[:, 0], uv[:, 1])
            su = np.where(uv[:, 0] + hu <= u1, hu, -hu)
            sv = np.where(uv[:, 1] + hv <= v1, hv, -hv)
            Su = (self(uv[:, 0] + su, uv[:, 1]) - S) / su[:, None]
            Sv = (self(uv[:, 0], uv[:, 1] + sv) - S) / sv[:, None]
            r = X - S
            a, b, c = (Su * Su).sum(1), (Su * Sv).sum(1), (Sv * Sv).sum(1)
            ru, rv = (Su * r).sum(1), (Sv * r).sum(1)
            det = a * c - b * b
            ok = np.abs(det) > 1e-30
            det = np.where(ok, det, 1.0)
            du = np.where(ok, (c * ru - b * rv) / det, 0.0)
            dv = np.where(ok, (a * rv - b * ru) / det, 0.0)
            uv[:, 0] = np.clip(uv[:, 0] + du, u0, u1)
            uv[:, 1] = np.clip(uv[:, 1] + dv, v0, v1)
        return uv


def _points_in_loops(pts, loops):
    """Even-odd test of pts (m, 2) against closed polygons (list of (k, 2))."""
    A = np.vstack([l for l in loops])
    B = np.vstack([np.roll(l, -1, axis=0) for l in loops])
    inside = np.zeros(len(pts), bool)
    step = max(1, int(3e6 // max(len(A), 1)))
    for c0 in range(0, len(pts), step):
        P = pts[c0:c0 + step]
        x, y = P[:, :1], P[:, 1:2]
        cond = (A[:, 1] > y) != (B[:, 1] > y)
        dy = np.where(B[:, 1] != A[:, 1], B[:, 1] - A[:, 1], 1.0)
        xi = A[:, 0] + (y - A[:, 1]) * (B[:, 0] - A[:, 0]) / dy
        inside[c0:c0 + step] = (cond & (x < xi)).sum(1) % 2 == 1
    return inside


_WARNED = set()


def triangulate_trimmed(srf, loops, h, max_n=600):
    """Triangles (nt, 3, 3) of the trimmed NURBS surface. ``loops``: closed
    (u, v) polygons (even-odd rule), already sampled at ~h in 3D."""
    u0, u1, v0, v1 = srf.dom
    cg = srf.grid(np.linspace(u0, u1, 17), np.linspace(v0, v1, 17))
    Lu = np.linalg.norm(np.diff(cg, axis=0), axis=-1).sum(0).max()
    Lv = np.linalg.norm(np.diff(cg, axis=1), axis=-1).sum(1).max()

    def params(which, L):
        a, b = (u0, u1) if which == 0 else (v0, v1)
        deg = srf.p if which == 0 else srf.q
        kn = srf.knots_in_domain(which)
        n = int(np.clip(math.ceil(L / h), 1, max_n))
        if deg > 1:
            n = max(n, 4 * (len(kn) - 1))
        t = np.unique(np.concatenate([np.linspace(a, b, n + 1), kn]))
        keep = np.concatenate([[True], np.diff(t) > 1e-9 * max(b - a, 1e-12)])
        return t[keep]

    us, vs = params(0, Lu), params(1, Lv)
    su = Lu / max(u1 - u0, 1e-300) if Lu > 0 else 1.0
    sv = Lv / max(v1 - v0, 1e-300) if Lv > 0 else 1.0
    scl = np.array([su, sv])
    bnd = np.vstack(loops)
    GU, GV = np.meshgrid(us, vs, indexing="ij")
    grid_uv = np.stack([GU.ravel(), GV.ravel()], 1)
    inside = _points_in_loops(grid_uv, loops)
    try:
        from scipy.spatial import Delaunay, cKDTree
    except ImportError:
        Delaunay = None
    if Delaunay is not None:
        cell = min(np.diff(us).min() * su, np.diff(vs).min() * sv, h)
        pin = grid_uv[inside]
        if len(pin):
            d, _ = cKDTree(bnd * scl).query(pin * scl)
            pin = pin[d > 0.4 * cell]
        pts = np.vstack([bnd, pin])
        q = np.round(pts * scl / (1e-9 * max(Lu, Lv, 1e-12))).astype(np.int64)
        _, first = np.unique(q, axis=0, return_index=True)
        pts = pts[np.sort(first)]
        if len(pts) < 3:
            return np.zeros((0, 3, 3))
        try:
            tri = Delaunay(pts * scl).simplices
        except Exception:
            tri = Delaunay(pts * scl, qhull_options="QJ").simplices
        cen = pts[tri].mean(1)
        tri = tri[_points_in_loops(cen, loops)]
        X = srf(pts[:, 0], pts[:, 1])
        return X[tri]
    if "scipy" not in _WARNED:
        _WARNED.add("scipy")
        print("NOTE: scipy not found: trimmed surfaces are meshed on a plain (u, v) grid, "
              "edges are staircase-like (pip install scipy for exact trimmed edges).")
    nu, nv = len(us), len(vs)
    X = srf.grid(us, vs).reshape(-1, 3)
    idx = np.arange(nu * nv).reshape(nu, nv)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    T = np.vstack([np.stack([a, b, c], 1), np.stack([a, c, d], 1)])
    cen = grid_uv[T].mean(1)
    return X[T[_points_in_loops(cen, loops)]]


# ##########################################################################
# IGES (native)
# ##########################################################################

_IGES_UNITS = {1: 0.0254, 2: 1e-3, 4: 0.3048, 5: 1609.344, 6: 1.0, 7: 1000.0,
               8: 2.54e-5, 9: 1e-6, 10: 0.01, 11: 2.54e-8}
_IGES_UNIT_NAMES = {"IN": 0.0254, "INCH": 0.0254, "MM": 1e-3, "FT": 0.3048, "MI": 1609.344,
                    "M": 1.0, "KM": 1e3, "MIL": 2.54e-5, "UM": 1e-6, "CM": 0.01, "UIN": 2.54e-8}
_HOLL = re.compile(r"\s*(\d+)H")


def _iges_tokenize(s, delim, rdelim):
    toks, i, n = [], 0, len(s)
    while i < n:
        m = _HOLL.match(s, i)
        if m:
            k, j = int(m.group(1)), m.end()
            toks.append(s[j:j + k])
            i = j + k
            while i < n and s[i] not in (delim, rdelim):
                i += 1
        else:
            j = i
            while j < n and s[j] not in (delim, rdelim):
                j += 1
            toks.append(s[i:j].strip())
            i = j
        if i >= n or s[i] == rdelim:
            break
        i += 1
    return toks


def _f(t):
    t = str(t).strip()
    return float(t.replace("D", "E").replace("d", "e")) if t else 0.0


def _i(t):
    return int(round(_f(t)))


class _IgesFile:
    def __init__(self, path):
        with open(path, "r", encoding="latin-1", errors="replace") as f:
            lines = [l.rstrip("\r\n").ljust(80) for l in f]
        sec = defaultdict(list)
        for l in lines:
            sec[l[72]].append(l)
        if sec.get("C") or (not sec.get("D") and lines and lines[0][72] == "C"):
            raise SurfaceReadError(f"{path}: compressed IGES is not supported")
        if not sec.get("D"):
            raise SurfaceReadError(f"{path}: not an IGES file (no directory section)")
        g = "".join(l[:72] for l in sec["G"])
        delim, rdelim = ",", ";"
        rest = g
        if g.startswith("1H"):
            delim, rest = g[2], g[3:]
        rest = rest[1:] if rest[:1] == delim else rest
        if rest.startswith("1H"):
            rdelim = rest[2]
        self.delim, self.rdelim = delim, rdelim
        G = _iges_tokenize(g, delim, rdelim)
        self.scale = 1.0
        try:
            flag = _i(G[13]) if len(G) > 13 and G[13].strip() else 1
            uname = G[14].strip().upper() if len(G) > 14 else ""
            self.scale = _IGES_UNITS.get(flag) or _IGES_UNIT_NAMES.get(uname, 1.0)
            if flag == 3 and uname in _IGES_UNIT_NAMES:
                self.scale = _IGES_UNIT_NAMES[uname]
            model_scale = _f(G[12]) if len(G) > 12 and G[12].strip() else 1.0
            if model_scale > 0:
                self.scale /= model_scale
            self.unit_name = uname or f"flag {flag}"
        except (ValueError, IndexError):
            self.unit_name = "unknown (taken as m)"
        self.ents = {}
        D = sec["D"]
        for k in range(0, len(D) - 1, 2):
            l1, l2 = D[k], D[k + 1]
            f1 = [l1[i * 8:(i + 1) * 8] for i in range(9)]
            f2 = [l2[i * 8:(i + 1) * 8] for i in range(9)]
            try:
                de = int(l1[73:80])
            except ValueError:
                de = k + 1

            def iv(s):
                s = s.strip()
                try:
                    return int(s) if s else 0
                except ValueError:
                    return 0
            self.ents[de] = dict(type=iv(f1[0]), pptr=iv(f1[1]), level=iv(f1[4]),
                                 xform=iv(f1[6]), status=f1[8].strip().zfill(8),
                                 color=iv(f2[2]), form=iv(f2[4]), de=de)
        self._plines = defaultdict(list)
        for l in sec["P"]:
            try:
                self._plines[int(l[64:72])].append(l[:64])
            except ValueError:
                pass
        self._cache = {}
        self.levels = {}
        for de, e in self.ents.items():
            if e["type"] == 406 and e["form"] == 3:
                t = self.p(de)
                if len(t) >= 3:
                    self.levels[_i(t[1])] = t[2]

    def p(self, de):
        if de not in self._cache:
            toks = _iges_tokenize("".join(self._plines.get(de, [])), self.delim, self.rdelim)
            self._cache[de] = toks[1:]
        return self._cache[de]

    def layer_names(self, de):
        lv = self.ents[de]["level"]
        if lv > 0:
            return [self.levels.get(lv, str(lv))]
        if lv < 0 and -lv in self.ents:                    # 406 form 1: several levels
            t = self.p(-lv)
            return [self.levels.get(_i(x), str(_i(x))) for x in t[1:1 + _i(t[0])]]
        return ["0"]

    def transform(self, de, P):
        """Apply the entity's transformation matrix chain to P (k, 3)."""
        x = self.ents.get(de, {}).get("xform", 0)
        seen = set()
        while x and x in self.ents and x not in seen:
            seen.add(x)
            t = [_f(v) for v in self.p(x)[:12]]
            if len(t) == 12:
                R = np.array([t[0:3], t[4:7], t[8:11]])
                T = np.array([t[3], t[7], t[11]])
                P = P @ R.T + T
            x = self.ents[x]["xform"]
        return P

    # ---------------- entities ----------------
    def surface(self, de):
        e = self.ents[de]
        if e["type"] != 128:
            raise NotImplementedError(f"surface entity {e['type']}")
        t = self.p(de)
        K1, K2, M1, M2 = (_i(x) for x in t[:4])
        A, B = 1 + K1 - M1 + 2 * M1, 1 + K2 - M2 + 2 * M2
        i = 9
        U = [_f(x) for x in t[i:i + A + 1]]; i += A + 1
        V = [_f(x) for x in t[i:i + B + 1]]; i += B + 1
        nw = (K1 + 1) * (K2 + 1)
        W = np.array([_f(x) for x in t[i:i + nw]]).reshape(K2 + 1, K1 + 1).T; i += nw
        P = np.array([_f(x) for x in t[i:i + 3 * nw]]).reshape(K2 + 1, K1 + 1, 3)
        P = P.transpose(1, 0, 2); i += 3 * nw
        dom = [_f(x) for x in t[i:i + 4]]
        if len(dom) < 4:
            dom = [U[M1], U[-M1 - 1], V[M2], V[-M2 - 1]]
        P = self.transform(de, P.reshape(-1, 3)).reshape(P.shape) * self.scale
        if np.any(W <= 0):
            W = np.ones_like(W)
        return NurbsSurface(U, V, M1, M2, W, P, dom)

    def curve_pieces(self, de, model):
        """Dense polylines (list of (k, 3)) of a curve; corners at piece ends.
        model=True: model space, converted to m; False: parameter space."""
        e = self.ents[de]
        typ, t = e["type"], self.p(de)
        pieces = []
        if typ == 102:
            for c in t[1:1 + _i(t[0])]:
                pieces += self.curve_pieces(_i(c), model)
            return pieces                                   # children transformed already
        if typ == 110:
            a, b = np.array([_f(x) for x in t[0:3]]), np.array([_f(x) for x in t[3:6]])
            pieces = [np.linspace(a, b, 2 if model else 33)]
        elif typ == 100:
            zt, x1, y1, x2, y2, x3, y3 = (_f(x) for x in t[:7])
            c = np.array([x1, y1])
            a0 = math.atan2(y2 - y1, x2 - x1)
            a1 = math.atan2(y3 - y1, x3 - x1)
            if a1 <= a0 + 1e-12:
                a1 += 2 * math.pi
            r = math.hypot(x2 - x1, y2 - y1)
            ang = np.linspace(a0, a1, max(9, int(64 * (a1 - a0) / math.pi)))
            pieces = [np.stack([c[0] + r * np.cos(ang), c[1] + r * np.sin(ang),
                                np.full_like(ang, zt)], 1)]
        elif typ == 126:
            K, M = _i(t[0]), _i(t[1])
            i = 6
            U = [_f(x) for x in t[i:i + K + M + 2]]; i += K + M + 2
            W = np.array([_f(x) for x in t[i:i + K + 1]]); i += K + 1
            P = np.array([_f(x) for x in t[i:i + 3 * (K + 1)]]).reshape(-1, 3); i += 3 * (K + 1)
            t0, t1 = _f(t[i]), _f(t[i + 1])
            if np.any(W <= 0):
                W = np.ones_like(W)
            crv = NurbsCurve(U, M, W, P, t0, t1)
            ks = crv.spans()
            ns = 2 if (model and M == 1) else (33 if M == 1 else 65)
            pieces = [crv(np.linspace(a, b, ns)) for a, b in zip(ks[:-1], ks[1:]) if b > a]
        elif typ == 106:
            ip, n = _i(t[0]), _i(t[1])
            if e["form"] in (1, 11, 63) or ip == 1:
                zt = _f(t[2])
                xy = np.array([_f(x) for x in t[3:3 + 2 * n]]).reshape(-1, 2)
                P = np.column_stack([xy, np.full(len(xy), zt)])
            else:
                P = np.array([_f(x) for x in t[2:2 + 3 * n]]).reshape(-1, 3)
            if e["form"] == 63:
                P = np.vstack([P, P[:1]])
            pieces = [np.linspace(a, b, 2 if model else 9) for a, b in zip(P[:-1], P[1:])]
        elif typ == 112:
            n = _i(t[3])
            T = [_f(x) for x in t[4:5 + n]]
            co = np.array([_f(x) for x in t[5 + n:5 + n + 12 * n]]).reshape(n, 3, 4)
            for k in range(n):
                s = np.linspace(0, T[k + 1] - T[k], 33)[:, None]
                c = co[k]
                pieces.append(np.stack([c[d, 0] + c[d, 1] * s[:, 0] + c[d, 2] * s[:, 0] ** 2
                                        + c[d, 3] * s[:, 0] ** 3 for d in range(3)], 1))
        else:
            raise NotImplementedError(f"curve entity {typ}")
        pieces = [self.transform(de, p) for p in pieces]
        return [p * self.scale for p in pieces] if model else pieces


def _chain(pieces):
    """Order / flip pieces so that each one starts where the previous ended."""
    pieces = [p for p in pieces if len(p) >= 2]
    if not pieces:
        return []
    out, rest = [pieces[0]], pieces[1:]
    while rest:
        end = out[-1][-1]
        ds = [min(np.linalg.norm(p[0] - end), np.linalg.norm(p[-1] - end)) for p in rest]
        k = int(np.argmin(ds))
        p = rest.pop(k)
        out.append(p if np.linalg.norm(p[0] - end) <= np.linalg.norm(p[-1] - end) else p[::-1])
    return out


def _resample(dense3, h):
    """Indices / fractions to resample a dense polyline at ~h (3D lengths)."""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(dense3, axis=0), axis=1))])
    n = max(1, int(math.ceil(s[-1] / h)))
    return s, np.linspace(0.0, s[-1], n + 1)


def _loop_from_uv_pieces(srf, pieces, h):
    out = []
    for p in _chain([q[:, :2] for q in pieces]):
        X = srf(p[:, 0], p[:, 1])
        s, st = _resample(X, h)
        if s[-1] <= 0:
            continue
        uv = np.stack([np.interp(st, s, p[:, 0]), np.interp(st, s, p[:, 1])], 1)
        out.append(uv if not out else uv[1:])
    return _close(out)


def _loop_from_model_pieces(srf, pieces, h):
    out = []
    for p in _chain(pieces):
        s, st = _resample(p, h)
        if s[-1] <= 0:
            continue
        X = np.stack([np.interp(st, s, p[:, d]) for d in range(3)], 1)
        out.append(srf.project(X) if not out else srf.project(X)[1:])
    return _close(out)


def _close(parts):
    if not parts:
        return None
    L = np.vstack(parts)
    if len(L) > 1 and np.linalg.norm(L[0] - L[-1]) < 1e-12 * (1 + np.abs(L).max()):
        L = L[:-1]
    return L if len(L) >= 3 else None


def _domain_loop(srf, h):
    u0, u1, v0, v1 = srf.dom
    c = np.array([[u0, v0, 0], [u1, v0, 0], [u1, v1, 0], [u0, v1, 0], [u0, v0, 0]], float)
    return _loop_from_uv_pieces(srf, [np.linspace(a, b, 33) for a, b in zip(c[:-1], c[1:])], h)


def read_iges(path, mesh_size, region=None, layers=None, verbose=True, info_only=False):
    ig = _IgesFile(path)
    h = float(mesh_size)
    faces, used = [], set()                 # (surface DE, loops spec, owner DE)
    for de, e in ig.ents.items():
        if e["type"] == 144:
            t = ig.p(de)
            pts, n1, n2, pto = _i(t[0]), _i(t[1]), _i(t[2]), _i(t[3])
            curves = ([pto] if (n1 and pto) else []) + [_i(x) for x in t[4:4 + n2]]
            faces.append(("144", pts, curves, de, n1 == 0 or not pto))
            used.add(pts)
        elif e["type"] == 143:
            t = ig.p(de)
            sp, n = _i(t[1]), _i(t[2])
            faces.append(("143", sp, [_i(x) for x in t[3:3 + n]], de, n == 0))
            used.add(sp)
    for de, e in ig.ents.items():
        if e["type"] in (128, 108, 118, 120, 122, 140) and de not in used \
                and e["status"][2:4] == "00":
            faces.append(("plain", de, [], de, True))
    per_layer = Counter()
    for f in faces:
        for ln in ig.layer_names(f[3]):
            per_layer[ln] += 1
    lay_txt = ", ".join(f"{k} ({v})" for k, v in sorted(per_layer.items()))
    _say(verbose, f"IGES {os.path.basename(path)}: unit {ig.unit_name} "
                  f"(x {ig.scale:g} -> m), {len(faces)} surfaces; layers: {lay_txt or '-'}")
    if info_only:
        return np.zeros((0, 3, 3))
    tris, skipped, n_mesh, warn = [], Counter(), 0, Counter()
    for kind, sde, curves, owner, whole in faces:
        if layers and not any(_layer_ok(n, layers) for n in ig.layer_names(owner)):
            continue
        try:
            srf = ig.surface(sde)
        except NotImplementedError as err:
            skipped[str(err)] += 1
            continue
        if region is not None:
            lo, hi = srf.P.reshape(-1, 3).min(0), srf.P.reshape(-1, 3).max(0)
            if np.any(hi < region[0]) or np.any(lo > region[1]):
                continue
        loops = []
        if whole:
            loops.append(_domain_loop(srf, h))
        for c in curves:
            loop = None
            try:
                if kind == "144":
                    tc = ig.p(c)                      # 142: CRTN SPTR BPTR CPTR PREF
                    bptr, cptr = _i(tc[2]), _i(tc[3])
                    if bptr:
                        try:
                            loop = _loop_from_uv_pieces(srf, ig.curve_pieces(bptr, False), h)
                        except NotImplementedError:
                            loop = None
                    if loop is None and cptr:
                        loop = _loop_from_model_pieces(srf, ig.curve_pieces(cptr, True), h)
                else:                                 # 141 boundary
                    tc = ig.p(c)
                    n, i, uvp, mp = _i(tc[3]), 4, [], []
                    for _ in range(n):
                        crv, sense, k = _i(tc[i]), _i(tc[i + 1]), _i(tc[i + 2])
                        ps = [_i(x) for x in tc[i + 3:i + 3 + k]]
                        i += 3 + k
                        if ps:
                            for q in ps:
                                uvp += ig.curve_pieces(q, False)
                        else:
                            mp += ig.curve_pieces(crv, True)
                    loop = (_loop_from_uv_pieces(srf, uvp, h) if uvp else
                            _loop_from_model_pieces(srf, mp, h))
            except NotImplementedError as err:
                warn[str(err)] += 1
            if loop is not None:
                loops.append(loop)
        loops = [l for l in loops if l is not None]
        if not loops:
            loops = [_domain_loop(srf, h)]
            warn["trimming curves unreadable: surface used untrimmed"] += 1
        t = triangulate_trimmed(srf, loops, h)
        if len(t):
            tris.append(t)
            n_mesh += 1
    for k, v in skipped.items():
        print(f"WARNING: {v} surface(s) skipped ({k} not supported natively; "
              "use reader='gmsh' or export them as NURBS from GiD)")
    for k, v in warn.items():
        print(f"WARNING: {v} trimming loop(s): {k}")
    out = np.concatenate(tris) if tris else np.zeros((0, 3, 3))
    _say(verbose, f"IGES: meshed {n_mesh} surfaces"
                  + (" near the beam" if region is not None else "")
                  + f" at ~{h:g} m -> {len(out)} triangles")
    return out


# ##########################################################################
# gmsh (STEP, BREP, IGES on request, gmsh .msh)
# ##########################################################################

def read_with_gmsh(path, mesh_size, region=None, verbose=True, cad=True):
    try:
        import gmsh
    except ImportError:
        raise SurfaceReadError(
            f"reading {os.path.basename(path)} needs gmsh (pip install gmsh). Without it, "
            "export the geometry from GiD as IGES (read natively), or as GiD mesh / STL.")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        if not cad:
            gmsh.open(path)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            X = np.asarray(coords, float).reshape(-1, 3)
            pos = np.zeros(int(tags.max()) + 1, np.int64)
            pos[tags.astype(np.int64)] = np.arange(len(tags))
            blocks = []
            for dim in (2, 3):
                for et, en in zip(*gmsh.model.mesh.getElements(dim)[::2]):
                    name, _, _, nn, _, nprim = gmsh.model.mesh.getElementProperties(et)
                    kind = {2: {3: "tri", 4: "quad"}, 3: {4: "tet", 8: "hex", 6: "wedge",
                                                          5: "pyr"}}[dim].get(nprim)
                    if kind:
                        blocks.append((kind, pos[np.asarray(en, np.int64)].reshape(-1, nn)))
            return cells_to_triangles(X, blocks)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.add("surfaces")
        gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        surfs = gmsh.model.getEntities(2)
        keep = []
        for dim, tag in surfs:
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, tag)
            if region is None or (x1 >= region[0, 0] and x0 <= region[1, 0] and
                                  y1 >= region[0, 1] and y0 <= region[1, 1] and
                                  z1 >= region[0, 2] and z0 <= region[1, 2]):
                keep.append((dim, tag))
        _say(verbose, f"gmsh: {os.path.basename(path)}: {len(surfs)} surfaces, meshing "
                      f"{len(keep)} (size {mesh_size:g} m) ...")
        if not keep:
            return np.zeros((0, 3, 3))
        if len(keep) < len(surfs):
            gmsh.model.setVisibility(gmsh.model.getEntities(), 0, recursive=True)
            gmsh.model.setVisibility(keep, 1, recursive=True)
            gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size / 20)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        pos = {int(t): i for i, t in enumerate(tags)}
        X = np.asarray(coords, float).reshape(-1, 3)
        tris = []
        for dim, tag in keep:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, tag)
            for et, en in zip(etypes, enodes):
                if et == 2:
                    idx = np.array([pos[int(k)] for k in en]).reshape(-1, 3)
                    tris.append(X[idx])
        return np.concatenate(tris) if tris else np.zeros((0, 3, 3))
    finally:
        gmsh.finalize()


# ##########################################################################
# GiD ASCII mesh
# ##########################################################################

_GID_KIND = {"triangle": "tri", "quadrilateral": "quad", "tetrahedra": "tet",
             "tetrahedron": "tet", "hexahedra": "hex", "hexahedron": "hex",
             "prism": "wedge", "pyramid": "pyr"}


def read_gid_mesh(path, layers=None, verbose=True, info_only=False):
    coords, blocks = {}, []
    cur, mode = None, None
    rows = []
    with open(path, "r", encoding="latin-1", errors="replace") as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("mesh"):
                m = re.search(r'"([^"]*)"', line)
                et = re.search(r"elemtype\s*=?\s*(\w+)", low)
                nn = re.search(r"nnode\s*=?\s*(\d+)", low)
                cur = dict(name=m.group(1) if m else f"mesh{len(blocks) + 1}",
                           etype=et.group(1) if et else "", nnode=int(nn.group(1)) if nn else 0,
                           elems=[])
                blocks.append(cur)
            elif low.startswith("end coordinates") or low.startswith("end elements"):
                if mode == "c" and rows:
                    a = np.array(rows, float)
                    for r in a:
                        coords[int(r[0])] = r[1:4] if len(r) >= 4 else np.r_[r[1:3], 0.0]
                elif mode == "e" and cur is not None and rows:
                    cur["elems"] = rows
                mode, rows = None, []
            elif low.startswith("coordinates"):
                mode, rows = "c", []
            elif low.startswith("elements"):
                mode, rows = "e", []
            elif mode == "c":
                rows.append([float(x) for x in line.split()[:4]])
            elif mode == "e":
                rows.append([int(x) for x in line.split()])
    if not coords:
        raise SurfaceReadError(f"{path}: no GiD 'Coordinates' block found")
    ids = np.array(sorted(coords))
    X = np.array([coords[i] for i in ids])
    lut = {int(v): k for k, v in enumerate(ids)}
    desc = ", ".join(f"'{b['name']}' {b['etype']} ({len(b['elems'])})" for b in blocks)
    _say(verbose, f"GiD mesh {os.path.basename(path)}: {len(X)} nodes; MESH blocks: {desc}")
    if info_only:
        return np.zeros((0, 3, 3))
    cells = []
    for b in blocks:
        kind = _GID_KIND.get(b["etype"].lower())
        if not kind or not b["elems"] or not _layer_ok(b["name"], layers):
            continue
        nc = _NCORNER[kind]
        conn = np.array([[lut[v] for v in r[1:1 + nc]] for r in b["elems"] if len(r) > nc])
        cells.append((kind, conn))
    return cells_to_triangles(X, cells)


# ##########################################################################
# VTK legacy
# ##########################################################################

_VTK_KIND = {5: ("tri", None), 9: ("quad", None), 22: ("tri", None), 23: ("quad", None),
             28: ("quad", None), 10: ("tet", None), 24: ("tet", None), 12: ("hex", None),
             25: ("hex", None), 29: ("hex", None), 13: ("wedge", None), 26: ("wedge", None),
             14: ("pyr", None), 27: ("pyr", None), 8: ("quad", [0, 1, 3, 2]),
             11: ("hex", [0, 1, 3, 2, 4, 5, 7, 6])}
_VTK_DT = {"float": ">f4", "double": ">f8", "int": ">i4", "unsigned_int": ">u4",
           "long": ">i8", "unsigned_long": ">u8", "short": ">i2", "unsigned_short": ">u2",
           "char": ">i1", "unsigned_char": ">u1", "vtktypeint64": ">i8", "vtktypeint32": ">i4",
           "vtkidtype": ">i8"}


class _VtkReader:
    def __init__(self, data, binary):
        self.d, self.pos, self.bin = data, 0, binary

    def line(self):
        while self.pos < len(self.d):
            j = self.d.find(b"\n", self.pos)
            j = len(self.d) if j < 0 else j
            s = self.d[self.pos:j].decode("latin-1").strip()
            self.pos = j + 1
            if s:
                return s
        return None

    def array(self, n, dtype):
        if self.bin:
            dt = np.dtype(_VTK_DT.get(dtype.lower(), ">f4"))
            a = np.frombuffer(self.d, dt, count=n, offset=self.pos)
            self.pos += n * dt.itemsize
            return a.astype(float if dt.kind == "f" else np.int64)
        out = []
        while len(out) < n:
            s = self.line()
            if s is None:
                break
            out += s.split()
        return np.array(out[:n], float)


def read_vtk_legacy(path, verbose=True):
    with open(path, "rb") as f:
        data = f.read()
    if data.lstrip()[:5] == b"<?xml":
        return read_with_meshio(path, verbose)
    r = _VtkReader(data, False)
    head = [r.line(), r.line(), r.line()]
    if not head[0] or "vtk" not in head[0].lower():
        raise SurfaceReadError(f"{path}: not a legacy VTK file")
    r.bin = head[2].upper().startswith("BINARY")
    X, cells, types, polys = None, None, None, []

    def read_cellarray(tok):
        n1, n2 = int(tok[1]), int(tok[2])
        save = r.pos
        nxt = r.line()
        if nxt and nxt.upper().startswith("OFFSETS"):           # VTK 5.x layout
            off = r.array(n1, nxt.split()[1]).astype(np.int64)
            c = r.line()
            conn = r.array(n2, c.split()[1]).astype(np.int64)
            return [conn[off[i]:off[i + 1]] for i in range(len(off) - 1)]
        r.pos = save
        a = r.array(n2, "int").astype(np.int64)
        out, i = [], 0
        while i < len(a):
            k = int(a[i])
            out.append(a[i + 1:i + 1 + k])
            i += 1 + k
        return out

    while True:
        s = r.line()
        if s is None:
            break
        tok = s.split()
        key = tok[0].upper()
        if key == "DATASET":
            if tok[1].upper() not in ("POLYDATA", "UNSTRUCTURED_GRID"):
                raise SurfaceReadError(f"{path}: VTK dataset {tok[1]} has no surface cells")
        elif key == "POINTS":
            X = r.array(3 * int(tok[1]), tok[2]).reshape(-1, 3)
        elif key == "CELLS":
            cells = read_cellarray(tok)
        elif key == "CELL_TYPES":
            types = r.array(int(tok[1]), "int").astype(int)
        elif key == "POLYGONS":
            polys += [list(c) for c in read_cellarray(tok)]
        elif key == "TRIANGLE_STRIPS":
            for c in read_cellarray(tok):
                polys += [[c[i], c[i + 1], c[i + 2]] for i in range(len(c) - 2)]
        elif key in ("LINES", "VERTICES"):
            read_cellarray(tok)
        elif key in ("POINT_DATA", "CELL_DATA", "FIELD"):
            break
    if X is None:
        raise SurfaceReadError(f"{path}: no POINTS in the VTK file")
    blocks = [("poly", polys)] if polys else []
    if cells is not None and types is not None:
        by = defaultdict(list)
        for c, ty in zip(cells, types):
            if ty in _VTK_KIND:
                kind, order = _VTK_KIND[ty]
                c = np.asarray(c)[: _NCORNER[kind]] if order is None else np.asarray(c)[order]
                by[kind].append(c)
            elif ty in (6, 7):
                blocks.append(("poly", [list(c)] if ty == 7 else
                               [[c[i], c[i + 1], c[i + 2]] for i in range(len(c) - 2)]))
        blocks += [(k, np.array(v)) for k, v in by.items()]
    _say(verbose, f"VTK {os.path.basename(path)}: {len(X)} points")
    return cells_to_triangles(X, blocks)


# ##########################################################################
# DXF
# ##########################################################################

def read_dxf(path, layers=None, verbose=True, info_only=False):
    with open(path, "r", encoding="latin-1", errors="replace") as f:
        raw = [l.strip() for l in f]
    if len(raw) % 2:
        raw.append("")
    pairs = [(int(raw[i]) if raw[i].lstrip("-").isdigit() else -999, raw[i + 1])
             for i in range(0, len(raw) - 1, 2)]
    units = None
    for k, (c, v) in enumerate(pairs):
        if c == 9 and v == "$INSUNITS" and k + 1 < len(pairs):
            units = {1: 0.0254, 2: 0.3048, 4: 1e-3, 5: 0.01, 6: 1.0}.get(_i(pairs[k + 1][1]))
            break
    # entities: list of (type, [(code, value), ...])
    ents, in_ent, cur = [], False, None
    for c, v in pairs:
        if c == 0:
            if v == "SECTION":
                in_ent = False
            if cur is not None:
                ents.append(cur)
                cur = None
            if in_ent and v not in ("ENDSEC", "EOF"):
                cur = (v, [])
            continue
        if c == 2 and v == "ENTITIES":
            in_ent = True
            continue
        if cur is not None:
            cur[1].append((c, v))
    if cur is not None:
        ents.append(cur)
    pts, faces, acis, lay_count = [], [], 0, Counter()

    def gv(data, code, default=0.0):
        for c, v in data:
            if c == code:
                return v
        return default

    i = 0
    while i < len(ents):
        typ, data = ents[i]
        layer = gv(data, 8, "0")
        if typ in ("3DFACE", "SOLID", "POLYLINE", "MESH", "3DSOLID", "BODY", "REGION", "SURFACE"):
            lay_count[layer] += 1
        ok = _layer_ok(layer, layers)
        if typ in ("3DFACE", "SOLID") and ok:
            P = [[_f(gv(data, 10 + k)), _f(gv(data, 20 + k)), _f(gv(data, 30 + k))] for k in range(4)]
            if typ == "SOLID":
                P = [P[0], P[1], P[3], P[2]]
            b = len(pts)
            pts += P
            if np.allclose(P[2], P[3]):
                faces.append([b, b + 1, b + 2])
            else:
                faces.append([b, b + 1, b + 2, b + 3])
        elif typ == "POLYLINE":
            flags = _i(gv(data, 70, 0))
            m, n = _i(gv(data, 71, 0)), _i(gv(data, 72, 0))
            verts, frecs = [], []
            i += 1
            while i < len(ents) and ents[i][0] == "VERTEX":
                d = ents[i][1]
                vf = _i(gv(d, 70, 0))
                if vf & 128 and not vf & 64:
                    frecs.append([abs(_i(gv(d, 71 + k, 0))) for k in range(4)])
                else:
                    verts.append([_f(gv(d, 10)), _f(gv(d, 20)), _f(gv(d, 30))])
                i += 1
            if ok and verts:
                b = len(pts)
                pts += verts
                if flags & 64:
                    for fr in frecs:
                        fr = [b + x - 1 for x in fr if x]
                        if len(fr) >= 3:
                            faces.append(fr)
                elif flags & 16 and m * n <= len(verts):
                    cm, cn = bool(flags & 1), bool(flags & 32)
                    idx = np.arange(m * n).reshape(m, n) + b
                    for a in range(m - 1 + cm):
                        for c in range(n - 1 + cn):
                            a1, c1 = (a + 1) % m, (c + 1) % n
                            faces.append([idx[a, c], idx[a1, c], idx[a1, c1], idx[a, c1]])
            continue
        elif typ == "MESH" and ok:
            vals = data
            k = 0
            verts, flist = [], []
            while k < len(vals):
                c, v = vals[k]
                if c == 92:
                    nv = _i(v)
                    xs = [(cc, vv) for cc, vv in vals[k + 1:] if cc in (10, 20, 30)][:3 * nv]
                    verts = np.array([_f(vv) for _, vv in xs]).reshape(-1, 3)
                elif c == 93:
                    nl = _i(v)
                    flist = [_i(vv) for cc, vv in vals[k + 1:k + 1 + nl] if cc == 90]
                    break
                k += 1
            b = len(pts)
            pts += [list(p) for p in verts]
            j = 0
            while j < len(flist):
                nn = flist[j]
                faces.append([b + x for x in flist[j + 1:j + 1 + nn]])
                j += 1 + nn
        elif typ in ("3DSOLID", "BODY", "REGION", "SURFACE"):
            acis += 1
        i += 1
    _say(verbose, f"DXF {os.path.basename(path)}: layers: "
                  + ", ".join(f"{k} ({v})" for k, v in sorted(lay_count.items())))
    if acis:
        print(f"WARNING: {acis} ACIS solid/surface entities in the DXF cannot be read "
              "(export IGES/STEP for those).")
    if info_only:
        return np.zeros((0, 3, 3)), units
    if not faces:
        raise SurfaceReadError(f"{path}: no faces (3DFACE / polyface / MESH) in the DXF. "
                               "GiD writes surfaces to DXF only as their mesh: mesh them "
                               "first, or export IGES instead.")
    return cells_to_triangles(np.array(pts, float), [("poly", faces)]), units


# ##########################################################################
# OBJ / OFF / PLY
# ##########################################################################

def read_obj(path):
    V, F = [], []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == "v":
                V.append([float(x) for x in t[1:4]])
            elif t[0] == "f":
                idx = [int(x.split("/")[0]) for x in t[1:]]
                F.append([i - 1 if i > 0 else len(V) + i for i in idx])
    return cells_to_triangles(np.array(V), [("poly", F)])


def read_off(path):
    with open(path, "r", encoding="latin-1") as f:
        toks = [l.split("#", 1)[0] for l in f]
    toks = " ".join(toks).split()
    if not toks[0].upper().endswith("OFF"):
        raise SurfaceReadError(f"{path}: not an OFF file")
    nv, nf = int(toks[1]), int(toks[2])
    X = np.array(toks[4:4 + 3 * nv], float).reshape(-1, 3)
    i, F = 4 + 3 * nv, []
    for _ in range(nf):
        k = int(toks[i])
        F.append([int(x) for x in toks[i + 1:i + 1 + k]])
        i += 1 + k
    return cells_to_triangles(X, [("poly", F)])


_PLY_DT = {"char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1", "short": "i2",
           "int16": "i2", "ushort": "u2", "uint16": "u2", "int": "i4", "int32": "i4",
           "uint": "u4", "uint32": "u4", "float": "f4", "float32": "f4", "double": "f8",
           "float64": "f8"}


def read_ply(path):
    with open(path, "rb") as f:
        data = f.read()
    end = data.find(b"end_header")
    head = data[:end].decode("latin-1").splitlines()
    pos = data.find(b"\n", end) + 1
    fmt, elems = "ascii", []
    for l in head:
        t = l.split()
        if not t:
            continue
        if t[0] == "format":
            fmt = t[1]
        elif t[0] == "element":
            elems.append([t[1], int(t[2]), []])
        elif t[0] == "property":
            elems[-1][2].append(t[1:])
    end_c = "<" if fmt == "binary_little_endian" else ">"
    X, F = None, []
    if fmt == "ascii":
        rows = data[pos:].decode("latin-1").split("\n")
        r = 0
        for name, n, props in elems:
            for _ in range(n):
                vals = rows[r].split()
                r += 1
                if name == "vertex":
                    if X is None:
                        X = []
                        names = [p[-1] for p in props]
                        ix = [names.index(c) for c in ("x", "y", "z")]
                    X.append([float(vals[k]) for k in ix])
                elif name == "face":
                    k = int(vals[0])
                    F.append([int(v) for v in vals[1:1 + k]])
        X = np.array(X)
    else:
        for name, n, props in elems:
            if all(p[0] != "list" for p in props):
                dt = np.dtype([(p[-1], end_c + _PLY_DT[p[0]]) for p in props])
                a = np.frombuffer(data, dt, count=n, offset=pos)
                pos += n * dt.itemsize
                if name == "vertex":
                    X = np.stack([a["x"], a["y"], a["z"]], 1).astype(float)
            else:
                for _ in range(n):
                    face = None
                    for p in props:
                        if p[0] == "list":
                            ct, it = np.dtype(end_c + _PLY_DT[p[1]]), np.dtype(end_c + _PLY_DT[p[2]])
                            k = int(np.frombuffer(data, ct, 1, pos)[0])
                            pos += ct.itemsize
                            v = np.frombuffer(data, it, k, pos)
                            pos += k * it.itemsize
                            if face is None:
                                face = v.tolist()
                        else:
                            pos += np.dtype(_PLY_DT[p[0]]).itemsize
                    if name == "face":
                        F.append(face)
    return cells_to_triangles(X, [("poly", F)])


# ##########################################################################
# meshio (generic), Exodus II, HDF5, JSON, Rhino
# ##########################################################################

_MESHIO_KIND = [("triangle", "tri"), ("quad", "quad"), ("tetra", "tet"), ("hexahedron", "hex"),
                ("wedge", "wedge"), ("pyramid", "pyr")]


def read_with_meshio(path, verbose=True):
    try:
        import meshio
    except ImportError:
        raise SurfaceReadError(f"reading {os.path.basename(path)} needs meshio "
                               "(pip install meshio), or export it from GiD as GiD mesh, "
                               "VTK (legacy), STL or IGES.")
    m = meshio.read(path)
    blocks = []
    for cb in m.cells:
        typ = cb.type
        for pre, kind in _MESHIO_KIND:
            if typ.startswith(pre):
                blocks.append((kind, np.asarray(cb.data)))
                break
        else:
            if typ == "polygon":
                blocks.append(("poly", [list(c) for c in cb.data]))
    X = np.asarray(m.points, float)
    if X.shape[1] == 2:
        X = np.column_stack([X, np.zeros(len(X))])
    _say(verbose, f"meshio {os.path.basename(path)}: {len(X)} points, cell blocks "
                  + ", ".join(f"{cb.type} ({len(cb.data)})" for cb in m.cells))
    return cells_to_triangles(X, blocks)


def _exo_kind(et):
    et = et.upper()
    for pre, kind in (("TRI", "tri"), ("QUAD", "quad"), ("SHELL", "quad"), ("TET", "tet"),
                      ("HEX", "hex"), ("WEDGE", "wedge"), ("PYRAMID", "pyr")):
        if et.startswith(pre):
            return kind
    return None


def read_exodus(path, layers=None, verbose=True, info_only=False):
    with open(path, "rb") as f:
        magic = f.read(8)
    vars_, attrs = {}, {}
    if magic[:3] == b"CDF":
        from scipy.io import netcdf_file
        with netcdf_file(path, "r", mmap=False) as nc:
            for k, v in nc.variables.items():
                vars_[k] = np.array(v.data)
                attrs[k] = {a: getattr(v, a) for a in v._attributes}
    elif magic[:4] == b"\x89HDF":
        try:
            import h5py
        except ImportError:
            raise SurfaceReadError("this Exodus file is netCDF-4 (HDF5): pip install h5py")
        with h5py.File(path, "r") as h:
            for k, v in h.items():
                if isinstance(v, h5py.Dataset):
                    vars_[k] = np.array(v)
                    attrs[k] = dict(v.attrs)
    else:
        raise SurfaceReadError(f"{path}: not an Exodus II (netCDF) file")
    if "coordx" in vars_:
        X = np.stack([vars_["coordx"], vars_["coordy"],
                      vars_.get("coordz", np.zeros_like(vars_["coordx"]))], 1).astype(float)
    elif "coord" in vars_:
        c = vars_["coord"]
        X = (c.T if c.shape[0] <= 3 else c).astype(float)
        if X.shape[1] == 2:
            X = np.column_stack([X, np.zeros(len(X))])
    else:
        raise SurfaceReadError(f"{path}: no node coordinates (coordx/coord) found")
    names = []
    if "eb_names" in vars_:
        for row in np.atleast_2d(vars_["eb_names"]):
            names.append(b"".join(np.asarray(row).astype("S1").tolist()).decode(
                "latin-1", "ignore").strip("\x00 "))
    blocks, desc = [], []
    keys = sorted((k for k in vars_ if re.fullmatch(r"connect\d+", k)), key=lambda s: int(s[7:]))
    for n, k in enumerate(keys):
        et = attrs[k].get("elem_type", b"")
        et = et.decode() if isinstance(et, bytes) else str(et)
        bname = names[n] if n < len(names) and names[n] else f"block{n + 1}"
        desc.append(f"'{bname}' {et} ({len(vars_[k])})")
        kind = _exo_kind(et)
        if kind and _layer_ok(bname, layers):
            blocks.append((kind, vars_[k].astype(np.int64) - 1))
    _say(verbose, f"Exodus {os.path.basename(path)}: {len(X)} nodes; blocks: " + ", ".join(desc))
    if info_only:
        return np.zeros((0, 3, 3))
    return cells_to_triangles(X, blocks)


def _guess_mesh(coords_cands, conn_cands, where):
    """Pick coordinates (n, 3) and connectivity arrays that fit them."""
    if not coords_cands:
        raise SurfaceReadError(f"{where}: no (n, 3) coordinate array found")
    name, X = max(coords_cands, key=lambda c: len(c[1]))
    n = len(X)
    blocks, used = [], []
    for cname, C in conn_cands:
        C = np.asarray(C, np.int64)
        if C.ndim != 2 or C.shape[1] not in (3, 4, 6, 8) or not len(C):
            continue
        base = 1 if (C.min() >= 1 and C.max() == n) else 0
        if C.min() - base < 0 or C.max() - base >= n:
            continue
        kind = {3: "tri", 4: "quad", 6: "tri", 8: "quad"}[C.shape[1]]
        blocks.append((kind, C - base))
        used.append(f"{cname} {C.shape}")
    if not blocks:
        raise SurfaceReadError(f"{where}: coordinates found ({name}, {n} points) but no "
                               "triangle/quad connectivity array that fits them")
    print(f"NOTE: {where}: using coordinates '{name}' ({n}) and faces " + ", ".join(used))
    return cells_to_triangles(X, blocks)


def read_hdf5(path, verbose=True):
    try:
        import h5py
    except ImportError:
        raise SurfaceReadError("reading HDF5 needs h5py (pip install h5py)")
    coords, conns, tree = [], [], []
    with h5py.File(path, "r") as h:
        if "coordx" in h or "connect1" in h:
            return read_exodus(path, verbose=verbose)

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                tree.append(f"  {name} {obj.shape} {obj.dtype}")
                a = obj
                if a.ndim == 2 and a.dtype.kind == "f" and 3 in a.shape:
                    arr = np.array(a)
                    coords.append((name, arr if arr.shape[1] == 3 else arr.T))
                elif a.ndim == 2 and a.dtype.kind in "iu":
                    conns.append((name, np.array(a)))
        h.visititems(visit)
    try:
        return _guess_mesh(coords, conns, os.path.basename(path))
    except SurfaceReadError as err:
        raise SurfaceReadError(f"{err}. Datasets in the file:\n" + "\n".join(tree[:60]) +
                               "\n(send this listing to have the format added, or export "
                               "GiD mesh / VTK instead)")


def read_json(path, verbose=True):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        doc = json.load(f)
    coords, conns, keys = [], [], set()

    def walk(o, p):
        if isinstance(o, dict):
            for k, v in o.items():
                keys.add(k)
                walk(v, f"{p}/{k}")
        elif isinstance(o, list) and o:
            try:
                a = np.array(o)
            except Exception:
                a = None
            if a is not None and a.dtype.kind in "iuf" and a.ndim == 2:
                if a.dtype.kind == "f" and a.shape[1] == 3:
                    coords.append((p, a))
                elif a.dtype.kind in "iu":
                    conns.append((p, a))
                    if a.shape[1] == 3 and re.search(r"(coord|vert|point|node)", p, re.I):
                        coords.append((p, a.astype(float)))
            elif a is not None and a.dtype.kind in "iuf" and a.ndim == 1 and len(a) % 3 == 0 \
                    and re.search(r"(coord|vert|point|node|position)", p, re.I):
                coords.append((p, a.reshape(-1, 3).astype(float)))
            elif a is not None and a.dtype.kind in "iu" and a.ndim == 1 and len(a) % 3 == 0 \
                    and re.search(r"(tri|face|index|indices|conn|elem)", p, re.I):
                conns.append((p, a.reshape(-1, 3)))
            else:
                for k, v in enumerate(o[:100000]):
                    if isinstance(v, (dict, list)):
                        walk(v, f"{p}[{k}]")
    walk(doc, "")
    coords = [c for c in coords if not any(c[1] is cc[1] for cc in conns)] or coords
    try:
        return _guess_mesh(coords, conns, os.path.basename(path))
    except SurfaceReadError as err:
        raise SurfaceReadError(f"{err}. Keys in the JSON: {sorted(keys)[:40]} (send a sample "
                               "to have this format added, or export IGES / GiD mesh instead)")


def read_3dm(path, mesh_size, verbose=True):
    try:
        import rhino3dm as r3
    except ImportError:
        raise SurfaceReadError("reading Rhino .3dm needs rhino3dm (pip install rhino3dm), "
                               "or export IGES from GiD instead")
    model = r3.File3dm.Read(path)
    if model is None:
        raise SurfaceReadError(f"{path}: rhino3dm could not read the file")
    tris, untrimmed = [], 0

    def add_mesh(m):
        V = np.array([[m.Vertices[i].X, m.Vertices[i].Y, m.Vertices[i].Z]
                      for i in range(len(m.Vertices))])
        F = []
        for i in range(len(m.Faces)):
            f = m.Faces[i]
            F.append([f[0], f[1], f[2]] if f[2] == f[3] else [f[0], f[1], f[2], f[3]])
        if len(V) and F:
            tris.append(cells_to_triangles(V, [("poly", F)]))

    for obj in model.Objects:
        g = obj.Geometry
        if isinstance(g, r3.Mesh):
            add_mesh(g)
        elif isinstance(g, (r3.Brep, r3.Extrusion)):
            got = False
            faces = g.Faces if isinstance(g, r3.Brep) else []
            for fi in range(len(faces)):
                m = faces[fi].GetMesh(r3.MeshType.Any)
                if m is not None and len(m.Faces):
                    add_mesh(m)
                    got = True
            if not got:
                srfs = [faces[fi] for fi in range(len(faces))] or [g]
                for s in srfs:
                    try:
                        du, dv = s.Domain(0), s.Domain(1)
                        n = 40
                        us, vs = np.linspace(du.T0, du.T1, n), np.linspace(dv.T0, dv.T1, n)
                        P = np.array([[[*(lambda q: (q.X, q.Y, q.Z))(s.PointAt(u, v))]
                                       for v in vs] for u in us])
                        idx = np.arange(n * n).reshape(n, n)
                        a, b = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()
                        c, d = idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
                        T = np.vstack([np.stack([a, b, c], 1), np.stack([a, c, d], 1)])
                        tris.append(P.reshape(-1, 3)[T])
                        untrimmed += 1
                    except Exception:
                        pass
    if untrimmed:
        print(f"WARNING: {untrimmed} Rhino surface(s) without render mesh were sampled "
              "UNTRIMMED (export IGES from GiD for exact trimmed surfaces).")
    return np.concatenate(tris) if tris else np.zeros((0, 3, 3))


# ##########################################################################
# Dispatcher
# ##########################################################################

def _sniff(path, n=4096):
    with open(path, "rb") as f:
        return f.read(n)


def read_surface_file(path, mesh_size=0.02, region=None, units="m", cad_reader="auto",
                      layers=None, verbose=True, info_only=False):
    """Read any supported file -> (nt, 3, 3) triangles in METRES.

    mesh_size  : target triangle size [m] for CAD (IGES / STEP / Rhino) surfaces
    region     : ((xmin,ymin,zmin),(xmax,ymax,zmax)) in m: CAD surfaces outside
                 are not meshed (speed); None = everything
    units      : unit of files that carry no unit (all mesh formats, DXF without
                 $INSUNITS); IGES and STEP are converted from their own unit
    cad_reader : "auto" | "native" | "gmsh"  (IGES: native unless "gmsh")
    layers     : list of layer / block names to keep (IGES layers, DXF layers,
                 GiD-mesh MESH names, Exodus block names); None = all
    """
    if not os.path.isfile(path):
        raise SurfaceReadError(f"file not found: {path!r}")
    if isinstance(layers, str):
        layers = [layers]
    region = None if region is None else np.asarray(region, float).reshape(2, 3)
    ext = os.path.splitext(path)[1].lower()
    fs = UNIT_TO_M[units] if isinstance(units, str) else float(units)
    head = _sniff(path)
    htxt = head.decode("latin-1", "ignore")
    is_iges = len(head) > 72 and head[72:73] == b"S"

    def done(tris, scale):
        return np.asarray(tris, float) * scale

    if ext in (".igs", ".iges") or (is_iges and ext not in (".stl",)):
        if cad_reader == "gmsh":
            return done(read_with_gmsh(path, mesh_size, region, verbose), 1.0)
        return done(read_iges(path, mesh_size, region, layers, verbose, info_only), 1.0)
    if ext in (".stp", ".step", ".brep", ".brp"):
        return done(read_with_gmsh(path, mesh_size, region, verbose), 1.0)
    if ext in (".sat", ".sab", ".asat", ".asab"):
        raise SurfaceReadError(
            "ACIS (.sat) cannot be read: there is no free ACIS reader. In GiD use "
            "Files > Export > IGES (read natively, keeps the trimmed NURBS and layers), "
            "STEP (needs gmsh) or GiD mesh / STL of the meshed surfaces.")
    if ext == ".3dm":
        return done(read_3dm(path, mesh_size, verbose), fs)
    if ext == ".stl" or htxt.lstrip().lower().startswith("solid"):
        return done(read_stl(path), fs)
    if ext == ".dxf":
        tris, u = read_dxf(path, layers, verbose, info_only)
        return done(tris, u if u else fs)
    if "$MeshFormat" in htxt:
        if cad_reader != "native":
            try:
                return done(read_with_gmsh(path, mesh_size, None, verbose, cad=False), fs)
            except SurfaceReadError:
                pass
        return done(read_with_meshio(path, verbose), fs)
    if re.search(r"^\s*MESH\b", htxt, re.I | re.M) and re.search(r"elemtype", htxt, re.I):
        return done(read_gid_mesh(path, layers, verbose, info_only), fs)
    if ext == ".vtk" or htxt.startswith("# vtk"):
        return done(read_vtk_legacy(path, verbose), fs)
    if ext in (".exo", ".e", ".ex2", ".exii", ".g", ".gen") or \
            (head[:3] == b"CDF" and ext not in (".nc",)):
        return done(read_exodus(path, layers, verbose, info_only), fs)
    if ext in (".h5", ".hdf5", ".hdf", ".he5") or head[:4] == b"\x89HDF":
        return done(read_hdf5(path, verbose), fs)
    if ext == ".obj":
        return done(read_obj(path), fs)
    if ext == ".off":
        return done(read_off(path), fs)
    if ext == ".ply":
        return done(read_ply(path), fs)
    if ext == ".json":
        return done(read_json(path, verbose), fs)
    if ext in (".msh",):
        return done(read_gid_mesh(path, layers, verbose, info_only), fs)
    try:
        return done(read_with_meshio(path, verbose), fs)
    except SurfaceReadError:
        raise
    except Exception as err:
        raise SurfaceReadError(
            f"{os.path.basename(path)}: format not recognised ({err}). Readable: IGES, STEP "
            "(gmsh), DXF, Rhino .3dm, VTK, GiD mesh, Exodus II, HDF5, OBJ/OFF/PLY, STL and "
            "anything meshio reads. GiD's 'Calculation file', 'Text data report' and "
            "'Nodes->Raster' are not geometry.")


# ##########################################################################
# Command line: convert / inspect
# ##########################################################################

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("-o", "--output", help="write the triangles to this STL (in --out-units)")
    ap.add_argument("--mesh-size", type=float, default=0.02, help="CAD mesh size [m]")
    ap.add_argument("--units", default="m", choices=list(UNIT_TO_M),
                    help="unit of mesh files without unit information")
    ap.add_argument("--out-units", default="m", choices=list(UNIT_TO_M))
    ap.add_argument("--reader", default="auto", choices=["auto", "native", "gmsh"])
    ap.add_argument("--layers", nargs="*", default=None, help="keep only these layers / blocks")
    ap.add_argument("--info", action="store_true", help="only list layers / blocks")
    a = ap.parse_args(argv)
    try:
        tris = read_surface_file(a.file, a.mesh_size, None, a.units, a.reader, a.layers,
                                 info_only=a.info)
    except SurfaceReadError as err:
        sys.exit(f"ERROR: {err}")
    if a.info:
        return
    if len(tris):
        lo, hi = tris.reshape(-1, 3).min(0), tris.reshape(-1, 3).max(0)
        print(f"{len(tris)} triangles, box [m]: {np.round(lo, 4)} .. {np.round(hi, 4)}")
    if a.output:
        write_stl(a.output, tris / UNIT_TO_M[a.out_units], name="surfaces")
        print(f"Wrote {a.output} ({a.out_units})")


if __name__ == "__main__":
    main()
