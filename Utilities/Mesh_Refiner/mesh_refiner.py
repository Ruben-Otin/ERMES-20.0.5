#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mesh_refiner.py - Uniform sub-meshing of ERMES 20.0 problems (.dat files)
==============================================================================

Reads the mesh that GiD writes with the ERMES.gid templates
(Problem-1.dat: nodes, Problem-4.dat: tetrahedra), builds a conforming refined
mesh and rewrites every mesh-dependent .dat file with the same format:

  file     content (template)                          treatment
  -------  ------------------------------------------  ---------------------------
  -1.dat   No[i] = p(x,y,z);                           originals kept, new appended
  -2.dat   No[i].V.Fix(v); / No[i].V.FixC([m,p]);      propagated to new nodes
  -3.dat   No[i].Sg.Fix(k);                            propagated to new nodes
  -4.dat   VE(n1,n2,n3,n4,mat);                        tetrahedra refined
  -5.dat   JE(n1,n2,n3,n4,mat);                        tetrahedra refined
  -6.dat   PEC/PMC/TEC = n([n1,n2,n3]);                triangles refined
  -7.dat   FF(..) PFF/GRC/RWP/COP(n1,n2,n3,mat);       triangles refined
  -8.dat   PBC(n1,n2,n3,mat,face);                     triangles refined
  -9.dat   CE = n([n1,n2,n3,n4,n5,n6]);                prisms refined (n1n2n3 <-> n4n5n6)
  -10.dat  PRWP/PCOP/PSIE(n1,n2,n3,mat);               triangles refined
           PVIE(n1,n2,n3,n4,mat);                      tetrahedra refined
  -19.dat  n1 n2 n3 id                                 triangles refined
  others   (materials, plasma, frequency, solver...)  untouched

Every child element inherits all the trailing values of its parent (material,
condition id, PBC face...). Tetrahedra keep their orientation and triangles
their winding (normal direction).

Refinement modes
----------------
  'levels'  Each edge is divided in LEVEL parts and every tetrahedron is filled
            with the regular barycentric lattice: new nodes on edges (L>=2), on
            faces (L>=3) and inside the element (L>=4). The lattice is cut in
            L^3 tetrahedra: up/down tetrahedra (similar to the parent) plus
            octahedra, each split in 4 along the diagonal that maximises the
            minimum dihedral angle for the actual element geometry. Hence there
            are at most 3 shape classes per parent and quality does NOT degrade
            when L grows. Boundary triangles are split in L^2 triangles.
  'ps'      Powell-Sabin / Worsey-Farin 3D split: edge midpoints, face centroids
            and element centroid are added; tetrahedron -> 24, triangle -> 6.

Edge, face and interior nodes are numbered through global tables, so the mesh
is conforming everywhere (also across material interfaces and between element
lists that repeat the same tetrahedra, e.g. VE / JE / PVIE).

Parallelism
-----------
  The refinement + text formatting of the elements and of the new nodes (by far
  the most expensive part) is split in many tasks of equal size (same number of
  CHILD elements/nodes), at least ~4 tasks per process, dispatched dynamically
  (largest first) so that all processes finish at about the same time. Each
  task writes its own part file; parts are then concatenated in order.
    
  - Linux / macOS / WSL : 'fork'  -> workers share the mesh arrays for free.
  - Windows (native)    : 'spawn' -> the mesh arrays are placed once in shared
                          memory and attached by every worker (no copies).
                          
  Reading, topology tables and the final concatenation are serial (vectorised
  NumPy / OS-level file copies).

How to run
----------
  1) Edit the USER SETTINGS block below and run the script (e.g. F5 in Spyder).
  2) Or from a terminal / IPython console; command-line options override the
     USER SETTINGS (settings not given on the command line keep their value):

     >> python mesh_refiner.py PROBLEM_DIR --level 4 --np 16
     >> %run   mesh_refiner.py PROBLEM_DIR --mode ps --stats
     >> python mesh_refiner.py --help        (list of options + defaults)

  On native Windows, parallel runs work from a terminal ("python script.py").
  Inside Spyder/IPython on Windows, 'spawn' workers may fail to start: the
  script then falls back to serial mode automatically (or set Spyder's Run
  configuration to "Execute in an external system terminal").

Typical use in the GiD -> ERMES workflow (see RefineExecute.sh)
---------------------------------------------------------------
  GiD writes <name>.dat, <name>-1.dat ... ; this script refines them IN PLACE
  (default, BACKUP = False) and ERMES is run afterwards on the refined mesh.
  !! The refinement is NOT idempotent: running it twice on the same files
  !! refines the already refined mesh (tets x L^3 each time). Always let GiD
  !! rewrite the .dat files first (or use --out to keep the input untouched).
  Under WSL with the problem on /mnt/c or /mnt/d add '--tmp /tmp' (much faster).

Command-line options  (option -> USER SETTING it overrides)
--------------------
  PROBLEM_DIR  (positional, optional)               -> PROBLEM_DIR
      Folder with the ERMES .dat files, e.g. ./static.gid . If omitted, the
      value of PROBLEM_DIR in USER SETTINGS is used.

  --name NAME                                       -> PROBLEM_NAME
      Problem name = prefix of the .dat files (NAME-1.dat, NAME-4.dat ...).
      Default: detected automatically from the folder.

  --mode {levels,ps}                                -> MODE
      levels : divide every edge in LEVEL parts (tets x L^3, triangles x L^2)
      ps     : Powell-Sabin split (tets x 24, triangles x 6)

  --level L                                         -> LEVEL
      Edge divisions for --mode levels (integer >= 2). Ignored with ps.

  --np N                                            -> NPROC
      Number of parallel processes. Default: all cores. 1 = serial.

  --start {auto,fork,spawn}                         -> START_METHOD
      How worker processes are created. auto = fork on Linux/macOS/WSL,
      spawn on native Windows (shared memory). Normally leave 'auto'.

  --out DIR                                         -> OUT_DIR
      Write the refined problem to DIR (complete copy) instead of overwriting
      the .dat files in PROBLEM_DIR.

  --tmp DIR                                         -> TMP_DIR
      Folder for temporary part files (default: the output folder). Under WSL
      with the problem on /mnt/c or /mnt/d, '--tmp /tmp' is much faster.

  --backup / --no-backup                            -> BACKUP (True / False)
      Keep / do not keep a copy of the original files in
      PROBLEM_DIR/mesh_backup_original (only when overwriting).

  --nodecond-rule {surface,loose}                   -> NODECOND_RULE
      How voltages (-2.dat) and singular nodes (-3.dat) are propagated to the
      new nodes (see USER SETTINGS for details).

  --stats / --no-stats                              -> STATS (True / False)
      Print mesh-quality statistics (min dihedral angles, volume check).

  --check / --no-check                              -> CHECK (True / False)
      Read the written mesh back and verify conformity (memory heavy).

  --task-size N                                     -> TASK_SIZE
      Child elements per parallel task. 0 = automatic (recommended).

  -h, --help
      Show the option list with the current defaults and exit.

  Examples
     >> python mesh_refiner.py ./static.gid --level 3 --np 12 --no-backup
     >> python mesh_refiner.py ./static.gid --mode ps --out ./static_ps.gid
     >> python mesh_refiner.py ./static.gid --level 2 --tmp /tmp --stats --check
"""

#==============================================================================
#==============================  USER SETTINGS  ===============================
#==============================================================================

# Folder containing the ERMES .dat files (Problem.dat, Problem-1.dat, ...).
# A relative path is searched from the current working directory first and
# then from the folder of this script.
# [command line: positional argument PROBLEM_DIR]
PROBLEM_DIR = './ProblemFolderName.gid'

# Problem name = prefix of the .dat files. None -> detected automatically
# (the folder must then contain exactly one "<name>-1.dat" + "<name>-4.dat").
# [command line: --name]
PROBLEM_NAME = None

# Refinement mode: 'levels' (lattice subdivision) or 'ps' (Powell-Sabin, x24)
# [command line: --mode]
MODE = 'levels'

# Number of divisions per edge for MODE = 'levels' (>= 2). Element growth:
#   tetrahedra x LEVEL^3 , surface triangles x LEVEL^2
# [command line: --level]
LEVEL = 2

# Number of parallel processes. None -> all available cores. 1 -> serial.
# [command line: --np]
NPROC = None

# Process start method: 'auto' (fork where available, else spawn), 'fork',
# 'spawn' (uses shared memory; the only option on native Windows).
# [command line: --start]
START_METHOD = 'auto'

# Output folder. None -> overwrite the .dat files inside PROBLEM_DIR.
# If a folder is given, a complete copy of the problem is written there.
# [command line: --out]
OUT_DIR = None

# Folder for the temporary part files. None -> inside the output folder.
# Tip: under WSL with the problem on /mnt/c or /mnt/d, use a Linux folder
# (e.g. '/tmp') - much faster I/O. It needs about the size of the output.
# [command line: --tmp]
TMP_DIR = None

# When overwriting, copy the original files to PROBLEM_DIR/mesh_backup_original
# (an existing backup is never overwritten, so it always holds the first mesh).
# Use False when the script runs from the GiD .bat (GiD rewrites the .dat files
# on every "Calculate").
# [command line: --backup / --no-backup]
BACKUP = False

# Propagation of node conditions (-2.dat voltages, -3.dat singular nodes):
#   'surface' : a new node gets the condition only if it lies on an edge/face of
#               a boundary face, a material interface or a listed surface
#               element, and all parent nodes carry the same condition (safe
#               default, avoids edges that cut through corner elements).
#   'loose'   : any edge/face whose parent nodes carry the same condition
#               (needed e.g. for voltages on lines embedded inside a volume).
# [command line: --nodecond-rule]
NODECOND_RULE = 'surface'

# Print mesh-quality statistics (min dihedral angle histogram, volume check)
# [command line: --stats / --no-stats]
STATS = False

# Read the written mesh back and check conformity (memory heavy, roughly
# 1 GB per 5e6 tetrahedra; use it only on moderate sizes)
# [command line: --check / --no-check]
CHECK = False

# Target number of child elements per parallel task. 0 -> automatic
# (total / (4 x NPROC), clipped to [5e4, 2e6]).
# [command line: --task-size]
TASK_SIZE = 0

#==============================================================================
#============================  END OF USER SETTINGS  ==========================
#==============================================================================

import argparse
import itertools
import os
import re
import shutil
import sys
import tempfile
import time

import numpy as np

# -----------------------------------------------------------------------------
#  Element keyword catalogue:  keyword -> (line style, element kind, node fields)
#    style 'call'  : KW(a,b,c,...);        style 'nlist' : KW = n([a,b,c,...]);
#    style 'plain' : a b c id              (only file -19)
# -----------------------------------------------------------------------------
KEYWORDS = {
    'VE':   ('call', 'tet', 4), 'JE':   ('call', 'tet', 4), 'PVIE': ('call', 'tet', 4),
    'FF':   ('call', 'tri', 3), 'PFF':  ('call', 'tri', 3), 'GRC':  ('call', 'tri', 3),
    'RWP':  ('call', 'tri', 3), 'COP':  ('call', 'tri', 3), 'PBC':  ('call', 'tri', 3),
    'PRWP': ('call', 'tri', 3), 'PCOP': ('call', 'tri', 3), 'PSIE': ('call', 'tri', 3),
    'PEC':  ('nlist', 'tri', 3), 'PMC': ('nlist', 'tri', 3), 'TEC': ('nlist', 'tri', 3),
    'CE':   ('nlist', 'prism', 6),
    'PLAIN': ('plain', 'tri', 3),
}

# What may appear in each numbered file (anything else is copied verbatim)
FILE_CONTENT = {
    2: 'nodecond', 3: 'nodecond',
    4: 'elem', 5: 'elem', 6: 'elem', 7: 'elem', 8: 'elem', 9: 'elem', 10: 'elem',
    19: 'plain',
}

RE_NODE = re.compile(r'^\s*No\[\s*\d+\s*\]\s*=\s*p\(')
RE_NODECOND = re.compile(r'^\s*No\[\s*(\d+)\s*\]\s*(\..*?)\s*$')

# byte translation tables for fast number parsing with np.fromstring
_KEEP_INT = set(b'0123456789-')
_KEEP_FLT = set(b'0123456789-+.eE')
TR_INT = bytes(c if c in _KEEP_INT else 32 for c in range(256))
TR_FLT = bytes(c if c in _KEEP_FLT else 32 for c in range(256))


def parse_numbers(lines, ncols, dtype, table):
    """Fast parsing of a block of lines with exactly ncols numbers per line."""
    raw = ' '.join(lines).encode('latin-1').translate(table)
    vals = np.fromstring(raw, dtype=dtype, sep=' ')
    if vals.size != len(lines) * ncols:
        raise ValueError('unexpected number of values (%d, expected %d x %d)'
                         % (vals.size, len(lines), ncols))
    return vals.reshape(len(lines), ncols)


# =============================================================================
#  1. Reference templates
# =============================================================================
class Template:
    """Refinement pattern of the reference tetrahedron and triangle.

    Points are given by integer barycentric weights that sum to D:
      W4 : (P4 x 4) weights of the tetrahedron points
      W3 : (P3 x 3) weights of the triangle points
      tet_variants : list of (C x 4) child tetrahedra (indices into W4); the
                     variants differ only in the octahedron diagonals
      tri : (C3 x 3) child triangles (indices into W3)
    All children are positively oriented with respect to the parent ordering.
    """

    def __init__(self, D, W4, tet_variants, W3, tri, n_fixed=0):
        self.D = D
        self.W4 = np.asarray(W4, dtype=np.int64)
        self.W3 = np.asarray(W3, dtype=np.int64)
        self.tet_variants = tet_variants
        self.tri = tri
        self.n_fixed = n_fixed                       # cells common to all variants
        self.choose_diagonal = len(tet_variants) > 1
        self.n_tet_children = tet_variants[0].shape[0]
        self.n_tri_children = tri.shape[0]


def _orient(cells, ref_xyz):
    """Swap two indices of every cell with negative orientation in reference coords."""
    cells = np.array(cells, dtype=np.int64)
    P = ref_xyz[cells]
    if cells.shape[1] == 4:
        det = np.einsum('ij,ij->i', np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), P[:, 3] - P[:, 0])
        a, b = 2, 3
    else:
        u, v = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
        det = u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
        a, b = 1, 2
    if np.any(np.abs(det) < 1e-14):
        raise RuntimeError('degenerate cell in refinement template')
    neg = det < 0
    cells[neg, a], cells[neg, b] = cells[neg, b].copy(), cells[neg, a].copy()
    return cells


def template_levels(L):
    """Lattice subdivision with L divisions per edge (L^3 tets, L^2 triangles)."""
    # --- tetrahedron lattice; lattice coordinates x = (j,k,l) = weights of v1,v2,v3
    pts, idx = [], {}
    for j in range(L + 1):
        for k in range(L + 1 - j):
            for l in range(L + 1 - j - k):
                idx[(j, k, l)] = len(pts)
                pts.append((L - j - k - l, j, k, l))
    W4 = np.array(pts, dtype=np.int64)
    ref = W4[:, 1:] / L                              # reference tet v0=0, v1=e1, v2=e2, v3=e3

    def at(x, *shifts):
        return idx[tuple(x[i] + sum(s[i] for s in shifts) for i in range(3))]

    e1, e2, e3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    fixed, octs = [], []
    for j in range(L):
        for k in range(L - j):
            for l in range(L - j - k):
                x, s = (j, k, l), j + k + l
                fixed.append([at(x), at(x, e1), at(x, e2), at(x, e3)])           # up tet
                if s <= L - 2:                                                     # octahedron
                    octs.append({'1': at(x, e1), '2': at(x, e2), '3': at(x, e3),
                                 '12': at(x, e1, e2), '13': at(x, e1, e3), '23': at(x, e2, e3)})
                if s <= L - 3:                                                     # down tet
                    fixed.append([at(x, e1, e2), at(x, e1, e3), at(x, e2, e3), at(x, e1, e2, e3)])

    # Octahedron split along one of its 3 diagonals -> 4 tets around the diagonal.
    # The equator ring is (P, Q, opposite(P), opposite(Q)).
    opposite = {'1': '23', '23': '1', '2': '13', '13': '2', '3': '12', '12': '3'}
    diagonals = [('1', '23', '2', '3'), ('2', '13', '1', '3'), ('3', '12', '1', '2')]
    variants = []
    for A, B, P, Q in diagonals:
        cells = list(fixed)
        for o in octs:
            ring = [o[P], o[Q], o[opposite[P]], o[opposite[Q]]]
            cells += [[o[A], o[B], ring[r], ring[(r + 1) % 4]] for r in range(4)]
        variants.append(_orient(cells, ref))
    if not octs:                                     # L = 1: nothing to choose
        variants = variants[:1]

    # --- triangle lattice
    tpts, tidx = [], {}
    for j in range(L + 1):
        for k in range(L + 1 - j):
            tidx[(j, k)] = len(tpts)
            tpts.append((L - j - k, j, k))
    W3 = np.array(tpts, dtype=np.int64)
    tri = []
    for j in range(L):
        for k in range(L - j):
            tri.append([tidx[(j, k)], tidx[(j + 1, k)], tidx[(j, k + 1)]])
            if j + k <= L - 2:
                tri.append([tidx[(j + 1, k)], tidx[(j + 1, k + 1)], tidx[(j, k + 1)]])
    tri = _orient(tri, W3[:, 1:] / L)
    return Template(L, W4, variants, W3, tri, n_fixed=len(fixed))


def template_ps():
    """Powell-Sabin / Worsey-Farin split (weights over D = 12)."""
    D = 12
    pts = [(12, 0, 0, 0), (0, 12, 0, 0), (0, 0, 12, 0), (0, 0, 0, 12)]
    mid, fc = {}, {}
    for a, b in itertools.combinations(range(4), 2):                  # edge midpoints
        w = [0] * 4
        w[a] = w[b] = 6
        mid[(a, b)] = mid[(b, a)] = len(pts)
        pts.append(tuple(w))
    for f in itertools.combinations(range(4), 3):                     # face centroids
        fc[f] = len(pts)
        pts.append(tuple(4 if v in f else 0 for v in range(4)))
    C = len(pts)                                                      # element centroid
    pts.append((3, 3, 3, 3))
    W4 = np.array(pts, dtype=np.int64)
    cells = []
    for f in itertools.combinations(range(4), 3):
        i, j, l = f
        for a, b in ((i, j), (j, l), (l, i)):
            cells += [[a, mid[(a, b)], fc[f], C], [mid[(a, b)], b, fc[f], C]]
    W3 = np.array([(12, 0, 0), (0, 12, 0), (0, 0, 12), (6, 6, 0), (0, 6, 6), (6, 0, 6), (4, 4, 4)])
    tri = []
    for (a, b), m in {(0, 1): 3, (1, 2): 4, (2, 0): 5}.items():
        tri += [[a, m, 6], [m, b, 6]]
    return Template(D, W4, [_orient(cells, W4[:, 1:] / D)], W3, _orient(tri, W3[:, 1:] / D))


# =============================================================================
#  2. Unique integer-tuple tables (edges, faces, tetrahedra)
# =============================================================================
class RowTable:
    """Set of unique sorted integer rows with vectorised lookup row -> index.

    Rows are packed into big-endian byte strings, whose lexicographic order
    equals the integer order, so np.unique / np.searchsorted work for any id
    size (no overflow of packed integer keys)."""

    def __init__(self, rows=None, keys=None, width=None):
        if keys is not None:                         # rebuilt in a worker (lookup only)
            self.keys, self.width, self.rows = keys, width, None
        else:
            rows = np.ascontiguousarray(rows, dtype=np.int64)
            self.width = rows.shape[1]
            self.keys, first = np.unique(self._pack(rows), return_index=True)
            self.rows = rows[first]
        self.n = len(self.keys)

    def _pack(self, rows):
        return np.ascontiguousarray(rows.astype('>i8')).view('V%d' % (8 * self.width)).ravel()

    def lookup(self, rows):
        k = self._pack(np.ascontiguousarray(rows, dtype=np.int64))
        i = np.minimum(np.searchsorted(self.keys, k), self.n - 1)
        if not np.all(self.keys[i] == k):
            raise RuntimeError('entity not found in table: inconsistent mesh')
        return i


# =============================================================================
#  3. Reading the .dat files
# =============================================================================
def read_text(path):
    """Return file text and its line terminator (GiD on Windows writes CRLF)."""
    with open(path, 'rb') as f:
        raw = f.read()
    return raw.decode('latin-1'), ('\r\n' if b'\r\n' in raw else '\n')


def read_nodes(path):
    """Parse -1.dat. Returns (non-node lines, node lines verbatim, ids, xyz, newline)."""
    text, nl = read_text(path)
    header, body = [], []
    for ln in text.split(nl):
        if ln.lstrip().startswith('No['):
            body.append(ln)
        elif ln.strip():
            header.append(ln)
    vals = parse_numbers([ln.replace('No[', ' ', 1) for ln in body], 4, float, TR_FLT)
    return header, body, vals[:, 0].astype(np.int64), vals[:, 1:], nl


def classify(line, content):
    """Return the block key of a line, or None for text copied verbatim.
    (fast string tests; the numbers are validated when the block is parsed)"""
    s = line.strip()
    if not s:
        return None
    if content == 'elem':
        if s.endswith(');') and '(' in s:
            head = s[:s.index('(')].rstrip()
            if head.endswith('= n'):
                kw, style = head[:-3].strip(), 'nlist'
            else:
                kw, style = head, 'call'
            if KEYWORDS.get(kw, ('',))[0] == style:
                return kw
    elif content == 'plain':
        if s[0].isdigit() and len(s.split()) == 4:
            return 'PLAIN'
    elif content == 'nodecond':
        if s.startswith('No[') and not RE_NODE.match(s):
            return 'NODECOND'
    return None


def parse_dat(path, content):
    """Split a .dat file into a list of segments:
         {'type':'text', 'line':...}                       verbatim line
         {'type':'elem', 'kw', 'kind', 'nodes', 'extra'}    block of elements
         {'type':'nodecond', 'ids', 'rest', 'lines'}        block of node conditions
    """
    text, nl = read_text(path)
    lines = text.split(nl)
    if lines and lines[-1] == '':
        lines.pop()
    segs, block, key = [], [], None

    def flush():
        if not block:
            return
        if key == 'NODECOND':
            ms = [RE_NODECOND.match(ln.strip()) for ln in block]
            if not all(ms):
                raise ValueError('%s: unexpected node condition line' % path)
            segs.append({'type': 'nodecond', 'lines': list(block),
                         'ids': np.array([int(m.group(1)) for m in ms], dtype=np.int64),
                         'rest': [m.group(2) for m in ms]})
        else:
            style, kind, nn = KEYWORDS[key]
            first = block[0]
            ncols = len(re.findall(r'-?\d+', first.split('(', 1)[1] if style != 'plain' else first))
            try:
                vals = parse_numbers(block, ncols, np.int64, TR_INT)
            except ValueError as err:
                raise ValueError('%s, block %s: %s' % (os.path.basename(path), key, err))
            segs.append({'type': 'elem', 'kw': key, 'style': style, 'kind': kind,
                         'nodes': vals[:, :nn].copy(), 'extra': vals[:, nn:].copy(),
                         'trail': ' ' if (style == 'plain' and first.endswith(' ')) else ''})
        block.clear()

    for ln in lines:
        k = classify(ln, content)
        if k is None:
            flush()
            key = None
            segs.append({'type': 'text', 'line': ln})
        else:
            if k != key:
                flush()
                key = k
            block.append(ln)
    flush()
    return segs, nl


# =============================================================================
#  4. Refinement kernels (run inside worker processes; state in global S)
# =============================================================================
S = {}          # shared state of the workers
_SHM = []       # shared-memory handles kept alive in 'spawn' workers


def point_ids(nodes, W, owner=None):
    """Global node ids of template points.

    nodes : (E x k) parent connectivity (k = 3 triangle, 4 tetrahedron)
    W     : (P x k) template barycentric weights
    owner : (E,) unique-tet index (only needed for interior points; parents
            must then be locally sorted by global id)
    Node numbering:  original ids | edge points | face points | interior points
    """
    out = np.empty((len(nodes), len(W)), dtype=np.int64)
    cache = {}
    for p, w in enumerate(W):
        sup = tuple(np.nonzero(w)[0])                       # local vertices involved
        if len(sup) == 1:                                   # vertex
            out[:, p] = nodes[:, sup[0]]
        elif len(sup) == 2:                                 # edge point
            if sup not in cache:
                na, nb = nodes[:, sup[0]], nodes[:, sup[1]]
                e = S['edges'].lookup(np.stack([np.minimum(na, nb), np.maximum(na, nb)], 1))
                cache[sup] = (e, na < nb)
            e, a_is_low = cache[sup]
            w_low = np.where(a_is_low, w[sup[0]], w[sup[1]])          # weight of lower id
            out[:, p] = S['base_e'] + e * S['n_ep'] + S['e_lut'][w_low]
        elif len(sup) == 3:                                 # face point
            if sup not in cache:
                sub = nodes[:, list(sup)]
                order = np.argsort(sub, axis=1, kind='stable')
                f = S['faces'].lookup(np.take_along_axis(sub, order, 1))
                cache[sup] = (f, order)
            f, order = cache[sup]
            ws = np.asarray(w)[list(sup)][order]                       # canonical order
            out[:, p] = S['base_f'] + f * S['n_fp'] + S['f_lut'][ws[:, 0], ws[:, 1]]
        else:                                               # interior point
            out[:, p] = S['base_i'] + owner * S['n_ip'] + S['i_lut'][w[0], w[1], w[2]]
    return out


def tet_quality(P):
    """Minimum dihedral angle [deg] and signed volume of tetrahedra P (E x 4 x 3)."""
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    vol = np.einsum('ij,ij->i', np.cross(b - a, c - a), d - a) / 6.0
    normals = []
    for p, q, r, o in ((b, c, d, a), (a, c, d, b), (a, b, d, c), (a, b, c, d)):
        n = np.cross(q - p, r - p)
        n *= np.sign(np.einsum('ij,ij->i', n, p - o))[:, None]       # outward
        normals.append(n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300))
    mind = np.full(len(P), 180.0)
    for i in range(4):
        for j in range(i + 1, 4):
            cos = np.clip(-np.einsum('ij,ij->i', normals[i], normals[j]), -1.0, 1.0)
            mind = np.minimum(mind, np.degrees(np.arccos(cos)))
    return mind, vol


def refine_tets(nodes, extra):
    """Refine tetrahedra. Returns (children connectivity, children extra columns)."""
    T = S['T']
    # Canonical local order (sorted ids): identical children for repeated tets
    # (VE/JE/PVIE); the permutation parity is used to restore the orientation.
    # (1) sort the 4 vertex ids of every tet -> children depend only on the
    #     geometry/ids, not on the order in which the element was written;
    # (2) refine the sorted tet; (3) if the sorting permutation was odd, swap two
    #     children vertices (see the [0, 1, 3, 2] swap below) so every child keeps
    #     the orientation (sign of the volume) of its parent.
    perm = np.argsort(nodes, axis=1, kind='stable')
    srt = np.take_along_axis(nodes, perm, 1)
    inversions = sum((perm[:, i] > perm[:, j]).astype(int) for i in range(4) for j in range(i + 1, 4))
    odd = inversions % 2 == 1
    owner = S['tets'].lookup(srt) if S['n_ip'] else None
    ids = point_ids(srt, T.W4, owner)

    choice = np.zeros(len(nodes), dtype=np.int64)
    if T.choose_diagonal:
        # All octahedra of a parent are translates of each other: evaluate the
        # 3 diagonal choices on the first one, keep the best minimum angle.
        V = S['X'][srt]                                              # E x 4 x 3
        Wf = T.W4 / T.D
        score = []
        for cells in T.tet_variants:
            q = np.full(len(nodes), 180.0)
            for cell in cells[T.n_fixed:T.n_fixed + 4]:
                q = np.minimum(q, tet_quality(np.einsum('pk,ekx->epx', Wf[cell], V))[0])
            score.append(q)
        choice = np.argmax(np.stack(score, 1), axis=1)

    sub = np.empty((len(nodes), T.n_tet_children, 4), dtype=np.int64)
    for c, cells in enumerate(T.tet_variants):
        m = choice == c
        if m.any():
            sub[m] = ids[m][:, cells]
    sub[odd] = sub[odd][:, :, [0, 1, 3, 2]]
    return sub.reshape(-1, 4), np.repeat(extra, T.n_tet_children, axis=0)


def refine_tris(nodes, extra):
    T = S['T']
    sub = point_ids(nodes, T.W3)[:, T.tri]
    return sub.reshape(-1, 3), np.repeat(extra, T.n_tri_children, axis=0)


def refine_prisms(nodes, extra):
    """Contact prism: both triangles refined with the same local pattern, so
    child bottom/top nodes keep the n1<->n4, n2<->n5, n3<->n6 pairing."""
    T = S['T']
    bot = point_ids(nodes[:, :3], T.W3)[:, T.tri]
    top = point_ids(nodes[:, 3:], T.W3)[:, T.tri]
    return np.concatenate([bot, top], 2).reshape(-1, 6), np.repeat(extra, T.n_tri_children, axis=0)


REFINE = {'tet': refine_tets, 'tri': refine_tris, 'prism': refine_prisms}


def format_rows(meta, rows, nl):
    """ERMES text for a block of element rows (same syntax as the .bas files)."""
    if len(rows) == 0:
        return ''
    ints = ','.join(['%d'] * rows.shape[1])
    fmt = {'call': meta['kw'] + '(' + ints + ');',
           'nlist': meta['kw'] + ' = n([' + ints + ']);',
           'plain': ' '.join(['%d'] * rows.shape[1]) + meta['trail']}[meta['style']] + nl
    out, step = [], 200000
    for s in range(0, len(rows), step):
        r = rows[s:s + step]
        out.append((fmt * len(r)) % tuple(r.ravel().tolist()))
    return ''.join(out)


def element_task(task):
    """Worker: refine one chunk of one element block and write it to a part file."""
    _, gid, s, e, path = task
    meta = S['meta'][gid]
    nodes, extra = S['b%d_nodes' % gid][s:e], S['b%d_extra' % gid][s:e]
    sub, ext = REFINE[meta['kind']](nodes, extra)
    rows = np.concatenate([sub, ext], axis=1) if ext.shape[1] else sub
    with open(path, 'w', newline='', encoding='latin-1') as f:
        f.write(format_rows(meta, rows, meta['nl']))
    info = {'kw': meta['kw'], 'n': len(sub)}
    if S['stats'] and meta['kw'] == 'VE' and e > s:
        mind, vol = tet_quality(S['X'][sub])
        mpar, vpar = tet_quality(S['X'][nodes])
        nch = len(sub) // (e - s)
        info.update(hist=np.histogram(mind, bins=18, range=(0, 90))[0],
                    hist_par=np.histogram(mpar, bins=18, range=(0, 90))[0],
                    min=float(mind.min()), min_par=float(mpar.min()),
                    vol=float(np.abs(vol).sum()), vol_par=float(np.abs(vpar).sum()),
                    flips=int(np.sum(np.sign(vol) != np.repeat(np.sign(vpar), nch))))
    return info


def node_task(task):
    """Worker: write coordinates of new nodes [s, e) to a part file."""
    _, s, e, path = task
    fmt = 'No[%d] = p(%.18f,%.18f,%.18f);' + S['nl_nodes']
    X = S['new_xyz']
    with open(path, 'w', newline='', encoding='latin-1') as f:
        for a in range(s, e, 200000):
            b = min(a + 200000, e)
            vals = []
            for i, (x, y, z) in zip(range(S['base_e'] + a, S['base_e'] + b), X[a:b].tolist()):
                vals += (i, x, y, z)
            f.write((fmt * (b - a)) % tuple(vals))
    return {'kw': 'nodes', 'n': e - s}


def run_task(task):
    return element_task(task) if task[0] == 'elem' else node_task(task)


# ---------------------------- process pool ----------------------------------
def _attach_worker(small, spec):
    """'spawn' worker initialiser: attach the shared-memory arrays and rebuild
    the refinement template and the lookup tables."""
    from multiprocessing import shared_memory
    S.clear()
    S.update(small)
    for name, (shm_name, shape, dtype) in spec.items():
        try:
            shm = shared_memory.SharedMemory(name=shm_name, track=False)      # Python >= 3.13
        except TypeError:
            shm = shared_memory.SharedMemory(name=shm_name)   # same tracker as the parent
        _SHM.append(shm)
        S[name] = np.ndarray(shape, dtype=np.dtype(dtype), buffer=shm.buf)
    mode, level = S['T_args']
    S['T'] = template_levels(level) if mode == 'levels' else template_ps()
    for tab in ('edges', 'faces', 'tets'):
        S[tab] = RowTable(keys=S[tab + '_keys'], width=S[tab + '_width']) if tab + '_keys' in S else None


def _ping():
    return os.getpid()


def _importable_self():
    """This script as a regular importable module (needed by 'spawn' workers,
    which cannot see functions defined in an interactive __main__)."""
    if __name__ != '__main__':
        return sys.modules[__name__]
    path = globals().get('__file__')
    if not path:
        return None
    folder, stem = os.path.split(os.path.abspath(path))
    stem = os.path.splitext(stem)[0]
    if not stem.isidentifier():
        return None
    if folder not in sys.path:
        sys.path.insert(0, folder)
    import importlib
    try:
        return importlib.import_module(stem)
    except Exception:
        return None


def run_parallel(tasks, nproc, method):
    """Run tasks dynamically (largest first). Returns the list of task infos.

    Largest-first + imap_unordered(chunksize=1) is the classic way to balance a
    pool: a worker that finishes early immediately picks the next pending task,
    and the small tasks at the end fill the gaps. Results are order-independent
    because every task writes its own part file (concatenated later in order).
    """
    order = sorted(range(len(tasks)), key=lambda i: -tasks[i][-1])          # weight
    tasks_sorted = [tasks[i][:-1] for i in order]
    if nproc <= 1 or len(tasks) <= 1:
        return [run_task(t) for t in tasks_sorted]
    import multiprocessing as mp
    if method == 'auto':
        method = 'fork' if 'fork' in mp.get_all_start_methods() else 'spawn'
    nproc = min(nproc, len(tasks))
    if method == 'fork':
        with mp.get_context('fork').Pool(nproc) as pool:
            return list(pool.imap_unordered(run_task, tasks_sorted, chunksize=1))

    # ---- spawn: big arrays in shared memory, small objects pickled once ----
    mod = _importable_self()
    if mod is None:
        print('  WARNING: script not importable for spawn workers; running serially')
        return [run_task(t) for t in tasks_sorted]
    from multiprocessing import shared_memory
    big, small = {}, {}
    for k, v in S.items():
        if k == 'T':
            continue                                   # rebuilt in the workers
        if isinstance(v, RowTable):
            big[k + '_keys'] = v.keys
            small[k + '_width'] = v.width
        elif isinstance(v, np.ndarray) and v.nbytes > 1 << 16:
            big[k] = v
        elif v is not None:
            small[k] = v
    handles, spec = [], {}
    try:
        for k, a in big.items():
            a = np.ascontiguousarray(a)
            shm = shared_memory.SharedMemory(create=True, size=max(a.nbytes, 1))
            np.ndarray(a.shape, dtype=a.dtype, buffer=shm.buf)[...] = a
            handles.append(shm)
            spec[k] = (shm.name, a.shape, a.dtype.str)
        pool = mp.get_context('spawn').Pool(nproc, initializer=mod._attach_worker,
                                            initargs=(small, spec))
        try:
            pool.apply_async(mod._ping).get(timeout=300)      # workers really started?
        except Exception as err:
            pool.terminate()
            print('  WARNING: spawn workers failed to start (%s); running serially'
                  % type(err).__name__)
            return [run_task(t) for t in tasks_sorted]
        with pool:
            return list(pool.imap_unordered(mod.run_task, tasks_sorted, chunksize=1))
    finally:
        for shm in handles:
            shm.close()
            shm.unlink()


# =============================================================================
#  5. Main steps
# =============================================================================
def find_problem_name(directory, name):
    if name:
        return name
    cands = sorted(f[:-6] for f in os.listdir(directory) if f.endswith('-1.dat')
                   and os.path.exists(os.path.join(directory, f[:-6] + '-4.dat')))
    if len(cands) != 1:
        sys.exit('Cannot determine the problem name in %s (found %s). Set PROBLEM_NAME / --name.'
                 % (directory, cands))
    return cands[0]


def load_problem(d, name):
    """Read nodes and all mesh-dependent files."""
    header, node_lines, ids, xyz, nl_nodes = read_nodes(os.path.join(d, name + '-1.dat'))
    if len(np.unique(ids)) != len(ids):
        sys.exit('Duplicated node ids in %s-1.dat' % name)
    max_id = int(ids.max())
    X = np.full((max_id + 1, 3), np.nan)
    X[ids] = xyz
    files, blocks = {}, []
    for fn in sorted(os.listdir(d), key=lambda f: (len(f), f)):
        m = re.match(re.escape(name) + r'-(\d+)\.dat$', fn)
        if not m or int(m.group(1)) not in FILE_CONTENT:
            continue
        segs, nl = parse_dat(os.path.join(d, fn), FILE_CONTENT[int(m.group(1))])
        if not any(sg['type'] != 'text' for sg in segs):
            continue                                  # only comments: nothing to do
        files[fn] = (segs, nl)
        for sg in segs:
            if sg['type'] == 'elem':
                sg['gid'], sg['nl'], sg['file'] = len(blocks), nl, fn
                n = sg['nodes']
                if n.size and (n.min() < 1 or n.max() > max_id or np.isnan(X[n]).any()):
                    sys.exit('%s: block %s references undefined nodes' % (fn, sg['kw']))
                blocks.append(sg)
    return dict(header=header, node_lines=node_lines, n_nodes=len(ids), max_id=max_id,
                X=X, nl_nodes=nl_nodes, files=files, blocks=blocks)


def build_tables(blocks):
    """Global tables of unique tetrahedra, faces and edges (all sorted rows)."""
    tet_rows, tri_rows = [], []
    for sg in blocks:
        if sg['kind'] == 'tet':
            tet_rows.append(np.sort(sg['nodes'], axis=1))
        elif sg['kind'] == 'tri':
            tri_rows.append(np.sort(sg['nodes'], axis=1))
        else:
            tri_rows += [np.sort(sg['nodes'][:, :3], axis=1), np.sort(sg['nodes'][:, 3:], axis=1)]
    if not tet_rows:
        sys.exit('No tetrahedra found (VE elements in -4.dat)')
    tets = RowTable(np.concatenate(tet_rows))
    t = tets.rows
    # surface triangles are included in case some are not faces of a tetrahedron
    faces = RowTable(np.concatenate([t[:, [0, 1, 2]], t[:, [0, 1, 3]], t[:, [0, 2, 3]], t[:, [1, 2, 3]]]
                                    + tri_rows))
    f = faces.rows
    edges = RowTable(np.concatenate([f[:, [0, 1]], f[:, [0, 2]], f[:, [1, 2]]]))
    return tets, faces, edges, tri_rows


def build_new_nodes(T, X, max_id, tets, faces, edges):
    """Pattern lookup tables, id offsets and coordinates of the new nodes."""
    D = T.D
    e_pat, f_pat = set(), set()
    for w in [tuple(w) for w in T.W4] + [tuple(w) for w in T.W3]:
        nz = [x for x in w if x]
        if len(nz) == 2:
            e_pat.update(nz)                          # pattern = weight of the lower id
        elif len(nz) == 3:
            f_pat.update(itertools.permutations(nz))  # pattern = weights in sorted order
    e_pat, f_pat = sorted(e_pat), sorted(f_pat)
    i_pat = [tuple(w) for w in T.W4 if np.count_nonzero(w) == 4]

    e_lut = np.full(D + 1, -1, dtype=np.int64)
    f_lut = np.full((D + 1, D + 1), -1, dtype=np.int64)
    i_lut = np.full((D + 1, D + 1, D + 1), -1, dtype=np.int64)
    for i, w in enumerate(e_pat):
        e_lut[w] = i
    for i, w in enumerate(f_pat):
        f_lut[w[0], w[1]] = i
    for i, w in enumerate(i_pat):
        i_lut[w[0], w[1], w[2]] = i

    # New node numbering (continues after the highest original id), so the
    # original nodes never change their number:
    #   [ original ids | edge nodes | face nodes | interior nodes ]
    #   id(edge node) = base_e + edge_index * n_ep + pattern
    #   id(face node) = base_f + face_index * n_fp + pattern
    #   id(tet  node) = base_i + tet_index  * n_ip + pattern
    # Edge/face/tet indices come from the global tables, which is what makes the
    # refined mesh conforming (two elements sharing an edge/face share its nodes).
    n_ep, n_fp, n_ip = len(e_pat), len(f_pat), len(i_pat)
    base_e = max_id + 1
    base_f = base_e + edges.n * n_ep
    base_i = base_f + faces.n * n_fp

    # coordinates, in the same order as the ids:  base + entity*n_pat + pattern
    parts = []
    er, fr, tr = edges.rows, faces.rows, tets.rows
    if n_ep:
        parts.append(np.stack([(w * X[er[:, 0]] + (D - w) * X[er[:, 1]]) / D for w in e_pat], 1))
    if n_fp:
        parts.append(np.stack([(a * X[fr[:, 0]] + b * X[fr[:, 1]] + c * X[fr[:, 2]]) / D
                               for a, b, c in f_pat], 1))
    if n_ip:
        parts.append(np.stack([sum(w[i] * X[tr[:, i]] for i in range(4)) / D for w in i_pat], 1))
    new_xyz = np.concatenate([p.reshape(-1, 3) for p in parts]) if parts else np.empty((0, 3))
    return dict(e_lut=e_lut, f_lut=f_lut, i_lut=i_lut, n_ep=n_ep, n_fp=n_fp, n_ip=n_ip,
                base_e=base_e, base_f=base_f, base_i=base_i, new_xyz=new_xyz,
                n_new=(edges.n * n_ep, faces.n * n_fp, tets.n * n_ip))


def propagate_node_conditions(files, blocks, faces, edges, tri_rows, ids, rule):
    """New lines for -2.dat / -3.dat. Returns {(file, segment index): [lines]}.

    A new node (on an edge or a face) inherits the condition of its parent nodes
    only when ALL of them carry the same condition text (e.g. 'V.Fix(0.0)').
    With rule 'surface' it must also lie on a surface face/edge (domain boundary,
    material interface or listed surface element) so that, for instance, a
    corner element does not get a fixed voltage through its interior edge.
    """
    out = {}
    nodeconds = [(fn, k, sg) for fn, (segs, _) in files.items()
                 for k, sg in enumerate(segs) if sg['type'] == 'nodecond']
    if not nodeconds:
        return out
    fr, er = faces.rows, edges.rows
    face_ok = edge_ok = None
    if rule == 'surface':
        # surface faces = mesh boundary + material interfaces + listed surface elements
        face_ok = np.zeros(faces.n, bool)
        vb = [sg for sg in blocks if sg['kw'] == 'VE']
        if vb:
            ve = np.concatenate([np.sort(sg['nodes'], 1) for sg in vb])
            mat = np.concatenate([sg['extra'][:, -1] for sg in vb])
            u = np.unique(np.column_stack([ve, mat]), axis=0)       # drop repeated tets
            ve, mat = u[:, :4], u[:, 4]
            fidx = np.concatenate([faces.lookup(ve[:, c]) for c in
                                   ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])])
            fmat = np.tile(mat, 4)
            cnt = np.bincount(fidx, minlength=faces.n)
            mmin = np.full(faces.n, np.iinfo(np.int64).max)
            mmax = np.full(faces.n, np.iinfo(np.int64).min)
            np.minimum.at(mmin, fidx, fmat)
            np.maximum.at(mmax, fidx, fmat)
            face_ok |= (cnt == 1) | ((cnt >= 2) & (mmin != mmax))
        for t in tri_rows:
            face_ok[faces.lookup(t)] = True
        sf = fr[face_ok]
        edge_ok = np.zeros(edges.n, bool)
        for c in ([0, 1], [0, 2], [1, 2]):
            edge_ok[edges.lookup(sf[:, c])] = True

    for fn, k, sg in nodeconds:
        # condition id per node (identical condition strings -> same id)
        cid = np.full(ids['max_id'] + 1, -1, dtype=np.int64)
        codes = {}
        for i, rest in zip(sg['ids'], sg['rest']):
            cid[i] = codes.setdefault(rest, len(codes))
        text_of = {v: r for r, v in codes.items()}
        new_ids, new_c = [], []
        if ids['n_ep']:
            ca = cid[er[:, 0]]
            m = (ca >= 0) & (ca == cid[er[:, 1]])
            if edge_ok is not None:
                m &= edge_ok
            sel = np.nonzero(m)[0]
            for p in range(ids['n_ep']):
                new_ids.append(ids['base_e'] + sel * ids['n_ep'] + p)
                new_c.append(ca[sel])
        if ids['n_fp']:
            ca = cid[fr[:, 0]]
            m = (ca >= 0) & (ca == cid[fr[:, 1]]) & (ca == cid[fr[:, 2]])
            if face_ok is not None:
                m &= face_ok
            sel = np.nonzero(m)[0]
            for p in range(ids['n_fp']):
                new_ids.append(ids['base_f'] + sel * ids['n_fp'] + p)
                new_c.append(ca[sel])
        if new_ids:
            ni, nc = np.concatenate(new_ids), np.concatenate(new_c)
            o = np.argsort(ni)
            out[(fn, k)] = ['No[%d]%s' % (i, text_of[c]) for i, c in zip(ni[o].tolist(), nc[o].tolist())]
        print('  %-34s %d node conditions -> +%d on new nodes'
              % (fn, len(sg['ids']), len(out.get((fn, k), []))))
    return out


def append_file(fo, path):
    """Append file 'path' to the open binary file 'fo', then delete it.
    Uses an in-kernel copy (Linux copy_file_range) when available."""
    fo.flush()
    with open(path, 'rb') as fi:
        size, done = os.fstat(fi.fileno()).st_size, 0
        if hasattr(os, 'copy_file_range'):
            try:
                while done < size:
                    n = os.copy_file_range(fi.fileno(), fo.fileno(), size - done)
                    if n <= 0:
                        break
                    done += n
            except OSError:
                pass
        fo.seek(0, os.SEEK_END)
        if done < size:
            fi.seek(done)
            shutil.copyfileobj(fi, fo, 16 << 20)
    os.remove(path)


def write_outputs(P, name, outdir, part_of, node_parts, extra_nodecond):
    """Assemble part files into the final .dat files (atomic replace).

    Every file is first written as <file>.tmp and then moved over the original
    with os.replace, so an interrupted run never leaves a half-written .dat file
    (but, when overwriting in place, files already replaced stay refined: use
    --out or keep a backup if you need to be able to roll back).
    Untouched lines (comments, other keywords) are copied verbatim; the new
    node conditions are appended at the end of their original block.
    """
    nl = P['nl_nodes']
    dst = os.path.join(outdir, name + '-1.dat')
    with open(dst + '.tmp', 'wb') as fo:
        fo.write(''.join(ln + nl for ln in P['header'] + P['node_lines']).encode('latin-1'))
        for p in node_parts:
            append_file(fo, p)
    os.replace(dst + '.tmp', dst)

    for fn, (segs, nl) in P['files'].items():
        dst = os.path.join(outdir, fn)
        with open(dst + '.tmp', 'wb') as fo:
            for k, sg in enumerate(segs):
                if sg['type'] == 'text':
                    fo.write((sg['line'] + nl).encode('latin-1'))
                elif sg['type'] == 'elem':
                    for p in part_of[sg['gid']]:
                        append_file(fo, p)
                else:
                    lines = sg['lines'] + extra_nodecond.get((fn, k), [])
                    fo.write(''.join(ln + nl for ln in lines).encode('latin-1'))
        os.replace(dst + '.tmp', dst)


def print_stats(infos):
    st = [i for i in infos if 'hist' in i]
    if not st:
        return
    h, hp = sum(i['hist'] for i in st), sum(i['hist_par'] for i in st)
    print('  Minimum dihedral angle: original %.2f deg -> refined %.2f deg'
          % (min(i['min_par'] for i in st), min(i['min'] for i in st)))
    print('  Histogram of element min dihedral angle (% of elements):')
    print('      bin [deg]   original   refined')
    for b in range(18):
        if hp[b] or h[b]:
            print('      %2d - %2d    %7.2f   %7.2f' % (5 * b, 5 * b + 5, 100. * hp[b] / hp.sum(),
                                                       100. * h[b] / h.sum()))
    print('  Total volume: original %.12e  refined %.12e ; orientation flips: %d'
          % (sum(i['vol_par'] for i in st), sum(i['vol'] for i in st), sum(i['flips'] for i in st)))


def check_mesh(outdir, name):
    """Read the written mesh back and verify conformity."""
    print('  CHECK: reading the refined mesh back ...')
    _, _, ids, xyz, _ = read_nodes(os.path.join(outdir, name + '-1.dat'))
    print('    duplicated node ids : %d' % (len(ids) - len(np.unique(ids))))
    q = np.round(xyz / (np.abs(xyz).max() + 1e-300) * 1e12).astype('>i8')
    print('    coincident nodes    : %d' % (len(q) - len(np.unique(np.ascontiguousarray(q).view('V24').ravel()))))
    segs, _ = parse_dat(os.path.join(outdir, name + '-4.dat'), 'elem')
    ve = np.concatenate([s['nodes'] for s in segs if s['type'] == 'elem' and s['kw'] == 'VE'])
    unused = len(ids) - len(np.unique(ve))
    f = np.sort(np.concatenate([ve[:, [0, 1, 2]], ve[:, [0, 1, 3]], ve[:, [0, 2, 3]], ve[:, [1, 2, 3]]]), 1)
    del ve
    ft = RowTable(f)
    cnt = np.bincount(ft.lookup(f))
    del f
    print('    nodes not used by VE: %d' % unused)
    print('    faces shared by 2 = %d, boundary = %d, non-manifold (>2) = %d'
          % ((cnt == 2).sum(), (cnt == 1).sum(), (cnt > 2).sum()))
    for fn in sorted(os.listdir(outdir), key=lambda x: (len(x), x)):
        m = re.match(re.escape(name) + r'-(\d+)\.dat$', fn)
        if not m or FILE_CONTENT.get(int(m.group(1))) not in ('elem', 'plain'):
            continue
        for s in parse_dat(os.path.join(outdir, fn), FILE_CONTENT[int(m.group(1))])[0]:
            if s['type'] == 'elem' and s['kind'] in ('tri', 'prism'):
                tris = [s['nodes'][:, :3]] + ([s['nodes'][:, 3:]] if s['kind'] == 'prism' else [])
                for t in tris:
                    try:
                        ft.lookup(np.sort(t, 1))
                        msg = 'all are faces of tetrahedra'
                    except RuntimeError:
                        msg = 'ERROR: some are not faces of tetrahedra'
                    print('    %-30s %-5s %9d triangles: %s' % (fn, s['kw'], len(t), msg))


# =============================================================================
#  6. Driver
# =============================================================================
def get_options(argv=None):
    """Command-line options. The defaults are the USER SETTINGS at the top of
    this file; any option given on the command line overrides its setting.
    The [SETTING] tag in each help text names the corresponding variable."""
    prog = os.path.basename(sys.argv[0]) if sys.argv and sys.argv[0] else 'mesh_refiner.py'
    ap = argparse.ArgumentParser(
        prog=prog,                       # whatever name the script has (mesh_refiner_v4.py ...)
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Uniform refinement of ERMES 20.0 meshes (.dat files).\n'
                    'Defaults come from the USER SETTINGS block of this script; options given\n'
                    'here override them. The name in [ ] is the corresponding USER SETTING.',
        epilog='examples:\n'
               '  python ' + prog + ' ./static.gid --level 3 --np 12 --no-backup\n'
               '  python ' + prog + ' ./static.gid --mode ps --out ./static_ps.gid\n'
               '  python ' + prog + ' ./static.gid --level 2 --tmp /tmp --stats --check\n'
               '  %run ' + prog + ' ./static.gid --level 2        (IPython / Spyder)')
    ap.add_argument('directory', nargs='?', default=PROBLEM_DIR, metavar='PROBLEM_DIR',
                    help='folder with the ERMES .dat files [PROBLEM_DIR] (default: %(default)s)')
    ap.add_argument('--name', default=PROBLEM_NAME,
                    help='problem name = prefix of the .dat files [PROBLEM_NAME] '
                         '(default: auto-detect)')
    ap.add_argument('--mode', choices=['levels', 'ps'], default=MODE,
                    help="'levels': edges divided in LEVEL parts; 'ps': Powell-Sabin x24 "
                         '[MODE] (default: %(default)s)')
    ap.add_argument('--level', type=int, default=LEVEL, metavar='L',
                    help='edge divisions for --mode levels, >= 2; tets x L^3, triangles x L^2 '
                         '[LEVEL] (default: %(default)s)')
    ap.add_argument('--np', type=int, default=NPROC or os.cpu_count() or 1, metavar='N',
                    help='number of parallel processes, 1 = serial [NPROC] (default: %(default)s)')
    ap.add_argument('--start', choices=['auto', 'fork', 'spawn'], default=START_METHOD,
                    help='process start method; auto = fork on Linux/WSL, spawn on Windows '
                         '[START_METHOD] (default: %(default)s)')
    ap.add_argument('--out', default=OUT_DIR, metavar='DIR',
                    help='write a refined copy of the problem to DIR instead of overwriting '
                         '[OUT_DIR] (default: overwrite)')
    ap.add_argument('--tmp', default=TMP_DIR, metavar='DIR',
                    help="folder for temporary part files, e.g. /tmp under WSL "
                         '[TMP_DIR] (default: output folder)')
    ap.add_argument('--backup', dest='backup', action='store_true', default=BACKUP,
                    help='keep the original files in mesh_backup_original [BACKUP=True] '
                         '(default: %(default)s)')
    ap.add_argument('--no-backup', dest='backup', action='store_false',
                    help='do not keep the original files [BACKUP=False]')
    ap.add_argument('--nodecond-rule', choices=['surface', 'loose'], default=NODECOND_RULE,
                    help='propagation of voltages (-2.dat) / singular nodes (-3.dat) to new '
                         'nodes [NODECOND_RULE] (default: %(default)s)')
    ap.add_argument('--stats', dest='stats', action='store_true', default=STATS,
                    help='print mesh-quality statistics [STATS=True] (default: %(default)s)')
    ap.add_argument('--no-stats', dest='stats', action='store_false',
                    help='do not print statistics [STATS=False]')
    ap.add_argument('--check', dest='check', action='store_true', default=CHECK,
                    help='read the result back and check conformity, memory heavy '
                         '[CHECK=True] (default: %(default)s)')
    ap.add_argument('--no-check', dest='check', action='store_false',
                    help='skip the conformity check [CHECK=False]')
    ap.add_argument('--task-size', type=int, default=TASK_SIZE, metavar='N',
                    help='child elements per parallel task, 0 = automatic '
                         '[TASK_SIZE] (default: %(default)s)')
    opt = ap.parse_args(argv)
    # basic validation (fail early with a clear message instead of deep in a worker)
    if opt.np < 1:
        ap.error('--np must be >= 1')
    if opt.task_size < 0:
        ap.error('--task-size must be >= 0 (0 = automatic)')
    return opt


def is_wsl():
    """True when running inside Windows Subsystem for Linux."""
    try:
        return 'microsoft' in os.uname().release.lower()
    except AttributeError:                       # os.uname does not exist on native Windows
        return False


def main(argv=None):
    # Line-buffered stdout: the progress lines appear immediately even when the
    # output is redirected to a file or a pipe (e.g. from RefineExecute.sh).
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):         # old Python / replaced stdout (IPython)
        pass
    opt = get_options(argv)
    t0 = time.time()

    def clock():
        return '[%6.1fs]' % (time.time() - t0)

    d = os.path.abspath(opt.directory)
    if not os.path.isdir(d) and not os.path.isabs(opt.directory):
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), opt.directory)
    if not os.path.isdir(d):
        sys.exit('Problem folder not found: %s' % d)
    name = find_problem_name(d, opt.name)
    if opt.mode == 'levels':
        if opt.level < 2:
            sys.exit('LEVEL must be >= 2 (LEVEL = 1 leaves the mesh unchanged)')
        T = template_levels(opt.level)
        label = 'levels, L = %d' % opt.level
    else:
        T = template_ps()
        label = 'Powell-Sabin (x24)'
    outdir = os.path.abspath(opt.out) if opt.out else d
    if is_wsl() and not opt.tmp and (d.startswith('/mnt/') or outdir.startswith('/mnt/')):
        print('  HINT: the problem is on a Windows drive (/mnt/...) under WSL; '
              'adding "--tmp /tmp" makes the run much faster.')
    print('ERMES refine | problem "%s" | %s | %d process(es), start=%s'
          % (name, label, opt.np, opt.start))
    print('  input : %s\n  output: %s%s' % (d, outdir, '  (overwrite)' if outdir == d else ''))

    # ---- read ---------------------------------------------------------------
    P = load_problem(d, name)
    print('  nodes: %d (max id %d) %s' % (P['n_nodes'], P['max_id'], clock()))
    for sg in P['blocks']:
        print('  %-34s %-5s %10d %s' % (sg['file'], sg['kw'], len(sg['nodes']), sg['kind']))

    # ---- topology and new nodes ---------------------------------------------
    tets, faces, edges, tri_rows = build_tables(P['blocks'])
    print('  unique tetrahedra %d, faces %d, edges %d %s' % (tets.n, faces.n, edges.n, clock()))
    ids = build_new_nodes(T, P['X'], P['max_id'], tets, faces, edges)
    ids['max_id'] = P['max_id']
    ne, nf, ni = ids['n_new']
    n_new = ne + nf + ni
    print('  new nodes: %d on edges + %d on faces + %d interior -> %d nodes in total %s'
          % (ne, nf, ni, P['n_nodes'] + n_new, clock()))

    # ---- shared state for workers -------------------------------------------
    # X: original nodes suffice to choose diagonals; all nodes only for STATS
    # S is the module-level dict that the workers read. With 'fork' it is simply
    # inherited by the child processes (copy-on-write, no pickling); with 'spawn'
    # run_parallel() ships it through shared memory (see _attach_worker).
    S.clear()
    S.update({k: v for k, v in ids.items() if k != 'n_new'})
    S.update(T=T, T_args=(opt.mode, opt.level), edges=edges, faces=faces, tets=tets if ids['n_ip'] else None,
             stats=opt.stats, nl_nodes=P['nl_nodes'],
             X=np.vstack([P['X'], ids['new_xyz']]) if opt.stats else P['X'],
             meta=[{k: sg[k] for k in ('kw', 'style', 'kind', 'trail', 'nl')} for sg in P['blocks']])
    for sg in P['blocks']:
        S['b%d_nodes' % sg['gid']] = sg['nodes']
        S['b%d_extra' % sg['gid']] = sg['extra']

    os.makedirs(outdir, exist_ok=True)
    tmp_base = os.path.abspath(opt.tmp) if opt.tmp else outdir
    os.makedirs(tmp_base, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix='ermes_refine_', dir=tmp_base)
    try:
        # ---- balanced tasks: equal number of children per task ----------------
        total = sum(len(sg['nodes']) * (T.n_tet_children if sg['kind'] == 'tet' else T.n_tri_children)
                    for sg in P['blocks']) + n_new
        target = opt.task_size if opt.task_size > 0 else \
            int(min(2e6, max(5e4, total / (4.0 * max(opt.np, 1)))))
        tasks, part_of = [], {}
        for sg in P['blocks']:
            n = len(sg['nodes'])
            nch = T.n_tet_children if sg['kind'] == 'tet' else T.n_tri_children
            chunk = max(1, target // nch)
            part_of[sg['gid']] = []
            for s in range(0, max(n, 1), chunk):
                e = min(s + chunk, n)
                p = os.path.join(tmpdir, 'blk%04d_%010d.part' % (sg['gid'], s))
                tasks.append(('elem', sg['gid'], s, e, p, (e - s) * nch))
                part_of[sg['gid']].append(p)
        node_parts = []
        nstep = max(1, target)
        for s in range(0, n_new, nstep):
            e = min(s + nstep, n_new)
            p = os.path.join(tmpdir, 'nodes_%010d.part' % s)
            tasks.append(('node', s, e, p, 2 * (e - s)))        # ~2x cost per node line
            node_parts.append(p)
        print('  %d parallel tasks (~%d children each) %s' % (len(tasks), target, clock()))
        infos = run_parallel(tasks, opt.np, opt.start)
        print('  elements and nodes refined/formatted %s' % clock())

        # ---- node conditions --------------------------------------------------
        extra_nc = propagate_node_conditions(P['files'], P['blocks'], faces, edges, tri_rows,
                                             ids, opt.nodecond_rule)

        # ---- backup and write -------------------------------------------------
        if outdir == d and opt.backup:
            bdir = os.path.join(d, 'mesh_backup_original')
            os.makedirs(bdir, exist_ok=True)
            for fn in [name + '-1.dat'] + list(P['files']):
                if not os.path.exists(os.path.join(bdir, fn)):
                    shutil.copy2(os.path.join(d, fn), os.path.join(bdir, fn))
            print('  original files saved in %s' % bdir)
        write_outputs(P, name, outdir, part_of, node_parts, extra_nc)
        if outdir != d:            # complete the problem folder with untouched files
            for fn in os.listdir(d):
                if fn.startswith(name) and fn.endswith('.dat') and not os.path.exists(os.path.join(outdir, fn)):
                    shutil.copy2(os.path.join(d, fn), os.path.join(outdir, fn))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ---- summary --------------------------------------------------------------
    counts = {}
    for i in infos:
        if i['kw'] != 'nodes':
            counts[i['kw']] = counts.get(i['kw'], 0) + i['n']
    print('  files written %s' % clock())
    print('  refined elements: ' + ', '.join('%s %d' % (k, counts[k]) for k in
                                             dict.fromkeys(sg['kw'] for sg in P['blocks'])))
    if opt.stats:
        print_stats(infos)
    if opt.check:
        check_mesh(outdir, name)
    print('Done in %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
