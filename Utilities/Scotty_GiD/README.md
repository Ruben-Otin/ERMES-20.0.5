# Scotty / Gaussian beam to GiD scripts (28 Oct 2026)

Tools that turn a microwave Gaussian beam (traced by Scotty in a plasma, or
propagated in vacuum) into geometry that GiD can import, plus a closed FEM
volume around the beam for ERMES, the ERMES launch data, an O-X conversion
estimate and an ERMES cost estimate.

Keep all scripts in the same folder: they import each other (plain names
first; `e_` names are also accepted).

| Script | What it does |
|---|---|
| `scotty_to_gid.py` | Runs Scotty (or reads its `.h5`), exports the beam surface (STL / NURBS IGES), FEM envelope, central ray, E/k fields, ERMES launch data, O-X diagnostics and ERMES FEM cost. Also holds the shared functions. |
| `beam_to_gid.py` | Paraxial Gaussian beam in vacuum (no Scotty, no plasma), same exports; straight-line O-X estimate. Its 3D plot shows the FEM envelope too. |
| `reflect_beam_to_gid.py` | Beam from either script above, reflected on PEC planes / spheres / cylinders / GiD surfaces, with the FEM envelope around all legs. Beam settings come from section 1 of the source script. |
| `gid_surface_readers.py` | Reads GiD exports (IGES, STEP, DXF, VTK, GiD mesh, Exodus, STL, ...) as triangles; used by `reflect_beam_to_gid.py` for mesh / CAD surfaces. |
| `eqdsk_to_scotty.py` | G-EQDSK + profiles -> Scotty input (`topfile`, `topfile.json`, `ne` / `Te` files). Its files are also read by the O-X diagnostics, and its limiter by `reflect_beam_to_gid.py` for `eqdsk:` contours. |

## Requirements

| Needed for | Packages |
|---|---|
| everything | Python >= 3.9, numpy |
| plots | matplotlib |
| O-X diagnostics, `eqdsk_to_scotty.py`, clean trimmed CAD edges | scipy |
| `scotty_to_gid.py` (and `reflect_beam_to_gid.py --beam-source scotty`) | h5py, plus [Scotty](https://github.com/beam-tracing/Scotty) to run it (not needed to export an existing `.h5`) |
| STEP files / IGES through gmsh | gmsh |
| `.vtu` and other meshio formats | meshio |
| Rhino `.3dm` | rhino3dm |

`beam_to_gid.py` and `reflect_beam_to_gid.py --beam-source vacuum` need
neither Scotty nor h5py.

## Folder layout

```
<folder>/
  scotty_to_gid.py  beam_to_gid.py  reflect_beam_to_gid.py
  gid_surface_readers.py  eqdsk_to_scotty.py  README.md
  Plasma/            input: G-EQDSK and profiles (sample files included)
  Scotty_input/      written by eqdsk_to_scotty.py, read by Scotty and the O-X diagnostics
  Scotty_results/    Scotty .h5 output
  GiD_input/         every file written for GiD / ERMES (out_dir)
```

Run the scripts from this folder: the default paths are relative to it.

## Parameters

All parameters are in section 1 ("USER PARAMETERS") of each script, with
their possible values in a comment. Any of them can be changed for one run
from the command line with `--name value` (same name, `-` instead of `_`;
`--name / --no-name` for True/False):

```
python beam_to_gid.py --show-params        # values in use
python beam_to_gid.py --help               # all options with defaults
python reflect_beam_to_gid.py --source-help   # options of its beam script
```

`reflect_beam_to_gid.py` passes any option it does not know to the beam
script (`scotty_to_gid.py` or `beam_to_gid.py`). Its `--plane`, `--sphere`
and `--cylinder` surfaces are ADDED to the `surfaces` list of its section 1;
there is no command-line option to remove those, so set `surfaces = []` in
section 1 to reflect only on the surfaces given on the command line.

The plots show the beam surface through all cross-sections
(`plot_points_around` = points around each one). The window must be closed
for the script to continue; use `--no-show-plot` to skip it.

## Typical workflow

1. **Plasma input** (only for Scotty or the O-X diagnostics):
   `python eqdsk_to_scotty.py` (edit its USER PARAMETERS block: input folder,
   EQDSK and profile files, `SUFFIX`). Writes `Scotty_input/topfile<suffix>`
   and `ne<suffix>.dat` (plus check plots `*_check.png`), and `Te<suffix>.dat`
   only if `TE_FILE` is set. Set `WRITE_JSON = True` if you use
   `find_B_method = "omfit"`.
2. **Beam**:
   * in plasma: `python scotty_to_gid.py` (runs Scotty, then exports), or
     `python scotty_to_gid.py --h5 Scotty_results/<file>.h5` (export only);
   * in vacuum: `python beam_to_gid.py`.
3. **Reflections** (optional): `python reflect_beam_to_gid.py` with the
   surfaces in its section 1 or on the command line.
4. **GiD**: Files > Import > IGES (`*_envelope.igs`, `*_beam_surface.igs`)
   or STL. Assign the ERMES port to `ENVFIRST`, PEC to `ENVPEC<k>` and the
   absorbing / open boundary to the other faces as your model needs.

### Examples

```
# vacuum beam of section 1, no plot window
python beam_to_gid.py --no-show-plot

# vacuum beam at 34.8 GHz, 1.4 m long (crosses the O cut-off of the sample
# plasma at about 0.75 m; at 175 GHz there is no cut-off)
python beam_to_gid.py --frequency 34.8e9 --length 1.4

# the same beam defined by the ERMES data (identical result)
python beam_to_gid.py --ermes-input

# export an existing Scotty run, 0.6 m vacuum extension after the last point
python scotty_to_gid.py --h5 ./Scotty_results/<scotty_output_file>.h5 --extend-last 0.6

# The next three ADD their surface to the section 1 'surfaces'. With the
# current section 1 (cylinder of radius 1 m) the beam hits that cylinder
# first (at 1.175 m); set surfaces = [] in section 1 to see only the new one.

# vacuum beam reflected on a flat mirror at X = 0.8 m facing +X
python reflect_beam_to_gid.py --beam-source vacuum --length 2 --plane 0.8 0 0 1 0 0

# reflection on a sphere of radius 1 m around the axis (convex, seen from outside), no plot window
python reflect_beam_to_gid.py --beam-source vacuum --length 2 --sphere 0 0 0.03 1.0 --no-show-plot

# reflection on geometry exported from GiD as IGES (no extra package; Inner.igs is not included)
python reflect_beam_to_gid.py --surface-files ./Inner.igs --surface-mesh-size 0.01
```

## Output files (in `out_dir`, default `./GiD_input/`)

`<stem>` is set by `file_stem`; by default it is built from the `.h5` name
(Scotty) or `vacuum_beam_<freq>GHz_L<length>m` (vacuum), plus `_refl` for
`reflect_beam_to_gid.py`.

| File | Content |
|---|---|
| `<stem>_beam_surface.stl / .igs` | Beam surface (1/e amplitude radius x `envelope_factor`), STL triangles or exact NURBS (IGES 128). `reflect_beam_to_gid.py`: one per leg, `<stem>_leg<k>_beam_surface.*` |
| `<stem>_envelope.stl / .igs` | FEM envelope: closed polyhedron with planar faces, at least `envelope_separation` from the beam surface |
| `<stem>_central_ray.dat` | Central ray, x y z per line (all legs) |
| `Efield_output_data.dat` (vacuum) / `E_final.dat` (Scotty, reflect) | E and k at both beam ends (name set by `efield_file`) |
| `<stem>_Efield_along_ray.dat`, `<stem>_Kvector_along_ray.dat` | E and k at every ray point |
| `<stem>_ermes_launch.txt` | Beam start in ERMES form (point, k0, Epol, \|Epar\|, \|Eper\|, phases), compared with the ERMES polarisation when `ermes_input = True` |
| `<stem>_scotty_launch.txt` | (vacuum) Scotty launch equivalent to the beam: R, zeta, Z, TORBEAM angles, width, curvature |
| `<stem>_beam_parameters.txt` | (vacuum) widths, curvatures, Psi, Gouy phase, waists at both ends |
| `<stem>_OX_along_ray.dat` | N_par, N_y, Mjolhus T along the ray (only with the plasma files) |
| `<stem>_ermes_fem_estimate.txt` / `<stem>_ermes_estimate.txt` (reflect) | Tetrahedra, unknowns, memory and solve time for the envelope |
| `<stem>_reflections.txt` | (reflect) hit points, normals, incidence angles, curvatures |
| `<stem>_surfaces.stl` | (reflect) the reflecting surfaces near the beam, to check them in GiD |

### FEM envelope face labels

| Label | Face | Typical ERMES condition |
|---|---|---|
| `ENVFIRST` | plane of the first beam cross-section (antenna / port) | port (launch the beam here) |
| `ENVLAST` | plane of the last beam cross-section | absorbing |
| `ENVPEC<k>` | on mirror k (`reflect_beam_to_gid.py`) | PEC |
| `ENVSIDE` | all other faces | absorbing |

A flat mirror can give more than one `ENVPEC<k>` face (one per leg
footprint); they lie in the same plane and the volume is still closed.

## Coordinates and conventions

* Cartesian X, Y, Z with Z = tokamak axis, origin on the axis; X, Z ->
  tokamak R, Z in the zeta = 0 plane, Y toroidal.
* Cylindrical (R, zeta, Z): X = R cos(zeta), Y = R sin(zeta).
* TORBEAM angles: (0, 0) = towards -R; poloidal > 0 tilts towards -Z,
  toroidal > 0 towards -zeta.
* Beam widths are 1/e AMPLITUDE radii (= 1/e^2 intensity radii).
  `envelope_factor` 1.0 / 1.5 / 2.0 holds 86.5 / 98.9 / 99.97 % of the power.
* Fields ~ exp(i k s + i w.Psi.w/2), time convention exp(-i w t) (Scotty);
  `time_convention = "+jwt"` only changes the handedness of `rcp` / `lcp`.
* ERMES `Epol` is perpendicular to k0, `Epar` is along Epol and `Eper` along
  k0 x Epol. A negative magnitude is a phase of pi:
  `(-0.78865, -1.5708)` = `(0.78865, +1.5708)`.

## Current configuration (section 1)

* Launch (all scripts): w0 = 17.92 mm at the launch point (curvature 0),
  ERMES data point X,Y,Z = (2.0725, 0.07381, 0.03249) m,
  k0 = (-0.961603, 0.222113, 0.161202), Epol = (0.025167, -0.513531, 0.857702),
  Epar = (0.614842, 0), Eper = (-0.78865, -1.5708). Scotty-style equivalent:
  (R, zeta, Z) = (2.073813918, 0.035598947, 0.03249), TORBEAM angles
  -9.276669 / -15.045861 deg, Jones vector
  (0.53434105+0.3901395j, -0.30415541+0.68539062j) in (x_hat, y_hat).
  `ermes_input = True` in `beam_to_gid.py` and `scotty_to_gid.py`; both
  definitions give the same beam.
* **Frequencies differ**: `beam_to_gid.py` (and so `reflect_beam_to_gid.py`
  with `beam_source = "vacuum"`) uses **175 GHz**, `scotty_to_gid.py` uses
  **34.8 GHz**, O mode (`mode_flag = +1`).
  * At 175 GHz the O cut-off density is ~3.8e20 m^-3, far above the sample
    plasma, so the O-X diagnostics print "the ray never gets near the O
    cut-off" (expected). Rayleigh length 0.59 m (0.117 m at 34.8 GHz).
  * The ERMES FEM estimate at 175 GHz for the 3 m envelope is ~1.9e11
    tetrahedra (~85 TB matrix): not a practical ERMES model at that size.
* `beam_to_gid.py`: `propagation_length = 3.0` m, `envelope_separation = 0.05`,
  `fem_solver = "direct"`.
* `scotty_to_gid.py`: `extend_first_from = "launch"`, `extend_first = -0.3`
  (geometry starts 0.3 m ahead of the antenna, waist stays at the antenna),
  `extend_last = 0.3`, `envelope_separation = 0.1`, `fem_solver = "iterative"`.
* `reflect_beam_to_gid.py`: `beam_source = "vacuum"`, one surface `Column`:
  cylinder of radius 1 m on the Z axis from Z = -2 to Z = 2 m. The central
  ray hits it at (0.942, 0.335, 0.222) m after 1.175 m (incidence 33.7 deg);
  the output has 2 legs and an `ENVPEC1` face.

## Sample data (`Plasma/`)

`plasma_equilibrium.eqdsk` (pyEquilibrium G-EQDSK, 65 x 65 grid,
R = 0.06 .. 2.00 m, Z = -2.2 .. 2.2 m, magnetic axis R = 0.964 m),
`rho_plasma.dat` (rho_pol = sqrt(psi_N), 43 points, 0.019 .. 1.247),
`ne_plasma.dat` (n_e in m^-3), `te_plasma.dat` (T_e in eV). All profile
files are one column. `eqdsk_to_scotty.py` is set up for them
(`SUFFIX = "_sample"`, `NE_UNITS = "m-3"`); `TE_FILE = None`, so
`te_plasma.dat` is not used (set `TE_FILE = "te_plasma.dat"`,
`TE_UNITS = "eV"` if Scotty runs with `relativistic_flag = True`).

* The EQDSK header current is a placeholder (1 A); the current from the psi
  map is -7.37e5 A. `eqdsk_to_scotty.py` warns about it ("no simple factor
  explains the mismatch"). B and psi are not affected.
* The file has a limiter (98 points, R = 0.26 .. 1.73 m, Z = -2.07 .. 2.06 m),
  so `type="revolved", contour="eqdsk:./Plasma/plasma_equilibrium.eqdsk"`
  works in `reflect_beam_to_gid.py` (with the current launch the beam hits it
  at R = 1.66 m after 0.437 m). Note that file names are case-sensitive on
  Linux.
* The launch at R = 2.074 m is outside the EQDSK grid (R <= 2.00 m). The
  straight-line O-X diagnostics handle it; check that Scotty does.

## Notes

* Scotty only handles launches at zeta = 0: the scripts launch it at
  zeta = 0 with the same local angles and rotate its output back by zeta.
* The beam in front of the plasma is built from the launch definition
  (`extend_first_from = "launch"`). Scotty's Psi at its first point is on the
  plasma side of the Psi_BC jump and must not be propagated back through
  vacuum (`"entry"` = old behaviour, gives a spurious elliptical beam).
* `extend_first` < 0 starts the geometry that far AHEAD of the antenna
  (shorter FEM domain); the waist stays at the antenna. The ERMES port is then
  the envelope's first face, with w0 = 17.92 mm and the waist |extend_first|
  behind it (see `<stem>_ermes_launch.txt`).
* For O-X conversion use `extend_last = 0`: past the turning point the beam
  is the reflected wave, and beam tracing is not valid at the cut-off.
* When re-exporting an existing `.h5`, the launch width / curvature must
  match the run (checked; falls back to `"entry"` with a warning otherwise).
* O-X diagnostics read `topfile<suffix>` (TORBEAM ASCII) or
  `topfile<suffix>.json`, and `ne<suffix>.dat`, from `magnetic_data_path` /
  `ne_data_path`. `beam_to_gid.py` follows the STRAIGHT line (no
  refraction): run Scotty for the refracted N_par / N_y.
* Reflections: a hit inside the part traced by Scotty (in the plasma) is
  still propagated in VACUUM after the reflection (a warning is printed).
* Curved mirrors: the `ENVPEC` face follows the real surface
  (`envelope_pec_conform = True`) only if the mirror covers the whole face;
  otherwise it stays on the tangent plane with a warning (reduce
  `envelope_separation`).
* The FEM cost estimate is order-of-magnitude: calibrate the fill / flop
  constants with one real ERMES run.

## Messages you may see

| Message | Meaning / action |
|---|---|
| `NOTE: ... skipped: scotty_to_gid.py not found` | `beam_to_gid.py` needs `scotty_to_gid.py` in the same folder for the ERMES launch file, O-X diagnostics and FEM envelope. |
| `WARNING: O-X diagnostics skipped: no topfile<suffix>(.json) in ...` | Run `eqdsk_to_scotty.py` first, or set `--no-ox-diagnostics`. Harmless otherwise. |
| `WARNING: ... the BEAM crosses surface '...' although the central ray does not hit it` | Part of the beam (1/e x `envelope_factor`) touches a surface the central ray misses; that part is not reflected. Physics warning; shorten the beam or move the surface if it matters. |
| `WARNING: no FEM envelope written: ...` | The envelope could not be built (e.g. too tight a bend); try another `envelope_separation` / `envelope_max_cell_length`. |
| `Ip from Ampere's law ... WARNING: no simple factor explains the mismatch` | `eqdsk_to_scotty.py`: header current inconsistent with psi (see Sample data). |

## Changes (28 Oct 2026)

* Sample data renamed: `plasma_equilibrium.eqdsk`, `rho_plasma.dat`,
  `ne_plasma.dat`, `te_plasma.dat` (new equilibrium, with limiter).
* `scotty_to_gid.py` now uses the same ERMES launch as `beam_to_gid.py`
  (the old R = 1.80 m test launch is gone).
* `scotty_to_gid.py`: beam-averaged O-X transmission was `nan` (with a
  RuntimeWarning) when the ray never reached any density gradient
  (L_n = inf); it is now 0.
* `reflect_beam_to_gid.py`: commented `eqdsk:` example points to
  `plasma_equilibrium.eqdsk` (was `Plasma_Equilibrium.eqdsk`, which fails on
  Linux). `eqdsk_to_scotty.py`: `SUFFIX` comment updated.

## Changes (27 Oct 2026)

* `beam_to_gid.py`: the 3D plot shows the FEM envelope (green first face,
  red last face, orange sides).
* `beam_to_gid.py` and `scotty_to_gid.py`: start point / launch_position
  given to full precision, (2.073813918, 0.035598947, 0.03249); it was
  rounded to (2.074, 0.035, 0.032), i.e. 1.2 mm / 0.5 mm away from the ERMES
  point and with a slightly different direction. `beam_to_gid.py`
  polarisation given to full precision. Now `ermes_input = True / False`
  really give the same beam.
* `scotty_to_gid.py`: `<stem>_ermes_launch.txt` of a vacuum beam no longer
  says "from Scotty ... at the plasma entry".
* Usage text: script names updated (`vacuum_beam_to_gid.py` ->
  `beam_to_gid.py`, `scotty_beam_to_gid.py` -> `scotty_to_gid.py`).
