# Plasma_Scripts: EQDSK → Scotty and EQDSK → ERMES

Scripts to turn an EFIT G-EQDSK equilibrium plus kinetic profiles into input
for the Scotty beam tracer and into ERMES cold-plasma models (MAST-U DBS
benchmark), and to build and run the ERMES model with GiD.

Part of the ERMES 20.0.5 Utilities (section 9 of `Utilities/README.txt`).
Requires ERMES 20.0.5 (problemtype `ERMES_20.0.5`): the GiD batch files use
its solver options and the scripts rely on its exit status and input checks
(see [ERMES 20.0.5](#ermes-2005)).

## Files

| File | Purpose |
|---|---|
| `eqdsk_to_scotty.py` | G-EQDSK + profiles → Scotty input (`topfile<sfx>`, `ne<sfx>.dat`, `Te<sfx>.dat`). Also the **shared EQDSK reader and field calculation** (`read_geqdsk`, `compute_fields`, `load_profile`) used by every other script. |
| `eqdsk_reader.py` | Plots an EQDSK and the profiles (ψ, ψ_N, B_p, B_t, J_p, J_φ, n_e, T_e) to check the inputs. |
| `eqdsk_to_ermes_simp.py` | Minimal EQDSK → ERMES plasma files (`ne_<name>.dat`, `Be_<name>.dat`), no turbulence. |
| `eqdsk_to_ermes_turb.py` | EQDSK → ERMES plasma files with optional turbulence and blobs. Front end of `torbeam_to_ermes.py`; also provides `Generate_Plasma_Files`, `Get_Point_Equilibrium_Data` and `Get_Pcrtc_Equilibrium_Data` for `problem_generator.py`. |
| `torbeam_to_ermes.py` | Engine: mesh reading, fields at the nodes, turbulence, blobs/holes, check plot. Standalone it reads a Scotty/TORBEAM topfile (default) or an EQDSK; every parameter can be set on the command line. |
| `problem_generator.py` | Builds the 2D ERMES model: beam direction and O-mode polarisation, GiD batch file, GiD meshing, plasma files, ERMES run. |
| `utilities.py` | Geometry, GiD batch-file writer (including the ERMES solver choice), O-mode polarisation, GiD node reader, (R, φ, Z) → (x, y, z) transform. |
| `permittivity_tensor.py` | Cold-plasma tensor K at every node → `Tensor_K_CPGE.bin` for the ERMES "Tensor" option. |
| `batch_script.sh` | SLURM job: `problem_generator.py` (GiD + plasma files) then ERMES. |
| `EQDSK_Samples/` | `mast-u-sample.eqdsk`; `rho_flux.txt` (ρ_pol = √ψ_N, 0 → 1.3), `rho_ne.txt` (n_e, m⁻³), `rho_te.txt` (T_e, eV), all 1-column with 2000 points. |

### Dependencies between scripts

```
eqdsk_to_scotty.py  ◄── eqdsk_reader.py
        ▲           ◄── eqdsk_to_ermes_simp.py ──► utilities.py
        │
torbeam_to_ermes.py ◄── eqdsk_to_ermes_turb.py ◄── problem_generator.py ──► utilities.py
```

All scripts must be in the same folder. Requirements: numpy, scipy,
matplotlib (pandas optional, speeds up writing large files).

## Workflows

Run every script from this folder: the default paths are relative to it.

**1. Check the inputs.** `python eqdsk_reader.py`.

**2. Scotty input.** Edit the USER PARAMETERS of `eqdsk_to_scotty.py` and run
it. It writes `Scotty_input/topfile_sample`, `ne_sample.dat` and
`Te_sample.dat`, prints sanity checks (vacuum field, ψ_N on axis, plasma
current from Ampère's law) and saves check plots. Settings for Scotty's
`beam_me_up()` are in the script docstring.

**3. Full ERMES model.** Set the beam, geometry and `Plasma_Options` in
`problem_generator.py` and run it. Steps:

1. `k0` from the toroidal/poloidal angles (`gu.Calculate_WaveVector`);
2. B at the critical density along the beam (`e2e.Get_Pcrtc_Equilibrium_Data`)
   → O-mode polarisation (`gu.Calculate_O_Mode`);
3. GiD batch file `BatchFiles/<Problem_Name>.bch` (`gu.Write_GiD_bch_File_2D`)
   → GiD creates, meshes and writes `<Problem_Name>.gid/`;
4. plasma files `Plasma/ne_<id>.dat`, `Plasma/Be_<id>.dat` at the mesh nodes
   (`e2e.Generate_Plasma_Files`);
5. ERMES, output in `<Problem_Name>.gid/<Problem_Name>.info`.

Names: `Problem_Name = MASTU-2D-t<tor>_p<pol>`, plasma id `2D-t<tor>_p<pol>`.
The GiD/ERMES executables are chosen by platform (Windows or Linux) at the
top of the script; the environment variables `GID_EXE`, `ERMES_EXE` override
them and `RUN_GID=0` / `RUN_ERMES=0` skip those steps. If ERMES stops with an
error (non-zero exit status in 20.0.5) the script raises an exception; the
error message is in the `.info` file.

Solver (written in the GiD problem data by `gu.Write_GiD_bch_File_2D`):

```python
PSettings.Solver   = 'Direct_solver'             # ERMES internal sparse direct solver
# PSettings.Solver = 'External_solver'           # then also:
PSettings.SolvName = 'bash'                      # Windows: 'wsl' with SolvFile = 'bash <path>'
PSettings.SolvFile = '../PETSc/ERMES2PETSc.sh'   # PETSc_Solvers/ERMES2PETSc.sh
```

`Direct_solver` (new in ERMES 20.0.5) needs nothing else and suits the 2D
models of this script; it prints the predicted memory and stops if the
problem does not fit. For large models use `External_solver` with
`Utilities/PETSc_Solvers/ERMES2PETSc.sh` (adapt the relative path to where
that script is).

The same `Plasma_Options` are used for the polarisation and for the plasma
files, so both always see the same field, e.g.

```python
Plasma_Options = dict(toroidal_convention="right-handed",
                      turbulence=True, turb_Lr=0.015, turb_Lb=0.03,
                      blobs=True, blob_packing_fraction=0.2, blob_amplitude=1.5)
```

Any parameter of `torbeam_to_ermes.py` can be given here.

**4. Plasma files for an existing mesh.** Either
`python eqdsk_to_ermes_turb.py` (USER PARAMETERS at the top) or the minimal
`python eqdsk_to_ermes_simp.py`. From a Scotty topfile instead of the EQDSK:
`python torbeam_to_ermes.py --mesh-file <problem>.gid/<problem>-1.dat`
(`--help` lists all parameters, `--show-params` the current values).

**5. Permittivity tensor (optional).** Set the file names, frequency,
damping and ion species in `permittivity_tensor.py` and run it.

**6. Cluster.** `sbatch batch_script.sh`. It runs `problem_generator.py` with
`RUN_ERMES=0` and the cluster GiD path, then runs ERMES itself and stops with
an error message if ERMES fails (so it can be chained with
`sbatch --dependency=afterok`). `PROBLEM_NAME` in the script must match
`Problem_Name`.

## Output files

| File | Content |
|---|---|
| `ne_<name>.dat` | `node  n_e [m⁻³]`, sorted by node id |
| `Be_<name>.dat` | `node  Bx By Bz [T]`, Cartesian, mesh frame |
| `dn_<name>.dat` | `node  δn/n` (turbulence or blobs only) |
| `check_<name>.png` | cut through the mesh (n_e, δn/n, \|B\|) and profiles |

`torbeam_to_ermes.py` standalone writes `ne_3D.dat`, `Be_3D.dat`,
`dn_3D.dat`, `plasma_check.png`.

## Conventions

**Frame.** z is the tokamak axis in all scripts: R = √(x² + y²),
φ = atan2(y, x), (R, φ, Z) right-handed. Lengths in metres. The previous
Plasma_Scripts (ERMES 20.0.4) used y as the tokamak axis; the change is the rotation
(x, y, z)_old → (x, −z, y): the mirror moved from (−1.818389, 0, 0.2846) m to
(−1.818389, −0.2846, 0) m and the GiD geometry is the old one rotated. GiD
problems and plasma files made with the y-axis scripts must be regenerated.

**Equilibrium.** ψ_N = (ψ − ψ_axis)/(ψ_LCFS − ψ_axis); B_R = −(1/R) ∂ψ/∂Z,
B_Z = (1/R) ∂ψ/∂R, B_t = F(ψ)/R with the vacuum F = `fpol[-1]` outside the
LCFS. ψ in Wb/rad (set `PSI_IN_WEBER = True` if the file is in Wb);
`FLIP_BP = True` for the opposite poloidal sign (other COCOS).

**Toroidal field direction.** `toroidal_convention = "right-handed"`
(default) follows the EQDSK sign (B0 = −0.600 T in the MAST-U sample).
`"legacy"` reverses B_t, reproducing the old y-axis runs (whose transform was
left-handed). For the 20°/20° launch the angle between k0 and B at the
critical density is 85° (right-handed) vs 130° (legacy).

**Density.** n_e(ψ_N) from the profile; beyond its last point (ρ = 1.3 in the
samples) `NE_OUTSIDE = "constant"` keeps the last value (default), `"zero"`
sets 0. In `eqdsk_to_ermes_turb.py`, points with ψ_N < 1 outside the LCFS
polygon (private-flux region below the X-point) are treated as outside the
plasma; nodes outside the EQDSK grid get n_e = 0.

**Agreement.** `eqdsk_to_ermes_simp.py` and `eqdsk_to_ermes_turb.py` agree
to ~1e-3 in n_e and ~5e-3 T in B on the MAST-U sample (spline orders and
private-flux handling). `torbeam_to_ermes.py` with the topfile written by
`eqdsk_to_scotty.py` reproduces the EQDSK route to ~1e-8.

## Turbulence and blobs

`turbulence=True` adds a field-aligned random density field: sum of random
Fourier modes in local (radial, binormal, parallel) coordinates, Gaussian
correlation with lengths `turb_Lr`, `turb_Lb`, `turb_Lpar`, rms δn/n vs ρ
from `turb_amplitude_profile`, lognormal or Gaussian statistics, seed
`turb_seed`.

`blobs=True` adds field-aligned filaments: blobs (density excess, born near
the separatrix and carried into the SOL) and holes (relative deficits just
inside the separatrix), with random positions and log-normal sizes and
amplitudes. Main knobs: `blob_packing_fraction`, `blob_amplitude` (× n_ref at
the separatrix), `blob_width_r` / `blob_width_b`, `blob_decay_length`,
`blob_rho_range`, `hole_fraction`, `hole_amplitude`, `blob_seed`.

Full description of both models: docstring of `torbeam_to_ermes.py`.

## ERMES 20.0.5

Changes in ERMES 20.0.5 that affect these scripts (full list:
`RELEASE_NOTES_20.0.5.txt`):

* **Gaussian beam amplitude.** The beam normalization is now the correct one
  for a 3D beam (20.0.4 overestimated the field away from the waist by
  √(w/w0)). Results of `problem_generator.py` models change with respect to
  20.0.4. For a true 2D beam in a thin slab, the waist radius of the invariant
  direction must be very large (e.g. 1e3), and the beam must travel in the
  plane of the slab (it does here: k0 ⟂ slab normal). See the note below.
* **Direct solver.** `PSettings.Solver = 'Direct_solver'` uses the new internal
  solver; `'External_solver'` keeps the PETSc route.
* **Exit status.** ERMES returns a non-zero exit status on error and also
  writes the messages to stderr; `problem_generator.py` (`check=True`) and
  `batch_script.sh` stop and point to the `.info` file.
* **Cold plasma.** With damping (`CPlasma.Damping`), ERMES replaces the
  Hermitian storages by `Full_matrix` automatically. Plasma file paths may
  contain spaces. The ne/Be files written here (sorted node ids, no header)
  pass the new input checks. CODATA 2018 constants and ion masses without the
  removed electrons, as in `permittivity_tensor.py`.
* **`permittivity_tensor.py`.** Off-diagonal sign of the Stix tensor corrected
  (K_xy = −iD; the old sign corresponded to a reversed B). Regenerate old
  `Tensor_K_CPGE.bin` files.

## Notes and known issues

* **2D beam waist.** `problem_generator.py` sets `GBeam.Rpar = GBeam.Rper =
  WaistSize`, i.e. a 3D circular beam evaluated in the slab. To model a 2D
  beam, the radius along the slab normal (`t0` in `Calculate_GiD_Points`,
  the extrusion direction) should be ~1e3. The beam radii follow the
  polarisation axes (Ep and k0 × Ep), and for the 20°/20° launch Ep is 79°
  from the slab normal, so neither radius is along the invariant direction:
  check the ERMES manual (equation 3.11) before changing the radii.
* `ProblemSettings` in `utilities.py` has no default for `Solver`: set
  `PSettings.Solver` (as `problem_generator.py` does) before calling
  `Write_GiD_bch_File_2D`.
* `mast-u-sample.eqdsk` has a corrupted plasma current in the header
  (7.4e11 A); Ampère's law on the boundary gives 9.6e5 A, so the current
  check in `eqdsk_to_scotty.py` prints a warning for this file.
* `utilities.Calculate_O_Mode`: the amplitude is fixed (0.5827, not
  normalised) and the in-plane component diverges as the angle between k0 and
  B approaches 90° (at 85° it is already 2.84).
* `utilities.Write_GiD_bch_File_2D`: the commented line
  `Geometry Create Line Join 1 5` is needed on Linux GiD to close the outline.

## Changes with respect to the previous Plasma_Scripts

* New scripts: `eqdsk_to_ermes_turb.py` (turbulence, blobs, private-flux
  region), `torbeam_to_ermes.py` (engine, topfile or EQDSK input),
  `eqdsk_to_ermes_simp.py`, and `eqdsk_to_scotty.py` (also in `Scotty_GiD`).
  `eqdsk2ermes_plasma_files.py` is replaced by `eqdsk_to_ermes_simp.py`.
* z is now the tokamak axis everywhere (see Conventions) and the toroidal
  field follows the EQDSK sign by default.
* One EQDSK reader (`eqdsk_to_scotty.read_geqdsk`) instead of several copies
  of a line-based reader; `eqdsk_reader.py` rewritten on top of it.
* Removed the plasma-file functions from `utilities.py`
  (`Read_Equilibrium_Files`, `Generate_Plasma_Files`,
  `Get_Point_Equilibrium_Data`, `Get_Pcrtc_Equilibrium_Data`): use
  `eqdsk_to_ermes_turb.py` (same signatures, equilibrium read once).
* `ProblemSettings.Freqency` → `Frequency`; new `ProblemSettings.Solver`
  written to the GiD problem data (`Direct_solver` or `External_solver`).
* `problem_generator.py`: platform-dependent executables,
  `RUN_GID` / `RUN_ERMES` / `GID_EXE` / `ERMES_EXE` environment variables,
  `~` expanded in Linux paths, stops if the mesh is missing or ERMES fails.
  `batch_script.sh` no longer runs ERMES twice.
* `permittivity_tensor.py`: parameters at the top, default files are the ones
  written by `problem_generator.py`.
