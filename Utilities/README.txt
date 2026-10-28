***************************************************************************************************
* ERMES 20.0.5 - Utilities
***************************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - October 2026
***************************************************************************************************

This folder contains auxiliary scripts that complement ERMES 20.0.5. They help to automate
simulations on Linux and HPC clusters, solve the ERMES matrices with external solvers,
generate files that can be imported into ERMES (plasma profiles, current sources, Robin
boundary conditions, element matrices, beam geometry) and reduce the size of the ERMES
results files. It includes the following folders:

- "Batch_Scripts"     : Examples of batch files to run GiD, ERMES and PETSc on Linux/SLURM.
- "Current_J_Scripts" : Generation of imported plasma current sources (J) from EQDSK files.
- "IRBC_Scripts"      : Generation of plane wave Imported Robin Boundary Conditions (IRBC).
- "IVEM_Scripts"      : Generation of imported element stiffness matrices (IVEM).
- "Mesh_Refiner"      : Parallel uniform refinement of the ERMES mesh (.dat files) and run.
- "NumPy_Solvers"     : Examples of Python NumPy/SciPy external solvers for ERMES.
- "PETSc_Python"      : Early Python interface to PETSc ( use for testing ).
- "PETSc_Solvers"     : Parallel C++ PETSc solver that reads the ERMES system files directly.
- "Plasma_Scripts"    : EQDSK readers, cold plasma files generation and problem automation.
- "Results_Filters"   : Keeps only the selected results of an ERMES/GiD results file (.res).
- "Scotty_GiD"        : Scotty beam tracing / vacuum Gaussian beams ( also reflected on PEC
                        surfaces ) exported as GiD geometry, with the FEM volume for ERMES.

All the scripts are examples intended to be copied and adapted to each specific problem.
Paths, file names, problem names and physical parameters are defined at the beginning of
each script and must be edited before use. Python scripts require Python 3 with NumPy and
SciPy (and Matplotlib for "Plasma_Scripts" and the "Scotty_GiD" plots). The "Scotty_GiD"
scripts that run beam tracing also require Scotty (https://github.com/beam-tracing/Scotty).
"Results_Filters" only uses the Python standard library.

***************************************************************************************************
1-) "Batch_Scripts"
***************************************************************************************************

Examples of shell scripts for running ERMES without the GiD graphical interface on Linux
workstations and HPC clusters.

- "ERMES_batch.sh"    : Runs the ERMES executable on an existing GiD problem and then calls
                        "ERMES2PETSc.sh" to solve the system with PETSc. Stops if ERMES or
                        the solver fail ( ERMES returns a non-zero exit status on error ).

- "GiD_batch.sh"      : Full workflow with GiD offscreen: meshes the problem, computes the
                        static currents, stores them in a "Currents" folder, computes the
                        fields, and solves with an external solver. It calls GiD batch
                        (.bch) macro files (Mesh.bch, Calculate_J.bch, Calculate_E.bch).

- "GiD_offscreen.bch" : Example of a GiD batch file. It creates a geometry, loads the ERMES
                        problemtype, defines Gaussian beam Robin coefficients and a cold
                        plasma material, sets the problem data (external PETSc solver),
                        assigns boundary conditions, meshes, and writes the calculation
                        file. Can also be run in GiD with "Import -> Batch file..." (Ctrl-b).

- "Linux_modules.sh"  : Example of Linux module configuration (gcc, cmake, python) for
                        ERMES. Adapt it to the modules available in your system.

- "PETSc_batch.sh"    : Example of a parallel PETSc solver call with petscmpiexec, including
                        solver options (e.g. MUMPS direct solver) and monitoring options.

- "SLURM_batch.sh"    : Example of SLURM job script for HPC clusters. It sets the PETSc
                        environment and runs ERMES in loops over several frequencies,
                        overwriting the ERMES input .dat files and collecting the results
                        into separate folders.

***************************************************************************************************
2-) "Current_J_Scripts"
***************************************************************************************************

Scripts that read an EQDSK plasma equilibrium file and compute the plasma current density
from the equilibrium profiles (p' and FF'). The result is written as an imported J source
file in the "Export_J_Sources" folder of the GiD problem. They require the EQDSK file and
the ERMES mesh files of the problem (nodes "*-1.dat" and volume elements "*-4.dat").

- "eqdsk2ermes_elevol_Jp.py" : J sources defined per tetrahedral volume element. The region
                               with current can be limited by coordinates or by material Id.

- "eqdsk2ermes_nodal_Jp.py"  : J sources defined per mesh node.

Note: by default these scripts read the EQDSK sample file from "./EQDSK_Samples/". Copy that
folder from "Plasma_Scripts" or edit the "EQDSK_File" path.

***************************************************************************************************
3-) "IRBC_Scripts"
***************************************************************************************************

- "IRBC_Plane_Wave_v2.py" : Generates a plane wave Imported Robin Boundary Condition for
                            EDG_1st and RME_1st elements. The user defines the frequency,
                            amplitude, phase, polarization, wave vector and material
                            properties of the plane wave. The script reads the mesh nodes
                            ("*-1.dat") and the IRBC surface elements ("*-19.dat") of the GiD
                            problem and writes the binary files "Vector_U_IRBC.bin" and
                            "Matrix_P_IRBC.bin" to be imported by ERMES.

***************************************************************************************************
4-) "IVEM_Scripts"
***************************************************************************************************

- "IVEM_K_Matrix.py" : Computes the element stiffness matrix of each tetrahedral element for
                       the RME_1st element type, given the material properties (frequency,
                       conductivity, permittivity, permeability). It reads the mesh files
                       ("*-1.dat" and "*-4.dat") and writes the binary file
                       "Matrix_K_IVEM.bin" to be imported by ERMES.

- "IVEM_Utils.py"    : Functions and objects used by "IVEM_K_Matrix.py" (mesh file reader,
                       tetrahedron shape function derivatives and volume, RME_1st element
                       matrix).

***************************************************************************************************
5-) "Mesh_Refiner"
***************************************************************************************************

Uniform refinement of an ERMES mesh directly on the .dat files written by GiD, without
re-meshing in GiD. It is useful to obtain very fine meshes (beyond what GiD can generate
comfortably) or to run convergence studies on the same geometry.

- "mesh_refiner.py" : Reads the nodes ("*-1.dat") and tetrahedra ("*-4.dat") of the problem,
                      builds a conforming refined mesh and rewrites every mesh-dependent
                      .dat file in the same format (nodal conditions "*-2.dat" and "*-3.dat",
                      volume and current elements "*-4.dat" and "*-5.dat", boundary
                      triangles "*-6.dat" to "*-8.dat", "*-10.dat" and "*-19.dat", prisms
                      "*-9.dat"). Original nodes keep their numbers, child elements inherit
                      the material, condition Id, etc. of their parent, and tetrahedra and
                      triangles keep their orientation. The other .dat files (materials,
                      plasma, frequency, solver...) are not modified. Two refinement modes:
                        'levels' : each edge divided in L parts (tetrahedron -> L^3,
                                   triangle -> L^2) without degrading the element quality.
                        'ps'     : Powell-Sabin / Worsey-Farin split (tetrahedron -> 24,
                                   triangle -> 6).
                      Runs in parallel on Linux, WSL and Windows. Settings are defined in the
                      USER SETTINGS block or on the command line, e.g.:
                        >> python mesh_refiner.py ./MyProblem.gid --level 2 --np 16
                        >> python mesh_refiner.py --help
                      IMPORTANT: by default the files are refined IN PLACE and the refinement
                      is NOT idempotent (running it twice refines the refined mesh). Let GiD
                      rewrite the .dat files first, or use "--out <dir>" to keep the input.

- "Refine2ERMES.sh" : Bash script (Linux or WSL) that runs the full pipeline inside the GiD
                      problem folder: removes previous results, converts the .dat files to
                      Unix line endings (dos2unix), refines the mesh with "mesh_refiner.py"
                      and runs ERMES on the refined mesh (log in the ERMES "*.info" file).
                      A stamp file (".refined_stamp") prevents refining the same mesh twice;
                      the refinement is re-enabled every time GiD rewrites the .dat files.
                      Problem name/folder, refiner and ERMES paths, mode, level, number of
                      processes and FORCE_REFINE (0 = auto, 1 = always, 2 = ERMES only) are
                      set at the top, or as arguments / environment variables, e.g.:
                        >> ./Refine2ERMES.sh MyProblem 3
                        >> LEVEL=4 NPROC=8 ./Refine2ERMES.sh
                        >> ./Refine2ERMES.sh -h

***************************************************************************************************
6-) "NumPy_Solvers"
***************************************************************************************************

Examples of Python scripts used as external solvers of ERMES. They read the matrix and
vector written by ERMES ("Matrix_A_cmplx.bin", "Matrix_A_int.bin", "Vector_B.bin"), solve
the system, and write the solution to "Vector_Xo.bin".

To use them, in the ERMES GiD interface select "Solving Parameters -> Solver -> Solver type
= External solver" and set: Solver executable = python, Solver parameters = path to the
script. The matrix storage format (Symmetric, Full-matrix, Hermitic-Symmetric or
Hermitic-Full) is detected automatically from the files, as in the PETSc solver: no editing
is needed. If the solver fails (e.g. BiCG does not converge) no "Vector_Xo.bin" is written
and ERMES stops with an error instead of using a wrong solution. Both scripts work with old
and current SciPy versions.

- "BiCG.py"    : Iterative BiConjugate Gradient solver with diagonal preconditioner.
- "SuperLU.py" : Direct sparse LU solver (SuperLU).

***************************************************************************************************
7-) "PETSc_Python"
***************************************************************************************************

Early version of the interface between ERMES and PETSc, based on a Python script. It is kept
for testing purposes ( e.g. quick checks of PETSc solvers and options, or comparison with
other solvers ). For performance, and for production runs on HPC clusters, the
"PETSc_Solvers" version is preferred.

The Python script reads the ERMES linear system ("Matrix_A_cmplx.bin", "Matrix_A_int.bin",
"Vector_B.bin"), converts it to PETSc binary format ("Matrix_A_PETSc.dat" and
"Vector_B_PETSc.dat"), calls the PETSc solver with petscmpiexec, and converts the PETSc
solution ("Vector_X_PETSc.dat") back to the ERMES format ("Vector_Xo.bin"). Compared with
"PETSc_Solvers", note that:

- The whole matrix is read and converted by a single Python process before calling PETSc, and
  intermediate files are written to disk. This limits the problem size and adds memory and
  time overhead.
- Only the Symmetric and Full-matrix storage formats are supported ( no Hermitic formats ),
  and the format is NOT detected automatically: "MatrixType" must be set in the script to
  match the ERMES storage. Remember that, with damping, ERMES replaces the Hermitian storages by
  Full_matrix automatically ( see section 12 ).
- It was written for ERMES 20.0: paths quoted in the comments of the files (e.g. "../PETSc/",
  "/External_Solvers/PETSc") refer to the old folder structure and must be adapted.

- "ERMES2PETSc.py" : Python script called by ERMES as external solver ( Solver executable =
                     python3, Solver parameters = path to the script ). Solver type and
                     options, matrix format, number of MPI processes, problem folder and path
                     to the solver executable are set in the "Interface with PETSc" section.
                     Requires NumPy and SciPy.

- "win2wsl.sh"     : Template to call "ERMES2PETSc.py" from Windows through WSL2 ( sets
                     PETSC_DIR, PETSC_ARCH and LD_LIBRARY_PATH ). Use "dos2unix win2wsl.sh" if
                     an end-of-line error appears.

- "Solver"         : C++ source code of the PETSc solver ("PETScSolver.cpp"), "makefile"
                     (">> make PETScSolver"), and the PETSc Python modules used by the script
                     to write the PETSc binary files ("PetscBinaryIO.py", "petsc_conf.py").
                     The precompiled "PETScSolver" executable included must be rebuilt for
                     your own PETSc installation ( complex scalars ).

- "README.txt"     : PETSc installation instructions and examples of solver calls.

***************************************************************************************************
8-) "PETSc_Solvers"
***************************************************************************************************

Direct interface between ERMES and PETSc. The C++ solver "ERMESPETScSolver" reads the ERMES
linear system in parallel (MPI-IO) directly from the ERMES binary files, solves it with any
PETSc solver (e.g. MUMPS direct solver or iterative KSP solvers), and writes the solution
"Vector_Xo.bin" read back by ERMES. No Python and no intermediate files are needed, and no
MPI process ever holds the whole matrix.

It reads "Matrix_A_cmplx.bin", "Matrix_A_int.bin" and "Vector_B.bin" and, for the Hermitic
formats, "Matrix_A_aux_cmplx.bin" and "Matrix_A_aux_int.bin". The matrix storage format
(Symmetric, Full-matrix, Hermitic-Symmetric or Hermitic-Full) is detected automatically and
printed in the ERMES "*.info" file. Requires PETSc >= 3.19 configured with complex scalars
and double precision.

- "ERMESPETScSolver.cpp" : Source code of the parallel PETSc solver. Accepts any PETSc option
                           plus "-ermes_folder <dir>" (folder with the ERMES files) and
                           "-ermes_monitor_every <n>" (print residual every n iterations).

- "makefile"             : Builds "ERMESPETScSolver" with the PETSc installation given by
                           PETSC_DIR and PETSC_ARCH (">> make ERMESPETScSolver").

- "ERMES2PETSc.sh"       : Script called by ERMES as external solver, on Linux
                           ("bash /path_to/PETSc_Solvers/ERMES2PETSc.sh") and on Windows
                           through WSL ("wsl bash /path_to/PETSc_Solvers/ERMES2PETSc.sh").
                           The PETSc configuration, solver type and options, number of MPI
                           processes, optional hybrid MPI+OpenMP and MPI launcher are set at
                           the top. Inside a SLURM job it uses all the nodes and tasks of the
                           job. All output and errors go to the ERMES "*.info" file. It can
                           also be run directly from a problem folder.

- "README.txt"           : Installation and usage instructions: PETSc download, Linux
                           modules, environment variables, example of PETSc configuration, 
                           compilation, use from ERMES, and solver options.

***************************************************************************************************
9-) "Plasma_Scripts"
***************************************************************************************************

Scripts for tokamak plasma problems (e.g. electron cyclotron waves in MAST-U): EQDSK reading,
Scotty input, ERMES cold plasma files ( with optional turbulence and blobs ) and automation of
the full GiD + ERMES workflow. Run them from the "Plasma_Scripts" folder ( default paths are
relative to it ); they import each other, so keep them together. Full description, workflows
and conventions in "Plasma_Scripts/README.md".

Frame used by all the scripts: z is the tokamak axis, R = sqrt(x^2 + y^2), (R, phi, Z)
right-handed, lengths in metres. The previous Plasma_Scripts used y as the tokamak axis: GiD
problems and plasma files made with them must be regenerated.

- "eqdsk_to_scotty.py"     : Converts an EQDSK file plus density and temperature profiles into
                             the Scotty input files ( as in "Scotty_GiD" ). It also contains
                             the EQDSK reader and the magnetic field calculation used by the
                             other scripts, so the Scotty beam and the ERMES plasma agree.

- "eqdsk_reader.py"        : Plots an EQDSK file and the profiles ( poloidal flux, normalized
                             flux, poloidal and toroidal magnetic field, current density,
                             electron density and temperature ) to check the inputs.

- "eqdsk_to_ermes_simp.py" : Minimal generator of the electron density (ne) and magnetic field
                             (Be) files of the ERMES cold plasma material from an EQDSK file and
                             a density profile. Requires the list of nodes of the problem
                             ("*-1.dat"). Replaces "eqdsk2ermes_plasma_files.py".

- "eqdsk_to_ermes_turb.py" : Same files, plus optional synthetic turbulence ( field-aligned
                             random density field ) and blobs/holes, and a check plot. Also
                             provides the functions used by "problem_generator.py" ( plasma
                             files, ne and B at a point, B at the critical density ).

- "torbeam_to_ermes.py"    : Engine of "eqdsk_to_ermes_turb.py". Standalone, it generates the
                             plasma files from a Scotty/TORBEAM topfile or an EQDSK file; all
                             parameters can be given on the command line ( "--help" ).

- "problem_generator.py"   : Automates a full simulation: Gaussian beam launch geometry and
                             O-mode polarization ( B at the critical density ), GiD batch file,
                             GiD offscreen to build and mesh the model, plasma files, and ERMES.
                             Windows or Linux executables are selected automatically; the
                             solver is "Direct_solver" ( new in 20.0.5 ) or "External_solver"
                             ( e.g. "PETSc_Solvers/ERMES2PETSc.sh" ). Raises an exception if
                             ERMES fails.

- "batch_script.sh"        : SLURM example that runs "problem_generator.py" ( GiD and plasma
                             files ) and then ERMES on the generated problem ( paths and
                             resources to adapt ).

- "utilities.py"           : Functions used by "problem_generator.py" (vector rotations, beam
                             and plasma geometry, GiD batch file writer, O-mode polarization,
                             GiD node reader, cylindrical to Cartesian transformation).

- "permittivity_tensor.py" : Computes the cold plasma (Stix) permittivity tensor from the ne and
                             Be files for a multi-species plasma and writes it to the binary
                             file "Tensor_K_CPGE.bin" in ERMES format ("Tensor" geometry
                             extrapolation). Same conventions and constants as ERMES: one
                             tensor per data line of the ne file, B matched by node Id,
                             text/comment lines skipped.

- "README.md"              : Description of the scripts, workflows, conventions and notes.

- "EQDSK_Samples"          : Sample input data for the scripts above:
                             - "mast-u-sample.eqdsk" : MAST-U EQDSK equilibrium file
                             - "rho_flux.txt"        : Normalized flux coordinate rho_pol
                             - "rho_ne.txt"          : Electron density profile [m^-3]
                             - "rho_te.txt"          : Electron temperature profile [eV]

***************************************************************************************************
10-) "Results_Filters"
***************************************************************************************************

ERMES writes every computed result to the GiD post-process file ( "<problem>.flavia.res" ),
which can be very large. This script writes a smaller copy with the same mesh ( coordinates,
elements, Gauss points ) and only the results you choose, to open in GiD quickly or to share.

- "filter_ermes_res.py" : Reads the list of results from the file itself ( rE, iE, mod(E),
                          Poynting_vector, Electron_density, Bext, Plasma_permittivity_*,
                          ... ) and keeps the selected ones. Works with binary GiDPostEx
                          files ( gzip, ERMES default ) and ASCII .res files; the output has
                          the same format as the input. Selection ( case-insensitive, comma
                          or space separated ): result numbers "3", ranges "2-5", exact names
                          "rE", groups "E" ( = rE, iE, mod(E); any field written as rX / iX /
                          mod(X) ), "S" or "P" ( = Poynting_vector ) and "all".
                          Settings in the USER SETTINGS block ( for Spyder, F5 ) or on the
                          command line, which overrides them, e.g.:
                            >> python filter_ermes_res.py case.flavia.res --list
                            >> python filter_ermes_res.py case.flavia.res E S Bext
                            >> python filter_ermes_res.py case.flavia.res "1-3, 6" -o out.res
                          Output file: "-o", else OUTPUT_FILE, else "<input>_<selection>.res"
                          in the input folder. Note that OUTPUT_FILE is also used by command
                          line runs without "-o": set OUTPUT_FILE = "" to get the automatic
                          name. The input file is never overwritten. On error the script
                          prints the message and the available results and stops ( it never
                          asks questions ). From the Spyder console, after one run:
                            >> filter_res(r"D:\path\case.flavia.res", "E, Bext")

***************************************************************************************************
11-) "Scotty_GiD"
***************************************************************************************************

Scripts that compute the path and shape of an electron cyclotron Gaussian beam and export it
as geometry that GiD can import (STL triangle mesh and/or exact NURBS surfaces in IGES). The
beam geometry can then be used to build and refine the ERMES model around the beam (e.g. the
MAST-U Doppler back-scattering launch). The beam surface is the 1/e amplitude contour of the
beam ( times an optional factor ). All scripts can also write an "FEM envelope": a closed
polyhedron with planar faces around the beam at a given distance, to be used as the ERMES
computational volume, with labelled faces for the port ( ENVFIRST ), the end ( ENVLAST ),
the mirrors ( ENVPEC<k> ) and the sides ( ENVSIDE ). Text files are written with the central
ray, E and k along the ray and at both ends, the beam start in ERMES form ( point, k0, Epol,
Epar, Eper ), a straight-line or Scotty O-X conversion estimate and an order-of-magnitude
estimate of the ERMES mesh and solver cost. Parameters are set in the "USER PARAMETERS" block
of each script, or on the command line with "--name value" ("--show-params" and "--help" list
them). Keep the scripts in the same folder ( they import each other ) and run them from it.
Full description, workflows, conventions and output files in "Scotty_GiD/README.md".

- "eqdsk_to_scotty.py"     : Converts an EQDSK equilibrium file plus electron density and
                             temperature profiles into the input files of the Scotty beam
                             tracer ("topfile<sfx>", "ne<sfx>.dat", "Te<sfx>.dat" and,
                             optionally, "topfile<sfx>.json"). It also indicates the Scotty
                             options to use with these files (find_B_method,
                             density_fit_method, etc.). The same files are read by the O-X
                             diagnostics of the other scripts. The copy in "Plasma_Scripts"
                             reads the samples of "EQDSK_Samples".

- "scotty_to_gid.py"       : Runs a Scotty beam tracing simulation (or reads an existing
                             "scotty_output*.h5" file) and writes the beam surface as STL and/or
                             NURBS (IGES), the FEM envelope and the text files above. The beam
                             can be extended through vacuum before its first and after its last
                             point. The launch can be given as Scotty data or as ERMES data
                             ( point, k0, Epol, ... ). Also holds the functions shared by the
                             other scripts. Requires h5py ( and Scotty to run it ).

- "beam_to_gid.py"         : Same outputs as "scotty_to_gid.py" for a Gaussian beam propagating
                             in vacuum, without Scotty and without plasma. The user defines the
                             frequency, start point, direction (vector or poloidal/toroidal
                             angles), beam profile (waist or launch width and curvature, also
                             elliptical/astigmatic beams), polarization and propagation length,
                             or the ERMES beam data. It also writes the equivalent Scotty launch.

- "reflect_beam_to_gid.py" : Beam of "beam_to_gid.py" or "scotty_to_gid.py" reflected on
                             perfectly conducting (PEC) surfaces: planes, spheres, cylinders,
                             surfaces of revolution ( e.g. the EQDSK limiter ) and geometry
                             exported from GiD. Writes one beam surface per leg, the FEM
                             envelope around all legs ( mirror faces ENVPEC<k> ), a report of
                             the reflections ( hit points, incidence angles, curvatures ) and
                             the reflecting surfaces as STL to check them in GiD.

- "gid_surface_readers.py" : Reads GiD exports as triangles ( IGES, STEP, DXF, VTK, GiD mesh,
                             Exodus, STL, OBJ, ... ) for "reflect_beam_to_gid.py". IGES, DXF,
                             VTK, GiD mesh and STL need no extra package; STEP needs gmsh.

- "README.md"              : Description of the scripts, requirements, workflow, examples,
                             output files, conventions, current configuration and notes.

- "Plasma"                 : Sample input data for "eqdsk_to_scotty.py":
                             - "plasma_equilibrium.eqdsk" : EQDSK equilibrium file ( with
                                                            limiter )
                             - "rho_plasma.dat"           : Radial coordinate rho_pol
                             - "ne_plasma.dat"            : Electron density profile [m^-3]
                             - "te_plasma.dat"            : Electron temperature profile [eV]

***************************************************************************************************

See "ERMES_20.0.5_Manual.pdf" for a description of the ERMES input/output files, imported
conditions, external solvers and batch execution used by these scripts.

ERMES 20.0.5 is licensed under the open-source 2-clause BSD license. Publications resulting
from the use of this software must cite the following article, which describes the program:

R. Otin, "ERMES 20.0: Open-source finite element tool for computational electromagnetics in
the frequency domain", Computer Physics Communications, Vol. 310, 109521, 2025.

For more information about ERMES 20.0.5, please contact: ruben.otin.bcn@gmail.com

***************************************************************************************************
12-) Compatibility with ERMES 20.0.5 ( summary of all changes: RELEASE_NOTES_20.0.5.txt )
***************************************************************************************************

The scripts were checked against the corrected ERMES source code and updated where needed:

- ERMES now returns a non-zero exit status when it stops because of an error ( it always
  returned 0 before ). "ERMES_batch.sh" and "GiD_batch.sh" stop at the first failure,
  "SLURM_batch.sh" reports the failed case and continues, "Refine2ERMES.sh" ( set -e ) stops
  with its error message, and "problem_generator.py" raises an exception ( check=True ).
  ERMES error messages are also written to stderr: the scripts redirect it to the ".info" file.
  "Plasma_Scripts/batch_script.sh" stops with a message if ERMES fails.

- "ERMES_batch.sh" and "GiD_batch.sh" called "ERMES2PETSc.py", which does not exist; they now
  call "ERMES2PETSc.sh". "GiD_batch.sh", "SLURM_batch.sh" and "Refine2ERMES.sh" had Windows
  line endings, which make bash fail on Linux; all shell scripts now have Unix line endings.
  "SLURM_batch.sh" starts with "#!/bin/bash" as SLURM requires.

- External solvers: ERMES deletes any old "Vector_Xo.bin" before calling the external solver
  and stops with an error if the solver does not write a complete solution ( before, it wrote
  a zero or stale solution as if it were valid ). It also deletes old auxiliary matrix files
  when the current storage has none, so the automatic format detection of "ERMESPETScSolver"
  and of the NumPy scripts cannot pick up files from a previous run. The PETSc solver needs no
  change.

- Cold plasma: with damping, the Hermitian storages ( Herm-Full, Herm-Symm ) are replaced by
  Full_matrix automatically ( they were wrong with damping ). External solvers detect it.
  Plasma file paths may now contain spaces. Profile files may contain "#" comment lines and
  lines starting with blanks or signs; 1D profiles must have ascending positions, and the B
  file may use a different grid from the density file ( the samples and the files written by
  the "Plasma_Scripts" generators already comply ).

- Gaussian beam: the amplitude normalization is now the correct one for a 3D beam, so models
  built by "problem_generator.py" give different fields than with 20.0.4. For a 2D beam in
  a thin slab the waist radius of the invariant direction must be very large ( e.g. 1e3 );
  "problem_generator.py" uses equal radii ( see "Plasma_Scripts/README.md", notes ).

- New "Direct_solver": "problem_generator.py" selects it with PSettings.Solver
  ( "utilities.py" writes it in the GiD problem data ); "External_solver" keeps the PETSc
  route through "PETSc_Solvers/ERMES2PETSc.sh".

- "Plasma_Scripts" reorganized: "eqdsk2ermes_plasma_files.py" replaced by
  "eqdsk_to_ermes_simp.py" / "eqdsk_to_ermes_turb.py", the plasma-file functions removed from
  "utilities.py", z as tokamak axis, and the scripts checked against 20.0.5 ( details in
  "Plasma_Scripts/README.md" ).

- "permittivity_tensor.py": the off-diagonal sign of the Stix tensor was opposite to the one
  used inside ERMES ( K_xy = +iD instead of -iD ), so "Tensor" mode results corresponded to a
  reversed B field. Corrected, together with the other points listed in the script header.
  Tensor files generated with the old script must be regenerated.

- "PETSc_Python" is the early ERMES 20.0 interface and was not updated for 20.0.5 ( see
  section 7 ). Use "PETSc_Solvers" for production runs.

- New folder "Results_Filters" ( "filter_ermes_res.py" ): filters the results files written by
  ERMES 20.0.5 ( binary GiDPostEx or ASCII ).

- "Scotty_GiD" extended ( not needed for 20.0.5 ): reflections on PEC surfaces
  ( "reflect_beam_to_gid.py", "gid_surface_readers.py" ), FEM envelope, ERMES launch data,
  O-X and ERMES cost estimates, new sample files in "Plasma" ( see "Scotty_GiD/README.md" ).

- Not affected: Current_J_Scripts, IRBC_Scripts, IVEM_Scripts, Mesh_Refiner ( mesh_refiner.py ),
  Scotty_GiD, ERMESPETScSolver.cpp.
