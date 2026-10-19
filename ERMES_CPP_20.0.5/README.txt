************************************************************************************
* ERMES 20.0 - C++ source code
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
************************************************************************************

ERMES is a finite element solver of the frequency-domain Maxwell equations
( full-wave, cold plasma and electrostatic formulations ), with GiD as pre- and
post-processor.

------------------------------------------------------------------------------------
1. Folder contents
------------------------------------------------------------------------------------

1) "ERMES" contains the finite element formulations, boundary conditions and source
   elements, one file per element order:

      1st, 2nd, 3sb  nodal ( Lagrange ) elements of 1st, 2nd and hierarchical order
      1eg, 2eg       edge ( Nedelec ) elements of 1st and 2nd order

   and the cold plasma model ( "ColdPlasma.cpp" ). The subfolder "Utils" contains
   Gauss integration tables, Lagrange basis functions, physical constants
   ( "PhysicalConstants.h", CODATA 2018 ) and functions to write and read matrices
   and vectors to/from files. The "Legacy_*" subfolders keep older elements that are
   not compiled.

2) "external_libraries" contains the libraries used by ERMES: gidpost ( GiD results
   files ), dirent ( directory functions, Windows only ), flex++bison++ ( headers of
   the input parser, plus the programs and templates to regenerate it, see 6 ) and
   metis + GKlib ( ordering of the direct solver, compiled automatically with ERMES,
   see external_libraries/metis/ERMES_README.txt ).
   Boost is no longer needed ( ERMES 20.0.5 uses std::shared_ptr ) and the folders
   "boost" and "cbp2make" can be deleted. If this source package only contains the
   folders metis and GKlib, copy them into your existing external_libraries folder
   ( next to gidpost, dirent and flex++bison++ ).

3) "includes" contains the C++ template definitions and header files of the kernel
   ( nodes, elements, properties, sparse matrix, model, modeler, kernel ).

4) "license" contains the details of the ERMES open-source license.

5) "linear_solvers" contains the internal linear solvers: iterative BiCG ( COCG )
   and QMR for complex symmetric matrices, their versions for the cold plasma
   matrix storages, and the internal direct solver ( sparse LU ). See
   linear_solvers/README.txt for how to choose one.

6) "sources" contains the main C++ objects of ERMES:

   - "modeler.cpp" implements the solving strategies: it builds the system matrix,
     treats discontinuities and singularities, computes the geometrical normals and
     calls the linear solver.
   - "kernel.cpp" connects the input data with the objects in "modeler.cpp".
   - "input_c_parser.cpp" and "input_c_scanner.cpp" read the input files generated
     by the GiD pre-processor. They are generated with flex++ and bison++ from the
     templates "input_c_parser.y" and "input_c_scanner.l" ( see
     "external_libraries/flex++bison++/README.txt" ). If they are regenerated, the
     virtual destructors added to the Statement and ValueStatement classes must be
     added to the template as well.
   - "gid_output.cpp" writes the results in GiD post-processor format using the
     gidpost library.
   - "update.cpp" is legacy code that is not compiled.

7) "tests" contains unit tests of the linear solvers that need no input files
   ( Linux: "make test" in the unix folder ).

8) "unix" contains the Makefile and instructions to compile ERMES for Linux.

9) "windows" contains the Visual Studio 2022 solution to compile ERMES for Windows.

"RELEASE_NOTES_20.0.5.txt" summarizes the changes of ERMES 20.0.5 with respect to
ERMES 20.0.4; "CHANGES.md" lists them in detail, file by file.

------------------------------------------------------------------------------------
2. Compilation
------------------------------------------------------------------------------------

Linux  : cd unix ; make -j$(nproc)     ( details in unix/README.txt )
Windows: open windows/ERMES.sln with Visual Studio 2022 and build "Release x64"
         ( details in windows/README.txt )

------------------------------------------------------------------------------------
3. Running ERMES
------------------------------------------------------------------------------------

ERMES is normally launched by GiD. From a terminal ( e.g. on a cluster ):

   ERMES <problem_name>

( "ERMES --version" prints the version; it is also the first line of the .info file. )

where <problem_name>.dat and the related files written by GiD are in the working
directory. Results are written in GiD format next to the input files.

Linear solver: selected in the GiD problem data. The environment variable

   ERMES_SOLVER=Direct

forces the internal direct solver without changing the GiD problemtype.

Exit status: 0 on success, non-zero when ERMES stops because of an error ( wrong
input, singular matrix, solver breakdown ). Batch dependencies such as SLURM
"--dependency=afterok" can rely on it. With the "External_solver" option ERMES
writes the system matrix and returns 0 after doing so.

------------------------------------------------------------------------------------
4. Optional libraries - METIS ( ordering library of the direct solver )
------------------------------------------------------------------------------------

METIS ( Apache-2.0 ) computes the nested dissection ordering of the internal
"Direct_solver". For large 3D meshes it reduces the work and memory of the factorization
a lot ( example: n = 165 856 unknowns, 11x fewer operations, 4x less memory, solution in
62 s instead of 582 s on 16 threads ).

METIS IS INCLUDED: its sources are bundled in external_libraries/metis and
external_libraries/GKlib and compiled together with ERMES. Nothing has to be installed:

   Linux  : cd unix ; make -j STATIC=1         ( METIS included by default )
   Windows: Visual Studio, "Release | x64"     ( METIS included in every configuration )

The executables contain METIS ( static ), so installing ERMES on other computers is still
"copy the problemtype folder".

Check: run a case with Direct_solver; the ".info" file shows

   Ordering: METIS nested dissection

( METIS is used for n >= 2000 unknowns ). ERMES_ORDERING=AMD switches it off at run time.

Options ( Linux Makefile ):

   make METIS=0                          build without METIS ( built-in minimum degree only )
   make METIS_DIR=/path METIS_LIB="..."  use an installed METIS instead of the bundled one

To update the bundled METIS, see external_libraries/README_METIS.txt.
