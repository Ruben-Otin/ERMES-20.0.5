************************************************************************************
* ERMES 20.0 - Windows compilation instructions
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
************************************************************************************

Requirements: Microsoft Visual Studio 2022 with the "Desktop development with C++"
workload. The "external_libraries" folder must be present next to this folder.

To compile ERMES for Windows:

1) Open "ERMES.sln" ( this folder ) with Visual Studio 2022.

2) Check that the toolbar shows "Release" and "x64" ( the default configuration of the
   solution ), then "Build > Build Solution". OpenMP is already enabled in the project.
   The project uses the Visual Studio 2022 toolset ( v143 ) and the latest installed
   Windows 10/11 SDK, so no "Retarget Projects" dialog should appear.
   Visual Studio remembers the last configuration used in the hidden ".vs" folder: if it
   opens in "Debug", select "Release" once, or delete ".vs" to return to the default.
   Do not use Debug builds for production: they are not optimized and much slower.

3) The executable "ERMES.exe" is written to "windows\x64\Release".

4) Rename "ERMES.exe" with the current version name, "ERMES_20.0.5.exe".

5) Copy "ERMES_20.0.5.exe" to the GiD problemtype folder "ERMES_20.0.5\ERMES.gid\execs\win",
   replacing the previous executable.

6) Check that "ERMES.win.bat" in "ERMES_20.0.5\ERMES.gid" calls "ERMES_20.0.5.exe" ( already set ):

   >> "%3\execs\win\ERMES_20.0.5.exe" %1 > "%2\%1.info" 2>&1

7) After compilation the following folders can be deleted to save disk space
   ( Visual Studio recreates them ):

   - "windows\x64"
   - "windows\Release"
   - "windows\.vs"

Notes:

- New source files must be added to the project ( "ERMES.vcxproj" ). This version adds
  the direct solver ( linear_solvers\direct_complex_solver.*, multifrontal_complex_solver.* ),
  ERMES\Utils\PhysicalConstants.h and the METIS / GKlib sources, already in the project.

- The internal direct solver is selected with "Direct_solver" in the GiD solver
  options ( when the problemtype offers it ) or by setting the environment
  variable ERMES_SOLVER=Direct before launching GiD or ERMES.

- METIS ( ordering of the direct solver ) is bundled in external_libraries\metis and
  external_libraries\GKlib and compiled with ERMES in every configuration: no installation
  ( vcpkg is no longer needed ).

- Boost is no longer needed: the "external_libraries\boost" folder can be deleted.

- "external_libraries\libc.lib" is no longer used and can be deleted: it was a 32-bit library
  ( from the lcc-win32 compiler ) that the 64-bit build ignored, and it is no longer referenced
  by the project.

- To compile ERMES for Linux, read "unix\README.txt".
