************************************************************************************
* ERMES 20.0 - Linux compilation instructions
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
************************************************************************************

Quick start ( GCC >= 5, any version up to current ones ):

   cd unix
   make -j$(nproc)                 ->  ./ERMES  ( -O3, OpenMP )

Options:

   make -j NATIVE=1                Tune for the CPU of the build machine ( -march=native ).
                                   Use it only if ERMES runs on the same CPU type
                                   ( e.g. build on a compute node of the cluster ).
   make -j STATIC=1                Fully static executable: runs on any x86-64 Linux, also with an
                                   older system C library than the build machine. Use it for the
                                   executable distributed with the GiD problemtype.
   make -j DEBUG=1                 -O0 -g with bounds-checked STL containers.
   make -j SANITIZE=1              Debug build with AddressSanitizer / UBSan.
   make -j METIS=0                 Build without the bundled METIS ordering library
                                   ( included by default, see INSTALL_OPTIONAL_LIBRARIES.txt ).
   make test                       Unit tests of the sparse products and the direct solver
                                   ( no input files needed ).
   make clean

On clusters it may be necessary to load a compiler module first, for instance:

   >> module load gcc

Number of threads: ERMES sets it from the problem data; OMP_NUM_THREADS is not needed.

Running from a terminal:  ./ERMES <problem_name>   ( <problem_name>.dat written by GiD )

Internal direct solver: "Direct_solver" in the GiD solver options, or

   export ERMES_SOLVER=Direct

to force it without changing the problemtype ( see linear_solvers/README.txt ).

Example SLURM job:

   #!/bin/bash
   #SBATCH --cpus-per-task=32
   #SBATCH --time=04:00:00
   module load gcc
   srun ./ERMES my_problem

Exit status: ERMES now returns a non-zero exit status when it stops because of an
error ( previous versions always returned 0 ). Batch dependencies such as SLURM
"--dependency=afterok" therefore work as expected.

The previous CodeBlocks-generated Makefile is kept as "Makefile_codeblocks"
( generated from the Visual Studio solution with cbp2make ). It is no longer needed:
the new Makefile picks up every .cpp file in ERMES/, ERMES/Utils/, linear_solvers/
and sources/ automatically, so adding a source file needs no Makefile change.
