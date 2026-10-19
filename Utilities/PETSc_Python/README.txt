************************************************************************************
* ERMES 20.0 - ERMES-PETSc interface
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - July 2024
************************************************************************************

This folder contains documentation about the PETSc library and the necessary files 
to connect ERMES 20.0 with the PETSc solvers.

This folder also contains the scripts that interface ERMES with PETSc (read the
comments on the scripts for more info).

The folder "/Solver" contains the PETSc executable and its C++ source code (read the
comments on the file "PETScSolver.cpp" for more info).

PETSc can be installed in Windows and Linux. We recommend its installation in Linux 
and, if required in Windows, call PETSc from the Windows Subsystem for Linux 2 
(WSL2). In this folder there is a bash file "win2wsl.sh" that can be used
as a template to call PETSc from Windows using WSL2. To install PETSc in Linux, 
follow these steps:

1) Download PETSc from its git repository to the local folder "/petsc": 

   >> git clone -b release https://gitlab.com/petsc/petsc.git petsc

2) Configure PETSc using the command "./configure" inside the folder "/petsc" (check
   PETSc manual for configuration parameters definitions). Note that multiple 
   configurations can be installed on the same machine, for instance:

   >> ./configure PETSC_ARCH=arch-complex-M1 --with-debugging=0 --download-mumps 
      --download-scalapack --download-parmetis --download-metis --download-ptscotch 
      --with-64-bit-indices=1 --with-scalar-type=complex --download-mpich 
      --download-cmake --with-openmp --download-hwloc --with-cc=gcc --with-cxx=g++ 
      --with-fc=gfortran --download-fblaslapack --download-bison --download-make
   
   >> ./configure PETSC_ARCH=arch-complex-S1 --with-debugging=0 --download-superlu_dist 
      --download-parmetis --download-metis --download-ptscotch --with-64-bit-indices=1 
      --with-scalar-type=complex --download-mpich --download-cmake --download-make 
      --download-bison --with-cc=gcc --with-cxx=g++ --with-fc=gfortran 
      --download-fblaslapack 
   
3) Follow the instructions provided during the configuration process. 

4) Add the selected active configuration PETSC_ARCH to ".bashrc":

   >> export PETSC_DIR=$HOME/petsc
   >> export PETSC_ARCH=arch-complex-M1
   >> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib
   
5) Run ".bashrc" or exit and log back in to your account:

   >> source $HOME/.bashrc
 
6) Compile "PETScSolver.cpp" for the selected configuration by typing inside its
   folder location (e.g. "/ERMES_20.0/External_Solvers/PETSc/Solver"):
   
   >> make PETScSolver

7) To run PETSc, use the Python scripts examples provided in this folder, or 
   the batch files examples located in "/ERMES_20.0/Documents/Batch_Scripts", or 
   call PETSc directly from the command line using, for instance:

   >> mpiexec -n 6 ./PETScSolver -ksp_type bcgs -pc_type asm 
      -ksp_converged_reason -ksp_monitor_true_residual 
      -log_view -memory_view
    
   >> mpiexec -n 8 ./PETScSolver -ksp_type richardson -pc_type lu 
      -pc_factor_mat_solver_type mumps -ksp_max_it 1 
      -log_view -ksp_error_if_not_converge
    
   >> mpiexec -n 4 ./PETScSolver -ksp_type gmres -pc_type lu 
      -pc_factor_mat_solver_type mumps -ksp_max_it 1 
      -ksp_monitor_true_residual
    
   >> mpiexec -n 8 ./PETScSolver -ksp_type richardson -pc_type lu 
      -pc_factor_mat_solver_type superlu_dist -ksp_max_it 1 
      -log_view -malloc_view
   
The above installation process can change depending on the Linux distribution. For 
further information, check the PETSc manual in this folder or visit the webpage: 
https://petsc.org.