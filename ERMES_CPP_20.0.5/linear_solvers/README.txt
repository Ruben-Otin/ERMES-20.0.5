************************************************************************************
* ERMES 20.0 - Linear solvers
************************************************************************************

All solvers derive from ComplexSolver ( complex_solver.h ) and are created in
Modeler::Generate_LinearSolver_FullWave / _ColdPlasma ( sources/modeler.cpp ) from
the solver name written by GiD.

------------------------------------------------------------------------------------
Matrix storages
------------------------------------------------------------------------------------

Full-wave        : complex symmetric matrix, upper triangle stored          -> Solve()
Cold plasma
  Full_matrix    : general matrix, all entries stored                       -> Solve()
  Herm-Symm      : A ( plasma, Hermitian, upper ) + A_aux ( symmetric, upper ) -> Solve_HermSymm()
  Herm-Full      : A ( plasma, Hermitian, upper ) + A_aux ( full )             -> Solve_HermFull()

Hermitian storages are only valid for a lossless plasma tensor; with damping ERMES
switches to Full_matrix automatically.

------------------------------------------------------------------------------------
Solvers
------------------------------------------------------------------------------------

GiD name                 Class                      Notes
-----------------------  -------------------------  --------------------------------
Bi_Conjugate_Gradient    BiCGComplexSolver          COCG for complex symmetric A.
                         BiCG_CP_ComplexSolver      BiCG for the plasma storages.
Quasi_Minimal_Residual   QMRComplexSolver           QMR for complex symmetric A
                         QMR_CP_ComplexSolver       ( default ). Smoother convergence.
Direct_solver            DirectComplexSolver        Sparse LU, every storage.
External_solver          ( none )                   ERMES writes A and b to files and
                                                    calls an external program
                                                    ( e.g. PETSc ), then reads x.

Iterative solvers: OpenMP parallel ( threads set from the GiD data ), optional
diagonal preconditioner. Convergence is checked on the recursive residual; the
true residual is printed at the end, with a warning when it is more than 10 times
the tolerance. A NaN/Inf breakdown stops ERMES with a failure exit status.

------------------------------------------------------------------------------------
Internal direct solver ( direct_complex_solver.cpp, multifrontal_complex_solver.cpp )
------------------------------------------------------------------------------------

- Supernodal multifrontal factorization ( the method used by MUMPS ):
  * analysis: approximate minimum degree ordering, elimination tree postorder,
    supernodes with relaxed amalgamation ( larger dense blocks );
  * complex symmetric matrices ( full-wave, electrostatic ): A = L D L^T, half the
    memory and work of LU; other storages ( cold plasma ): A = L U;
  * dense frontal matrices with cache-blocked kernels, parallel with OpenMP ( threads
    from the GiD data );
  * symmetric equilibration ( Ruiz scaling ) of the matrix before factorizing;
  * pivoting inside each front: partial pivoting ( LU ) or symmetric 1x1 pivoting
    ( LDL^T ) among the fully summed variables; pivots still smaller than 1e-8 x max|A|
    are perturbed ( static pivoting, as in PARDISO );
  * solution with GMRES preconditioned by the factorization: 1-2 iterations when the
    factorization is exact or nearly exact, and it recovers the accuracy when many
    pivots were perturbed ( strongly indefinite or saddle-point matrices );
  * if LDL^T breaks down or GMRES does not converge with it, ERMES refactorizes the
    matrix as a general one ( LU with partial pivoting ) automatically.
- Blocked front factorization, vectorized dense kernel ( runtime selection of generic, AVX2/FMA
  or AVX-512 code ) and two-level OpenMP parallelism ( independent subtrees in parallel, then
  parallel dense kernels for the large top fronts ).
- Measured on ONE core ( 3D grids ): n = 27 000 ( 7-point ) 4.6 s instead of 45 s with the first
  version; n = 10 600 ( 27-point, FEM-like ) 1.4 s instead of 6.9 s. Multi-core machines are
  faster still ( not measured ).
- Before factorizing, the predicted memory ( factors + work ), operations and a rough time are
  printed. If the memory exceeds 80% of the physical memory, ERMES stops with an error.
  Progress is printed every 10% of the predicted work.
- Recommended up to roughly 10^5 - 3x10^5 unknowns, depending on memory and cores. Larger
  problems: iterative solvers or the external solver ( PETSc + MUMPS, distributed memory ).
- Singular matrices ( missing boundary conditions, unconnected regions ) stop ERMES with a
  failure exit status.

Selection: "Direct_solver" in the GiD solver options ( requires adding this
option to the problemtype ), or the environment variable ERMES_SOLVER=Direct, which
overrides the GiD choice.

------------------------------------------------------------------------------------
Sparse products ( complex_matrix_operations.cpp )
------------------------------------------------------------------------------------

Symmetric / Hermitian products scatter the implicit lower triangle into reusable
per-thread buffers that are reduced in parallel. Vector operations are parallel
above ERMES_OMP_MIN_SIZE ( 20 000 ) entries. Tested by tests/test_matrix_operations.cpp.

------------------------------------------------------------------------------------
METIS ordering ( bundled )
------------------------------------------------------------------------------------

METIS ( Apache-2.0, https://github.com/KarypisLab/METIS ) provides the nested dissection
ordering: far less work and memory than minimum degree for large 3D meshes. Its sources are
bundled in external_libraries/metis and external_libraries/GKlib and compiled with ERMES
( Linux and Windows ), so it is always available. It is used for n >= 2000;
ERMES_ORDERING=AMD forces the built-in minimum degree ordering. Linux: make METIS=0 builds
without it.

Advanced ( not required ): the dense kernels can also call an external multithreaded BLAS
( make BLAS=openblas on Linux, see unix/Makefile ). On Windows this was measured 2x slower
than ERMES's own parallel kernel, so it is not used there. The ".info" file shows the
ordering and the dense kernel used ( "Ordering: ...,  dense kernels: ..." ).
