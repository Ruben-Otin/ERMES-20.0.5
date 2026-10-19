*******************************************************************************
* ERMES 20.0.5 - Examples
*******************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - October 2026
*******************************************************************************

This folder contains three ready-to-open GiD projects that show how a
complete ERMES problem is set up: problem data, materials, boundary
conditions, sources and requested results. They are intended as a starting
point for new users and as templates for new models.

- Eddy_Currents.gid    : Low-frequency eddy currents (A-V potentials, 50 Hz)
- Microwave_Filter.gid : Rectangular waveguide filter with TE10 ports (9 GHz)
- SAM-Head.gid         : Plane-wave illumination of a head phantom (1 GHz)

Each example uses a different FEM formulation, excitation and set of
boundary conditions, so together they cover the most common ways of driving
an ERMES simulation.

*******************************************************************************
1-) How to open and run an example
*******************************************************************************

1. Install GiD and the ERMES problem type (see the main ERMES documentation).
2. In GiD go to  File > Open  and select one of the *.gid folders.
   GiD treats the whole folder as a single project; do not open the
   individual files inside it.
3. Explore the set-up (see section 2).
4. Generate the mesh:  Mesh > Generate mesh.
   The examples are distributed WITHOUT a mesh to keep the folder small,
   so this step is required before calculating.
5. Run the analysis:  Calculate > Calculate.
6. When the calculation finishes, switch to post-processing
   (Files > Postprocess, or the pre/post toggle button) to visualise the
   results.

Tip: before modifying an example, save a copy with  File > Save as  so the
original project remains unchanged.

*******************************************************************************
2-) What to explore in each project
*******************************************************************************

- Problem data (Data > Problem data):
    Problem mode, symmetry, FEM formulation (element type, A-V potentials
    on/off), length units, frequency, solver settings and the list of
    results to be written. Every field has a help text explaining the
    available options.

- Materials (Data > Materials):
    IHL materials (conductivity, permittivity, permeability), cold plasma,
    port definitions (rectangular and coaxial), Robin coefficients (plane
    waves, Gaussian beams...), current-source definitions and integration
    regions. Use the "Draw" option in the materials window to see which
    volume each material is assigned to.

- Conditions (Data > Conditions):
    Dirichlet conditions (PEC, PMC, TEC, PBC, voltages), Robin conditions
    (far field, full-wave Robin, ports), field integrals and current
    sources. Use the "Draw" option in the conditions window to highlight
    the entities where each condition has been applied.

All values are given in SI units (length units: m).

*******************************************************************************
3-) "Eddy_Currents.gid"
*******************************************************************************

- Description:
    A conducting plate with a square hole and a coil placed above it,
    enclosed in an air box. A time-harmonic current in the coil induces
    eddy currents in the plate.

- Problem data:
    Problem mode ........ Full_wave, 3D
    Formulation ......... RME_1st (regularized, 1st-order nodal elements)
    Potentials .......... On  (A-V potentials formulation, recommended for
                               low-frequency / conductor problems)
    Frequency ........... 50 Hz
    Solver .............. Quasi-Minimal Residual, tolerance 1e-6

- Materials:
    Vacuum .............. surrounding air
    Coil ................ sigma = 0, eps_r = 1, mu_r = 1
    Plate ............... sigma = 1e6 S/m, eps_r = 1, mu_r = 1
    J_axisymm ........... azimuthal current density Ja = 1e6 A/m^2,
                          phase 0, axis along Z through (0,0,0)

- Conditions:
    PEC ................. on the six outer faces of the air box
    FD_current_density .. J_axisymm applied to the coil volume

- Results requested:
    E, B, J, Joule heating (frequency domain) and B(t), J(t) (time domain)

- Things to try:
    - Change the plate conductivity or the frequency and observe the skin
      depth and the distribution of the induced currents.
    - Switch Potentials to Off to compare with the E-field formulation.
    - Select Solver = Direct_solver: the internal direct solver ( new in
      20.0.5 ) is very robust for low-frequency A-V problems like this one.

*******************************************************************************
4-) "Microwave_Filter.gid"
*******************************************************************************

- Description
    A WR-90 rectangular waveguide (22.86 mm x 10.16 mm cross-section,
    160 mm long) with internal discontinuities forming a microwave filter.
    The structure is excited at one end with the TE10 mode and the
    transmitted/reflected fields are evaluated at both ports.

- Problem data
    Problem mode ........ Full_wave, 3D
    Formulation ......... RME_1st
    Potentials .......... Off (E-field formulation)
    Frequency ........... 9 GHz (X band). A frequency sweep can be enabled
                          in the Frequency section of Problem data.
    Solver .............. Quasi-Minimal Residual, tolerance 1e-6

- Materials
    Vacuum .............. inside the waveguide
    RWPort_1 ............ Port ID 1, incident field |E| = 1 V/m, phase 0
                          (excitation port, located at x = 0)
    RWPort_2 ............ Port ID 2, incident field |E| = 0
                          (matched receiving port, located at x = 0.16 m)
    The port geometry is defined by a corner point (RW_00) plus the
    "Height" and "Width" vectors of the waveguide cross-section.

- Conditions
    PEC ................. on all metallic walls of the waveguide and filter
    RW_port_TE10 ........ RWPort_1 and RWPort_2 on the two end faces
    Projection_RWTE10 ... on the same two faces, used to compute the
                          scattering parameters Sij

- Results requested
    E, H, Poynting vector (frequency domain) and E(t), H(t) (time domain)

- Things to try
    - Enable Sweep_frequency_mode to obtain the filter response (S11, S21)
      over a frequency band.
    - Excite the structure from port 2 instead of port 1.

*******************************************************************************
5-) "SAM-Head.gid"
*******************************************************************************

- Description:
    A SAM (Specific Anthropomorphic Mannequin) head phantom, modelled as a
    homogeneous lossy dielectric, illuminated by a plane wave inside a
    computational box truncated with Robin boundary conditions. A typical
    set-up for absorbed-power / SAR-type studies.

- Problem data:
    Problem mode ........ Full_wave, 3D
    Formulation ......... RME_1st
    Potentials .......... On  (A-V potentials formulation)
    Frequency ........... 1 GHz
    Solver .............. Quasi-Minimal Residual, tolerance 1e-4

- Materials:
    Vacuum .............. surrounding air
    Dielectric .......... head tissue: sigma = 0.1 S/m, eps_r = 2.0, mu_r = 1
    Robin_coefficients_Plane_waves
                          |E| = 1 V/m, phase 0,
                          polarisation along +Y,
                          propagation (wave vector) along -Z
    Volume_1 ............ region used for the volume field integral

- Conditions:
    Robin_condition_Full_wave .. plane-wave Robin coefficients on the six
                                 outer faces of the box (incident wave plus
                                 absorbing condition for the scattered field)
    Surface_Voltage_Full_wave .. V = 0 on the same six outer faces
                                 (fixes the scalar potential of the A-V
                                 formulation)
    Field_volume_integral ...... on Volume_1

- Results requested:
    E, H, Poynting vector, Joule heating (frequency domain) and E(t), H(t)
    (time domain)

- Things to try:
    - Change the polarisation or the direction of incidence of the plane
      wave (the polarisation must be perpendicular to the wave vector).
    - Modify the tissue properties or the frequency and compare the
      absorbed power given by the volume integral.

*******************************************************************************
6-) Notes
*******************************************************************************

- Mesh density: the element size must resolve the wavelength in every
  material (and the skin depth in conductors). If you increase the
  frequency, refine the mesh accordingly.

- Solver: the examples use the built-in iterative solver with 8 processors
  ("Processors" field in Problem data); adjust it to your machine. If the
  iterative solver converges slowly, select "Direct_solver" ( internal
  sparse direct solver, new in ERMES 20.0.5 ), which also uses the
  "Processors" threads. For models too large for the memory of one computer
  the "External_solver" option can be used to export the system and solve it
  with an external tool (see the Utilities folder of the ERMES distribution,
  e.g. the PETSc-based solver).

- The example projects were prepared with ERMES 20.0.4 and updated for
  20.0.5 ( new solver option and help texts ). If GiD asks to transform
  the project to the current problem type, accept.

- Help: hovering over / clicking the help icon next to any field in the
  Problem data, Materials or Conditions windows shows a description of the
  option and its allowed values.
