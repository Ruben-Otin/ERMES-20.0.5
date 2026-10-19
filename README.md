<div align="center">

<img src="Images/ERMESLogo.png" alt="ERMES logo" width="320"/>

# ERMES 20.0.5

**E**lectric **R**egularized **M**axwell's **E**quations with **S**ingularities

*Open-source finite element tool for computational electromagnetics in the frequency domain*

[![License: BSD-2-Clause](https://img.shields.io/badge/License-BSD_2--Clause-blue.svg)](#license)
![Version](https://img.shields.io/badge/version-20.0.5-brightgreen)
![Platforms](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-lightgrey)
![Language](https://img.shields.io/badge/C%2B%2B-source-00599C?logo=cplusplus&logoColor=white)
![GUI](https://img.shields.io/badge/GUI-GiD-orange)
![Solvers](https://img.shields.io/badge/solvers-Direct%20%7C%20Iterative%20%7C%20PETSc-8A2BE2)

[Overview](#overview) •
[What's new](#whats-new) •
[Features](#features) •
[Package contents](#package-contents) •
[Installation](#installation) •
[PETSc](#petsc) •
[Workflow](#workflow) •
[Gallery](#gallery) •
[Citation](#citation) •
[Previous versions](#previous-versions) •
[Contact](#contact)

</div>

---

<a id="overview"></a>

## 🔭 Overview

ERMES is an open-source C++ finite element (FEM) code that solves Maxwell's equations in the frequency domain. Version 20.0 was a major upgrade, adding new modules, FEM formulations and HPC capabilities aimed at the demanding electromagnetic problems found in the design and analysis of nuclear fusion reactors. Version **20.0.5** adds an internal parallel direct solver, corrects several bugs and makes ERMES easier to build and install ( see [What's new](#whats-new) ).

ERMES works across the **static, quasi-static and high-frequency** regimes. It has been applied to microwave engineering, bioelectromagnetics, electromagnetic compatibility and electromagnetic forming, and, through its electrostatic and cold-plasma modules, to fusion problems such as plasma disruption induced currents and forces, plasma control, electric-arc initiation, and wave–plasma–wall interaction.

The graphical user interface is fully integrated into the pre/post-processor **[GiD](https://www.gidsimulation.com)**, which handles geometry, data input, meshing, and visualisation. GiD is the recommended option for setting up simulations and post-processing results, but it is not required. Any pre-processor can be used, provided that it generates the ERMES input files in the format described in the manual.

---

<a id="whats-new"></a>

## 🆕 What's new in 20.0.5

- **Internal direct solver** — parallel ( OpenMP ) supernodal multifrontal factorization with **METIS** nested dissection ordering, pivoting and GMRES refinement. Robust where iterative solvers struggle ( resonances, low losses, high-contrast materials, low-frequency A–V formulations ). Example: 165 856 unknowns solved in about one minute on 16 threads.
- **METIS bundled** — compiled automatically with ERMES on Windows and Linux; nothing to install.
- **Ready-to-run executables** — a static Linux executable that runs on any x86-64 Linux ( workstations and clusters ) and a Windows executable; installation is still copy & paste.
- **Corrections that change results** — Gaussian-beam amplitude normalization ( see the note for 2D models below ), Gauss-point evaluation of Gaussian-beam and plane-wave sources in nodal elements, cold-plasma Hermitian storage with damping, 1D B-field profiles, permittivity-tensor sign in `permittivity_tensor.py`.
- **Robustness** — clear error messages and a non-zero exit status on errors ( useful for scripts and SLURM ), solver breakdown detection, input checks, safe external-solver workflow.
- **Faster** iterative solvers ( parallel matrix–vector and vector operations ) and **cleaner build** ( no Boost, simple `make` on Linux, Release build on Windows ).
- **Version check** — the first line of the `.info` file shows the ERMES version, and `ERMES --version` prints it.

> **Gaussian beams in 2D models** ( geometry invariant in one direction, e.g. thin slabs ): set the waist radius of the invariant direction to a very large value ( e.g. `1e3` ). The beam then becomes exactly the 2D Gaussian beam, with |E| at the waist equal to the *GB module* input.

The full list of changes is in `RELEASE_NOTES_20.0.5.txt` ( summary ) and `ERMES_CPP_20.0.5/CHANGES.md` ( detailed ).

---

<a id="features"></a>

## ✨ Features

### 🧮 FEM formulations

ERMES 20.0.5 solves the **frequency-domain, full-wave Maxwell's equations with the Finite Element Method (FEM)**. Three formulations are available: **the regularized Maxwell's equations with nodal elements, the double-curl Maxwell's equations with edge elements, and the local $L^2$ projection method with nodal and bubble elements**. Each can be solved in terms of the electric field **E** or in terms of the magnetic vector potential **A** together with the scalar electric potential V.

The double-curl and local $L^2$ projection formulations can be **stabilised with a Lagrange multiplier**, which can be switched on or off. This matters most in the low-frequency and quasi-static regimes, where these formulations tend to become ill-conditioned. The regularized nodal formulation includes dedicated treatment of **field discontinuities** at material interfaces and of **singularities** at sharp edges and corners. 

Having several formulations available means the most stable one can be chosen for each problem and produce the best-conditioned matrix, which can often be solved with low-memory iterative methods. 

| Formulation | Elements | A–V potentials | Lagrange-multiplier stabilisation |
|---|---|:---:|:---:|
| Regularized Maxwell's equations | Nodal | ✅ | — |
| Double-curl Maxwell's equations | Edge | ✅ | ✅ (switchable) |
| Local $L^2$ projection method | Nodal + bubble | ✅ | ✅ (switchable) |

### 🧱 Materials
- **IHL materials** — isotropic, homogeneous, linear media
- **Cold plasma** — configurable plasma geometry, composition and problem settings

### ⚡ Sources
- **Volumetric current densities J** — Cartesian, axisymmetric, combined Cartesian + axisymmetric, helicoidal, and plasma modes
- **Boundary sources** — rectangular waveguide ports (TE10), coaxial ports (TEM), plane waves, Gaussian beams, quasi-static and electrostatic current fluxes

### 🧭 Boundary conditions
- **Dirichlet** — PEC, PMC, TEC, periodic/cyclic (PBC), full-wave voltage, electrostatic voltage, singularity
- **Robin** — full-wave and cold-plasma far-field conditions (FW, RW, LW, PW waves), plane-wave / Gaussian-beam / quasi-static Robin coefficients, imported Robin data, electrostatic Robin, coaxial TEM and RW-port TE10

### 📊 Outputs
- Frequency-domain, time-domain, cold-plasma and electrostatic field visualisation
- Surface and volume integrals, scattering parameters, and fields-on-nodes files
- Results in GiD format, plus generic formats adaptable to other post-processors

### 🚀 Solvers & HPC
- Built-in **BiCG** and **QMR** iterative solvers with diagonal preconditioning and **OpenMP**, with a low memory footprint
- Built-in **direct solver** ( new in 20.0.5 ): parallel multifrontal LDLᵀ / LU with **METIS** ordering, for shared-memory computers when the iterative solvers struggle to converge
- **PETSc interface** that reads the ERMES system directly and in parallel, and solves it using direct or iterative solvers over **MPI** on workstations or SLURM clusters, with no Python or intermediate files.
- Interface to **Python NumPy**, with scripts to export the linear system in other formats
- Editable `.dat` input files and batch/Python scripts for parametric and cluster runs

---

<a id="package-contents"></a>

## 📦 Package contents

```
ERMES 20.0.5/
├── ERMES_20.0.5_Manual.pdf   📘 User manual
├── ERMES_CPP_20.0.5/         🧩 C++ source code
├── ERMES_20.0.5/             🖥️ GiD user interface (problemtype, with executables)
├── Examples/                 🧪 GiD usage examples
└── Utilities/                🛠️ Python tools, solvers and batch scripts
```

Download only what you need:

| Folder | Get it if you want to… |
|---|---|
| `ERMES_CPP_20.0.5` | work with or modify the C++ source code |
| `ERMES_20.0.5` | use ERMES through the GiD graphical interface ( ready-to-run executables included ) |
| `Examples` | explore ready-made usage examples |
| `Utilities` | customise boundary conditions, plasmas, sources and solvers |

<details>
<summary><b>🧩 Inside the C++ source tree</b></summary>

| Subfolder | Contents |
|---|---|
| `ERMES/` | FEM formulations, boundary conditions, source elements, Gaussian integration tables, file I/O |
| `external_libraries/` | gidpost, dirent, `flex++bison++`, and METIS + GKlib ( ordering of the direct solver, compiled automatically ) |
| `includes/` | Object template definitions and header files |
| `license/` | Open-source license details |
| `linear_solvers/` | Internal linear solvers: iterative ( BiCG, QMR ) and direct ( multifrontal ) |
| `sources/` | Core objects: `modeler.cpp` (solving strategies, system matrix, discontinuities/singularities, normals), `kernel.cpp`, the input parser/scanner and `gid_output.cpp` |
| `tests/` | Unit tests of the linear solvers ( `make test` ) |
| `unix/` | Linux `Makefile` ( `make`, or `make STATIC=1` for a portable executable ) |
| `windows/` | Visual Studio 2022 solution `ERMES.sln` ( Release x64 ) |

</details>

---

<a id="installation"></a>

## 🛠️ Installation

### 1. Install GiD
Download the latest GiD from [gidsimulation.com](https://www.gidsimulation.com), then register it via **Help → Register GiD** (local, USB, floating or named-user licence; a one-month free licence is available and renewable).

### 2. Install ERMES
Copy the entire `ERMES_20.0.5` folder into GiD's `problemtypes` directory (e.g. C:\MySoftware\GiD 16.0.2\problemtypes). Then, restart GiD and select **Data → Problem type → ERMES_20.0.5 → ERMES**. The ERMES logo and menu bar should appear. The same procedure applies to Windows and Linux.

The executables included in the problem type are ready to run, with nothing else to install: `ERMES.gid/execs/win/ERMES_20.0.5.exe` ( with `vcomp140.dll` in the same folder ) runs on any 64-bit Windows computer, and `ERMES.gid/execs/unix/ERMES_20.0.5` is statically linked and runs on any x86-64 Linux computer, including HPC clusters. The first line of the `.info` file of a run, `ERMES 20.0.5 analysis started...`, confirms the version used.

### Compiling ( only if needed )
Recompiling is only necessary if the source code is modified:
- **Windows** — open `ERMES_CPP_20.0.5/windows/ERMES.sln` with Visual Studio 2022 and build **Release | x64**.
- **Linux** — `cd ERMES_CPP_20.0.5/unix` and run `make -j8` ( or `make -j8 STATIC=1` for a portable, fully static executable ).

See the `README.txt` files in those folders for details.

### 3. Install PETSc
For high-performance computing and large problems, use the ERMES–PETSc interface in `Utilities/PETSc_Solvers`. See installation details in **[ERMES–PETSc interface](#petsc)** below.

---

<a id="petsc"></a>

## ⚡ ERMES–PETSc interface

`ERMESPETScSolver` reads the ERMES linear system **A X₀ = B** directly and in parallel, solves it with **[PETSc](https://petsc.org)**, and writes the solution `Vector_Xo.bin` back for ERMES. **No Python and no intermediate files** are needed.

📁 Source: [`Utilities/PETSc_Solvers`](Utilities/PETSc_Solvers)

| File | Purpose |
|---|---|
| `ERMESPETScSolver.cpp` | Solver source code |
| `makefile` | Builds the solver |
| `ERMES2PETSc.sh` | Script called by ERMES (all settings at the top) |

### 📥 Files read from the problem folder

| File | Content |
|---|---|
| `Matrix_A_cmplx.bin`, `Matrix_A_int.bin` | System matrix (coefficients, indices) |
| `Matrix_A_aux_cmplx.bin`, `Matrix_A_aux_int.bin` | Auxiliary matrix (Hermitic formats) |
| `Vector_B.bin` | Right-hand side |

### 🔍 Matrix formats (detected automatically)

| Format | `Matrix_A` | `Matrix_A_aux` |
|---|---|---|
| **Symmetric** | Diagonal + upper triangle of a complex symmetric matrix | — |
| **Full-matrix** | Entire matrix | — |
| **Hermitic-Symmetric** | Diagonal + upper triangle of the Hermitian volumetric contribution | Diagonal + upper triangle of the complex symmetric Robin boundary contribution |
| **Hermitic-Full** | As above | Entire Robin contribution, added as stored |

If `Matrix_A_aux` files are present the format is Hermitic, and the number of stored triangles selects the variant. The detected format is printed in the `*.info` file.

### 🛠️ Installation

```bash
# 1. Download PETSc
git clone -b release https://gitlab.com/petsc/petsc.git petsc

# 2. Environment (also needed before compiling the solver;
#    use the same PETSC_DIR and PETSC_ARCH in section 1 of ERMES2PETSc.sh)
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex-M1
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH

# 3. Configure PETSc (several configurations can coexist, one per PETSC_ARCH)
cd $PETSC_DIR
./configure PETSC_ARCH=arch-complex-M1 \
  --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
  --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices=1 \
  --download-mpich --download-hwloc --download-openblas \
  --download-scalapack --download-mumps \
  --download-metis --download-parmetis --download-ptscotch \
  --download-cmake --download-bison --download-make
#    ...then run the "make" commands printed at the end of configure

# 4. Compile the solver (inside Utilities/PETSc_Solvers)
make ERMESPETScSolver

# 5. Execution permissions
chmod +x *.sh ERMESPETScSolver
```

### ▶️ Use from ERMES

In ERMES, set **Solver type** to **External solver** and enter in **External solvers settings**:

| OS | Command |
|---|---|
| 🐧 Linux | `bash /path_to/PETSc_Solvers/ERMES2PETSc.sh` |
| 🪟 Windows | `wsl bash /path_to/PETSc_Solvers/ERMES2PETSc.sh` |

ERMES runs the script from the problem folder. Inside a **SLURM** job, the solver runs on all the nodes and tasks of the job. All solver output, including errors, goes to the ERMES `*.info` file.

The script can also be run without ERMES, from a folder that contains the matrix and vector files:

```bash
cd /path_to/problem_folder
bash /path_to/PETSc_Solvers/ERMES2PETSc.sh
```

For more information, see README.txt in `Utilities/PETSc_Solvers` or visit [petsc.org](https://petsc.org).

---

<a id="workflow"></a>

## 🔁 Workflow


```mermaid
%%{init: {"themeVariables": {"fontSize": "13px"}, "flowchart": {"rankSpacing": 10, "nodeSpacing": 1, "padding": 10, "diagramPadding": 1, "subGraphTitleMargin": {"top": 10, "bottom": 10}}}}%%
flowchart TB
    subgraph PRE["🧱 GiD pre-processor"]
        direction TB
        A["#nbsp;#nbsp;#nbsp;Create / import geometry#nbsp;#nbsp;#nbsp;"] --> B[Assign materials & BCs]
        B --> C[Set problem parameters]
        C --> D[Mesh]
        D --> E[Write .dat files]
    end
    subgraph SOLVE["⚙️ ERMES"]
        direction TB
        F[Read .dat input files] --> G[Assemble & solve]
        G --> H[Write results]
    end
    subgraph POST["📊 GiD post-processor"]
        direction TB
        I[Open results] --> J[Visualise & analyse]
    end
    PRE -----> SOLVE
    SOLVE -----> POST
```

GiD is recommended but not required. Any pre-processor can be used as long as it writes the `.dat` input files in the format defined by the `.bas` templates in `ERMES_20.0.5/ERMES.gid`. Those `.dat` files can also be edited directly for parametric batch runs.

---

<a id="gallery"></a>

## 🖼️ Gallery

<p align="center">
  <img src="Images/GiD-GeoRender.png"/><br/>
  <sub>ERMES GiD interface with CAD geometry of a coil and conductive plate (TEAM problem 7)</sub>
</p>

<p align="center">
  <img src="Images/CheckRWPort.png"/><br/>
  <sub>ERMES GiD interface showing the rectangular waveguide ports applied to a microwave filter</sub>
</p>

<p align="center">
  <img src="Images/JET-A2-E-3D.png"/><br/>
  <sub>Electric field generated by the JET A2 antenna</sub>
</p>

<p align="center">
  <img src="Images/JET-A2-E-FRT.png"/><br/>
  <sub>Cut planes of the electric field generated by the JET A2 antenna</sub>
</p>

<p align="center">
  <img src="Images/Disruption-ITER.png"/><br/>
  <sub>Eddy currents induced on ITER walls by plasma disruption</sub>
</p>

<p align="center">
  <img src="Images/STEP-EddyCurrents.png"/><br/>
  <sub>Eddy currents induced on STEP components by plasma and control coils</sub>
</p>

<p align="center">
  <img src="Images/ITER-Arc-J.png"/><br/>
  <sub>Electric currents following an electric arc strike on ITER components</sub>
</p>

<p align="center">
  <img src="Images/ITER-Quench-Voltage.png"/><br/>
  <sub>Voltage induced by a quench in a superconducting poloidal coil of ITER</sub>
</p>

<p align="center">
  <img src="Images/BWShield.png"/><br/>
  <sub>Fields leaking from a coaxial cable braided shield</sub>
</p>

<p align="center">
  <img src="Images/BioEM.png"/><br/>
  <sub>Electric field exposure simulation (EU NextGEM project)</sub>
</p>

---

<a id="documentation"></a>

## 📘 Documentation

The full user manual, `ERMES_20.0.5_Manual.pdf`, covers:

1. **Introduction** — user interface, C++ source code, license
2. **Installation** — GiD, ERMES, PETSc
3. **Pre-process** — geometry, materials, sources, boundary conditions, and problem settings
4. **Post-process** — output files and GiD visualisation
5. **Appendix A** — electromagnetic theory
6. **Appendix B** — finite element formulations

---

<a id="citation"></a>

## 📝 Citation

Publications resulting from the use of ERMES **must cite**:

R. Otin, "ERMES 20.0: Open-source finite element tool for computational electromagnetics in the frequency domain", *Computer Physics Communications*, Vol. 310, 109521, 2025.

```bibtex
@article{Otin2025ERMES,
  author  = {R. Otin},
  title   = {{ERMES} 20.0: Open-source finite element tool for computational
             electromagnetics in the frequency domain},
  journal = {Computer Physics Communications},
  volume  = {310},
  pages   = {109521},
  year    = {2025}
}
```

---

<a id="license"></a>

## ⚖️ License

ERMES 20.0.5 source code, binaries and interface are released under the **2-clause BSD license**. The executables also contain third-party open-source components ( METIS / GKlib under Apache-2.0, and the compiler runtime libraries ); their licenses are in `ERMES_20.0.5/ERMES.gid/execs/licenses`.

<details>
<summary>Full license text</summary>

```
Copyright (c) 2013 Ruben Otin, CIMNE (International Center for Numerical
Methods in Engineering). All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

</details>

---

<a id="previous-versions"></a>

## 🗂️ Previous versions

Earlier ERMES versions, including their executables, manuals and examples, remain available at:

- 🌐 [Ruben Otin software page ( blog )](https://ruben-otin.blogspot.com/2015/04/ruben-otin-software-ruben-otin-april-19.html)
- 📁 [ERMES versions archive ( Google Drive )](https://drive.google.com/drive/folders/1Gipqbgu3bHuv8TdSJd-eZ2W2nOFT75kS)

---

<a id="contact"></a>

## 📬 Contact

**Ruben Otin** — United Kingdom Atomic Energy Authority (UKAEA), Oxford, UK

- ✉️ [ruben.otin@ukaea.uk](mailto:ruben.otin@ukaea.uk)
- ✉️ [ruben.otin.bcn@gmail.com](mailto:ruben.otin.bcn@gmail.com)
- 🌐 [ruben-otin.blogspot.com](https://ruben-otin.blogspot.com)

<div align="center">
<sub>ERMES 20.0.5 · October 2026</sub>
</div>
