# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

DG-SWEM (Discontinuous Galerkin Shallow Water Equation Model) — a Fortran/C++ computational model for simulating shallow water flow (coastal flooding, storm surge, compound flooding). Source files use fixed-form Fortran `.F` with preprocessor macros and free-form `.F90`/`.f90`.

## Build System

DG-SWEM uses [Meson](https://mesonbuild.com/). Two pre-existing build directories are present: `build/` (serial) and `build-parallel/` (parallel with MPI).

**Serial (single-core):**
```bash
FC=gfortran CC=gcc CXX=g++ meson setup build [OPTIONS]
cd build && meson compile -v
# Produces: dgswem_serial, adcprep, adcpost
```

**Parallel (MPI):**
```bash
FC=gfortran CC=gcc CXX=g++ MPIFC=mpif90 meson setup build-parallel -Dparallel=true
cd build-parallel && meson compile -v
# Produces: dgswem, dgswem_serial, adcprep, adcpost
```

**Key build options** (`-D<opt>=<val>`):

| Option | Values | Default |
|---|---|---|
| `parallel` | `true/false` | `false` |
| `gpu` | `true/false` | `false` (NVHPC only) |
| `precision` | `single/double` | `double` |
| `netcdf` | `true/false` | `false` |

Wipe and reconfigure: `meson setup build [OPTIONS] --wipe`

View all options: `meson configure build`

**Conda environments** for reproducible builds are in the repo root:
- `environment_gcc.yml` — GCC + OpenMPI + netCDF
- `environment_intel.yml` — Intel `ifx`/`icx` + IntelMPI

CI uses `export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig` before `meson setup` when using conda-installed netCDF.

## Testing

Tests live in `tests/` and require `pytest`. Run from the `tests/` directory:

```bash
cd tests

# Run all tests against a build directory
pytest test_run.py --binpath ../build

# Run a single test
pytest test_run.py::test_quarter_annular --binpath ../build

# Run serial tests only (parallel skipped if dgswem not built)
pytest test_run.py -k "not parallel" --binpath ../build
```

Tests compare the last snapshot of `fort.63` against `fort.63.true`. The `--binpath` argument is required and must point to a Meson build directory containing the executables.

Performance profiling tests (`test_performance_*`) require Intel VTune; they are auto-skipped if `vtune` is not on `PATH`.

## Running

All executables must be run from the directory containing the ADCIRC-format input files (`fort.14`, `fort.15`, `fort.dg`, etc.).

**Serial:**
```bash
./dgswem_serial
```

**Parallel:**
```bash
./adcprep          # domain decomposition, prompts for num ranks + input file names
mpirun -np <N> ./dgswem
./adcpost          # merge subdomain outputs into global files
```

**Key input files:**
- `fort.14` — mesh (nodes + elements)
- `fort.15` — run parameters (ADCIRC format; see adcirc.org docs)
- `fort.dg` — DG-SWEM control parameters (keyword format; sample in `work/fort.dg`)

`fort.dg` controls solver-specific options not in the ADCIRC spec: polynomial order (`px`, `pl`, `ph`), flux type (`fluxtype`: 1=Roe, 2=LLF, 3=HLLC, 4=NCP), slope limiter (`slopeflag`), rainfall model (`rainfall` 0–4), p-enrichment (`padapt`), sediment, artificial diffusion.

Output files follow ADCIRC conventions (`fort.63` = water surface elevation time series, etc.).

## Code Architecture

### Source organization

```
src/          — all Fortran/C++ solver source
prep/         — adcprep (domain decomposition) and adcpost (output assembly) source
swan/         — SWAN wave model source (coupled via -DSWAN preprocessor flag)
subprojects/  — Meson wrap for METIS graph partitioner
tests/        — pytest integration tests + test case directories
work/         — sample input files and run scripts
```

### Key modules (in `src/`)

- `global.F` / `sizes.F` — global data structures and mesh arrays (node coordinates, element connectivity, state variables). Nearly everything else `USE GLOBAL`.
- `dgswem.F` — main program entry point; orchestrates cold/hot start → time loop → output.
- `dg.F` — DG module containing mesh-level DG data structures (modal coefficients, edge data).
- `read_input.F` / `read_fort_dg.F90` — mesh and parameter I/O.
- `DG_timestep.F` / `DG_hydro_timestep.F` — Runge-Kutta time integration.
- `rhs_dg_hydro.F` — right-hand side assembly (calls edge flux routines).
- `numerical_flux.F` — Riemann solvers (Roe, LLF, HLLC, NCP).
- `*_edge_hydro.F` — boundary condition flux handlers: `internal_edge`, `ocean_edge`, `land_edge`, `flow_edge`, `radiation_edge`, `ibarrier_edge`, `ebarrier_edge`.
- `slopelimiter.F` / `prep_slopelim.F` — slope limiting for shock capture.
- `wetdry.F` — wetting and drying treatment.
- `precipitation.F` / `owi_rain.F` — rainfall forcing.
- `messenger.F` / `messenger_elem.F` — MPI halo exchange (compiled in only with `-DCMPI`).
- `write_output.F` / `write_results.F` — output routines (ADCIRC-format ASCII; netCDF if `-Dnetcdf=true`).
- `orthogonal_basis_v1.F` / `quad_rules_general.F` — DG basis functions and quadrature rules.
- `p_enrichment.F` — hp-adaptive polynomial order selection.

### Preprocessor flags

| Flag | Purpose |
|---|---|
| `CMPI` | Enables MPI parallel code paths (set by `-Dparallel=true`) |
| `RKSSP` | Always defined; selects SSP Runge-Kutta timestepping |
| `LINUX` | Always defined; platform portability |
| `SLOPE5` | Always defined; enables slope limiter variant |
| `SWAN` | Couples to the SWAN wave model |
| `VF` | Intel Visual Fortran on Windows (not used in current CI) |

### Parallel workflow

Domain decomposition uses METIS (via `adcprep`). The parallel executable (`dgswem`) runs one MPI rank per subdomain; each rank reads from its `PE****` subdirectory. `adcpost` reassembles the per-rank output files.

## Code Style

- Fixed-form Fortran (`.F`) uses the 132-column line length limit.
- Free-form Fortran (`.F90`) is linted by `fortitude` and formatted by `fprettify` via pre-commit hooks.
- Pre-commit checks run on `.F90`, `.py`, and `.cpp` files only (not `.F`).

To install and run pre-commit locally:
```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

## TACC Systems

**Lonestar6 (GNU):** load `gcc/11.2.0`, `mvapich2/2.3.7`, `TACC`.

**Vista (NVIDIA GPU):** load `nvidia/24.7`, `openmpi/5.0.5`, `TACC`; build with `-Dgpu=true`; check unified memory support first with `nvidia-smi -q | grep -i 'addressing mode'`. Run with one MPI rank per GPU (`--tasks-per-node=1`).

On TACC, set `PKG_CONFIG_PATH` to the netCDF pkg-config path before `meson setup` (the CI script `build_gcc.sh` shows the conda-based approach).
