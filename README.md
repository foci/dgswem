dgswem release 11.13
=========
Discontinuous Galerkin Shallow Water Equation Model

Building the CPU version
--------
Compilation should be done in the `work/` directory. The program can be compiled
with several compilers based on the environmental variable `FC`; check the file
`cmplrflags.mk` for details. Floating point precision is set using the `sz` variable.
```
cd work/
export FC=ifx        # Intel compilers
export FC=nvfortran  # NVHPC compilers
export FC=           # Else use GNU compilers
make all sz=8 -j4    # Use double-precision (default)
make all sz=4 -j4    # Use single-precision
```
To build everything from scratch, use
```
make clobber
```
Building the GPU version
---
DG-SWEM also runs on NVIDIA gpus through OpenACC, but requires NVHPC and unified memory support. 
To find out if your system supports this, run
```
nvidia-smi -q | grep -i 'addressing mode'
```
If the output is `HMM` or `ATS`, then the GPU code can be built. Some systems like Grace Hopper support the latter by default. 
To enable `HMM` on a system with a consumer GPU like the RTX series, read https://forums.developer.nvidia.com/t/issue-activating-hmm-feature-on-nvidia-rtx-a4500-with-cuda-toolkit-12-4-on-debian-bookworm/285142/3 .

To compile, run
```
cd work/
export FC=nvfortran
make all sz=8 gpu=1 -j4  # Use double-precision (default)
make all sz=4 gpu=1 -j4  # Use single-precision
```
Testing
---
Several test cases are contained in `work/test`. This requires pytest. After having compiled the executables as above, run
```
make test
```
Running
---
For the most part, ADCIRC input (`fort.14`, `fort.15`, etc.) will work minus some features like netCDF support. 
Consult https://adcirc.org/home/documentation/users-manual-v53/input-file-descriptions/ for their formats.

An additional control file called
`fort.dg` is required; a sample can be found in `work/`. An important option is the `rainfall` flag.

- `rainfall = 0` : No rain
- `rainfall = 1` : Constant 1 inch/hour rain throughout the whole domain
- `rainfall = 2` : Use given rainfall data in OWI format
- `rainfall = 3` : Use the R-CLIPER parametric model to compute rain (only `NWS=20` is supported)
- `rainfall = 4`: Use the IPET parametric model to compute rain (only `NWS=20` is supported)

**Serial run**

To run on a single CPU or GPU, use
```
./dgswem_serial
```
in the same directory as the input files. Consult https://adcirc.org/home/documentation/users-manual-v53/output-file-descriptions/ for
the output file formats.

**Parallel run**

To run on multiple CPUs or GPUs, first perform domain decomposition

    ./adcprep
which will ask for the number of MPI ranks and the names of the input files.
This creates multiple `PE****` folders, one per rank.
Once finished, we can run

    mpirun -np <N> ./dgswem

Finally, to agglomerate the subdomain-specific output files into global ones we need to run

    ./adcpost

Specific instructions for TACC systems
---
**Lonestar6**

Modules used (for compiling with GNU):

- `gcc/11.2.0`
- `mvapich2/2.3.7`
- `TACC`

**Vista**

Modules used:

- `nvidia/24.7`
- `openmpi/5.0.5`
- `TACC`

To run on multiple GPU nodes, set `-N` to be the number of GPUs and `--tasks-per-node` to 1; this
assigns one CPU core per GPU. 


dgswem Software License
---    
DG-SWEM - The Discontinuous Galerkin Shallow Water Equation Model

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
