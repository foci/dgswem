DG-SWEM
=========
Discontinuous Galerkin Shallow Water Equation Model

For reference, see

- Kubatko, Ethan J., Joannes J. Westerink, and Clint Dawson. 2006. “Hp Discontinuous Galerkin Methods for Advection Dominated Problems in Shallow Water Flow.” Computer Methods in Applied Mechanics and Engineering 196 (1): 437–51.
- Dawson, Clint, Ethan J. Kubatko, Joannes J. Westerink, Corey Trahan, Christopher Mirabito, Craig Michoski, and Nishant Panda. 2011. “Discontinuous Galerkin Methods for Modeling Hurricane Storm Surge.” Advances in Water Resources 34 (9): 1165–76.
- Wichitrnithed, Chayanon, Eirik Valseth, Ethan J. Kubatko, Younghun Kang, Mackenzie Hudson, and Clint Dawson. 2024. “A Discontinuous Galerkin Finite Element Model for Compound Flood Simulations.” Computer Methods in Applied Mechanics and Engineering 420 (February): 116707.

Building
---
DG-SWEM uses the [Meson build system](https://mesonbuild.com/) for configuration and compilation. It can be installed through `pip`:

    pip install meson

As a first example, first enter the root directory of this repo and run

    FC=gfortran CC=gcc CXX=g++ meson setup build-gcc

to create a build directory called `build-gcc` and use GNU compilers to build single-core DG-SWEM.

To also build the parallel DG-SWEM code, specify the MPI wrapper and enable the parallel option

    FC=gfortran CC=gcc CXX=g++ MPIFC=mpif90 meson setup build-gcc -Dparallel=true

Other compilers currently supported are:

- Intel LLVM (`ifx`, `icx`, `icpx`, `mpiifx`)
- NVIDIA HPC (`nvfortran`, `nvc`, `nvc++`, `mpif90`)

For systems with multiple compilers, it may be best to specify the full path to the desired compilers.

Options are specified in the format `-D<option>=<value>`. To see the available options for the above build, run

    meson configure build-gcc

Options specific to DG-SWEM are:

| Option name | Description | Possible values | Default |
| --- | --- | --- | --- |
| gpu | Enable GPU build through OpenACC | true, false | false |
| parallel | Build parallel DG-SWEM, adcprep, and adcpost | true, false | false |
| precision | Floating point precision | single, double | double |

Once configured, we can compile the program in the build directory:

    cd build-gcc 
    meson compile -v   # show compilation commands

This will produce the executable `dgswem_serial`, and if `-Dparallel=true`, also `dgswem`, `adcprep`, and `adcpost`. 

To redo a clean build, say with a different configuration, go back to the root directory and run

    meson setup build-gcc [OPTIONS] --wipe

and repeat the process.

***Notes on building the GPU version***

DG-SWEM (both serial and parallel) can be configured to run on NVIDIA gpus through OpenACC using `-Dgpu=true`, but requires NVHPC and unified memory support. 
To find out if your system supports this, run
```
nvidia-smi -q | grep -i 'addressing mode'
```
If the output is `HMM` or `ATS`, then the GPU code can be built. Some systems like Grace Hopper support the latter by default. 
To enable `HMM` on a system with a consumer GPU like the RTX series, read https://forums.developer.nvidia.com/t/issue-activating-hmm-feature-on-nvidia-rtx-a4500-with-cuda-toolkit-12-4-on-debian-bookworm/285142/3 .


Testing
---
Several test cases are contained in `tests`. This requires `pytest`. After having compiled the executables as above, go inside `tests` and run
```
pytest test_run.py --binpath <build_dir>
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
