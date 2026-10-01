---
trigger: always_on
---

# Machine specifics
Use the modern intel library (ifx, mpiifx) where applicable. Use Spack to install dependencies, do not use any other package manager. Always use the 'dagswem' spack environment in any shell or subprocess, by activating with 'spacktivate dagswem'. For debugging, you can use the intel gdb debugger `gdb-oneapi`.

#Reference material
`/home/wukenton/starpu.pdf` has the complete StarPU reference handbook. 
`/home/wukenton/dagswem/examples` for some StarPU test codes
`/home/wukenton/starpu-1.4.12/examples/native_fortran` for some more examples from StarPU

#StarPU debugging
`/home/wukenton/starpu-1.4.12/tools/gdbinit` provides some gdb tools to help debug StarPU, which is activated with `(gdb) source gdbinit`. `(gdb) help starpu` to view a description of all the tools. 