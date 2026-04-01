export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
export FC=ifx CC=icx CXX=icpx MPIFC=mpiifx MPIF90=mpiifx
meson setup build_intel -Dparallel=true 
cd build_intel
meson compile


