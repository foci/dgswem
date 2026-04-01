export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
meson setup build_gcc -Dparallel=true -Dnetcdf=true
cd build_gcc
meson compile
