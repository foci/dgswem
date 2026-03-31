project('dg-swem', ['c', 'cpp', 'fortran'],
        version: '4.0.0',
        default_options: ['buildtype=debugoptimized'],
        meson_version: '>=1.5.0')

add_project_arguments('-DRKSSP',
                      '-DLINUX',
                      '-DSLOPE5',
                      language: 'fortran')

compiler_name = meson.get_compiler('fortran').get_id()
compiler = meson.get_compiler('fortran')
message('Using compiler', compiler_name)

if get_option('precision') == 'double'
  sz = 8
else
  sz = 4
endif 

extra_adcprep_args = []

if compiler_name == 'gcc'
  if get_option('gpu')
    warning(f'GPU not supported for @compiler_name@. Compiling for CPUs instead.')
  endif

  add_project_arguments('-ffixed-line-length-132', language: 'fortran')
  if sz == 8
    add_project_arguments('-fdefault-real-8', 
                          '-fdefault-double-8',
                          '-fallow-argument-mismatch',
                          '-march=native',
                          language: 'fortran')
  endif
  extra_adcprep_args = ['-Wno-incompatible-pointer-types']

elif compiler_name == 'nvidia_hpc'
  add_project_arguments('-Mextend',
                        f'-r@sz@',
                        language: 'fortran')
  if get_option('gpu')
    add_project_arguments('-cuda',
                          '-acc',
                          '-gpu=mem:unified,lineinfo',
                          '-Minfo=accel',
                          language: 'fortran')
    add_project_link_arguments('-cuda',
                               '-acc',
                               '-gpu=mem:unified',
                               language: 'fortran')
  endif
elif compiler_name == 'intel-llvm'
  if get_option('gpu')
    warning(f'GPU not supported for @compiler_name@. Compiling for CPUs instead.')
  endif

  add_project_arguments(f'-r@sz@',
                        '-xHost',
                        '-132',
                        language: 'fortran')
else
  error(f'Compiler @compiler_name@ not supported.')
endif

if get_option('netcdf')
  netcdf_dep = [dependency('netcdf', language: 'fortran'),
               dependency('netcdf', language: 'c')]
  compiler.run(
  '''
  PROGRAM netCDF3Test
    USE NETCDF
    IMPLICIT NONE
    INTEGER :: IERR, NCID
    IERR = NF90_OPEN('test.nc', NF90_NOWRITE, NCID)
  END PROGRAM ''' ,
  required: true,
  dependencies: netcdf_dep)
else
  netcdf_dep = []
endif
    


subdir('src')
subdir('prep')

# Serial DGSWEM Main Program
executable('dgswem_serial',
           sources: dgswem_src + c_interface_src,
           dependencies: [cpp_dep] + netcdf_dep,
           link_language: 'fortran')

# Parallel DGSWEM Main Program
if get_option('parallel')
  mpi_dep = dependency('mpi', language: 'fortran')
  executable('dgswem', sources: dgswem_src + mpi_handlers + c_interface_src,
             dependencies: [mpi_dep] + cpp_dep + netcdf_dep,
             fortran_args: '-DCMPI',
             link_language: 'fortran')
endif

executable('adcprep', sources: prep_src,
           dependencies: [metis_dep] + netcdf_dep,
           c_args: extra_adcprep_args,
           link_language: 'fortran',
        ) #have to specify link_language to avoid meson trying to link with C compiler

executable('adcpost',
           sources: post_src,
           dependencies: netcdf_dep,
          link_language: 'fortran')
