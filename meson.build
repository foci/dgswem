project('dg-swem', ['c', 'cpp', 'fortran'],
        version: '4.0.0',
        default_options: ['buildtype=debugoptimized'],
        meson_version: '>=1.5.0')



#####################################
# SET UP DEPENDENCIES AND ARGUMENTS #
#####################################

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

# Get hostname
hostname_cmd = run_command('hostname',  check: true)
hostname = hostname_cmd.stdout().strip()

if get_option('netcdf')
  if not hostname.contains('tacc.utexas.edu')
    netcdf_dep = [dependency('netcdf', language: 'fortran', required: true),
                 dependency('netcdf', language: 'c', required: true)]
  else
    # Problem: pkg-config --cflags netcdf-fortran returns blank -I flag
    # Fix: filter it out manually
    nf_flags_raw = run_command('nf-config', '--fflags', check: true).stdout().strip().split()
    nf_libs_raw = run_command('nf-config', '--flibs', check: true).stdout().strip().split()
    nc_flags_raw = run_command('nc-config', '--cflags', check: true).stdout().strip().split()
    nc_libs_raw = run_command('nc-config', '--libs', check: true).stdout().strip().split()

    # Filter out the empty '-I'
    clean_flags = []
    foreach flag : nf_flags_raw + nc_flags_raw
      if flag != '-I' and flag != '' and flag not in clean_flags
        clean_flags += flag
      endif
    endforeach

    netcdf_dep = declare_dependency(
      compile_args: clean_flags,
      link_args: nf_libs_raw + nc_libs_raw
    )
  endif

  # Temporary fix for netcdf pkg-config bug on TACC systems,
  # skip netcdf test program

  if not hostname.contains('tacc.utexas.edu')
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
  endif

else
  netcdf_dep = []
endif



#############################
# DEFINE SOURCE DIRECTORIES #
#############################

subdir('src')
subdir('prep')

# NetCDF test program
if get_option('netcdf')
    subdir('netcdf')
    subdir('nctest')
endif



######################
# DEFINE EXECUTABLES #
######################

# Serial DGSWEM Main Program
if get_option('netcdf')
    executable(
        'dgswem_serial',
        sources: dgswem_src + c_interface_src + ncdg_src + dgswem_netcdf_src,
        dependencies: [cpp_dep] + netcdf_dep,
        link_language: 'fortran',
        fortran_args: '-DNETCDF'
    )
else
    executable(
        'dgswem_serial',
        sources: dgswem_src + c_interface_src,
        dependencies: [cpp_dep] + netcdf_dep,
        link_language: 'fortran'
    )
endif

# Parallel DGSWEM Main Program
if get_option('parallel')
    mpi_dep = dependency('mpi', language: 'fortran')
    if get_option('netcdf')
        executable(
            'dgswem',
            sources: dgswem_src + mpi_handlers + c_interface_src + ncdg_src + dgswem_netcdf_src,
            dependencies: [mpi_dep] + cpp_dep + netcdf_dep,
            link_language: 'fortran',
            fortran_args: ['-DCMPI', '-DNETCDF']
        )
      else
        executable(
            'dgswem',
            sources: dgswem_src + mpi_handlers + c_interface_src,
            dependencies: [mpi_dep] + cpp_dep,
            fortran_args: '-DCMPI',
            link_language: 'fortran'
        )
    endif
endif

# Pre-processing program
if get_option('netcdf')
    executable(
        'adcprep',
        sources: prep_src + ncdg_src + prep_netcdf_src,
        dependencies: [metis_dep] + netcdf_dep,
        c_args: extra_adcprep_args,
        link_language: 'fortran', # have to specify link_language to avoid meson trying to link with C compiler
        fortran_args: ['-DNETCDF']
)
else
    executable(
        'adcprep',
        sources: prep_src,
        dependencies: [metis_dep],
        c_args: extra_adcprep_args,
        link_language: 'fortran', # have to specify link_language to avoid meson trying to link with C compiler
    )
endif

# Post-processing program
executable(
    'adcpost',
    sources: post_src,
    # dependencies: netcdf_dep,
    link_language: 'fortran'
)

# NetCDF test program
if get_option('netcdf')
    executable(
        'nctest.exe',
        sources: ncdg_src + nctest_src,
        dependencies: netcdf_dep,
        link_language: 'fortran'
    )

    if get_option('parallel')
        executable(
            'nctest_parallel.exe',
            sources: ncdg_src + nctest_parallel_src,
            dependencies: [mpi_dep] + netcdf_dep,
            link_language: 'fortran',
            fortran_args: '-DCMPI',
        )
    endif
endif
