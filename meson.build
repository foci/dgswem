project('dagswem',['cpp','fortran'], meson_version: '>=1.0.0')

# DEVELOPER NOTES 
  # Many common flags are implictly handled by Meson. Please review them here, before you
  # add custom logic. 
  # https://mesonbuild.com/Builtin-options.html

  # Future ideas (search WIP) for development to improve robustness: 
  # - Add options for single vs double precision builds (REAL4 vs REAL8)
  # - Add more comprehensive compiler detection and optimization logic
  # - Add testing framework (meson has built-in support)
  # - Add system adaptation to preprocessor forking
  # - Support legacy feature of also grabbing graphchk from metis build

# CONFIG LOGIC 
  compiler = meson.get_compiler('fortran')

  if compiler.has_argument('-132') 
    add_project_arguments('-132', language: 'fortran') 
  endif 
  if compiler.has_argument('-ffixed-line-length-132')
    add_project_arguments('-ffixed-line-length-132', language: 'fortran')
  endif 
  
  add_project_arguments('-DREAL8', language: 'fortran') #WIP Let's ignore single precision builds :) 
  add_project_arguments('-DLINUX', language: 'fortran') 
  add_project_arguments('-DCSCA', language: 'fortran')  

  # Custom optimization options
  if get_option('overfit')
    if get_option('optimization') != '3'
      error('Enable O3 optimization level before trying overfitting.')
    endif  
    # Intel 
    intel_overfit_flags = ['-xHost', '-ipo']
    foreach flag : intel_overfit_flags
      if compiler.has_argument(flag)
        add_project_arguments(flag, language: 'fortran')
      endif
    endforeach
    # GNU
    gnu_overfit_flags = ['-march=native', '-mtune=native']
    foreach flag : gnu_overfit_flags
      if compiler.has_argument(flag)
        add_project_arguments(flag, language: 'fortran')
      endif
    endforeach
  endif

  if get_option('reckless')
    if get_option('optimization') != '3'
      error('Enable O3 optimization level before trying reckless execution.')
    endif
    # Intel 
    intel_reckless_flags = ['-fp-model fast=2', '-no-prec-div']
    foreach flag : intel_reckless_flags
      if compiler.has_argument(flag)
        add_project_arguments(flag, language: 'fortran')
      endif
    endforeach
    # GNU
    gnu_reckless_flags = ['-ffast-math', '-Ofast']
    foreach flag : gnu_reckless_flags
      if compiler.has_argument(flag)
        add_project_arguments(flag, language: 'fortran')
      endif
    endforeach
  endif
    

# BUILD LOGIC 
  mpi_dep = dependency('mpi', required: true, fallback: ['mpi', 'mpi_dep'], language: 'fortran')
  metis_dep = dependency('metis', version: '>=4.0.1', required: true, fallback: ['metis', 'metis_dep'])
  
  #ADCPREP Preprocessor 
  subdir('src')

# Summarize
  summary({'bindir': get_option('bindir')}, section: 'Install locations')
  summary({'fortran compiler': compiler.get_id(), 'version': compiler.version()}, section: 'Compiler')
  summary({'c compiler': meson.get_compiler('c').get_id()}, section: 'Compiler')
  summary({'metis_used': metis_dep.type_name()}, section: 'Metis')
  summary({'mpi_used': mpi_dep.type_name()}, section: 'MPI')
  summary({'general optimization level': get_option('optimization'), 'hardware optimization': get_option('overfit')}, section: 'Optimizations')



