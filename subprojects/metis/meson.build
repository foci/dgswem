project('metis-subproject', 'c', version: '4.0.1')

# Retrieve the source list from subbuild file

subdir('Lib')

compiler = meson.get_compiler('c')
math_dep = compiler.find_library('m', required: true)
inc = include_directories('Lib') 

c_args = []
if compiler.get_id() == 'gcc'
    c_args = ['-Wno-incompatible-pointer-types']
endif

metis_lib = library('metis',
    sources,
    dependencies: math_dep, 
    include_directories: inc,
    install: false,
    c_args: c_args
)
metis_dep = declare_dependency(
    link_with: metis_lib,
    include_directories: inc
)
