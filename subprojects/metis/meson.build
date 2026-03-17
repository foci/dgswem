project('metis-subproject', 'c', version: '4.0.1')

# Retrieve the source list from subbuild file

subdir('lib')

math_dep = meson.get_compiler('c').find_library('m', required: true)
inc = include_directories('inc') 

metis_lib = library('metis',
    sources,
    dependencies: math_dep, 
    include_directories: inc,
    install: false
)
metis_dep = declare_dependency(
    link_with: metis_lib,
    include_directories: inc
)