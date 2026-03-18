prep_src = files([
'adcprep.F',
'decomp.F',
'interp.F',
'machdep.F',
'metis.F',
'metis2.F',
'picomsg.F',
'pre_global.F',
'prep.F',
'presizes.F',
'read_global.F'
])

inc = include_directories('inc')

executable('adcprep', sources: prep_src, dependencies: [metis_dep, mpi_dep], include_directories: inc, 
        link_language: 'fortran',
        ) #have to specify link_language to avoid meson trying to link with C compiler
