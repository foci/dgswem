if get_option('netcdf')
    ncdg_src = files([
        'ncutil.f90',
        'ncdg.f90',
        'ncdg_file_sub.f90',
        'ncdg_nodal_sub.f90',
        'ncdg_63_sub.f90',
        'ncdg_64_sub.f90',
        'ncdg_maxele_sub.f90',
        'ncdg_station_sub.f90'
    ])
    if get_option('parallel')
    piodg_src = files([
        'pioutil.f90'
    ])
    endif
endif
