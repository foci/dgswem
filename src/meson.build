# Kenton: Could do a glob here but I think it's better practice to be explicit (this isn't make!) 
# Also promotes better code coverage. 
dgswem_src = files([
'DG_hydro_timestep.F',
'DG_timestep.F',
'LDG_hydro.F',
'calc_normal.F',
'couple2swan.F',
'create_edge_data.F',
'cstart.F',
'detector.F',
'dg.F',
'dgswem.F',
'diff45_41.F',
'ebarrier_edge_hydro.F', 
'edge_int_hydro.F',
'flow_edge_hydro.F',
#'flow_edge_sed.F',         looks like deprecated file from sediment transport 
'fparser.F90',
'fparser2.F90',
'global.F',
'global_3dvs.F',
'globalio.F',
'harm.F',
'hstart.F',
'ibarrier_edge_hydro.F',
'ibarrier_fluxes.F',
'internal_edge_hydro.F',
'land_edge_hydro.F',
'met_forcing.F',
#'modal2nodal.F',           looks like deprecated file from modal to nodal basis conversion? 
'nodalattr.F',
'numerical_flux.F',
'ocean_edge_hydro.F',
#'orthobasis.F',            For some reason, this file is still here, but isn't picked up by the makefile? bizzare. 
'orthogonal_basis_v1.F',
'owi_rain.F',
'owiwind.F',
'p_enrichment.F',
'parameters.f90',
#'pdg_debug.F',             Looks like a debugging file that isn't used in the main code?
'precipitation.F',
'prep_DG.F',
'prep_slopelim.F',
'quad_rules_general.F',
'quadrature_v1.F',
'radiation_edge_hydro.F', 
'read_fort_dg.F90',
'read_input.F',
'rhs_dg_hydro.F',
'sizes.F',
'slopelimiter.F',
'sta_basis.F',
#'sta_location.F',          looks like a deprecated file for station locations? 
'tidal_potential.F',
'vortex.F',
'wetdry.F',
'wind.F',
'write_output.F',
'write_results.F'])

mpi_handlers = files(['messenger.F', 'messenger_elem.F'])
c_interface_src = files(['count_edges.cpp', 'c_interface.F90'])

cpp = meson.get_compiler('cpp')
# This should find the C++ standard library 90% of the time. 
# The Fortran cross-language linking doesn't link in the C++ standard library by default, 
# (only the C standard library) so we have to do it manually here.
cpp_dep = cpp.find_library('stdc++', required: true)
if not cpp_dep.found()
    error('Could not find C++ standard library using the name stdc++. You may need to specify the name using a native file.')
endif

