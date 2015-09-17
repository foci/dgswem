set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline ${ADDITIONAL_COMPILE_FLAGS}")

list(APPEND DGSWEM_FORTRAN 
  calc_normal.F90
  create_edge_data.F90
  cstart.F90
  cwrappers.F90
  dg.F90
  DG_hydro_timestep.F90
  dgswem_init.F90
  DG_timestep_advance.F90
  DG_timestep.F90
  diff45_41.F90
  ebarrier_edge_hydro.F90
  flow_edge_hydro.F90
  fort_dg_module.F90
  get_neighbors.F90
  global.F90
  hpx_read_n_domains.F90
  ibarrier_edge_hydro.F90
  init_fileunits.F90
  internal_edge_hydro.F90
  land_edge_hydro.F90
  LDG_hydro.F90
  lgd_yield_subdomain_coord.F90
  make_dirname.F90
  message_table_hpx.F90
  messenger_hpx.F90
  met_forcing.F90
  nodalattr.F90
  nodalattr_subs.F90
  numerical_flux.F90
  ocean_edge_hydro.F90
  orthobasis_area.F90
  orthobasis_edge.F90
  p_enrichment.F90
  prep_DG.F90
  prep_slopelim.F90
  quad_pts_area.F90
  quad_pts_edge.F90
  radiation_edge_hydro.F90
  read_fort_dg.F90
  read_input.F90
  rhs_dg_hydro.F90
  sizes.F90
  slopelimiter.F90
  sta_basis.F90
  tidal_potential.F90
  wetdry.F90
  wind.F90
  write_results.F90
)

add_executable(
  dgswem_cpp_serial
  dgswem_cpp.cpp
  fname.h
  fortran_declarations.hpp
  ${DGSWEM_FORTRAN}
)


add_executable(
  dgswem_lgd_serial
  dgswem_lgd.cpp
  fname.h
  fortran_declarations.hpp
  ${DGSWEM_FORTRAN}
)

#set_property(TARGET dgswem_cpp_serial PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries(dgswem_cpp_serial gfortran)
target_link_libraries(dgswem_lgd_serial gfortran)
target_link_libraries(dgswem_lgd_serial geodecomp)

#target_link_libraries(dgswem_lgd_serial libgeodecomp)