program nctest
    
    use netcdf_fort, only : fort_63_nc
    
    implicit none
    
    Type(fort_63_nc) :: my_fort_63_nc
        
    print *, "Creating fort.63.nc"
    call my_fort_63_nc%init("fort.63.nc")
    print *, "Setting metadata"
    call my_fort_63_nc%set_metadata( &
        time_dimsize = 1, &
        node_dimsize = 1, &
        nele_dimsize = 1, &
        nvertex_dimsize = 1, &
        nope_dimsize = 1, &
        neta_dimsize = 1, &
        max_nvdll_dimsize = 1, &
        nbou_dimsize = 1, &
        nvel_dimsize = 1, &
        max_nvell_dimsize = 1, &
        mesh_dimsize = 1 &
    )
    print *, "Closing fort.63.nc"
    call my_fort_63_nc%close()
    print *, "Program complete"
    
end program nctest
