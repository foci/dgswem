program nctest
    
    use ncdg, only : ncdg_63
    
    implicit none
    
    type(ncdg_63) :: my_63
        
    print *, "Creating fort.63.nc"
    call my_63%init("fort.63.nc")
    print *, "Setting metadata"
    call my_63%ncdg_63_set_metadata( &
        time_dimsize = 2, &
        node_dimsize = 3, &
        nele_dimsize = 4, &
        nvertex_dimsize = 5, &
        nope_dimsize = 6, &
        neta_dimsize = 7, &
        max_nvdll_dimsize = 8, &
        nbou_dimsize = 9, &
        nvel_dimsize = 10, &
        max_nvell_dimsize = 11, &
        mesh_dimsize = 12 &
    )
    print *, "Closing fort.63.nc"
    call my_63%close()
    print *, "Program complete"
    
end program nctest
