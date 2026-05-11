module nctest_util

    use ncdg, only: ncdg_63
    
    implicit none

    contains

    subroutine print_metadata(file)
        type(ncdg_63), intent(in) :: file

        print *, "Printing metadata:"

        print *, "DGSWEM Dimension IDs"
        print *, file%time_dimname, file%time_dimid
        
        print *, "Nodal Dimension IDs"
        print *, file%node_dimname, file%node_dimid
        print *, file%nele_dimname, file%nele_dimid
        print *, file%nvertex_dimname, file%nvertex_dimid
        print *, file%nope_dimname, file%nope_dimid
        print *, file%neta_dimname, file%neta_dimid
        print *, file%max_nvdll_dimname, file%max_nvdll_dimid
        print *, file%nbou_dimname, file%nbou_dimid
        print *, file%nvel_dimname, file%nvel_dimid
        print *, file%max_nvell_dimname, file%max_nvell_dimid
        print *, file%mesh_dimname, file%mesh_dimid
        
        print *, "DGSWEM Variable IDs"
        print *, file%time_varname, file%time_varid
        
        print *, "Nodal Variable IDs"
        print *, file%x_varname, file%x_varid
        print *, file%y_varname, file%y_varid
        print *, file%element_varname, file%element_varid
        print *, file%adcirc_mesh_varname, file%adcirc_mesh_varid
        print *, file%neta_varname, file%neta_varid
        print *, file%nvdll_varname, file%nvdll_varid
        print *, file%max_nvdll_varname, file%max_nvdll_varid
        print *, file%ibtypee_varname, file%ibtypee_varid
        print *, file%nbdv_varname, file%nbdv_varid
        print *, file%nvel_varname, file%nvel_varid
        print *, file%nvell_varname, file%nvell_varid
        print *, file%max_nvell_varname, file%max_nvell_varid
        print *, file%ibtype_varname, file%ibtype_varid
        print *, file%nbvv_varname, file%nbvv_varid
        print *, file%depth_varname, file%depth_varid

        print *, "fort.63.nc Variable IDs"
        print *, file%zeta_varname, file%zeta_varid
    end subroutine print_metadata
    
end module nctest_util



program nctest

    use nctest_util
    use ncdg, only : ncdg_63
    
    implicit none
    
    type(ncdg_63) :: my_63
        
    print *, "Creating fort.63.nc"
    ! call print_metadata(my_63)
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
    ! call print_metadata(my_63)
    print *, "Closing fort.63.nc"
    call my_63%close()
    ! call print_metadata(my_63)
    print *, "Opening fort.63.nc"
    call my_63%open()
    ! call print_metadata(my_63)
    print *, "Closing fort.63.nc"
    call my_63%close()
    print *, "Program complete"
    
end program nctest
