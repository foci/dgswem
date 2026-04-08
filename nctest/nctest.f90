program nctest
    
    use netcdf_fort, only : fort_63_nc
    
    implicit none
    
    Type(fort_63_nc) :: my_fort_63_nc
        
    print *, "Creating fort.63.nc"
    call my_fort_63_nc%init("fort.63.nc")
    print *, "Closing fort.63.nc"
    call my_fort_63_nc%close()
    print *, "Program complete"
    
end program nctest
