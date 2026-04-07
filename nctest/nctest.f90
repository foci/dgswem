program nctest
    
    use netcdf_fort, only fort_63_nc
    
    implicit none
    
    Type(fort_63_nc) :: my_fort_63_nc
        
    my_fort_63_nc%init("fort.63.nc")
    
end program nctest