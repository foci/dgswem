subroutine init_fileunits(s)
  use sizes
  implicit none
  
  type (sizes_type) :: s
  
  !      fort61unit=61*100+s%myproc
  ! multiply myproc by 1000, add file unit
  
  s%fort12unit = s%myproc*1000+12

  s%fort14unit = s%myproc*1000+14
  s%fort15unit = s%myproc*1000+15

  s%fortdgunit = s%myproc*1000+25

  s%fort16unit = s%myproc*1000+16

  s%fort19unit = s%myproc*1000+19
  s%fort20unit = s%myproc*1000+20
  s%fort22unit = s%myproc*1000+22
  s%fort23unit = s%myproc*1000+23
  s%fort24unit = s%myproc*1000+24
  

  s%fort80unit = s%myproc*1000+80
  s%fort81unit = s%myproc*1000+81
  s%fort82unit = s%myproc*1000+82

  s%fort71unit = s%myproc*1000+71
  s%fort72unit = s%myproc*1000+72
  s%fort73unit = s%myproc*1000+73
  s%fort74unit = s%myproc*1000+74

  s%fort83unit = s%myproc*1000+83
  s%fort84unit = s%myproc*1000+84
  s%fort85unit = s%myproc*1000+85

  s%fort88unit = s%myproc*1000+88
  s%fort89unit = s%myproc*1000+89

  s%fort94unit = s%myproc*1000+94
  s%fort96unit = s%myproc*1000+96

  s%fort4lunit = s%myproc*1000+4441

  s%rads64unit = s%myproc*1000+164

#ifdef HPX
  s%fortdgunit = s%myproc*1000+18
#endif

end subroutine init_fileunits
