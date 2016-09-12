subroutine hpx_read_n_domains(n_domains)
  implicit none
  ! This is to read the fort.80 file to get the number of domains
  ! at the beginning of the execution.

  integer n_domains

  !
  INTEGER IDUM80
  CHARACTER CDUM80
   

  !     Read in number of domains from fort.80 file
  !OPEN(80,FILE='fort.80')
  !READ(80,'(A)') CDUM80     !Skip global_here%RUNDES
  !READ(80,'(A)') CDUM80     !Skip global_here%RUNID
  !READ(80,'(A)') CDUM80     !Skip global_here%AGRID
  !READ(80,*) IDUM80         !Skip NELG & NNODG
  !READ(80,*) IDUM80         !Read in NPROC
  !CLOSE(80)
  !n_domains = IDUM80

  n_domains = 1

end subroutine hpx_read_n_domains
