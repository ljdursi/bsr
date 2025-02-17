!======================================================================
!     PROGRAM       B S R _ P O L                          version. 3                        
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!    INPUT  ARGUMENTS:
!
!     klsp  - the indexes of partial waves under consideration  
!
!    INPUT FILES:
!
!     bsr_par       -  description of partial waves
!     cfg.nnn       -  c-file for close-coupling expansion
!     bsr_mat.nnn   -  interaction matrix
!     dv.nnn        -  dipole matrix
!
!    OUTPUT FILES:
!   
!     bsr_pol.log   -  running information 
!     pol.nnn       -  bound-like solusions for polarized pseudo-states
!
!     Above, nnn indicates index of the partial wave
!
!=====================================================================
      Use bsr_pol
      
      Implicit none
      Real(8) :: t1,t2

      Call inf_bsr_pol

      Call CPU_time(t1)

! ... read data:

      Call Read_data
      
! ... read interaction matrix:

      Call Read_bsrmat

! ... read transition matrix:

      Call Read_dipmat

! ... additional orthogonality

      Call Read_nortb

! ... solve the dipole equation:

      Call Solv_mat

      Call CPU_time(t2)
      write(pri,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'
      write(*  ,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'

      End  ! program BSR_POL


