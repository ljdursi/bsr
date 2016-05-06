!=====================================================================
      MODULE param_LS
!=====================================================================
!     main parameter in library ZCONF_LS
!---------------------------------------------------------------------

      Implicit none

      Integer, parameter :: msh = 25     ! max. number of shells behind core
      Integer, parameter :: ibc = 2**15  ! packing parameter
      Integer, parameter :: noper=7      ! different operators
      Integer, parameter :: ksmax=61*61  ! max. set number,  see ASET in ELF4
       
      END MODULE param_LS


