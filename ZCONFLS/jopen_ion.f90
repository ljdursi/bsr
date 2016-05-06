!====================================================================
      Integer Function jopen_ion (e,ilsp)
!====================================================================
!     number of open channels
!--------------------------------------------------------------------

      USE target_ion
      USE channels_ion

      Implicit  none

      Real(8) :: e
      Integer :: ilsp,i

      jopen_ion = 0
      Do i = 1,nch(ilsp)
       if(e.lt.etarg(iptar(ilsp,i))) Exit
       jopen_ion = jopen_ion + 1
      End do

      End Function jopen_ion