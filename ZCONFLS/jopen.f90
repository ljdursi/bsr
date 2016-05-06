!====================================================================
      Integer Function jopen (e,ilsp)
!====================================================================
!     number of open channels in the jk-case
!--------------------------------------------------------------------
      USE target
      USE channels

      Implicit  none
      Real(8) :: e
      Integer :: ilsp,i

      jopen = 0
      Do i = 1,nch(ilsp)
       if(e.lt.etarg(iptar(ilsp,i))) Exit
       jopen = jopen + 1
      End do

      End Function jopen

