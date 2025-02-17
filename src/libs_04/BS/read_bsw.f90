!======================================================================
      Subroutine Read_bsw(nu)
!======================================================================
!     read B-spline orbitals from bsw-file (unit nu)
!     only for orbitals which are in the list "spline_orbitals"
!----------------------------------------------------------------------
      Use spline_param
      Use spline_atomic
      Use spline_orbitals

      Implicit none
      Character elw*4
      Real(8) :: x, zw,hw,hmw,rmw
      Integer :: nu, i,ii, ksw,nsw,mw, nw,lw,kw
      Integer, external :: Ifind_nlk  

      rewind(nu)
    1 read(nu,end=2) elw,zw,hw,hmw,rmw,ksw,nsw,mw
      if(zw.ne.z) Call Stop_mpi(0,0,' Read_bsw:  z <> zw')
      if(hw.ne.h) Call Stop_mpi(0,0,' Read_bsw:  h <> hw')
      if(abs(hmw-hmax).gt.1.d-15) Call Stop_mpi(0,0,' Read_bsw:  hmw <> hmax')
      if(rmw.ne.rmax) Call Stop_mpi(0,0,' Read_bsw:  rmw <> rmax')
      if(ksw.ne.ks)   Call Stop_mpi(0,0,' Read_bsw:  ksw <> ks')
      if(nsw.ne.ns)   Call Stop_mpi(0,0,' Read_bsw:  nsw <> ns')
      Call EL4_nlk(elw,nw,lw,kw)
      ii = Ifind_nlk(nw,lw,kw,0)

      if(ii.gt.0) then
        mbs(ii)  = mw
        nbs(ii)  = nw
        lbs(ii)  = lw
        kbs(ii)  = kw
        ebs(ii)  = elw
        read(nu) pbs(1:mw,ii)
        if(mw.lt.ns) pbs(mw+1:ns,ii) = 0.d0
      else
        read(nu) (x,i=1,mw)
      end if

      go to 1
    2 Close(nu)

      End Subroutine Read_bsw

