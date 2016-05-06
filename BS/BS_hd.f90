!====================================================================
    MODULE spline_atomic
!====================================================================
!   contains some atomic parameters
!--------------------------------------------------------------------
    
    IMPLICIT NONE

    REAL(KIND=8) :: z  = 1.d0   ! nuclear charge
    REAL(KIND=8) :: EC = 0.d0    ! core energy

    REAL(KIND=8) :: fine = 0.25D0/(137.036D0)**2

    LOGICAL :: rel = .FALSE.    ! relativistic corrections
    INTEGER :: irel  =  0       ! the same
    INTEGER :: kclosd = 0       ! closed shells
    INTEGER :: MASS =   0       ! mass-corrections
    INTEGER :: ioo =    0       ! orbit-orbit interaction

    END MODULE spline_atomic

!====================================================================
    MODULE spline_param
!====================================================================
!   contains basic spline parameters
!--------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(4) :: grid_type = 1   !  type of grid
    INTEGER(4) :: nug = 99 
    Character(40) :: AF_grid = 'knot.dat'

    INTEGER(4) :: ks = 8  !   order of B-splines
    INTEGER(4) :: ns = 0  !   number of splines
    INTEGER(4) :: nv = 0  !   number of intervals ( = ns-ks+1 )
    INTEGER(4) :: ml = 0  !   number of intervals from 0 to 1 (=1/h)
    INTEGER(4) :: me = 0  !   number of intervals in the exponential region

    REAL(8) :: h = 0.125    !   initial step in the knot sequence for z*r
    REAL(8) :: hmax = 1.00  !   maximum step, t(ns+1) - t(ns) 
    REAL(8) :: rmax = 40.00 !   border radius, t(ns+1)

    END MODULE spline_param


!====================================================================
     MODULE spline_grid
!====================================================================
!
!    defines the values of splines at the gaussian points for each
!    interval of a grid; included in the module is the gaussian data
!    for performing integrations on the grid
!
!--------------------------------------------------------------------

     IMPLICIT NONE

! .. knot sequence, t(1:ns+ks)

     REAL(8), DIMENSION(:), ALLOCATABLE:: t

! .. arrays for spline values in gausian points
!
!    bsp(1:nv+1,1:ks,1:ks), bspd(1:nv+1,1:ks,1:ks,2) 
!
!    bsp(i,m,ith)  - values of the i+ith-1 B-spline in interval i
!                    at gausian point m
!    bspd(i,m,ith,1|2)  - corresponding values of first and second
!                         derivatives 
!    bsp(nv+1,1,.) and bspd(nv+1,1,.,.) - corresponding values at
!                                         last knot point (rmax)

     REAL(8), DIMENSION(:,:,:), ALLOCATABLE::   bsp, bsq
     REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: bspd

! .. arrays for gaussian data
!
!    gr(1:nv;1:ks),  grm(1:nv;1:ks),  grw(1:nv;1:ks)  
!
!    gr(i,m)  - gaussian points m in the interval i
!    grm(i,m) - reciprocal value of gr(i,m)
!    grw(i,m) - gaussian weights at corresponding points

     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: gr, grm, grw


     CONTAINS


!====================================================================
      SUBROUTINE allocate_grid
!====================================================================
!
! ... allocates space of the arrays in MODULE spline_grid
!
!--------------------------------------------------------------------

      USE spline_param

      INTEGER(4) :: ierr
      
      if(Allocated(bsp)) Deallocate(bsp,bsq,bspd,gr,grm,grw)
      ALLOCATE(bsp(nv+1,ks,ks), bsq(nv+1,ks,ks), bspd(nv+1,ks,ks,2), &
               gr(nv,ks), grm(nv,ks), grw(nv,ks))
      bsp = 0.d0; bsq = 0.d0; bspd = 0.d0
	  gr = 0.d0; grm = 0.d0; grw = 0.d0

      END SUBROUTINE allocate_grid

    END MODULE spline_grid

!====================================================================
   MODULE spline_galerkin
!====================================================================
!
!  contains common arrays used in the application of splines and 
!  the Galerkin method 
!
!--------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: r1, rm1, rm2, rm3
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: sb, bs, db1, db2, bb

! ...  sb (1:ns,1:ks)  -  sb (i,j) -->  <B_i|B_i+j-k>
! ...  r1 (1:ns,1:ks)  -  r1 (i,j) -->  <B_i|r|B_i+j-k>
! ...  rm1(1:ns,1:ks)  -  rm1(i,j) -->  <B_i|1/r|B_i+j-k>
! ...  rm2(1:ns,1:ks)  -  rm2(i,j) -->  <B_i|1/r^2|B_i+j-k>
! ...  rm3(1:ns,1:ks)  -  rm2(i,j) -->  <B_i|1/r^3|B_i+j-k>
! ...  db1(1:ns,1:ks)  -  db1(i,j) -->  <B_i|B'_i+j-k>
! ...  db2(1:ns,1:ks)  -  db2(i,j) -->  <B_i|B'_i+j-k>
! ...  bs(1:ks,1:ns)   -  factorization of sb

! ...  all arrays (except bs) in the symmetric lower-columb storage
! ...  mode, db2 - 'almost' symmetruc, db1 - antisymmetric


    CONTAINS


!====================================================================
    SUBROUTINE allocate_galerkin
!====================================================================
!
!   allocates space for the arrays in the MODULE spline_galerkin
!
!--------------------------------------------------------------------

    USE spline_param

    INTEGER :: ierr

    if(Allocated(r1)) Deallocate(r1,rm1,rm2,rm3,sb,bs,db1,db2)
    ALLOCATE( r1(ns,ks), rm1(ns,ks), rm2(ns,ks), rm3(ns,ks), &
              sb(ns,ks), bs(ks,ns), db1(ns,ks), db2(ns,ks), bb(ns,ns) )
    r1 = 0.d0; rm1 = 0.d0; rm2 = 0.d0; sb = 0.d0; bs = 0.d0; bb = 0.d0
    db1 = 0.d0; db2 = 0.d0


    END SUBROUTINE allocate_galerkin

    END MODULE spline_galerkin
!=========================================================================
      SUBROUTINE define_grid(z)
!=========================================================================
!
!     gets input data for the grid and sets up the spline knots 
!
!     Requires file  'knot.dat' with grid parameters 
!
! -------------------------------------------------------------------------

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      Real(8) :: z
      Logical :: EX

! ... get input data for the grid

      Inquire(file=AF_grid,exist=EX)

      if(.not.EX) then
       CALL mkgrid1(z)
       Open(nug,file=AF_grid)
       Call Print_grid_01(z)
       Stop 'check knot.dat file for splines parameters ...'
      end if 

      Open(nug,file=AF_grid); Rewind(nug)

      grid_type = 0
      Call Read_ipar(nug,'grid_type',grid_type)      

      Select case(grid_type)
       
       Case(0);  Call define_grid_00(z) 
       Case(1);  Call define_grid_01(z) 
       Case default; Stop ' Unknown grid_type '

      End Select

      Close(nug)

      END SUBROUTINE define_grid


!=========================================================================
      SUBROUTINE define_grid_00(z)
!=========================================================================

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      Real(8) :: z

      rewind(nug)
      read(nug,*) ks
      read(nug,*) ns
      read(nug,*) z
      read(nug,*) h
      read(nug,*) hmax
      read(nug,*) rmax

! ... set up the knots for spline

      CALL mkgrid(z)

! ... print grid in the file 'knot.dat'

      rewind(nug)
      write(nug,'(i12,a)') ks,' ==>  order  of splines (ks)'
      write(nug,'(i12,a)') ns,' ==>  number of splines (ns)'
      write(nug,'(f12.5,a)') z,' ==>  nuclear charge (z)'
      write(nug,'(f12.5,a)') h,' ==>  step size from 0 to 1 (h for z*r, = (1/2)^n)'
      write(nug,'(f12.5,a)') hmax,' ==>  maximum step size (hmax for r)'
      write(nug,'(f12.5,a)') rmax,' ==>  maximum r (rmax)'
      write(nug,'(a)') '***'
      write(nug,'(5f12.5)') t(1:ns+ks)
      write(nug,'(a)') '***'
      write(nug,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
      write(nug,'(i12,a)') me,' ==>  number of knots in the exponential region '
      write(nug,'(i12,a)') nv,' ==>  number of intervals '
      write(nug,'(f12.5,a)') (1/hmax)**2,' ==>  max k^2 (Ry) '


      END SUBROUTINE define_grid_00


!=====================================================================
    SUBROUTINE mkgrid(z)
!=====================================================================
!
!   sets up the knots for spline (original version)
!
!---------------------------------------------------------------------

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      INTEGER(4) :: n, i, m, me1, me2, nt
      REAL(8):: z,hp1, h2, tmax, tx, h0

      ! .. determine ml, the number of distinct points from 0 to 1

       h0 = h
       ml = 1.d0/h + 0.1
       h = 1.d0/ml
       hp1 = 1.d0 + h

      ! .. determine tmax

      tmax = z*rmax
      hmax = z*hmax

      ! .. determine final point of "exponential" grid
      ! .. me: number of points from 1 to (1+h)**me
      ! .. m:  number of points from (1+h)**me to tmax

      me1 = MAX(0.0d0, LOG(hmax/h)/LOG(hp1)+1.d0)
      me2 = LOG(tmax)/LOG(hp1)+1

      IF ( me2 <= me1 ) THEN
        me = me2
        m = 0
      ELSE
        me = me1
        tx = hp1**me
        h2 = h*tx/hp1
        m = NINT((tmax-tx)/h2)
      END IF

      n = ml + me + m + ks -1
      ns = n
      nv = ns - (ks -1)
      nt = ns + ks

      ! .. establish the grid for z*r

      ALLOCATE (t(nt)); t = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h0
      END DO

      DO i = ks+ml+1, ks+me+ml
        t(i) = t(i-1)*hp1
      END DO

      DO i = ks+me+ml+1, n+1
        t(i) = t(i-1) + h2
      END DO
      t(n+2:nt) = t(n+1)

      ! .. scale the values to the R variable

      t = t/z
      hmax = hmax/z

    END SUBROUTINE mkgrid

!=========================================================================
      SUBROUTINE define_grid_01(z)
!=========================================================================

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      Real(8) :: z

      Call Read_ipar(nug,'ks',ks)      
      Call Read_ipar(nug,'ns',ns)      
      Call Read_rpar(nug,'z' ,z)      
      Call Read_rpar(nug,'h' ,h)      
      Call Read_rpar(nug,'hmax',hmax)      
      Call Read_rpar(nug,'rmax',rmax)      

! ... set up the knots for spline

      CALL mkgrid1(z)

! ... print grid in the file 'knot.dat'

      Call Print_grid_01(z)

      END SUBROUTINE define_grid_01


!=========================================================================
      SUBROUTINE Print_grid_01(z)
!=========================================================================

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      Real(8) :: z

      rewind(nug)
      write(nug,'(a,i3,T20,a)') 'grid_type  = ',grid_type,' !  type of grid '
      write(nug,*)
      write(nug,'(a,i10,T20,a)')   'ks   =',ks,' ! order of splines '
      write(nug,'(a,i10,T20,a)')   'ns   =',ns,' ! number of splines '
      write(nug,'(a,f10.5,T20,a)') 'z    =',z, ' ! nuclear charge'
      write(nug,'(a,f10.5,T20,a)') 'h    =',h, ' ! step size from 0 to 1, for z*r'
      write(nug,'(a,f10.5,T20,a)') 'hmax =',hmax,' ! maximum step size for r '
      write(nug,'(a,f10.5,T20,a)') 'rmax =',rmax,' ! maximum r '
      write(nug,'(a)') '***'
      write(nug,'(5f12.5)') t(1:ns+ks)
      write(nug,'(a)') '***'
      write(nug,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
      write(nug,'(i12,a)') me,' ==>  number of knots in the exponential region '
      write(nug,'(i12,a)') nv,' ==>  number of intervals '
      write(nug,'(f12.5,a)') (1/hmax)**2,' ==>  max k^2 (Ry) '

      END SUBROUTINE Print_grid_01


!=====================================================================
      SUBROUTINE mkgrid1(z)
!=====================================================================
!
!     sets up the knots for splines in semi-exponential grid:
!
!-----------------------------------------------------------------------

      USE spline_param
      USE spline_grid

      IMPLICIT NONE
      INTEGER(4):: i,nt
      Real(8) :: z 

      rmax = rmax * z
      hmax = hmax * z
      nt = (rmax / h + 1 + ks)*2
      if(Allocated(t)) Deallocate(t)
      ALLOCATE (t(nt))
      t = 0.d0
      ns = ks
	  
! ... determine ml, the number of intervals from 0 to 1
! ... first make sure that h = 1/n as recomended

      ml = 1.d0/h + 0.5;  h = 1.d0/ml
      Do i=1,ml 
       ns = ns + 1; t(ns) = t(ns-1) + h
      End do

! ... determine me, the number of intervals in "exponential" grid

      if(hmax.lt.h) hmax = h
      me = 0 
      Do 
       t(ns+1) = t(ns)*(1.d0+h)
       if(t(ns+1)-t(ns).gt.hmax) Exit
       ns = ns + 1
       if(t(ns).gt.rmax) Exit
       me = me + 1
      End do
        
! ... rest of interval with step = hmax

      if(t(ns).lt.rmax) then
       Do 
        t(ns+1) = t(ns) + hmax
        ns = ns + 1
        if(t(ns).ge.rmax) Exit
       End do
      end if

      if(t(ns).gt.rmax) then
       ns = ns - 1; me = me - 1
       t(ns+1:ns+ks) = rmax
       t(ns) = (t(ns+1)+t(ns-1))/2.d0
      end if

      nv = ns - ks + 1

! ... scale to the R variable

      t = t/z
      hmax = hmax / z
      rmax = rmax / z

      END SUBROUTINE mkgrid1




!======================================================================
      Subroutine DEF_BS
!======================================================================

      Use spline_atomic, only: z

      Implicit none

      CALL define_grid(z); Call define_spline   

      End Subroutine DEF_BS


!======================================================================
    SUBROUTINE define_spline
!======================================================================
!
!   initializes the values of the spline and its derivatives
!   and evaluates the spline basic arrays (elementary operators in
!   spline basis).
!
!   SUBROUTINE called:
!       gauss
!       allocate_grid
!       initvb
!       allocate_galerkin
!       initas
!
!   calling sequence:
!                        define_spline
!                   ----------------------
!                  / |          ||      ||
!                 /  |        initvb   initas
!                /   |           |     // \  \\
!           gauss    |        vbsplvd mdb mrm facsb
!                    |           ||
!       allocate_grid,galerkin vbsplvb
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin

    IMPLICIT NONE
    REAL(8), DIMENSION(ks) :: gx, gw

    ! .. initializes variables for gaussian integration

    CALL gauss(ks,gx,gw)

    ! .. initializes the values of the spline and its derivatives
    
    CALL allocate_grid
    
    CALL initvb
    
    ! .. initializes the spline array (operators in spline basis)

    Call allocate_galerkin

    CALL initas

    CONTAINS

!=======================================================================
    SUBROUTINE initvb
!=======================================================================
!
!   Sets (or Initializes) the arrays
!       gr      The gaussian points for the nint intervals of [0,Rmax]
!       grm     Reciprocals of the values of gr
!               (to avoid repeated division on the CRAY)
!       grw     Gaussian weights for each of the points in gr
!       bsp     array of B-spline values at the gaussian points
!       bspd    array of values of the first and second derivative of
!               the B-splines at the gaussian points
!
!   Calling sequence:
!        initvb
!          |
!       vbsplvd
!         ||
!       vbsplvb
!
!-----------------------------------------------------------------------
!   on entry
!   --------
!       gx      the gaussian points for the interval [0,1]
!       gw      the gaussian weights for the interval [0,1]
!
!   working arrays
!   --------------
!       dbiatx  working array of dimension (ns,ks,ks) which contains the
!               values and the second derivative values of of B-spline at
!               gaussian points in the interval [0,1]
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    REAL(8), Allocatable, DIMENSION(:,:,:) :: dbiatx
    INTEGER :: m, i
    Allocate(dbiatx(nv,ks,ks)); dbiatx = 0.d0 

    Do m=1,ks
      gr(1:nv,m)=(t(1+ks:nv+ks)-t(ks:nv+ks-1))*gx(m)+t(ks:nv+ks-1)
      grm(1:nv,m) = 1.d0/gr(1:nv,m)
      grw(1:nv,m)=(t(1+ks:nv+ks)-t(ks:nv+ks-1))*gw(m)
      Call vbsplvd(t,ks,nv,gr(1,m),3,dbiatx)
      bsp(1:nv,m,1:ks)    = dbiatx(1:nv,1:ks,1)
      bspd(1:nv,m,1:ks,1) = dbiatx(1:nv,1:ks,2)
      bspd(1:nv,m,1:ks,2) = dbiatx(1:nv,1:ks,3)
      Do i = 1,ks
       bsq(1:nv,m,i) = bspd(1:nv,m,i,1) - grm(1:nv,m)*bsp(1:nv,m,i)
      End do
    End do

    ! .. store also the values at the last knot

    call vbsplvd(t,ns,1,t(ns+1),3,dbiatx)

    bsp(nv+1,1,1:ks)    = dbiatx(1,1:ks,1)
    bspd(nv+1,1,1:ks,1) = dbiatx(1,1:ks,2)
    bspd(nv+1,1,1:ks,2) = dbiatx(1,1:ks,3)
    bsq(nv+1,1,1:ks)    = bspd(nv+1,1,1:ks,1)-bsp(nv+1,1,1:ks)/t(ns+1)

    Deallocate(dbiatx) 

    END SUBROUTINE initvb


!==================================================================
    SUBROUTINE initas
!==================================================================
!
!   Sets ( or Initializes ) the array in symmetric storage mode:
!
!       db1 --- matrix of integral <B_i,B'_j>
!       db2 --- matrix of integral <B_i,B"_j>
!       sb  --- matrix of integral <B_i,B_j>
!       r1  --- matrix of integral <B_i,r B_j>
!       rm1 --- matrix of integral <B_i,(1/r)B_j>
!       rm2 --- matrix of integral <B_i,(1/r^2)B_j>
!       rm3 --- matrix of integral <B_i,(1/r^3)B_j>
!               where i=1,..,ns, j=1,...ks
!
!   Calling sequence:
!       initas
!        /  \
!      mdb  mrm
!
!------------------------------------------------------------------


    IMPLICIT NONE


    CALL mdb1      ! .. sets db1 --- matrix of integral <B_i,B'_j>

    CALL mdb2      ! .. sets db2 --- matrix of integral <B_i,B"_j>

    CALL mrm(0, sb)              ! .. sets sb  ---  <B_i,B_j>

    Call Full_mat_sym(ns,ks,sb,bb,'l') 

    CALL facsb                   ! .. factorizes sb

    CALL mrm(1, r1)              ! .. sets r1  ---  <B_i,r B_j>

    CALL mrm(-1, rm1)            ! .. sets rm1 ---  <B_i,(1/r)B_j>

    CALL mrm(-2, rm2)            ! .. sets rm2 ---  <B_i,(1/r^2)B_j>

    CALL mrm(-3, rm3)            ! .. sets rm2 ---  <B_i,(1/r^3)B_j>





   END SUBROUTINE initas


!====================================================================
   SUBROUTINE mdb1
!====================================================================
!
!  Computes the matrix elements <B_i,B'_i+j-1> in the B-spline basis
!
!--------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables

    INTEGER :: i, irow, jcol, ith, jth

    db1 = 0.d0

    do ith = 1,ks
      do jth = 1,ith
        jcol = jth-ith+ks
        do i = 1,nv
          irow = i+ith-1
          db1(irow,jcol) = db1(irow,jcol) &
        + SUM(grw(i,:)*bsp(i,:,ith)*bspd(i,:,jth,1))
        end do
      end do
    end do

  END SUBROUTINE mdb1


!====================================================================
   SUBROUTINE mdb2
!====================================================================
!
!  Computes the matrix elements <B_i,B"_i+j-1> in the B-spline basis
!
!--------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables

    INTEGER :: i, irow, jcol, ith, jth

    db2 = 0.d0

    do ith = 1,ks
      do jth = 1,ith
        jcol = jth-ith+ks
        do i = 1,nv
          irow = i+ith-1
          db2(irow,jcol) = db2(irow,jcol) &
        + SUM(grw(i,:)*bsp(i,:,ith)*bspd(i,:,jth,2))
        end do
      end do
    end do

    db2(1,1) = db2(1,1) + SUM(grw(nv,:)*bsp(nv,:,ks-1)*bspd(nv,:,ks,2))

  END SUBROUTINE mdb2


  END SUBROUTINE define_spline
!====================================================================
    SUBROUTINE mvcv(l,vc)
!====================================================================
!
!   Computes the matrix elements for the mass-velocity correction
!   in the B-spline basis.  The lack of symmetry in the d^2/dr^2
!   operator is ignored.
!
!     VC(i,j) = INT [  (d^2/dr^2 - l(l+1)/r) B_i(r) *
!                      (d^2/dr^2 - l(l+1)/r) B_j(r)  ] dr
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       l    the angular momentum
!
!   on exit
!   -------
!       vc   the mass velocity correction in symmetric storage mode
!
!--------------------------------------------------------------------

    USE spline_param; USE spline_atomic;  USE spline_grid
    
    IMPLICIT NONE

    INTEGER(4), INTENT(in) :: l
    REAL(8), INTENT(inout), Dimension(ns,ks) :: vc

    INTEGER(4) :: m, ith, jth, i, irow, jcol
    REAL(8) :: fll, y1, y2, S, B
    Real(8), External :: AZL

    ! .. initialize the vc array

    vc = 0.d0;  fll = l*(l+1);  nv = ns-ks+1

    ! .. compute the matrix elements

    do m = 1,ks
      do i = 1,nv
        S = fll*grm(i,m)*grm(i,m)

! ... cutoff correction

!        B = gr(i,m)/(gr(i,m)+2*fine*Z);  B = B*B*B

        do ith = 1,ks
          irow = i+ith-1
          do jth = 1,ith
          jcol = jth-ith+ks

            y1 = bspd(i,m,ith,2) - S*bsp(i,m,ith)
            y2 = bspd(i,m,jth,2) - S*bsp(i,m,jth)
            vc(irow,jcol) = vc(irow,jcol) + grw(i,m)*y1*y2 !* B

          end do
        end do
      end do
    end do

    vc = vc * fine

! ... one-electron Darwin correction:

    if(l.eq.0) then
      S = azl(z,h,ks,l+1);  vc(2,ks) = vc(2,ks) - z*S*S*fine
    end if


    END SUBROUTINE mvcv


!=====================================================================
    FUNCTION bvmv(n,k,array,type,x,y)
!=====================================================================
!
!   Returns bvmv = <x,array y> for the band array.
!
!   The original size of array is n*n. The bandwidth of array is 2k-1
!   for the non-symmetric case, and is k for symmetric case.
!   The lower-band column storage mode is assumed.
!---------------------------------------------------------------------
!
!   on entry
!   --------
!       n      the leading dimension.
!       k      the band width.
!       array  a banded matrix in the lower-band column storage mode.
!              Its element (i,j) represents the element (i,i+j-k) of
!              the original matrix. For a symmetric matrix, only the
!              lower diagonals are needed.
!       type   'g' for general;
!              's' for symmetric;
!              'a' for "almost"  symmetric, as db2
!       x      vector with length n.
!       y      vector with length n.
!
!----------------------------------------------------------------------

    REAL(KIND=8) :: bvmv, ans
    INTEGER, INTENT(in) :: n, k
    CHARACTER(LEN=1), INTENT(in) :: type
    REAL(KIND=8), DIMENSION(n,*), INTENT(in) :: array
    REAL(KIND=8), DIMENSION(n), INTENT(in) :: x, y

    ! .. Local variables

    INTEGER :: i,j,jp,kp

    bvmv = 0.d0

    if ( type /= 'g') then

      ! .. symmetric

      do jp = 1,k-1
        do i = k+1-jp,n
          j = i+jp-k
          bvmv = bvmv + array(i,jp) * (x(i)*y(j) + x(j)*y(i))
        end do
      end do

      ! .. make correction for 'a' (almost symmetric)

      if (type == 'a')  &
          bvmv = bvmv + x(n-1)*y(n)*(array(1,1)-array(n,k-1))

    else

      do jp = 1,k-1
        do i = k-jp+1,n
          j = i+jp-k
          kp = 2*k-jp
          ans = ans + array(i,jp)*x(i)*y(j)
          ans = ans + array(j,kp)*x(j)*y(i)
        end do
      end do

    end if

    ! .. add central diagonal

    bvmv = bvmv + SUM(array(:,k)*x(:)*y(:))

  END FUNCTION bvmv

!=======================================================================
    SUBROUTINE mrm(mm,rm)
!=======================================================================
!
!   Computes the matrix representation of the operator
!
!             r^mm
!
!   in the B-spline basis.
!
!-----------------------------------------------------------------------
!
!   on entry
!   --------
!       mm      integer defining the power of r
!
!   on exit
!   -------
!       rm      <B_i,r^mm B_j> in symmetric storage mode
!
!-----------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ns,ks), INTENT(OUT) :: rm
    INTEGER, INTENT(IN) :: mm

    ! .. local variables

    INTEGER :: i, irow, jcol, ith, jth

    ! .. clear the rm array

    rm = 0.d0

    ! .. assemble the matrix elements

    if (mm == 0) then
      do  ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm == -1) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*grm(i,:)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm <= -2) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*grm(i,:)**(-mm)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm == 1) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
             + SUM(grw(i,:)*bsp(i,:,ith)*gr(i,:)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm >= 2) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
             + SUM(grw(i,:)*bsp(i,:,ith)*gr(i,:)**mm*bsp(i,:,jth))
          end do
        end do
      end do
    end if

  END SUBROUTINE mrm
!=======================================================================
    SUBROUTINE gauss(k,x,w)
!=======================================================================
!   Looks up the values of gaussian coordinates and guassian weights
!   for k-point gaussian quadrature over the interval [0, 1].
!
!   on entry:  k     the number of points in the quadrature
!   -------
!
!   on exit:  x(i)   Gaussian coordinates of the points
!   -------   w(i)   Gaussian weight to point x(i)
!
!   Restriction:  1<= k <= 15   FOR k POINT CASE
!   -----------
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(k), INTENT(OUT) :: x, w

    IF (k < 1 .OR. k > 15)   &
       Stop 'Error in GAUSS: number of points is out of range '

    SELECT CASE ( k )
    CASE (1)
      x(1) = .5d0
      w(1) = 1.d0
    CASE (2)
      x(1) = .211324865405187d0
      x(2) = .788675134594813d0
      w(1) = .5d0
      w(2) = .5d0
    CASE (3)
      x(1) = .112701665379258d0
      x(2) = .5d0
      x(3) = .887298334620742d0
      w(1) = .277777777777778d0
      w(2) = .444444444444444d0
      w(3) = .277777777777778d0
    CASE (4)
      x(1) = .0694318442029737d0
      x(2) = .330009478207572d0
      x(3) = .669990521792428d0
      x(4) = .930568155797026d0
      w(1) = .173927422568727d0
      w(2) = .326072577431273d0
      w(3) = .326072577431273d0
      w(4) = .173927422568727d0
    CASE (5)
      x(1) = .046910077030668d0
      x(2) = .230765344947158d0
      x(3) = .5d0
      x(4) = .769234655052842d0
      x(5) = .953089922969332d0
      w(1) = .118463442528095d0
      w(2) = .239314335249683d0
      w(3) = .284444444444444d0
      w(4) = .239314335249683d0
      w(5) = .118463442528095d0
    CASE (6)
      x(1) = .033765242898424d0
      x(2) = .169395306766868d0
      x(3) = .380690406958402d0
      x(4) = .619309593041598d0
      x(5) = .830604693233132d0
      x(6) = .966234757101576d0
      w(1) = .0856622461895852d0
      w(2) = .180380786524069d0
      w(3) = .233956967286346d0
      w(4) = .233956967286346d0
      w(5) = .180380786524069d0
      w(6) = .0856622461895852d0
    CASE (7)
      x(1) = .0254460438286207d0
      x(2) = .129234407200303d0
      x(3) = .297077424311301d0
      x(4) = .5d0
      x(5) = .702922575688699d0
      x(6) = .870765592799697d0
      x(7) = .974553956171379d0
      w(1) = .0647424830844348d0
      w(2) = .139852695744638d0
      w(3) = .19091502525256d0
      w(4) = .208979591836735d0
      w(5) = .19091502525256d0
      w(6) = .139852695744638d0
      w(7) = .0647424830844348d0
    CASE (8)
      x(1) = .0198550717512319d0
      x(2) = .101666761293187d0
      x(3) = .237233795041835d0
      x(4) = .408282678752175d0
      x(5) = .591717321247825d0
      x(6) = .762766204958164d0
      x(7) = .898333238706813d0
      x(8) = .980144928248768d0
      w(1) = .0506142681451881d0
      w(2) = .111190517226687d0
      w(3) = .156853322938944d0
      w(4) = .181341891689181d0
      w(5) = .181341891689181d0
      w(6) = .156853322938944d0
      w(7) = .111190517226687d0
      w(8) = .0506142681451881d0
    CASE (9)
      x(1) = .015919880246187d0
      x(2) = .0819844463366821d0
      x(3) = .193314283649705d0
      x(4) = .337873288298095d0
      x(5) = .5d0
      x(6) = .662126711701904d0
      x(7) = .806685716350295d0
      x(8) = .918015553663318d0
      x(9) = .984080119753813d0
      w(1) = .0406371941807872d0
      w(2) = .0903240803474287d0
      w(3) = .130305348201468d0
      w(4) = .156173538520001d0
      w(5) = .16511967750063d0
      w(6) = .156173538520001d0
      w(7) = .130305348201468d0
      w(8) = .0903240803474287d0
      w(9) = .0406371941807872d0
    CASE (10)
      x(1) = .0130467357414141d0
      x(2) = .0674683166555077d0
      x(3) = .160295215850488d0
      x(4) = .283302302935376d0
      x(5) = .425562830509184d0
      x(6) = .574437169490816d0
      x(7) = .716697697064624d0
      x(8) = .839704784149512d0
      x(9) = .932531683344492d0
      x(10)= .986953264258586d0
      w(1) = .0333356721543441d0
      w(2) = .0747256745752903d0
      w(3) = .109543181257991d0
      w(4) = .134633359654998d0
      w(5) = .147762112357376d0
      w(6) = .147762112357376d0
      w(7) = .134633359654998d0
      w(8) = .109543181257991d0
      w(9) = .0747256745752903d0
      w(10)= .0333356721543441d0
    CASE (11)
      x(1) = .0108856709269715d0
      x(2) = .0564687001159523d0
      x(3) = .134923997212975d0
      x(4) = .240451935396594d0
      x(5) = .365228422023827d0
      x(6) = .5d0
      x(7) = .634771577976172d0
      x(8) = .759548064603406d0
      x(9) = .865076002787025d0
      x(10)= .943531299884048d0
      x(11)= .989114329073028d0
      w(1) = .0278342835580868d0
      w(2) = .0627901847324523d0
      w(3) = .0931451054638672d0
      w(4) = .116596882295995d0
      w(5) = .131402272255123d0
      w(6) = .13646254338895d0
      w(7) = .131402272255123d0
      w(8) = .116596882295995d0
      w(9) = .0931451054638672d0
      w(10)= .0627901847324523d0
      w(11)= .0278342835580868d0
    CASE (12)
      x(1) = .00921968287664038d0
      x(2) = .0479413718147626d0
      x(3) = .115048662902848d0
      x(4) = .206341022856691d0
      x(5) = .31608425050091d0
      x(6) = .437383295744266d0
      x(7) = .562616704255734d0
      x(8) = .68391574949909d0
      x(9) = .793658977143309d0
      x(10)= .884951337097152d0
      x(11)= .952058628185237d0
      x(12)= .99078031712336d0
      w(1) = .0235876681932559d0
      w(2) = .0534696629976592d0
      w(3) = .0800391642716731d0
      w(4) = .101583713361533d0
      w(5) = .116746268269177d0
      w(6) = .124573522906701d0
      w(7) = .124573522906701d0
      w(8) = .116746268269177d0
      w(9) = .101583713361533d0
      w(10)= .0800391642716731d0
      w(11)= .0534696629976592d0
      w(12)= .0235876681932559d0
    CASE (13)
      x(1) = .00790847264070593d0
      x(2) = .041200800388511d0
      x(3) = .099210954633345d0
      x(4) = .17882533027983d0
      x(5) = .275753624481777d0
      x(6) = .384770842022433d0
      x(7) = .5d0
      x(8) = .615229157977567d0
      x(9) = .724246375518223d0
      x(10)= .82117466972017d0
      x(11)= .900789045366655d0
      x(12)= .958799199611489d0
      x(13)= .992091527359294d0
      w(1) = .0202420023826579d0
      w(2) = .0460607499188642d0
      w(3) = .0694367551098937d0
      w(4) = .0890729903809729d0
      w(5) = .103908023768444d0
      w(6) = .113141590131449d0
      w(7) = .116275776615437d0
      w(8) = .113141590131449d0
      w(9) = .103908023768444d0
      w(10)= .0890729903809729d0
      w(11)= .0694367551098937d0
      w(12)= .0460607499188642d0
      w(13)= .0202420023826579d0
    CASE (14)
      x(1) = .00685809565159383d0
      x(2) = .0357825581682132d0
      x(3) = .0863993424651175d0
      x(4) = .156353547594157d0
      x(5) = .242375681820923d0
      x(6) = .340443815536055d0
      x(7) = .445972525646328d0
      x(8) = .554027474353672d0
      x(9) = .659556184463945d0
      x(10)= .757624318179077d0
      x(11)= .843646452405843d0
      x(12)= .913600657534882d0
      x(13)= .964217441831787d0
      x(14)= .993141904348406d0
      w(1) = .0175597301658759d0
      w(2) = .0400790435798801d0
      w(3) = .0607592853439516d0
      w(4) = .0786015835790968d0
      w(5) = .092769198738969d0
      w(6) = .102599231860648d0
      w(7) = .107631926731579d0
      w(8) = .107631926731579d0
      w(9) = .102599231860648d0
      w(10)= .092769198738969d0
      w(11)= .0786015835790968d0
      w(12)= .0607592853439516d0
      w(13)= .0400790435798801d0
      w(14)= .0175597301658759d0
    CASE (15)
      x(1) = .00600374098975728d0
      x(2) = .031363303799647d0
      x(3) = .0758967082947864d0
      x(4) = .137791134319915d0
      x(5) = .214513913695731d0
      x(6) = .302924326461218d0
      x(7) = .399402953001283d0
      x(8) = .5d0
      x(9) = .600597046998717d0
      x(10)= .697075673538782d0
      x(11)= .785486086304269d0
      x(12)= .862208865680085d0
      x(13)= .924103291705214d0
      x(14)= .968636696200353d0
      x(15)= .993996259010243d0
      w(1) = .0153766209980586d0
      w(2) = .0351830237440541d0
      w(3) = .053579610233586d0
      w(4) = .0697853389630772d0
      w(5) = .083134602908497d0
      w(6) = .0930805000077812d0
      w(7) = .0992157426635559d0
      w(8) = .101289120962781d0
      w(9) = .0992157426635559d0
      w(10)= .0930805000077812d0
      w(11)= .083134602908497d0
      w(12)= .0697853389630772d0
      w(13)= .053579610233586d0
      w(14)= .0351830237440541d0
      w(15)= .0153766209980586d0
   END SELECT
  END SUBROUTINE gauss
!====================================================================
    SUBROUTINE facsb
!====================================================================
!
!   Factorizes bs matrix which is a transpose of overlap matrix sb,
!   <B_i,B_j>,  with the correct boundary condition at r=0 
!
!   SUBROUTINES called:  dpbfa (from LINPACK)
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin

    IMPLICIT NONE
    INTEGER :: m, ierr

    ! .. copy the array, converting to row oriented band storage mode

    bs = TRANSPOSE(sb)

    ! .. apply boundary condition at r=0

    do m = 1,ks-1
      bs(m,ks-m+1)=0.d0
    end do
    bs(ks,1) = 1.d0

!    CALL dpbfa(bs,ks,ns,ks-1,ierr)       
!    if (ierr /= 0 ) STOP 'facsb: dpbfa (LINPACK) failed'

    Call DPBTRF('U',ns,ks-1,bs,ks,ierr)
    if (ierr.ne.0)  Stop 'facsb: dpbtrf (LAPACK) failed'

  END SUBROUTINE facsb

!=======================================================================
   SUBROUTINE vbsplvd(t, kg, ni, x, nderiv, dbiatx)
!=======================================================================
!
!  This routine calculates the values of the B-splines and their deriva-
!  tives, of order up to nderiv, that do not vanish at x(i), i=1,..ni
!  There are ks such B-splines at each point.
!
!  This routine is a vector version of bsplvd written by C. de Boor,
!  ``A Practical Guide to Splines".
!
!  subroutine contained: vbsplvb
!
!  calling sequence:
!       vbsplvd
!          ||
!       vbsplvb
!
!-----------------------------------------------------------------------
!  on entry
!  --------
!  t     the knot array, of length nt >=nv+2ks-1.  It is assumed
!        that t(i) < t(i+1) for each interval containing an x(i)
!        Division by zero will result otherwise (in vbsplvb).
!
!  kg    gives the beginning interval from which the B-splines are
!        evaluated at the Gaussian points.
!
!  ni    the number of intervals in which B-splines are to be evaluated
!        at all Gaussian points, its uplimit is nt.
!
!  x     the point array at which these values are sought,
!        one per interval, of length ni.
!
!  nderiv   an integer indicating that values of B-splines and their
!        derivatives up to but not including the  nderiv-th  are asked
!        for.
!
!  working area
!  ------------
!  w31   a three dimensional array, w31(i,j,m) (j=1,..,ks m=1,..,ks) con-
!        tains B-coeff.s of the derivatives of a certain order of the
!        ks B-splines of interest at point x(i)
!
!  w1,w2      one dimensional arrays
!
!  on return
!  ---------
!  dbiatx     a three dimensional array. its entry (i,j,m) contains
!        value of  (m-1)st  derivative of  (l-ks+j)-th B-spline of
!        order ks at point x(i) for knot sequence t, i=1..ni,
!        j=m..ks; m=1..nderiv;and l=kg..kg+ni-1
!
!  method
!  ------
!  values at x(i) of all the relevant B-splines of order ks,ks-1,...,
!  ks+1-nderiv  are generated via vbsplvb and stored temporarily
!  in dbiatx. then, the B-coeffs of the required derivatives of the
!  B-splines of interest are generated by differencing, each from the
!  preceding one of lower order, and combined with the values of B-
!  splines of corresponding order in dbiatx to produce the desired
!  values.
!----------------------------------------------------------------------

    USE spline_param

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kg, ni, nderiv
    REAL(KIND=8), DIMENSION(ns+ks), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(nv,ks,ks), INTENT(INOUT):: dbiatx
    REAL(KIND=8), DIMENSION(ni), INTENT(IN):: x

    ! local variables
    REAL(KIND=8):: fkpimm
    REAL(KIND=8), DIMENSION(ni):: w1,w2
    REAL(KIND=8), DIMENSION(ni,ks,ks):: w31
    REAL(KIND=8), DIMENSION(ni,ks) :: deltar, deltal
    INTEGER:: i, j, n, m, mhigh, kp1, jhigh, ideriv
    INTEGER:: ldummy, kpimm, jlow, jpimid, il

    mhigh = MAX(MIN(nderiv,ks),1)   !mhigh is usually equal to nderiv.
    kp1 = ks+1
    jhigh = kp1-mhigh
    CALL vbsplvb(kg,ni,x,jhigh,1,dbiatx)
    IF(mhigh == 1) RETURN

    ! ..the first row of dbiatx always contains the B-spline values
    ! ..for the current order. these are stored in row ks+1-current
    ! ..order before vbsplvb is called to put values for the next
    ! ..higher order on top of it. Vbsplvb only uses the first two dimensions

    ideriv = mhigh
    DO m = 2, mhigh
      jpimid = 1
      DO j = ideriv,ks
       dbiatx(1:ni,j,ideriv) = dbiatx(1:ni,jpimid,1)
       jpimid = jpimid+1
      END DO
      ideriv = ideriv-1
      jhigh = kp1-ideriv	
      CALL vbsplvb(kg,ni,x,jhigh,2,dbiatx)
    END DO

    ! at this point,  b(.,n-ks+i, ks+1-j)(x) is in dbiatx(.,i,j) for
    ! n=kg..kg+ni-1,i=j..ks,j=1..mhigh('='nderiv).in particular,the
    ! first row of  dbiatx  is already in final form. to obtain cor-
    ! responding derivatives of B-splines in subsequent rows, gene-
    ! rate their B-repr. by differencing, then evaluate at x(.).

    jlow = 1
    DO i = 1,ks
      w31(1:ni,jlow:ks,i) = 0.d0	
      jlow = i
      w31(1:ni,i,i) = 1.d0
    END DO

    ! at this point, w31(.,.,j) contains the B-coeffs for the j-th of the
    ! ks B-splines of interest here.

    DO m = 2,mhigh
      kpimm = kp1-m
      fkpimm = kpimm
      i = ks
      il = 0

      ! for j=1,...,ks, construct B-coeffs of  (m-1)st  derivative of
      ! B-splines from those for preceding derivative by differencing
      ! and store again in  w31(.,.,j). the fact that w31(i,j) = 0  for
      ! i < j is used.

      DO ldummy = 1, kpimm
       DO n = kg-il,ni+kg-il-1
        w1(n-kg+il+1) = fkpimm/(t(n+kpimm)-t(n))
       END DO

        ! the assumption that t(n) < t(n+1) makes denominator
        ! in w1(1..ni) nonzero.

        DO j = 1,i
         w31(1:ni,i,j) = (w31(1:ni,i,j)-w31(1:ni,i-1,j))*w1(1:ni)
        END DO
        il = il+1
        i = i-1
      END DO

      ! for i=1,...,ks, combine B-coeffs a(.,.,i) with B-spline values
      ! stored in dbiatx(.,.,m) to get value of (m-1)st  derivative of
      ! i-th B-spline (of interest here) at x(.), and store in
      ! dbiatx(.,i,m). storage of this value over the value of a B-spline
      ! of order m there is safe since the remaining B-spline derivat-
      ! ive of the same order do not use this value due to the fact
      ! that  a(.,j,i) = 0  for j .lt. i .

      DO i = 1,ks
        w2(1:ni) = 0.d0
        jlow = MAX(i,m)
        DO j = jlow,ks
          w2(1:ni) = w2(1:ni) + w31(1:ni,j,i)*dbiatx(1:ni,j,m)
        END DO
        dbiatx(1:ni,i,m) = w2(1:ni)
      END DO

    END DO

    CONTAINS


    !===================================================================
      SUBROUTINE vbsplvb(kg, ni, x, jhigh, index, biatx)
    !=====================================================================
    !  This routine calculates the values of all possibly nonzero B-splines
    !  at x(i) (i=1,..ni) of order
    !               jout=max(jhigh,(j+1)*(index-1))
    !  with knot sequence  t .
    !
    !  This routine is a vector version of bsplvb written by C. de Boor,
    !  "A Practical Guide to Splines", Chapter X, page 135
    !
    !  on entry
    !  --------
    !  t    -  knot sequence, of length nt=ns+ks, assumed to be nonde-
    !          creasing, that is t(i) <= t(i+1)
    !
    !  jhigh-  choose jhigh=1 to get the B-spline values directly
    !            by calling vbsplvb.
    !
    !  kg   -  gives the beginning interval from which the B-splines
    !           are to be evaluated at Gaussin points.
    !
    !  ni   -  the number of intervals in which B-splines are
    !            evaluated at all Gaussian points, its uplimit is nv.
    !
    !  x    -  the points at which the B-splines are to be evaluated,
    !            its length is ni,
    !
    !  index-  integers which determine the order  jout = max(jhigh,
    !            (j+1)*(index-1))  of the B-splines whose values at x(i)
    !            are to be returned.  index is used to avoid recalcula-
    !            tions when several columns of the triangular array of
    !            B-spline values are needed (e.g., in vbsplvd ).
    !            More precisely,
    !                     if index = 1 ,
    !            the calculation starts from scratch and the entire
    !            triangular array of B-spline values of orders
    !            1,2,...,jhigh is generated, order by order ,
    !            i.e., column by column .
    !                     if  index = 2 ,
    !            only the B-spline values of order  j+1, j+2, ..., jout
    !            are generated, the assumption being that  biatx,j,
    !            deltal,deltar are, on  entry, as they were on exit at the
    !            previous call. In particular, if  jhigh = 0, then
    !            jout = j+1, i.e., just the next column of B-spline
    !            values is generated.
    !
    !  working area
    !  ------------
    !  deltal, deltar: two dimensional arrays
    !  term, saved:    one dimensional arrays.
    !
    !  on return
    !  ---------
    !  biatx.....two dimensional array, with biatx(j-k+1,i)(j=k..ni)
    !        containing the value at x(j-k+1) of the polynomial of order
    !        jout which agrees with the B-spline b(j-jout+1,jout,t) on
    !        the interval (t(j),t(j+1)).
    !
    !  method
    !  ------
    !  The recurrence relation
    !
    !                       x - t(i)              t(i+j+1) - x
    !     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
    !                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
    !
    !  is used (repeatedly) to generate the (j+1)-vector  b(l-j,j+1)(x),
    !  ...,b(l,j+1)(x)  from the j-vector  b(l-j+1,j)(x),...,
    !  b(l,j)(x), storing the new values in  biatx  over the old. the
    !  facts that
    !            b(i,1) = 1  if  t(i) <= x < t(i+1)
    !  and that
    !            b(i,j)(x) = 0  unless  t(i) <= x < t(i+j)
    !  are used. the particular organization of the calculations follows al-
    !  gorithm  (8)  in chapter x of the text.
    !-----------------------------------------------------------------------

        USE spline_param
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(nv,ks), INTENT(INOUT):: biatx
        REAL(KIND=8), DIMENSION(nv), INTENT(IN):: x
        INTEGER, INTENT(IN):: kg, ni, jhigh, index

        ! .. Local variables
        INTEGER:: i, jp1, m
        INTEGER, SAVE:: j=1
        REAL(KIND=8), DIMENSION(ni) :: term, saved

        IF(index == 1) THEN
          j=1
          biatx(1:ni,1)=1.d0
          IF (j >= jhigh)  RETURN
        END IF

        DO
          jp1=j+1
          saved(1:ni)=0.d0

          DO i=1,ni
            deltar(i,j)=t(i+kg-1+j)-x(i)
            deltal(i,j)=x(i)-t(i+kg-j)
          END DO 	

          DO m=1,j
            DO i=1,ni
              term(i)=biatx(i,m)/(deltar(i,m)+deltal(i,jp1-m))
              biatx(i,m)=saved(i)+deltar(i,m)*term(i)
              saved(i)=deltal(i,jp1-m)*term(i)
            END DO
          END DO

          biatx(1:ni,jp1)=saved(1:ni)
          j=jp1
          IF (j >= jhigh) EXIT
        END DO
      END SUBROUTINE vbsplvb

  END SUBROUTINE vbsplvd!=======================================================================
    Real(8) FUNCTION azl(z,h,ks,lp)
!=======================================================================
!
!   Value of B_(lp+1)/r^lp at r = 0 where lp = l+1
!----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ks,lp
    REAL(8), INTENT(in) :: z,h

    INTEGER(4) :: j
    REAL(8) :: c

    IF (lp < ks ) THEN
      azl = 1.d0
      c = z/h
      do j = 1,lp
        azl = (azl*c*(ks-j))/(j*j)
      end do
    ELSE
      azl = 0.d0
    END IF

   END FUNCTION azl
