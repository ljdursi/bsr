!======================================================================
      Function BHL(i,j)
!======================================================================
!     <P_i| L |P_j>  with inclusion of rel.shift if rel = .true.
!----------------------------------------------------------------------
      USE spline_param
      USE spline_orbitals
      USE spline_hl
      USE spline_atomic

      IMPLICIT NONE
      INTEGER, INTENT(in) :: i,j
      REAL(KIND=8) :: BHL
      REAL(KIND=8), EXTERNAL :: BVMV

      if (iabs(LBS(i)-LBS(j)) .NE. 0) Stop ' HL:  LI <> LJ'

      Call HLM (lbs(i))

      BHL = BVMV(ns,ks,hl,'a',pbs(1,i),pbs(1,j))

      End Function BHL

