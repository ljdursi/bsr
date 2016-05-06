!=======================================================================
      Subroutine ZRMAT(E,RA,nch,nhm,RMAT,VALUE,WMAT)
!=======================================================================
!     CALCULATE THE R-MATRIX FOR GIVEN ENERGY "E" BASED ON THE
!     SURFACE AMPLITUDES AND EIGENVALUES 
!-----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in)  :: E,RA
      Integer, intent(in)  :: nch,nhm
      Real(8), intent(in)  :: VALUE(nhm)
      Real(8), intent(in)  :: WMAT(nch,nhm)
      Real(8), intent(out) :: RMAT(nch,nch)
      Integer :: I,J,K
      Real(8) :: A,B,C

! ... CALCULATE R-MATRIX 

      RMAT = 0.d0
      A = 1.d0/RA
      DO K = 1,nhm
        B = A / (VALUE(K)-E)
        DO I = 1,nch
          C = WMAT(I,K)*B
          DO J = 1,nch
           RMAT(I,J) = RMAT(I,J) + WMAT(J,K)*C
          END DO
        END DO
      END DO

      End Subroutine ZRMAT
