!=======================================================================
  Subroutine Det_d5 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell d5 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)
 
  Integer, parameter :: iq_d5 =   5
  Integer, parameter :: kd_d5 = 252
 
  Integer :: Idet_d5 (iq_d5,kd_d5)
 
  Integer :: ML_d5 (kd_d5)
  Integer :: MS_d5 (kd_d5)
 
  if(id.le.0.or.id.gt.kd_d5) Stop "Det_d5: index id is out of range"
 
  ML = ML_d5 (id)
  MS = MS_d5 (id)
 
  Idet (1:iq_d5)= Idet_d5 (:,id)
 

  Data Idet_d5 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6  /

  Data Idet_d5 ( 2,:)/ &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7, 7  /

  Data Idet_d5 ( 3,:)/ &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, &
   5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, &
   6, 6, 6, 7, 7, 7, 8, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 5, 5, 5, 5, 5, 5, 5, &
   5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, &
   6, 6, 7, 7, 7, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8  /

  Data Idet_d5 ( 4,:)/ &
   4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 6, 6, 6, 6, 7, 7, &
   7, 8, 8, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 7, 7, 7, &
   8, 8, 9, 8, 8, 9, 9, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, &
   5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 6, 6, 6, 6, 7, 7, 7, &
   8, 8, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 7, 7, 7, 8, &
   8, 9, 8, 8, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 8, 8, 9, 9, 9, 9  /

  Data Idet_d5 ( 5,:)/ &
   5, 6, 7, 8, 9,10, 6, 7, 8, 9,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 6, 7, 8, 9,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 7, 8, 9,10, 8, 9, &
  10, 9,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 6, 7, 8, 9,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 8, 9,10, &
   9,10,10, 9,10,10,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, &
   6, 7, 8, 9,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 7, 8, 9,10, 8, 9,10, &
   9,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 7, 8, 9,10, 8, 9,10, 9,10,10, 8, 9,10, 9, &
  10,10, 9,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 9,10,10,10,10,10  /

  Data ML_d5 / &
 -1, -1, -7, -7,  1,  1,  3, -3, -3,  5,  5, -3, -3,  5,  5, -9, -1, -1, -1, -1,  7,  3, -3, -3,  5,  5, -3, -3,  5,  5, -9, -1, &
 -1, -1, -1,  7,  1,  1,  9,  9, -5,  3,  3,  3,  3, 11, -5,  3,  3,  3,  3, 11, -3, -3,  5,  5,  1, -5, -5,  3,  3, -5, -5,  3, &
  3,-11, -3, -3, -3, -3,  5, -1, -1,  7,  7, -7,  1,  1,  1,  1,  9, -7,  1,  1,  1,  1,  9, -5, -5,  3,  3, -1, -1,  7,  7, -7, &
  1,  1,  1,  1,  9, -7,  1,  1,  1,  1,  9, -5, -5,  3,  3, -3,  5,  5,  5,  5, 13, -1, -1,  7,  7, -1, -1,  7,  7,  1,  1, -5, &
 -5,  3,  3, -5, -5,  3,  3,-11, -3, -3, -3, -3,  5, -1, -1,  7,  7, -7,  1,  1,  1,  1,  9, -7,  1,  1,  1,  1,  9, -5, -5,  3, &
  3, -1, -1,  7,  7, -7,  1,  1,  1,  1,  9, -7,  1,  1,  1,  1,  9, -5, -5,  3,  3, -3,  5,  5,  5,  5, 13, -1, -1,  7,  7, -1, &
 -1,  7,  7,  1, -3, -3,  5,  5, -9, -1, -1, -1, -1,  7, -9, -1, -1, -1, -1,  7, -7, -7,  1,  1, -5,  3,  3,  3,  3, 11, -3, -3, &
  5,  5, -3, -3,  5,  5, -1, -5,  3,  3,  3,  3, 11, -3, -3,  5,  5, -3, -3,  5,  5, -1,  1,  1,  9,  9,  3,  3  /

  Data MS_d5 / &
  0,  2,  0,  2,  0,  2,  0, -2,  0, -2,  0,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  0,  2,  2,  4,  2,  4,  2,  0, &
  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0, -2,  0,  0,  2,  0, &
  2,  0, -2,  0,  0,  2,  0, -2,  0, -2,  0, -2, -4, -2, -2,  0, -2,  0, -2,  0,  0,  2,  0, -2,  0, -2,  0,  0,  2,  0,  2,  0, &
 -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0, -2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  0, &
  2,  0,  2,  2,  4,  2,  4,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0, &
  2,  2,  4,  2,  4,  2,  0,  2,  2,  4,  2,  4,  2,  4,  4,  6,  4,  2,  4,  2,  4,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  2, &
  4,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0, -2,  0, &
 -2,  0,  0,  2,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  2,  4,  2,  4,  2,  0,  2,  0,  2,  0,  2  /
 
  End Subroutine Det_d5 
