!=======================================================================
! Checks SET_TEMPERING_WEIGHTS against thesis Table 2.1 (Eq. 2.12).
! Built and run by tests/test_weights.sh.
!=======================================================================
PROGRAM test_weights
  USE common
  USE common_namelist
  IMPLICIT NONE
  INTEGER :: nfail = 0

  !Thesis Table 2.1 (6 decimals)
  CALL check(2, 0.0d0, (/0.5d0, 0.5d0/), 5.0d-7)
  CALL check(2, 1.0d0, (/0.182426d0, 0.817574d0/), 5.0d-7)
  CALL check(2, 2.0d0, (/0.047426d0, 0.952574d0/), 5.0d-7)
  CALL check(2, 3.0d0, (/0.010987d0, 0.989013d0/), 5.0d-7)
  CALL check(3, 0.0d0, (/1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0/), 5.0d-7)
  CALL check(3, 1.0d0, (/0.043896d0, 0.324352d0, 0.631752d0/), 5.0d-7)
  CALL check(3, 2.0d0, (/0.003806d0, 0.207815d0, 0.788379d0/), 5.0d-7)
  CALL check(3, 3.0d0, (/0.000295d0, 0.119168d0, 0.880537d0/), 5.0d-7)
  !alpha_s = 1, Nt = 4 and 5 (6 decimals, computed independently in Python)
  CALL check(4, 1.0d0, (/0.011942d0, 0.145488d0, 0.334765d0, 0.507804d0/), 5.0d-7)
  CALL check(5, 1.0d0, (/0.003482d0, 0.069933d0, 0.190098d0, 0.313418d0, 0.423070d0/), 5.0d-7)

  !ntemp = 1 must give exactly 1, whatever alpha_s is
  ntemp = 1 ; alpha_s = 2.5d0 ; alpha_list = -1.0d0
  CALL SET_TEMPERING_WEIGHTS()
  IF( temp_alpha(1) /= 1.0d0 )THEN
    WRITE(*,*) 'FAIL ntemp=1: temp_alpha(1) = ', temp_alpha(1)
    nfail = nfail + 1
  ENDIF

  !Explicit list, normalized
  ntemp = 2 ; alpha_list = -1.0d0 ; alpha_list(1:2) = (/0.3d0, 0.7d0/)
  CALL SET_TEMPERING_WEIGHTS()
  IF( ABS(temp_alpha(1)-0.3d0) > 1.0d-12 .OR. ABS(temp_alpha(2)-0.7d0) > 1.0d-12 )THEN
    WRITE(*,*) 'FAIL alpha_list: ', temp_alpha(1:2)
    nfail = nfail + 1
  ENDIF

  IF( nfail == 0 )THEN
    WRITE(*,*) 'test_weights: all checks passed'
  ELSE
    WRITE(*,*) 'test_weights: ', nfail, ' checks failed'
    STOP 1
  ENDIF

CONTAINS

  SUBROUTINE check(nt, as, ref, tol)
    INTEGER, INTENT(IN)      :: nt
    REAL(r_size), INTENT(IN) :: as, ref(:), tol
    ntemp = nt ; alpha_s = as ; alpha_list = -1.0d0
    CALL SET_TEMPERING_WEIGHTS()
    IF( MAXVAL( ABS( temp_alpha(1:nt) - ref ) ) > tol .OR. &
        ABS( SUM( temp_alpha(1:nt) ) - 1.0d0 ) > 1.0d-12 )THEN
      WRITE(*,'(A,I2,A,F4.1,A,6F10.6)') 'FAIL Nt=', nt, ' alpha_s=', as, ': ', temp_alpha(1:nt)
      nfail = nfail + 1
    ELSE
      WRITE(*,'(A,I2,A,F4.1,A,6F10.6)') 'ok   Nt=', nt, ' alpha_s=', as, ': ', temp_alpha(1:nt)
    ENDIF
  END SUBROUTINE check

END PROGRAM test_weights
