!=======================================================================
! Reads ./letkf.namelist and prints the tempering settings.
! Built and run by tests/test_weights.sh.
!=======================================================================
PROGRAM test_namelist_read
  USE common
  USE common_namelist
  IMPLICIT NONE
  CALL READ_NAMELIST()
  WRITE(*,'(A,I3,A,6F10.6)') 'ntemp=', ntemp, ' temp_alpha=', temp_alpha(1:ntemp)
END PROGRAM test_namelist_read
