MODULE User
  USE Params
  USE Grids
  USE Levels
  IMPLICIT NONE
CONTAINS
  SUBROUTINE init_user
    STOP 'No user-supplied initialization'
  END SUBROUTINE init_user
END MODULE User
