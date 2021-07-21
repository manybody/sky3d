MODULE Params
  IMPLICIT NONE
  !**********************************************************************
  !     data type definition                                            *
  !**********************************************************************
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  !**********************************************************************
  !     useful constants                                                *
  !**********************************************************************
  REAL(db),PARAMETER :: pi=3.14159265358979D0
  REAL(db),PARAMETER :: hbc=197.32164D0
  REAL(db),PARAMETER :: e2=1.43989D0
  !**********************************************************************
  !  names of files and units to be used                                *
  !**********************************************************************
  CHARACTER(LEN=80) :: wffile='none'
  CHARACTER(LEN=80) :: converfile='conver.res'
  CHARACTER(LEN=80) :: monopolesfile='monopoles.res'
  CHARACTER(LEN=80) :: dipolesfile='dipoles.res'
  CHARACTER(LEN=80) :: momentafile='momenta.res'
  CHARACTER(LEN=80) :: energiesfile='energies.res'
  CHARACTER(LEN=80) :: diffenergiesfile='diffenergies.res'
  CHARACTER(LEN=80) :: quadrupolesfile='quadrupoles.res'
  CHARACTER(LEN=80) :: spinfile='spin.res'
  CHARACTER(LEN=80) :: extfieldfile='extfield.res'
  INTEGER,PARAMETER :: scratch=11, scratch2=12
  !**********************************************************************
  !     basic parameters controlling the job                            *
  !**********************************************************************
  LOGICAL :: tcoul=.TRUE.
  LOGICAL :: tstatic
  LOGICAL :: tdynamic
  LOGICAL :: tfft=.TRUE.
  LOGICAL :: trestart=.FALSE.
  LOGICAL,PARAMETER ::  toddls=.TRUE.    
  LOGICAL,PARAMETER ::  todd=.TRUE.   
  !**********************************************************************
  ! parameters controlling printout frequency etc.
  !**********************************************************************
  INTEGER :: mprint=100
  INTEGER :: mplot=0
  INTEGER :: mrest=0
  INTEGER :: iter
  REAL(db) :: time
  LOGICAL :: wflag,printnow
  INTEGER,PARAMETER :: nselect=10
  CHARACTER(LEN=nselect) :: writeselect='r'
  LOGICAL :: write_isospin=.FALSE.
  INTEGER,PARAMETER ::  mnof=4   ! maximum number of fragments
  INTEGER :: nof                 ! real number of fragments
  REAL(db) :: r0=1.2D0
  INTEGER :: mlocalize=0
  LOGICAL :: tlocalize=.FALSE.
END MODULE Params
