PROGRAM spectral_analysis
!
! This program computes the spectral distribution of quadrupole 
! oscillations produced as output on file 'quadrupoles.res'. 
! Isoscalar and isovector quadrupole moments are composed from the
! information on x**2, y**2, and z**2 in 'quadrupoles.res'.
! A limiting time profile of type COS**n is multiplied to the signal
! to guarantee that the signal approaches zero at the end of the
! time interval (called time filtering, or windowing).
! It is assumed that the quadrupole motion was initialized by an
! instantaneous isoscalar quadrupole boost. The spectral strength is 
! then simply the imaginary part of the Fourier transform of the 
! isoscalar signal. The code also produces the Fourier transform
! of the isovector signal. This has to be taken with care because
! the excitation operator and the analzing operator are not the same
! (what they should be for spectral analysis).
!
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  REAL(db),PARAMETER :: PI=3.14159265358979d0
  REAL(db),PARAMETER :: hbc=197.32D0        ! to convert 1/fm to MeV
  INTEGER,PARAMETER :: nfilter=4       ! order of damping COS**nfilter
!
  REAL(db), ALLOCATABLE :: qext(:)
  REAL(db), ALLOCATABLE :: qfilter(:)
  REAL(db) :: time,timax,timin,deltime,timold,filter
  REAL(db) :: delomega,omega
  COMPLEX(db) :: cfac,cprod,cssum,cvsum,ctsum
  INTEGER :: ndim
  INTEGER :: nt,iomega


!
! A first quick reading to check the length and integrity of 
! the data file.
!
  OPEN(1,file='extfield.res')
  READ(1,'(1x)')
  DO nt=1,999999
    READ(1,*,ERR=19,END=19) timax
    IF(nt==1) timin = timax
    IF(nt==2) deltime=timax-timold
    IF(nt>2 .AND. ABS(timax-timold-deltime)>1D-10) THEN
      WRITE(*,*) 'problem at line ',nt+1
      STOP "time steps not equidistant in 'extfield.res'"
    END IF
    ndim = nt
    timold = timax
  END DO
19 CONTINUE
  REWIND(1)
  WRITE(*,'(a,i10,f12.2)') 'nr. of time steps, max.time=',ndim,timax
  ALLOCATE(qext(ndim),qfilter(ndim))

!
! Read data and compose quadrupole moments (with time filter).
! The coupling to isovector moment assumes an N=Z nucleus.
!
  OPEN(2,file='extfield_filt.res')
  WRITE(2,'(a)') '# external field signal in time domain:'
  WRITE(2,'(a)') '# time[fm/c]     raw      filtered '
  READ(1,'(1x)')
  DO nt=1,ndim
    READ(1,*,ERR=29,END=29) timax,qext(nt)
    filter = COS((nt-1)*PI/(2D0*(ndim-1)))**nfilter    
    qfilter(nt) = qext(nt)*filter
    WRITE(2,'(f10.2,4(1pg13.5))') (nt-1)*deltime,qfilter(nt)
  END DO
  CLOSE(1)
  CLOSE(2)
  WRITE(*,*) "filtered signal written on 'extfield_filt.res''"

!
! Perform Fourier transformation and print.
!
  OPEN(2,file='extfield_spectrum.res')
  WRITE(2,'(a)') '# spectral distributions of external-field strenght and power:'
  WRITE(2,'(a)') '# omega[MeV]  REAL(FT signal)  strength     power'
  deltime   = (timax-timin)/(ndim-1)
  delomega  = (PI+PI)/(ndim*deltime)
  DO iomega=1,ndim/4
    omega  = iomega*delomega-delomega
    cfac   = CEXP(CMPLX(0D0,omega*deltime))
    cprod  = CMPLX(1D0,0D0)/cfac
    cssum  = CMPLX(0D0,0D0)
    DO nt=1,ndim
      cprod = cprod*cfac
      cssum = cprod*CMPLX(-qfilter(nt),0D0) + cssum
    END DO
    cssum = cssum*deltime
    WRITE(2,'(f11.5,6(1pg13.5))') omega*hbc,cssum,ABS(cssum)**2
  END DO
  CLOSE(2)
  WRITE(*,*) "spectra written on 'extfield_spectrum.res''"

STOP "regular end"

29 STOP "error in reading data"

END PROGRAM spectral_analysis
