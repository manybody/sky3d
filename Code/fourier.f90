MODULE Fourier
  USE params, ONLY: db,wflag
  USE Grids, ONLY: nx,ny,nz
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(C_LONG),SAVE :: pforward,pbackward,xforward,xbackward, &
       yforward,ybackward,zforward,zbackward
CONTAINS
  SUBROUTINE init_fft
    INCLUDE 'fftw3.f'
    COMPLEX(db),ALLOCATABLE :: p(:,:,:,:,:)
    INTEGER,SAVE :: FFTW_planflag             
! set option for FFTW setup here
!    FFTW_planflag=FFTW_ESTIMATE
!    FFTW_planflag=FFTW_MEASURE
    FFTW_planflag=FFTW_PATIENT
!    FFTW_planflag=FFTW_EXHAUSTIVE
!    FFTW_planflag=FFTW_MEASURE+FFTW_UNALIGNED
    ALLOCATE(p(nx,ny,nz,2,2))
    CALL dfftw_plan_dft_3d(pforward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_FORWARD, FFTW_planflag)
    CALL dfftw_plan_dft_3d(pbackward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_BACKWARD, FFTW_planflag)
    CALL dfftw_plan_many_dft(xforward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,2),0,1,nx, &
         FFTW_FORWARD, FFTW_planflag)
    CALL dfftw_plan_many_dft(xbackward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,1),0,1,nx, &
         FFTW_BACKWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(yforward,1,(/ny/),nx, &
       p(:,:,:,:,1),0,nx,1,p(:,:,:,:,2),0,nx,1, &
       FFTW_FORWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(ybackward,1,(/ny/),nx, &
       p(:,:,:,:,1),0,nx,1,p(:,:,:,:,1),0,nx,1, &
       FFTW_BACKWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(zforward,1,(/nz/),nx*ny, &
       p(:,:,:,:,1),0,nx*ny,1,p(:,:,:,:,2),0,nx*ny,1, &
       FFTW_FORWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(zbackward,1,(/nz/),nx*ny, &
       p(:,:,:,:,1),0,nx*ny,1,p(:,:,:,:,1),0,nx*ny,1, &
       FFTW_BACKWARD, FFTW_planflag)
    IF(wflag) WRITE(*,*) '***** FFTW3 plans established *****'
    DEALLOCATE(p)
  END SUBROUTINE init_fft
END MODULE Fourier
