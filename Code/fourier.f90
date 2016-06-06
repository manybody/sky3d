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
    ALLOCATE(p(nx,ny,nz,2,2))
    CALL dfftw_plan_dft_3d(pforward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_FORWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_dft_3d(pbackward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_BACKWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_many_dft(xforward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,2),0,1,nx, &
         FFTW_FORWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_many_dft(xbackward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,1),0,1,nx, &
         FFTW_BACKWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_guru_dft(yforward,1,ny,nx,nx, &
         3,(/nx,nz,2/),(/1,nx*ny,nx*ny*nz/),(/1,nx*ny,nx*ny*nz/), &
         p(:,:,:,:,1),p(:,:,:,:,2), &
         FFTW_FORWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_guru_dft(ybackward,1,ny,nx,nx, &
         3,(/nx,nz,2/),(/1,nx*ny,nx*ny*nz/),(/1,nx*ny,nx*ny*nz/), &
         p(:,:,:,:,1),p(:,:,:,:,1), &
         FFTW_BACKWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_guru_dft(zforward,1,nz,nx*ny,nx*ny, &
         3,(/nx,ny,2/),(/1,nx,nx*ny*nz/),(/1,nx,nx*ny*nz/), &
         p(:,:,:,:,1),p(:,:,:,:,2), &
         FFTW_FORWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    CALL dfftw_plan_guru_dft(zbackward,1,nz,nx*ny,nx*ny, &
         3,(/nx,ny,2/),(/1,nx,nx*ny*nz/),(/1,nx,nx*ny*nz/), &
         p(:,:,:,:,1),p(:,:,:,:,1), &
         FFTW_BACKWARD, FFTW_MEASURE+FFTW_UNALIGNED)
    IF(wflag) WRITE(*,*) '***** FFTW3 plans established *****'
    DEALLOCATE(p)
  END SUBROUTINE init_fft
END MODULE Fourier
