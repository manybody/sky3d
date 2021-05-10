MODULE Fourier
  USE params, ONLY: db,wflag
  USE Grids, ONLY: nx,ny,nz
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(C_LONG),SAVE :: pforward,pbackward,xforward,xbackward, &
       yforward,ybackward,zforward,zbackward
CONTAINS
  SUBROUTINE init_fft
#ifdef CUDA
    USE cufft_m
    ! Create plans with cufft
    ! Use cufftExecZ2Z(plan,in,out,direction) to execute double-complex FFTs
    ! cufftPlan3d(plan,x,y,z,type)
    CALL cufftPlan3d(pforward,nx,ny,nz,CUFFT_Z2Z)
    CALL cufftPlan3d(pbackward,nx,ny,nz,CUFFT_Z2Z)

    ! cufftPlanMany(plan,rank,size,
    !    inembed,istride,idist,onembed,ostride,odist,type,batch
    CALL cufftPlanMany(xforward,1,(/nx/), &
         0,1,nx,0,1,nx,CUFFT_Z2Z,2*ny*nz)
    CALL cufftPlanMany(xbackward,1,(/nx/), &
         0,1,nx,0,1,nx,CUFFT_Z2Z,2*ny*nz)
    CALL cufftPlanMany(yforward,1,(/ny/), &
         0,nx,1,0,nx,1,CUFFT_Z2Z,nx)
    CALL cufftPlanMany(ybackward,1,(/ny/), &
         0,nx,1,0,nx,1,CUFFT_Z2Z,nx)
    CALL cufftPlanMany(zforward,1,(/nz/), &
         0,nx*ny,1,0,nx*ny,1,CUFFT_Z2Z,nx*ny)
    CALL cufftPlanMany(zbackward,1,(/nz/), &
         0,nx*ny,1,0,nx*ny,1,CUFFT_Z2Z,nx*ny)
#else
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
    DEALLOCATE(p)
#endif
    IF(wflag) WRITE(*,*) '***** FFT plans established *****'
  END SUBROUTINE init_fft
END MODULE Fourier
