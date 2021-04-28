MODULE Fourier
  USE params, ONLY: db,wflag
  USE Grids, ONLY: nx,ny,nz
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(C_LONG),SAVE :: pforward,pbackward,xforward,xbackward, &
       yforward,ybackward,zforward,zbackward
CONTAINS
  SUBROUTINE init_fft
    INCLUDE 'cufftw.h'
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
    IF(wflag) WRITE(*,*) '***** FFT plans established *****'
    DEALLOCATE(p)
  END SUBROUTINE init_fft
END MODULE Fourier
