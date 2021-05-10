MODULE Coulomb
  USE Params, ONLY: db,pi,e2
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz,wxyz,periodic
  USE Densities, ONLY: rho
#ifdef CUDA
  USE cufft_m
#endif
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER,PRIVATE :: nx2,ny2,nz2
  INTEGER(C_LONG),PRIVATE,SAVE :: coulplan1,coulplan2
  REAL(db),ALLOCATABLE,SAVE :: wcoul(:,:,:)
  COMPLEX(db),PRIVATE,ALLOCATABLE,SAVE :: q(:,:,:)
#ifdef CUDA
  COMPLEX(db),PRIVATE,DEVICE :: q_d(:,:,:)
#endif
  PUBLIC :: poisson,coulinit,wcoul
  PRIVATE :: initiq
CONTAINS
  !***************************************************
  SUBROUTINE poisson
    COMPLEX(db),ALLOCATABLE :: rho2(:,:,:)
#ifdef CUDA
    COMPLEX(db),DEVICE :: rho2_d(:,:,:)
#endif
    ALLOCATE(rho2(nx2,ny2,nz2))
    ! put proton density into array of same or double size, zeroing rest
    IF(.NOT.periodic) rho2=(0.D0,0.D0)
    rho2(1:nx,1:ny,1:nz)=rho(:,:,:,2)
    ! transform into momentum space
#ifdef CUDA
    ! Copy to device, perform FFT, copy back to host
    rho2_d = rho2
    CALL cufftExecZ2Z(coulplan1,rho2_d,rho2_d,CUFFT_FORWARD)
    rho2 = rho2_d
#else
    CALL dfftw_execute_dft(coulplan1,rho2,rho2)
#endif
    ! add charge factor and geometric factors
    ! note that in the periodic case q has only a real part
    IF(periodic) THEN
       rho2=4.D0*pi*e2*REAL(q)*rho2
    ELSE
       rho2=e2*wxyz*q*rho2
    END IF
    ! transform back to coordinate space and return in wcoul
#ifdef CUDA
    rho2_d = rho2
    CALL cufftExecZ2Z(coulplan2,rho2_d,rho2_d,CUFFT_INVERSE)
    rho2 = rho2_d
#else
    CALL dfftw_execute_dft(coulplan2,rho2,rho2)
#endif
    wcoul=REAL(rho2(1:nx,1:ny,1:nz))/(nx2*ny2*nz2)
    DEALLOCATE(rho2)
  END SUBROUTINE poisson
  !***************************************************
  SUBROUTINE coulinit
#ifndef CUDA
    INCLUDE 'fftw3.f'
#endif
    REAL(db),ALLOCATABLE :: iqx(:),iqy(:),iqz(:)
    INTEGER :: i,j,k
    IF(ALLOCATED(q)) RETURN ! has been initialized already
    ! dimensions will be doubled for isolated distribution
    IF(periodic) THEN
       nx2=nx; ny2=ny; nz2=nz
    ELSE
       nx2=2*nx; ny2=2*ny; nz2=2*nz
    END IF
    ! allocated helper arrays
    ALLOCATE (wcoul(nx,ny,nz),q(nx2,ny2,nz2),iqx(nx2),iqy(ny2),iqz(nz2))
    ! set up FFTW plans
#ifdef CUDA
    CALL cufftPlan3d(coulplan1,nx2,ny2,nz2,CUFFT_Z2Z)
    CALL cufftPlan3d(coulplan2,nx2,ny2,nz2,CUFFT_Z2Z)
#else
    CALL dfftw_plan_dft_3d(coulplan1,nx2,ny2,nz2,q,q, &
         FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
    CALL dfftw_plan_dft_3d(coulplan2,nx2,ny2,nz2,q,q, &
         FFTW_BACKWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
#endif
    ! calculate coordinate contributions to 1/k^2 or 1/r^2, respectively
    CALL initiq(nx2,dx,iqx)
    CALL initiq(ny2,dy,iqy)
    CALL initiq(nz2,dz,iqz)
    ! combine proper values, using dummy value at (1,1,1) to avoid div0
    FORALL(k=1:nz2,j=1:ny2,i=1:nx2) q(i,j,k)=iqx(i)+iqy(j)+iqz(k)
    q(1,1,1)=1.D0
    IF(periodic) THEN
       q=1.D0/REAL(q)
       q(1,1,1)=(0.D0,0.D0)
    ELSE
       q=1.D0/SQRT(REAL(q))
       q(1,1,1)=2.84D0/(dx*dy*dz)**(1.D0/3.D0)
#ifdef CUDA
       ! Copy to device, perform FFT, copy back to host
       q_d = q
       CALL cufftExecZ2Z(coulplan1,q_d,q_d,CUFFT_FORWARD)
       q = q_d
#else
       CALL dfftw_execute_dft(coulplan1,q,q)
#endif
    END IF
    DEALLOCATE(iqx,iqy,iqz)
  END SUBROUTINE coulinit
  !***************************************************
  SUBROUTINE initiq(n,d,iq)
    INTEGER,INTENT(IN) :: n
    REAL(db),INTENT(IN) :: d
    REAL(db),INTENT(OUT) :: iq(:)
    INTEGER :: i,ii
    DO i=1,n
       IF(i<=n/2) THEN
          ii=i-1
       ELSE
          ii=i-n-1
       ENDIF
       IF(periodic) THEN
          iq(i)=(2.D0*pi*ii/(n*d))**2
       ELSE
          iq(i)=(d*ii)**2
       END IF
    ENDDO
  END SUBROUTINE initiq
END MODULE Coulomb
