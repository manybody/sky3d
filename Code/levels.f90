MODULE Levels
  USE Params, ONLY: db,pi
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz
  USE Fourier
#ifdef CUDA
  USE cudafor
  USE cufft_m
#endif
  IMPLICIT NONE
  SAVE
  INTEGER :: nstmax,nstloc,nneut,nprot,npmin(2),npsi(2)
  REAL(db) :: charge_number,mass_number
  COMPLEX(db),POINTER :: psi(:,:,:,:,:)  
  COMPLEX(db), ALLOCATABLE :: hmatr(:,:)
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_energy,sp_efluct1, &
       sp_kinetic,sp_norm,sp_efluct2,sp_parity,wocc
  REAL(db),ALLOCATABLE :: sp_orbital(:,:),sp_spin(:,:)
  INTEGER, ALLOCATABLE :: isospin(:)
  REAL(db) :: start_levels, finish_levels, timer
  COMMON /Timer/ timer
CONTAINS
  !************************************************************
  SUBROUTINE alloc_levels
    ALLOCATE(psi(nx,ny,nz,2,nstloc), &
         sp_energy(nstmax),sp_efluct1(nstmax),sp_kinetic(nstmax),& 
         sp_norm(nstmax),sp_efluct2(nstmax),sp_parity(nstmax), &
         sp_orbital(3,nstmax),sp_spin(3,nstmax),wocc(nstmax), &
         isospin(nstmax))
    isospin(1:npsi(1))=1
    isospin(npmin(2):npsi(2))=2
  END SUBROUTINE alloc_levels
  !************************************************************
  SUBROUTINE cdervx(psin,d1psout,d2psout)
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)
#ifdef CUDA
    INTEGER :: istat, psin_size, d1psout_size, d2psout_size
    COMPLEX(db), ALLOCATABLE, DEVICE :: psin_d(:,:,:,:), d1psout_d(:,:,:,:), d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: ix
    kfac=(PI+PI)/(dx*nx)
    CALL cpu_time(start_levels)
#ifdef CUDA
    ALLOCATE(psin_d, MOLD=psin)
    ALLOCATE(d1psout_d, MOLD=d1psout)
    psin_size = SIZE(psin)
    d1psout_size = SIZE(d1psout)
    ! Copy to device, perform FFT, copy back to host
    istat = cudaMemcpy(psin_d(1,1,1,1), psin(1,1,1,1), psin_size)
    istat = cudaMemcpy(d1psout_d(1,1,1,1), d1psout(1,1,1,1), d1psout_size)
    CALL cufftExecZ2Z(xforward,psin_d,d1psout_d,CUFFT_FORWARD)
    istat = cudaMemcpy(psin(1,1,1,1), psin_d(1,1,1,1), psin_size)
    istat = cudaMemcpy(d1psout(1,1,1,1), d1psout_d(1,1,1,1), d1psout_size)
    istat = cudaDeviceSynchronize()
#else
    CALL dfftw_execute_dft(xforward,psin,d1psout)
#endif
    CALL cpu_time(finish_levels)
    timer=timer+finish_levels-start_levels
    WRITE(*,*) "FFT_DEBUG: levels.1 = ", finish_levels-start_levels
    WRITE(*,*) "FFT_DEBUG: timer = ", timer
    IF(PRESENT(d2psout)) THEN
       DO ix=1,nx/2
          d2psout(ix,:,:,:)=-((ix-1)*kfac)**2*d1psout(ix,:,:,:)/REAL(nx)
          d2psout(nx-ix+1,:,:,:)=-(ix*kfac)**2*d1psout(nx-ix+1,:,:,:)/REAL(nx)
       ENDDO
       CALL cpu_time(start_levels)
#ifdef CUDA
       ALLOCATE(d2psout_d, MOLD=d2psout)
       d2psout_size = SIZE(d2psout)

       istat = cudaMemcpy(d2psout_d(1,1,1,1), d2psout(1,1,1,1), d2psout_size)
       CALL cufftExecZ2Z(xbackward,d2psout_d,d2psout_d,CUFFT_INVERSE)
       istat = cudaMemcpy(d2psout(1,1,1,1), d2psout_d(1,1,1,1), d2psout_size)
       istat = cudaDeviceSynchronize()
#else
       CALL dfftw_execute_dft(xbackward,d2psout,d2psout)
#endif
       CALL cpu_time(finish_levels)
       timer=timer+finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: levels.2 = ", finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: timer = ", timer
    ENDIF
    d1psout(1,:,:,:)=(0.D0,0.D0)
    DO ix=2,nx/2
       d1psout(ix,:,:,:)=CMPLX(0.0D0,((ix-1)*kfac),db)*d1psout(ix,:,:,:)/REAL(nx)
    ENDDO
    DO ix=1,nx/2-1
       d1psout(nx-ix+1,:,:,:)=CMPLX(0.0D0,-(ix*kfac),db)*d1psout(nx-ix+1,:,:,:) &
            /REAL(nx)
    ENDDO
    d1psout(nx/2+1,:,:,:)=(0.D0,0.D0)
    CALL cpu_time(start_levels)
#ifdef CUDA
    istat = cudaMemcpy(d1psout_d(1,1,1,1), d1psout(1,1,1,1), d1psout_size)
    CALL cufftExecZ2Z(xbackward,d1psout_d,d1psout_d,CUFFT_INVERSE)
    istat = cudaMemcpy(d1psout(1,1,1,1), d1psout_d(1,1,1,1), d1psout_size)
    istat = cudaDeviceSynchronize()

    DEALLOCATE(psin_d, d1psout_d, d2psout_d)
#else
    CALL dfftw_execute_dft(xbackward,d1psout,d1psout)
#endif
    CALL cpu_time(finish_levels)
    timer=timer+finish_levels-start_levels
    WRITE(*,*) "FFT_DEBUG: levels.3 = ", finish_levels-start_levels
    WRITE(*,*) "FFT_DEBUG: timer = ", timer
  END SUBROUTINE cdervx
  !************************************************************
  SUBROUTINE cdervy(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)
#ifdef CUDA
    INTEGER :: istat, psin_size, d1psout_size, d2psout_size
    COMPLEX(db), ALLOCATABLE, DEVICE :: psin_d(:,:,:,:), d1psout_d(:,:,:,:), d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iy,is,k
    kfac=(PI+PI)/(dy*ny)
#ifdef CUDA
    ALLOCATE(psin_d, MOLD=psin)
    ALLOCATE(d1psout_d, MOLD=d1psout)
    ALLOCATE(d2psout_d, MOLD=d2psout)
    psin_size = SIZE(psin, 1) + SIZE(psin, 2)
    d1psout_size = SIZE(d1psout, 1) + SIZE(d1psout, 2)
    d2psout_size = SIZE(d2psout, 1) + SIZE(d2psout, 2)
#endif
    DO is=1,2
       DO k=1,nz
          CALL cpu_time(start_levels)
#ifdef CUDA
          ! Copy to device, perform FFT, copy back to host
          istat = cudaMemcpy(psin_d(1,1,k,is), psin(1,1,k,is), psin_size)
          istat = cudaMemcpy(d1psout_d(1,1,k,is), d1psout(1,1,k,is), d1psout_size)
          CALL cufftExecZ2Z(yforward,psin_d(:,:,k,is),d1psout_d(:,:,k,is),CUFFT_FORWARD)
          istat = cudaMemcpy(psin(1,1,k,is), psin_d(1,1,k,is), psin_size)
          istat = cudaMemcpy(d1psout(1,1,k,is), d1psout_d(1,1,k,is), d1psout_size)
          istat = cudaDeviceSynchronize()
#else
          CALL dfftw_execute_dft(yforward,psin(:,:,k,is),d1psout(:,:,k,is))
#endif
          CALL cpu_time(finish_levels)
          timer=timer+finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: levels.4 = ", finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: timer = ", timer
       END DO
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iy=1,ny/2
          d2psout(:,iy,:,:)=-((iy-1)*kfac)**2*d1psout(:,iy,:,:)/REAL(ny)
          d2psout(:,ny-iy+1,:,:)=-(iy*kfac)**2*d1psout(:,ny-iy+1,:,:)/REAL(ny)
       ENDDO
       DO is=1,2
          DO k=1,nz
             CALL cpu_time(start_levels)
#ifdef CUDA
             istat = cudaMemcpy(d2psout_d(1,1,k,is), d2psout(1,1,k,is), d2psout_size)
             CALL cufftExecZ2Z(ybackward,d2psout_d(:,:,k,is),d2psout_d(:,:,k,is),CUFFT_INVERSE)
             istat = cudaMemcpy(d2psout(1,1,k,is), d2psout_d(1,1,k,is), d2psout_size)
             istat = cudaDeviceSynchronize()
#else
             CALL dfftw_execute_dft(ybackward,d2psout(:,:,k,is),d2psout(:,:,k,is))
#endif
             CALL cpu_time(finish_levels)
             timer=timer+finish_levels-start_levels
             WRITE(*,*) "FFT_DEBUG: levels.5 = ", finish_levels-start_levels
             WRITE(*,*) "FFT_DEBUG: timer = ", timer
          END DO
       END DO
    ENDIF
    d1psout(:,1,:,:)=(0.D0,0.D0)
    DO iy=2,ny/2
       d1psout(:,iy,:,:)=CMPLX(0.0D0,((iy-1)*kfac),db)*d1psout(:,iy,:,:)/REAL(ny)
    ENDDO
    DO iy=1,ny/2-1
       d1psout(:,ny-iy+1,:,:)=CMPLX(0.0D0,-(iy*kfac),db)*d1psout(:,ny-iy+1,:,:) &
            /REAL(ny)
    ENDDO
    d1psout(:,ny/2+1,:,:)=(0.D0,0.D0)
    DO is=1,2
       DO k=1,nz
          CALL cpu_time(start_levels)
#ifdef CUDA
          istat = cudaMemcpy(d1psout_d(1,1,k,is), d1psout(1,1,k,is), d1psout_size)
          CALL cufftExecZ2Z(ybackward,d1psout_d(:,:,k,is),d1psout_d(:,:,k,is),CUFFT_INVERSE)
          istat = cudaMemcpy(d1psout(1,1,k,is), d1psout_d(1,1,k,is), d1psout_size)
          istat = cudaDeviceSynchronize()
#else
          CALL dfftw_execute_dft(ybackward,d1psout(:,:,k,is),d1psout(:,:,k,is))
#endif
          CALL cpu_time(finish_levels)
          timer=timer+finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: levels.6 = ", finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: timer = ", timer
       END DO
    END DO
#ifdef CUDA
    DEALLOCATE(psin_d, d1psout_d, d2psout_d)
#endif
  END SUBROUTINE cdervy
  !************************************************************
  SUBROUTINE cdervz(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)
#ifdef CUDA
    INTEGER :: istat, psin_size, d1psout_size, d2psout_size
    COMPLEX(db), ALLOCATABLE, DEVICE :: psin_d(:,:,:,:), d1psout_d(:,:,:,:), d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iz,is
    kfac=(PI+PI)/(dz*nz)
#ifdef CUDA
    ALLOCATE(psin_d, MOLD=psin)
    ALLOCATE(d1psout_d, MOLD=d1psout)
    ALLOCATE(d2psout_d, MOLD=d2psout)
    psin_size = SIZE(psin, 1) + SIZE(psin, 2) + SIZE(psin, 3)
    d1psout_size = SIZE(d1psout, 1) + SIZE(d1psout, 2) + SIZE(d1psout, 3)
    d2psout_size = SIZE(d2psout, 1) + SIZE(d2psout, 2) + SIZE(d2psout, 3)
#endif
    DO is=1,2
       CALL cpu_time(start_levels)
#ifdef CUDA
       ! Copy to device, perform FFT, copy back to host
       istat = cudaMemcpy(psin_d(1,1,1,is), psin(1,1,1,is), psin_size)
       istat = cudaMemcpy(d1psout_d(1,1,1,is), d1psout(1,1,1,is), d1psout_size)
       CALL cufftExecZ2Z(zforward,psin_d(:,:,:,is),d1psout_d(:,:,:,is),CUFFT_FORWARD)
       istat = cudaMemcpy(psin(1,1,1,is), psin_d(1,1,1,is), psin_size)
       istat = cudaMemcpy(d1psout(1,1,1,is), d1psout_d(1,1,1,is), d1psout_size)
       istat = cudaDeviceSynchronize()
#else
       CALL dfftw_execute_dft(zforward,psin(:,:,:,is),d1psout(:,:,:,is))
#endif
       CALL cpu_time(finish_levels)
       timer=timer+finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: levels.7 = ", finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: timer = ", timer
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iz=1,nz/2
          d2psout(:,:,iz,:)=-((iz-1)*kfac)**2*d1psout(:,:,iz,:)/REAL(nz)
          d2psout(:,:,nz-iz+1,:)=-(iz*kfac)**2*d1psout(:,:,nz-iz+1,:)/REAL(nz)
       ENDDO
       DO is=1,2
          CALL cpu_time(start_levels)
#ifdef CUDA
          istat = cudaMemcpy(d2psout_d(1,1,1,is), d2psout(1,1,1,is), d2psout_size)
          CALL cufftExecZ2Z(zbackward,d2psout_d(:,:,:,is),d2psout_d(:,:,:,is),CUFFT_INVERSE)
          istat = cudaMemcpy(d2psout(1,1,1,is), d2psout_d(1,1,1,is), d2psout_size)
          istat = cudaDeviceSynchronize()
#else
          CALL dfftw_execute_dft(zbackward,d2psout(:,:,:,is),d2psout(:,:,:,is))
#endif
          CALL cpu_time(finish_levels)
          timer=timer+finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: levels.8 = ", finish_levels-start_levels
          WRITE(*,*) "FFT_DEBUG: timer = ", timer
       END DO
    ENDIF
    d1psout(:,:,1,:)=(0.D0,0.D0)
    DO iz=2,nz/2
       d1psout(:,:,iz,:)=CMPLX(0.0D0,((iz-1)*kfac),db)*d1psout(:,:,iz,:)/REAL(nz)
    ENDDO
    DO iz=1,nz/2-1
       d1psout(:,:,nz-iz+1,:)=CMPLX(0.0D0,-(iz*kfac),db)*d1psout(:,:,nz-iz+1,:) &
            /REAL(nz)
    ENDDO
    d1psout(:,:,nz/2+1,:)=(0.D0,0.D0)
    DO is=1,2
       CALL cpu_time(start_levels)
#ifdef CUDA
       istat = cudaMemcpy(d1psout_d(1,1,1,is), d1psout(1,1,1,is), d1psout_size)
       CALL cufftExecZ2Z(zbackward,d1psout_d(:,:,:,is),d1psout_d(:,:,:,is),CUFFT_INVERSE)
       istat = cudaMemcpy(d1psout(1,1,1,is), d1psout_d(1,1,1,is), d1psout_size)
       istat = cudaDeviceSynchronize()
#else
       CALL dfftw_execute_dft(zbackward,d1psout(:,:,:,is),d1psout(:,:,:,is))
#endif
       CALL cpu_time(finish_levels)
       timer=timer+finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: levels.9 = ", finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: timer = ", timer
    END DO
#ifdef CUDA
    DEALLOCATE(psin_d, d1psout_d, d2psout_d)
#endif
  END SUBROUTINE cdervz
  !************************************************************
  SUBROUTINE laplace(psin,psout,e0inv)  
    USE Forces, ONLY: h2ma
    USE Grids, ONLY: dx,dy,dz
    COMPLEX(db), INTENT(IN)   :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT)  :: psout(:,:,:,:)
#ifdef CUDA
    INTEGER :: istat, psout_size
    COMPLEX(db), ALLOCATABLE, DEVICE :: psout_d(:,:,:,:)
#endif
    REAL(db), INTENT(IN), OPTIONAL :: e0inv
    REAL(db) :: kfacx, kfacy, kfacz
    REAL(db) :: k2facx(nx),k2facy(ny),k2facz(nz)
    INTEGER :: ix, iy, iz, is
    kfacz=(PI+PI)/(dz*nz)
    kfacy=(PI+PI)/(dy*ny)
    kfacx=(PI+PI)/(dx*nx)
    DO iz=1,nz/2
       k2facz(iz)=-((iz-1)*kfacz)**2
       k2facz(nz-iz+1)=-(iz*kfacz)**2
    ENDDO
    DO iy=1,ny/2
       k2facy(iy)=-((iy-1)*kfacy)**2
       k2facy(ny-iy+1)=-(iy*kfacy)**2
    ENDDO
    DO ix=1,nx/2
       k2facx(ix)=-((ix-1)*kfacx)**2
       k2facx(nx-ix+1)=-(ix*kfacx)**2
    ENDDO
    psout=psin
#ifdef CUDA
    ALLOCATE(psout_d, MOLD=psout)
    psout_size = SIZE(psout, 1) + SIZE(psout, 2) + SIZE(psout, 3)
#endif
    DO is=1,2
       CALL cpu_time(start_levels)
#ifdef CUDA
       ! Copy to device, perform FFT, copy back to host
       istat = cudaMemcpy(psout_d(1,1,1,is), psout(1,1,1,is), psout_size)
       CALL cufftExecZ2Z(pforward,psout_d(:,:,:,is),psout_d(:,:,:,is),CUFFT_FORWARD)
       istat = cudaMemcpy(psout(1,1,1,is), psout_d(1,1,1,is), psout_size)
       istat = cudaDeviceSynchronize()
#else
       CALL dfftw_execute_dft(pforward,psout(:,:,:,is),psout(:,:,:,is))
#endif
       CALL cpu_time(finish_levels)
       timer=timer+finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: levels.10 = ", finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: timer = ", timer
    END DO
    IF(PRESENT(e0inv)) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz,is=1:2)
          psout(ix,iy,iz,is)=psout(ix,iy,iz,is)  &
               /(e0inv-h2ma*(k2facx(ix)+k2facy(iy)+k2facz(iz))) &
               / DBLE(nx*ny*nz)
       END FORALL
    ELSE
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz,is=1:2)
          psout(ix,iy,iz,is)=psout(ix,iy,iz,is)  &
               *(k2facx(ix)+k2facy(iy)+k2facz(iz)) &
               /DBLE(nx*ny*nz)
       END FORALL
    ENDIF
    DO is=1,2
       CALL cpu_time(start_levels)
#ifdef CUDA
       istat = cudaMemcpy(psout_d(1,1,1,is), psout(1,1,1,is), psout_size)
       CALL cufftExecZ2Z(pbackward,psout_d(:,:,:,is),psout_d(:,:,:,is),CUFFT_INVERSE)
       istat = cudaMemcpy(psout(1,1,1,is), psout_d(1,1,1,is), psout_size)
       istat = cudaDeviceSynchronize()
#else
       CALL dfftw_execute_dft(pbackward,psout(:,:,:,is),psout(:,:,:,is))
#endif
       CALL cpu_time(finish_levels)
       timer=timer+finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: levels.11 = ", finish_levels-start_levels
       WRITE(*,*) "FFT_DEBUG: timer = ", timer
    END DO
#ifdef CUDA
    DEALLOCATE(psout_d)
#endif
  END SUBROUTINE laplace
  !************************************************************
  SUBROUTINE schmid
    USE Trivial, ONLY: rpsnorm, overlap
    COMPLEX(db) :: oij,omptemp(nx,ny,nz,2)
    INTEGER :: nst,j,iq
    DO iq=1,2
       ! gram-schmidt procedure - neutrons
       IF(npsi(iq)-npmin(iq)+1 /= 1) THEN
          DO nst=npmin(iq),npsi(iq)
             omptemp=0
             !$OMP PARALLEL DO PRIVATE(j,oij) REDUCTION(+:omptemp) 
             DO j=npmin(iq),nst-1
                oij=overlap(psi(:,:,:,:,j),psi(:,:,:,:,nst))
                omptemp(:,:,:,:)=omptemp(:,:,:,:)+oij*psi(:,:,:,:,j)
             END DO
             !$OMP END PARALLEL DO
             psi(:,:,:,:,nst)=psi(:,:,:,:,nst)-omptemp
             sp_norm(nst)=rpsnorm(psi(:,:,:,:,nst))
             psi(:,:,:,:,nst)=psi(:,:,:,:,nst)/SQRT(sp_norm(nst))
          END DO
       END IF
    END DO
  END SUBROUTINE schmid
END MODULE Levels
