MODULE Levels
  USE Params, ONLY: db,pi
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz
  USE Fourier
#ifdef CUDA
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
    COMPLEX(db), DEVICE ::  psin_d(:,:,:,:)
    COMPLEX(db), DEVICE:: d1psout_d(:,:,:,:)
    COMPLEX(db), DEVICE :: d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: ix
    kfac=(PI+PI)/(dx*nx)
#ifdef CUDA
    ! Copy to device, perform FFT, copy back to host
    psin_d = psin
    d1psout_d = d1psout
    CALL cufftExecZ2Z(xforward,psin_d,d1psout_d,CUFFT_FORWARD)
    psin = psin_d
    d1psout = d1psout_d
#else
    CALL dfftw_execute_dft(xforward,psin,d1psout)
#endif
    IF(PRESENT(d2psout)) THEN
       DO ix=1,nx/2
          d2psout(ix,:,:,:)=-((ix-1)*kfac)**2*d1psout(ix,:,:,:)/REAL(nx)
          d2psout(nx-ix+1,:,:,:)=-(ix*kfac)**2*d1psout(nx-ix+1,:,:,:)/REAL(nx)
       ENDDO
#ifdef CUDA
       d2psout_d = d2psout
       CALL cufftExecZ2Z(xbackward,d2psout_d,d2psout_d,CUFFT_INVERSE)
       d2psout = d2psout_d
#else
       CALL dfftw_execute_dft(xbackward,d2psout,d2psout)
#endif
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
#ifdef CUDA
    d1psout_d = d1psout
    CALL cufftExecZ2Z(xbackward,d1psout_d,d1psout_d,CUFFT_INVERSE)
    d1psout = d1psout_d
#else
    CALL dfftw_execute_dft(xbackward,d1psout,d1psout)
#endif
  END SUBROUTINE cdervx
  !************************************************************
  SUBROUTINE cdervy(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)
#ifdef CUDA
    COMPLEX(db), DEVICE :: psin_d(:,:,:,:)
    COMPLEX(db), DEVICE :: d1psout_d(:,:,:,:)
    COMPLEX(db), DEVICE :: d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iy,is,k
    kfac=(PI+PI)/(dy*ny)
    DO is=1,2
       DO k=1,nz
#ifdef CUDA
          ! Copy to device, perform FFT, copy back to host
          psin_d = psin
          d1psout_d = d1psout
          CALL cufftExecZ2Z(yforward,psin_d(:,:,k,is),d1psout_d(:,:,k,is),CUFFT_FORWARD)
          psin = psin_d
          d1psout = d1psout_d
#else
          CALL dfftw_execute_dft(yforward,psin(:,:,k,is),d1psout(:,:,k,is))
#endif
       END DO
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iy=1,ny/2
          d2psout(:,iy,:,:)=-((iy-1)*kfac)**2*d1psout(:,iy,:,:)/REAL(ny)
          d2psout(:,ny-iy+1,:,:)=-(iy*kfac)**2*d1psout(:,ny-iy+1,:,:)/REAL(ny)
       ENDDO
       DO is=1,2
          DO k=1,nz
#ifdef CUDA
             d2psout_d = d2psout
             CALL cufftExecZ2Z(ybackward,d2psout_d(:,:,k,is),d2psout_d(:,:,k,is),CUFFT_INVERSE)
             d2psout = d2psout_d
#else
             CALL dfftw_execute_dft(ybackward,d2psout(:,:,k,is),d2psout(:,:,k,is))
#endif
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
#ifdef CUDA
          d1psout_d = d1psout
          CALL cufftExecZ2Z(ybackward,d1psout_d(:,:,k,is),d1psout_d(:,:,k,is),CUFFT_INVERSE)
          d1psout = d1psout_d
#else
          CALL dfftw_execute_dft(ybackward,d1psout(:,:,k,is),d1psout(:,:,k,is))
#endif
       END DO
    END DO
  END SUBROUTINE cdervy
  !************************************************************
  SUBROUTINE cdervz(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)
#ifdef CUDA
    COMPLEX(db), DEVICE :: psin_d(:,:,:,:)
    COMPLEX(db), DEVICE :: d1psout_d(:,:,:,:)
    COMPLEX(db), DEVICE :: d2psout_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iz,is
    kfac=(PI+PI)/(dz*nz)
    DO is=1,2
#ifdef CUDA
       ! Copy to device, perform FFT, copy back to host
       psin_d = psin
       d1psout_d = d1psout
       CALL cufftExecZ2Z(zforward,psin_d(:,:,:,is),d1psout_d(:,:,:,is),CUFFT_FORWARD)
       psin = psin_d
       d1psout = d1psout_d
#else
       CALL dfftw_execute_dft(zforward,psin(:,:,:,is),d1psout(:,:,:,is))
#endif
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iz=1,nz/2
          d2psout(:,:,iz,:)=-((iz-1)*kfac)**2*d1psout(:,:,iz,:)/REAL(nz)
          d2psout(:,:,nz-iz+1,:)=-(iz*kfac)**2*d1psout(:,:,nz-iz+1,:)/REAL(nz)
       ENDDO
       DO is=1,2
#ifdef CUDA
          d2psout_d = d2psout
          CALL cufftExecZ2Z(zbackward,d2psout_d(:,:,:,is),d2psout_d(:,:,:,is),CUFFT_INVERSE)
          d2psout = d2psout_d
#else
          CALL dfftw_execute_dft(zbackward,d2psout(:,:,:,is),d2psout(:,:,:,is))
#endif
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
#ifdef CUDA
       d1psout_d = d1psout
       CALL cufftExecZ2Z(zbackward,d1psout_d(:,:,:,is),d1psout_d(:,:,:,is),CUFFT_INVERSE)
       d1psout = d1psout_d
#else
       CALL dfftw_execute_dft(zbackward,d1psout(:,:,:,is),d1psout(:,:,:,is))
#endif
    END DO
  END SUBROUTINE cdervz
  !************************************************************
  SUBROUTINE laplace(psin,psout,e0inv)  
    USE Forces, ONLY: h2ma
    USE Grids, ONLY: dx,dy,dz
    COMPLEX(db), INTENT(IN)   :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT)  :: psout(:,:,:,:)
#ifdef CUDA
    COMPLEX(db), DEVICE :: psout_d(:,:,:,:)
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
    DO is=1,2
#ifdef CUDA
       ! Copy to device, perform FFT, copy back to host
       psout_d = psout
       CALL cufftExecZ2Z(pforward,psout_d(:,:,:,is),psout_d(:,:,:,is),CUFFT_FORWARD)
       psout = psout_d
#else
       CALL dfftw_execute_dft(pforward,psout(:,:,:,is),psout(:,:,:,is))
#endif
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
#ifdef CUDA
       psout_d = psout
       CALL cufftExecZ2Z(pbackward,psout_d(:,:,:,is),psout_d(:,:,:,is),CUFFT_INVERSE)
       psout = psout_d
#else
       CALL dfftw_execute_dft(pbackward,psout(:,:,:,is),psout(:,:,:,is))
#endif
    END DO
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
