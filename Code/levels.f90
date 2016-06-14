MODULE Levels
  USE Params, ONLY: db,pi
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz,bangx,bangy,bangz
  USE Fourier
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
  SUBROUTINE cdervx0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: ix
    kfac=(PI+PI)/(dx*nx)
    CALL dfftw_execute_dft(xforward,psin,d1psout)
    IF(PRESENT(d2psout)) THEN
       DO ix=1,nx/2
          d2psout(ix,:,:,:)=-((ix-1)*kfac)**2*d1psout(ix,:,:,:)/REAL(nx)
          d2psout(nx-ix+1,:,:,:)=-(ix*kfac)**2*d1psout(nx-ix+1,:,:,:)/REAL(nx)
       ENDDO
       CALL dfftw_execute_dft(xbackward,d2psout,d2psout)
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
    CALL dfftw_execute_dft(xbackward,d1psout,d1psout)
  END SUBROUTINE cdervx0
  !*************************************************************
  SUBROUTINE cdervx(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangx)>0.00001) THEN
       DO i=1,nx
          p_temp(i,:,:,:)=psin(i,:,:,:)*exp(CMPLX(0.0d0,-bangx/nx*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervx0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervx0(p_temp,d1psout)
    END IF
    IF(abs(bangx)>0.00001) THEN
       DO i=1,nx
          d1psout(i,:,:,:)=d1psout(i,:,:,:)*exp(CMPLX(0.0d0,bangx/nx*i,db))&
               +CMPLX(0.0d0,bangx/(nx*dx),db)*psin(i,:,:,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,nx
             d2psout(i,:,:,:)=d2psout(i,:,:,:)*exp(CMPLX(0.0d0,bangx/nx*i,db))&
                  +(bangx/(nx*dx))**2*psin(i,:,:,:)&
                  +CMPLX(0.0d0,2.0d0*bangx/(nx*dx),db)*d1psout(i,:,:,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervx
  !************************************************************
  SUBROUTINE cdervy0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: iy
    kfac=(PI+PI)/(dy*ny)
    CALL dfftw_execute_dft(yforward,psin,d1psout)
    IF(PRESENT(d2psout)) THEN
       DO iy=1,ny/2
          d2psout(:,iy,:,:)=-((iy-1)*kfac)**2*d1psout(:,iy,:,:)/REAL(ny)
          d2psout(:,ny-iy+1,:,:)=-(iy*kfac)**2*d1psout(:,ny-iy+1,:,:)/REAL(ny)
       ENDDO
       CALL dfftw_execute_dft(ybackward,d2psout,d2psout)
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
    CALL dfftw_execute_dft(ybackward,d1psout,d1psout)
  END SUBROUTINE cdervy0
  !************************************************************
  SUBROUTINE cdervy(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangy)>0.00001) THEN
       DO i=1,ny
          p_temp(:,i,:,:)=psin(:,i,:,:)*exp(CMPLX(0.0d0,-bangy/ny*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervy0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervy0(p_temp,d1psout)
    END IF
    IF(abs(bangy)>0.00001) THEN
       DO i=1,ny
          d1psout(:,i,:,:)=d1psout(:,i,:,:)*exp(CMPLX(0.0d0,bangy/ny*i,db))&
               +CMPLX(0.0d0,bangy/(ny*dy),db)*psin(:,i,:,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,ny
             d2psout(:,i,:,:)=d2psout(:,i,:,:)*exp(CMPLX(0.0d0,bangy/ny*i,db))&
                  +(bangy/(ny*dy))**2*psin(:,i,:,:)&
                  +CMPLX(0.0d0,2.0d0*bangy/(ny*dy),db)*d1psout(:,i,:,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervy
  !************************************************************
  SUBROUTINE cdervz0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: iz
    kfac=(PI+PI)/(dz*nz)
    CALL dfftw_execute_dft(zforward,psin,d1psout)
    IF(PRESENT(d2psout)) THEN
       DO iz=1,nz/2
          d2psout(:,:,iz,:)=-((iz-1)*kfac)**2*d1psout(:,:,iz,:)/REAL(nz)
          d2psout(:,:,nz-iz+1,:)=-(iz*kfac)**2*d1psout(:,:,nz-iz+1,:)/REAL(nz)
       ENDDO
       CALL dfftw_execute_dft(zbackward,d2psout,d2psout)
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
    CALL dfftw_execute_dft(zbackward,d1psout,d1psout)
  END SUBROUTINE cdervz0
  !************************************************************
  SUBROUTINE cdervz(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangz)>0.00001) THEN
       DO i=1,nz
          p_temp(:,:,i,:)=psin(:,:,i,:)*exp(CMPLX(0.0d0,-bangz/nz*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervz0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervz0(p_temp,d1psout)
    END IF
    IF(abs(bangz)>0.00001) THEN
       DO i=1,nz
          d1psout(:,:,i,:)=d1psout(:,:,i,:)*exp(CMPLX(0.0d0,bangz/nz*i,db))&
               +CMPLX(0.0d0,bangz/(nz*dz),db)*psin(:,:,i,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,nz
             d2psout(:,:,i,:)=d2psout(:,:,i,:)*exp(CMPLX(0.0d0,bangz/nz*i,db))&
                  +(bangz/(nz*dz))**2*psin(:,:,i,:)&
                  +CMPLX(0.0d0,2.0d0*bangz/(nz*dz),db)*d1psout(:,:,i,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervz
  !************************************************************
  SUBROUTINE laplace(psin,psout,e0inv)  
    USE Forces, ONLY: h2ma
    USE Grids, ONLY: dx,dy,dz
    COMPLEX(db), INTENT(IN)   :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT)  :: psout(:,:,:,:)
    REAL(db), INTENT(IN), OPTIONAL :: e0inv
    REAL(db) :: kfacx, kfacy, kfacz
    REAL(db) :: k2facx(nx),k2facy(ny),k2facz(nz)
    INTEGER :: ix, iy, iz, is
    IF(.NOT.PRESENT(e0inv).AND.(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001))&
      STOP 'Laplace does not work with Bloch boundaries'
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
    IF(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          psout(ix,iy,iz,:)=psin(ix,iy,iz,:)*exp(CMPLX(0.0d0,-bangx/nx*ix,db))*exp(CMPLX(0.0d0,-bangy/ny*iy,db))*&
               exp(CMPLX(0.0d0,-bangz/nz*iz,db))
       END FORALL
    ELSE
       psout=psin
    END IF
    DO is=1,2
       CALL dfftw_execute_dft(pforward,psout(:,:,:,is),psout(:,:,:,is))
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
       CALL dfftw_execute_dft(pbackward,psout(:,:,:,is),psout(:,:,:,is))
    END DO
    IF(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          psout(ix,iy,iz,:)=psout(ix,iy,iz,:)*exp(CMPLX(0.0d0,bangx/nx*ix,db))*exp(CMPLX(0.0d0,bangy/ny*iy,db))*&
               exp(CMPLX(0.0d0,bangz/nz*iz,db))
       END FORALL
    END IF
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
             omptemp=0.0d0
             DO j=npmin(iq),nst-1
                oij=overlap(psi(:,:,:,:,j),psi(:,:,:,:,nst))
                omptemp(:,:,:,:)=omptemp(:,:,:,:)+oij*psi(:,:,:,:,j)
             END DO
             psi(:,:,:,:,nst)=psi(:,:,:,:,nst)-omptemp
             sp_norm(nst)=rpsnorm(psi(:,:,:,:,nst))
             psi(:,:,:,:,nst)=psi(:,:,:,:,nst)/SQRT(sp_norm(nst))
          END DO
       END IF
    END DO
  END SUBROUTINE schmid
END MODULE Levels
