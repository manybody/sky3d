MODULE User
  USE Params
  USE Grids
  USE Levels
  USE Trivial, ONLY: rpsnorm
  IMPLICIT NONE
CONTAINS
  SUBROUTINE init_user
  IMPLICIT NONE
  INTEGER  :: nst, iq
  !
  INTEGER :: i,j,l,ii,jj,kk
  INTEGER :: kf(3,npsi(2)),temp_k(3)
  LOGICAL :: check
  INTEGER :: ki(3,8*7**3),ki_t(3,7**3)
  REAL(db) :: temp_e,temp_energies(8*7**3)
  WRITE(*,*)
  WRITE(*,*)'*****init plane waves:*****'
  psi=(0.d0,0.d0)
  !***********************************************************************
  !                                                                      *
  !                           calculate all k                            *
  !                                                                      *
  !***********************************************************************
    j=0
    ii=0
    jj=0
    kk=0
    DO i=1,7**3
      IF (ii==7) THEN
        ii=0
        jj=jj+1
      END IF
      IF (jj==7) THEN
        jj=0
        kk=kk+1
      END IF
      ki(1,i)=ii
      ki(2,i)=jj
      ki(3,i)=kk
      ii=ii+1
    END DO
    ki_t(:,1:7**3)=ki(:,1:7**3)
    l=1
    DO i=1,7**3
      DO j=1,8
        ki(:,l)=ki_t(:,i)
        IF(j==2.OR.j==4.OR.j==6.OR.j==8) THEN 
          ki(1,l)=-ki_t(1,i)
        END IF
        IF(j==3.OR.j==4.OR.j==7.OR.j==8) THEN
          ki(2,l)=-ki_t(2,i)
        END IF
        IF(j==5.OR.j==6.OR.j==7.OR.j==8) THEN
          ki(3,l)=-ki_t(3,i)
        END IF
        temp_energies(l)=epw(ki(1,l),ki(2,l),ki(3,l))
        l=l+1
      END DO
    END DO
    !insertion_sort
    DO i=2,8*7**3
      temp_e=temp_energies(i)
      temp_k(:)=ki(:,i)
      j=i
      DO WHILE (j>1 .AND. temp_energies(j-1)>temp_e)
        temp_energies(j)=temp_energies(j-1)
        ki(:,j)=ki(:,j-1)
        j=j-1
      END DO
      temp_energies(j)=temp_e
      ki(:,j)=temp_k(:)
    END DO
    nst = 1
  DO iq = 1,2  
    i=1 !counts nobs/2 (spin)
    j=1 !counts "-"-signs
    l=1 !counts ki
    DO WHILE(nst<=npsi(iq))
      kf(:,i)=ki(:,l)
      CALL check_kf(kf,i,check)
      IF(check) THEN 
        CALL background(nst,kf(1,i),kf(2,i),kf(3,i),1)
        nst=nst+1
        CALL background(nst,kf(1,i),kf(2,i),kf(3,i),-1)
        nst=nst+1
        i=i+1
      END IF
      l=l+1
    END DO
  END DO
  END SUBROUTINE init_user
  !
REAL(db) FUNCTION epw(kx,ky,kz) RESULT(e)
  USE FORCES, ONLY: nucleon_mass
  INTEGER,INTENT(IN) :: kx,ky,kz
  REAL(db) :: dx,dy,dz
  dx=x(2)-x(1)
  dy=y(2)-y(1)
  dz=z(2)-z(1)
  e=(hbc**2)/(2*nucleon_mass)*(((2*pi*kx+bangx)/(nx*dx))**2&
  +((2*pi*ky+bangy)/(ny*dy))**2+((2*pi*kz+bangz)/(nz*dz))**2)
END FUNCTION

SUBROUTINE background(nst,kx,ky,kz,s)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nst,kx,ky,kz,s
  INTEGER :: ix,iy,iz
  REAL(db) :: facx,facy,facz,norm
  COMPLEX(db) :: fy,fz
  wocc(nst)=0.0d0
  IF(nst>=npmin(1) .AND. (nst<npmin(1)+nneut)) wocc(nst)=1.0d0
  IF(nst>=npmin(2) .AND. (nst<npmin(2)+nprot)) wocc(nst)=1.0d0
  DO iz = 1,nz  
     facz=REAL(iz-1)*((2.D0*pi*REAL(kz)+bangz)/FLOAT(nz))
     fz=CMPLX(COS(facz),SIN(facz),db)
     DO iy=1,ny
        facy=REAL(iy-1)*((2.D0*pi*REAL(ky)+bangy)/FLOAT(ny))
        fy=CMPLX(COS(facy),SIN(facy),db)
        DO ix=1,nx
           facx=REAL(ix-1)*((2.D0*pi*REAL(kx)+bangx)/FLOAT(nx))
           IF(s>0) THEN
              psi(ix,iy,iz,1,nst)=fz*fy*CMPLX(COS(facx),SIN(facx),db)
              psi(ix,iy,iz,2,nst)=0.D0
           ELSE
              psi(ix,iy,iz,2,nst)=fz*fy*CMPLX(COS(facx),SIN(facx),db)
              psi(ix,iy,iz,1,nst)=0.D0
           END IF
        ENDDO
     ENDDO
  ENDDO
  norm=SQRT(rpsnorm(psi(:,:,:,:,nst)))
  psi(:,:,:,:,nst)=psi(:,:,:,:,nst)/norm
  WRITE(*,'(A14,3I2,A7,I2,A12,I4,A8,F9.5,A10,F6.3)')'state with k=(',kx,ky,kz,&
  '), spin= ',s,' at position',nst,' energy ',epw(kx,ky,kz),' and wocc=',wocc(nst)
END SUBROUTINE background
!
SUBROUTINE check_kf(k,i,check)
  INTEGER,INTENT(IN) :: k(3,npsi(2)),i
  LOGICAL,INTENT(OUT) :: check
  INTEGER :: j
  check=.TRUE.
  IF(i==1) RETURN
  DO j=1,i-1
    IF(k(1,j)==k(1,i).AND.k(2,j)==k(2,i).AND.k(3,j)==k(3,i)) check=.FALSE.
  END DO
END SUBROUTINE
END MODULE User
