MODULE Trivial
  USE Params, ONLY: db
  USE Grids, ONLY: nx,ny,nz,wxyz
  IMPLICIT NONE  
CONTAINS
  !***********************************************************************
  ! Multiply wave function by matrix in x-direction
  !***********************************************************************
  PURE SUBROUTINE cmulx(xmat,pinn,pout,ifadd)  
    REAL(db) :: xmat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: xmat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz
    IF(ifadd==0) pout=0.0D0
    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(xmat(ix,:)*pinn(:,iy,iz,is))
    END FORALL
  END SUBROUTINE cmulx
  !***********************************************************************
  ! Multiply wave function by matrix in y-direction
  !***********************************************************************
  PURE SUBROUTINE cmuly(ymat,pinn,pout,ifadd)  
    REAL(db) :: ymat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: ymat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz
    IF(ifadd==0) pout=0.0D0
    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(ymat(iy,:)*pinn(ix,:,iz,is))
    END FORALL
  END SUBROUTINE cmuly
  !***********************************************************************
  ! Multiply wave function by matrix in z-direction
  !***********************************************************************
  PURE SUBROUTINE cmulz(zmat,pinn,pout,ifadd)  
    REAL(db) :: zmat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: zmat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz,izz
    IF(ifadd==0) pout=0.0D0
    !
    !    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
    !       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(zmat(iz,:)*pinn(ix,iy,:,is))
    !    END FORALL
    DO is=1,2
       DO iz=1,nz
          DO izz=1,nz
             DO iy=1,ny
                DO ix=1,nx
                   pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+ &
                        zmat(iz,izz)*pinn(ix,iy,izz,is)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE cmulz
  !***********************************************************************
  ! Calculate total norm of wave function
  !***********************************************************************
  PURE FUNCTION rpsnorm(ps) RESULT(r)
    COMPLEX(db),INTENT(IN) :: ps(:,:,:,:)
    REAL(db) :: r
    r=wxyz*SUM(REAL(CONJG(ps)*ps))
  END FUNCTION rpsnorm
  !***********************************************************************
  ! Calculate overlap of two wave functions
  !***********************************************************************
  PURE FUNCTION overlap(pl,pr)  RESULT(c)
    COMPLEX(db) :: c,pl(:,:,:,:),pr(:,:,:,:)
    INTENT(IN) :: pl,pr
    c=wxyz*SUM(CONJG(pl)*pr)
  END FUNCTION overlap
  !***********************************************************************
  ! Multiply real field with matrix in x-direction
  !***********************************************************************
  PURE SUBROUTINE rmulx(xmat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: xmat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: xmat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(xmat(ix,:)*finn(:,iy,iz))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(xmat(ix,:)*finn(:,iy,iz))
       END FORALL
    ENDIF
  END SUBROUTINE rmulx
  !***********************************************************************
  ! Multiply real field with matrix in y-direction
  !***********************************************************************
  PURE SUBROUTINE rmuly(ymat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: ymat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: ymat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(ymat(iy,:)*finn(ix,:,iz))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(ymat(iy,:)*finn(ix,:,iz))
       END FORALL
    ENDIF
  END SUBROUTINE rmuly
  !***********************************************************************
  ! Multiply real field with matrix in z-direction
  !***********************************************************************
  PURE SUBROUTINE rmulz(zmat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: zmat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: zmat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(zmat(iz,:)*finn(ix,iy,:))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(zmat(iz,:)*finn(ix,iy,:))
       END FORALL
    ENDIF
    RETURN  
  END SUBROUTINE rmulz
  !***********************************************************************
END MODULE Trivial
