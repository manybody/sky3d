MODULE Densities
  USE Params, ONLY: db,tfft
  USE Grids, ONLY: nx,ny,nz,der1x,der1y,der1z
  USE Levels, ONLY: cdervx,cdervy,cdervz
  USE Trivial, ONLY: cmulx,cmuly,cmulz
  IMPLICIT NONE
  SAVE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: rho,tau
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: current,sdens,sodens
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_densities
    ALLOCATE(rho(nx,ny,nz,2),tau(nx,ny,nz,2),current(nx,ny,nz,3,2), &
         sdens(nx,ny,nz,3,2),sodens(nx,ny,nz,3,2))
  END SUBROUTINE alloc_densities
  !***********************************************************************
  SUBROUTINE add_density(iq,weight,psin,lrho,ltau,lcurrent,lsdens,lsodens)  
    COMPLEX(db),INTENT(INOUT) :: psin(nx,ny,nz,2)
    REAL(db),DIMENSION(:,:,:,:),INTENT(INOUT) :: lrho,ltau
    REAL(db),DIMENSION(:,:,:,:,:),INTENT(INOUT) :: lcurrent,lsdens,lsodens
    INTEGER,INTENT(IN) :: iq
    REAL(db),INTENT(IN) :: weight
    COMPLEX(db) :: ps1(nx,ny,nz,2)  
    INTEGER :: ix,iy,iz
    IF(weight<=0.D0) RETURN
    !***********************************************************************
    ! non-derivative terms
    !***********************************************************************
    lrho(:,:,:,iq)=lrho(:,:,:,iq)+weight* &
         (psin(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         psin(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(iz=1:nz,iy=1:ny,ix=1:nx)
      lsdens(ix,iy,iz,1,iq)=lsdens(ix,iy,iz,1,iq)+2.D0*weight &
           *REAL(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,2))
      lsdens(ix,iy,iz,2,iq)=lsdens(ix,iy,iz,2,iq)+2.D0*weight &
           *AIMAG(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,2))
      lsdens(ix,iy,iz,3,iq)=lsdens(ix,iy,iz,3,iq)+weight &
           *(REAL(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,1)) &
           -REAL(CONJG(psin(ix,iy,iz,2))*psin(ix,iy,iz,2)))
    END FORALL
    !***********************************************************************
    ! x-derivatives
    !***********************************************************************
    IF(TFFT) THEN
      CALL cdervx(psin,ps1)  
    ELSE
      CALL cmulx(der1x,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,1,iq)=lcurrent(:,:,:,1,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
      lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)-weight &
           *(  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
           - AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2))) )
      lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
           *(  REAL(psin(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,2))) &
           - REAL(psin(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1))) )
    END FORALL
    !***********************************************************************
    ! y-derivatives
    !***********************************************************************
    IF(TFFT) THEN
      CALL cdervy(psin,ps1)  
    ELSE
      CALL cmuly(der1y,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,2,iq)=lcurrent(:,:,:,2,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
      lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
           *AIMAG( ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1)) &
           -ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)) )
      lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
           *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
           +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
    END FORALL
    !***********************************************************************
    ! z-derivatives
    !***********************************************************************
    IF(TFFT) THEN
      CALL cdervz(psin,ps1)  
    ELSE
      CALL cmulz(der1z,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,3,iq)=lcurrent(:,:,:,3,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
      lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
           *REAL( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
           -ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
      lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)+weight &
           *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
           +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
    END FORALL
  END SUBROUTINE add_density
END MODULE Densities
