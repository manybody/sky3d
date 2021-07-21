MODULE Densities
  USE Params, ONLY: db,tfft
  USE Grids, ONLY: nx,ny,nz,der1x,der1y,der1z
  USE Levels, ONLY: cdervx,cdervy,cdervz
  USE Trivial, ONLY: cmulx,cmuly,cmulz
  USE Forces
  IMPLICIT NONE
  SAVE
   REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: rho,tau
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: current,sdens,sodens,tdens, &
  scurrentx,scurrenty,scurrentz, fdens,nablarho
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_densities
    ALLOCATE(rho(nx,ny,nz,2),tau(nx,ny,nz,2),current(nx,ny,nz,3,2), &
           sdens(nx,ny,nz,3,2),sodens(nx,ny,nz,3,2), tdens(nx,ny,nz,3,2), &
         scurrentx(nx,ny,nz,3,2),scurrenty(nx,ny,nz,3,2),scurrentz(nx,ny,nz,3,2), &
         fdens(nx,ny,nz,3,2))
         IF(mlocalize/=0)  ALLOCATE(nablarho(nx,ny,nz,3,2))
  END SUBROUTINE alloc_densities
  !***********************************************************************
  SUBROUTINE add_density(iq,weight,psin,lrho,ltau,lcurrent,lsdens,lsodens,&
              ltdens,lfdens,lscurrentx,lscurrenty,lscurrentz)
    COMPLEX(db),INTENT(INOUT) :: psin(nx,ny,nz,2)
    COMPLEX(db) :: psx(nx,ny,nz,2),psy(nx,ny,nz,2),psz(nx,ny,nz,2)
    REAL(db),INTENT(INOUT) :: lrho(nx,ny,nz,2),ltau(nx,ny,nz,2)
    REAL(db),INTENT(INOUT) :: lcurrent(nx,ny,nz,3,2),lsdens(nx,ny,nz,3,2)
    REAL(db),INTENT(INOUT) :: lsodens(nx,ny,nz,3,2)
    REAL(db),INTENT(INOUT) :: lfdens(nx,ny,nz,3,2),ltdens(nx,ny,nz,3,2)
    REAL(db),INTENT(INOUT) :: lscurrentx(nx,ny,nz,3,2),lscurrenty(nx,ny,nz,3,2)
    REAL(db),INTENT(INOUT) :: lscurrentz(nx,ny,nz,3,2)
    REAL(db) :: lnablarho(nx,ny,nz,3,2)
    INTEGER,INTENT(IN) :: iq
    REAL(db),INTENT(IN) :: weight
    COMPLEX(db) :: ps1(nx,ny,nz,2)  
    INTEGER :: ix,iy,iz
    ALLOCATE(ps1(nx,ny,nz,2))
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
     IF(tlocalize) lnablarho(:,:,:,1,iq)=lnablarho(:,:,:,1,iq)+weight* &
         REAL(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
     lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)-weight &
          *(  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
          - AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2))) )
     lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
          *(  REAL(psin(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,2))) &
          - REAL(psin(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1))) )
    END FORALL

     ! T_x components        
     IF(ston) THEN
     FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
     ltdens(ix,iy,iz,1,iq) = ltdens(ix,iy,iz,1,iq)+2.d0*weight &
          *(  REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,2,iq) = ltdens(ix,iy,iz,2,iq)+2.d0*weight &
          *(  AIMAG(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,3,iq) = ltdens(ix,iy,iz,3,iq)+weight &
          *(  REAL(ps1(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,1))) &
          -   REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,2)))  )
     END FORALL
     ENDIF

     ! x components J_xj 
     IF(jfon .OR. j2on) THEN
     FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
       lscurrentx(ix,iy,iz,1,iq) = lscurrentx(ix,iy,iz,1,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1))) &
          +  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)))  )
       lscurrentx(ix,iy,iz,2,iq) = lscurrentx(ix,iy,iz,2,iq)+weight* &
          (  REAL(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2))) &
          -  REAL(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)))  )
       lscurrentx(ix,iy,iz,3,iq) = lscurrentx(ix,iy,iz,3,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
          -  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)))  )
     END FORALL
     ENDIF
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
     IF(tlocalize) lnablarho(:,:,:,2,iq)=lnablarho(:,:,:,2,iq)+weight* &
         REAL(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
          *AIMAG( ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1)) &
          -ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)) )
     lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
          *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
          +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
     END FORALL

     ! T_y components
     IF(ston) THEN
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     ltdens(ix,iy,iz,1,iq) = ltdens(ix,iy,iz,1,iq)+2.d0*weight &
          *(  REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,2,iq) = ltdens(ix,iy,iz,2,iq)+2.d0*weight &
          *(  AIMAG(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,3,iq) = ltdens(ix,iy,iz,3,iq)+weight &
          *(  REAL(ps1(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,1))) &
          -   REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,2)))  )
     END FORALL
     ENDIF

     ! y components J_yj
     IF(jfon .OR. j2on) THEN
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     lscurrenty(ix,iy,iz,1,iq) = lscurrenty(ix,iy,iz,1,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1))) &
          +  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)))  )
     lscurrenty(ix,iy,iz,2,iq) = lscurrenty(ix,iy,iz,2,iq)+weight* &
          (  REAL(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2))) &
          -  REAL(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)))  )
     lscurrenty(ix,iy,iz,3,iq) = lscurrenty(ix,iy,iz,3,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
          -  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)))  )
     END FORALL
     ENDIF
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
     IF(tlocalize) lnablarho(:,:,:,3,iq)=lnablarho(:,:,:,3,iq)+weight* &
         REAL(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
          *REAL( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
          -ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
     lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)+weight &
          *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
          +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
     END FORALL

     ! T_z components
     IF(ston) THEN
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     ltdens(ix,iy,iz,1,iq) = ltdens(ix,iy,iz,1,iq)+2.d0*weight &
          *(  REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,2,iq) = ltdens(ix,iy,iz,2,iq)+2.d0*weight &
          *(  AIMAG(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1)))  )
     ltdens(ix,iy,iz,3,iq) = ltdens(ix,iy,iz,3,iq)+weight &
          *(  REAL(ps1(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,1))) &
          -   REAL(ps1(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,2)))  )
     END FORALL
     ENDIF

     ! z components J_zj
     IF(jfon .OR. j2on) THEN
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     lscurrentz(ix,iy,iz,1,iq) = lscurrentz(ix,iy,iz,1,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1))) &
          +  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)))  )
     lscurrentz(ix,iy,iz,2,iq) = lscurrentz(ix,iy,iz,2,iq)+weight* &
          (  REAL(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2))) &
          -  REAL(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)))  )
     lscurrentz(ix,iy,iz,3,iq) = lscurrentz(ix,iy,iz,3,iq)+weight* &
          (  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
          -  AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)))  )
     END FORALL
     ENDIF
    !***********************************************************************
    ! mixed-derivatives
    !***********************************************************************
     IF(sfon) THEN
     IF(TFFT) THEN
       CALL cdervx(psin,psx)  
       CALL cdervy(psin,psy)  
       CALL cdervz(psin,psz)  
     ELSE
       CALL cmulx(der1x,psin,psx,0)  
       CALL cmuly(der1y,psin,psy,0)  
       CALL cmulz(der1z,psin,psz,0)  
     ENDIF
     FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
     ! F_x components
     lfdens(ix,iy,iz,1,iq) = lfdens(ix,iy,iz,1,iq)+2.0d0*weight &
          *(  REAL(psx(ix,iy,iz,2)*CONJG(psx(ix,iy,iz,1)))  )
     lfdens(ix,iy,iz,1,iq) = lfdens(ix,iy,iz,1,iq)-weight &
          *(  AIMAG(psx(ix,iy,iz,1)*CONJG(psy(ix,iy,iz,2))) &
          +   AIMAG(psy(ix,iy,iz,1)*CONJG(psx(ix,iy,iz,2)))  )
     lfdens(ix,iy,iz,1,iq) = lfdens(ix,iy,iz,1,iq)+weight &
          *(  REAL(psx(ix,iy,iz,1)*CONJG(psz(ix,iy,iz,1))) &
          -   REAL(psx(ix,iy,iz,2)*CONJG(psz(ix,iy,iz,2)))  )
     ! F_y components 
     lfdens(ix,iy,iz,2,iq) = lfdens(ix,iy,iz,2,iq)+weight &
          *(  REAL(psy(ix,iy,iz,2)*CONJG(psx(ix,iy,iz,1))) &
          +   REAL(psy(ix,iy,iz,1)*CONJG(psx(ix,iy,iz,2)))  )
     lfdens(ix,iy,iz,2,iq) = lfdens(ix,iy,iz,2,iq)-2.0d0*weight &
          *(  AIMAG(psy(ix,iy,iz,1)*CONJG(psy(ix,iy,iz,2))) )
     lfdens(ix,iy,iz,2,iq) = lfdens(ix,iy,iz,2,iq)+weight &
          *(  REAL(psy(ix,iy,iz,1)*CONJG(psz(ix,iy,iz,1))) &
          -   REAL(psy(ix,iy,iz,2)*CONJG(psz(ix,iy,iz,2)))  )
     ! F_z components 
     lfdens(ix,iy,iz,3,iq) = lfdens(ix,iy,iz,3,iq)+weight &
          *(  REAL(psz(ix,iy,iz,2)*CONJG(psx(ix,iy,iz,1))) &
           +   REAL(psz(ix,iy,iz,1)*CONJG(psx(ix,iy,iz,2)))  )
     lfdens(ix,iy,iz,3,iq) = lfdens(ix,iy,iz,3,iq)-weight &
          *(  AIMAG(psz(ix,iy,iz,1)*CONJG(psy(ix,iy,iz,2))) &
          +   AIMAG(psy(ix,iy,iz,1)*CONJG(psz(ix,iy,iz,2)))  )
     lfdens(ix,iy,iz,3,iq) = lfdens(ix,iy,iz,3,iq)+weight &
          *(  REAL(psz(ix,iy,iz,1)*CONJG(psz(ix,iy,iz,1))) &
          -   REAL(psz(ix,iy,iz,2)*CONJG(psz(ix,iy,iz,2)))  )
    END FORALL
    ENDIF
    DEALLOCATE(ps1)
  END SUBROUTINE add_density
END MODULE Densities
