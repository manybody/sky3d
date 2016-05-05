Module Meanfield
  USE Params, ONLY: db,tcoul
  USE Densities
  USE Forces 
  USE Grids, ONLY: nx,ny,nz,der1x,der2x,der1y,der2y,der1z,der2z
  USE Coulomb, ONLY: poisson,wcoul
  IMPLICIT NONE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: upot,bmass,divaq
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: aq,spot,wlspot,dbmass
  PRIVATE :: divaq,aq,wlspot,dbmass
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_fields
    ALLOCATE(upot(nx,ny,nz,2),bmass(nx,ny,nz,2),divaq(nx,ny,nz,2), &
         aq(nx,ny,nz,3,2),spot(nx,ny,nz,3,2),wlspot(nx,ny,nz,3,2), &
         dbmass(nx,ny,nz,3,2))
    upot=0.D0
    bmass=0.D0
    wlspot=0.D0
    aq=0.D0
    divaq=0.D0
    dbmass=0.D0
  END SUBROUTINE alloc_fields
  !***********************************************************************
  SUBROUTINE skyrme 
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    REAL(db),PARAMETER :: epsilon=1.0d-25  
    REAL(db) :: rotspp,rotspn
    REAL(db),ALLOCATABLE :: workden(:,:,:,:),workvec(:,:,:,:,:)
    INTEGER :: ix,iy,iz,ic,iq,icomp
    ALLOCATE(workden(nx,ny,nz,2),workvec(nx,ny,nz,3,2))
    !  Step 1: 3-body contribution to upot.
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=(rho(:,:,:,1)+rho(:,:,:,2))**f%power * &
            ((b3*(f%power+2.D0)/3.D0-2.D0*b3p/3.D0)*rho(:,:,:,iq) &
            +b3*(f%power+2.D0)/3.D0*rho(:,:,:,ic) &
            -(b3p*f%power/3.D0)*(rho(:,:,:,1)**2+rho(:,:,:,2)**2)/ &
            (rho(:,:,:,1)+rho(:,:,:,2)+epsilon))
    ENDDO
    ! Step 2: add divergence of spin-orbit current to upot
    DO iq=1,2
       CALL rmulx(der1x,sodens(:,:,:,1,iq),workden(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),workden(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),workden(:,:,:,iq),1)
    ENDDO
    DO iq=1,2  
       ic=3-iq 
       upot(:,:,:,iq)=upot(:,:,:,iq) &
            -(b4+b4p)*workden(:,:,:,iq)-b4*workden(:,:,:,ic)
    ENDDO
    ! Step 3: Coulomb potential
    IF(tcoul) THEN
       CALL poisson
       upot(:,:,:,2)=upot(:,:,:,2)+wcoul
       IF(f%ex/=0) &
            upot(:,:,:,2)=upot(:,:,:,2)-slate*rho(:,:,:,2)**(1.0D0/3.0D0)
    ENDIF
    ! Step 4: remaining terms of upot
    DO iq=1,2
       CALL rmulx(der2x,rho(:,:,:,iq),workden(:,:,:,iq),0)  
       CALL rmuly(der2y,rho(:,:,:,iq),workden(:,:,:,iq),1)  
       CALL rmulz(der2z,rho(:,:,:,iq),workden(:,:,:,iq),1)
    ENDDO
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=upot(:,:,:,iq)+(b0-b0p)*rho(:,:,:,iq)+b0*rho(:,:,:,ic) &
                                ! t1,t2, and tau-dependent part      !
            +(b1-b1p)*tau(:,:,:,iq)+b1*tau(:,:,:,ic) &
                                ! two-body laplacian*rho-dependent part
            -(b2-b2p)*workden(:,:,:,iq)-b2*workden(:,:,:,ic)
       ! Step 5: effective mass
       bmass(:,:,:,iq)=f%h2m(iq)+(b1-b1p)*rho(:,:,:,iq)+b1*rho(:,:,:,ic)
       ! Step 6: calculate grad(rho) and wlspot
       CALL rmulx(der1x,rho(:,:,:,iq),workvec(:,:,:,1,iq),0)
       CALL rmuly(der1y,rho(:,:,:,iq),workvec(:,:,:,2,iq),0)
       CALL rmulz(der1z,rho(:,:,:,iq),workvec(:,:,:,3,iq),0)
    ENDDO
    DO iq=1,2
       ic=3-iq
       wlspot(:,:,:,:,iq)= &
            (b4+b4p)*workvec(:,:,:,:,iq)+b4*workvec(:,:,:,:,ic)
    END DO
    ! Step 7: calculate curl of spin density vector, store in workvec
    DO iq=1,2  
       CALL rmuly(der1y,sdens(:,:,:,3,iq),workvec(:,:,:,1,iq),0)
       CALL rmulz(der1z,sdens(:,:,:,2,iq),workvec(:,:,:,1,iq),-1)
       CALL rmulz(der1z,sdens(:,:,:,1,iq),workvec(:,:,:,2,iq),0)
       CALL rmulx(der1x,sdens(:,:,:,3,iq),workvec(:,:,:,2,iq),-1)
       CALL rmulx(der1x,sdens(:,:,:,2,iq),workvec(:,:,:,3,iq),0)
       CALL rmuly(der1y,sdens(:,:,:,1,iq),workvec(:,:,:,3,iq),-1)
    ENDDO
    ! Step 8: calculate A_q vector
    DO iq=1,2
       ic=3-iq
       aq(:,:,:,:,iq)=-2.0D0*(b1-b1p)*current(:,:,:,:,iq) &
            -2.0D0*b1*current(:,:,:,:,ic) &
            -(b4+b4p)*workvec(:,:,:,:,iq)-b4*workvec(:,:,:,:,ic)
    ENDDO
    ! Step 9: calculate the curl of the current density, stopr in spot
    DO iq=1,2  
       CALL rmuly(der1y,current(:,:,:,3,iq),spot(:,:,:,1,iq),0)
       CALL rmulz(der1z,current(:,:,:,2,iq),spot(:,:,:,1,iq),-1)
       CALL rmulz(der1z,current(:,:,:,1,iq),spot(:,:,:,2,iq),0)
       CALL rmulx(der1x,current(:,:,:,3,iq),spot(:,:,:,2,iq),-1)
       CALL rmulx(der1x,current(:,:,:,2,iq),spot(:,:,:,3,iq),0)
       CALL rmuly(der1y,current(:,:,:,1,iq),spot(:,:,:,3,iq),-1)
    ENDDO
    ! Step 10: combine isospin contributions
    DO icomp=1,3  
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx  
                rotspp=spot(ix,iy,iz,icomp,1)  
                rotspn=spot(ix,iy,iz,icomp,2)  
                spot(ix,iy,iz,icomp,1)=-(b4+b4p)*rotspp-b4*rotspn  
                spot(ix,iy,iz,icomp,2)=-(b4+b4p)*rotspn-b4*rotspp  
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! Step 11: calculate divergence of aq in divaq 
    DO iq=1,2  
       CALL rmulx(der1x,aq(:,:,:,1,iq),divaq(:,:,:,iq),0)
       CALL rmuly(der1y,aq(:,:,:,2,iq),divaq(:,:,:,iq),1)
       CALL rmulz(der1z,aq(:,:,:,3,iq),divaq(:,:,:,iq),1)
    ENDDO
    ! Step 12: calculate the gradient of the effective mass in dbmass
    DO iq=1,2  
       CALL rmulx(der1x,bmass(:,:,:,iq),dbmass(:,:,:,1,iq),0)
       CALL rmuly(der1y,bmass(:,:,:,iq),dbmass(:,:,:,2,iq),0)
       CALL rmulz(der1z,bmass(:,:,:,iq),dbmass(:,:,:,3,iq),0)
    ENDDO
    DEALLOCATE(workden,workvec)
  END SUBROUTINE skyrme
  !***********************************************************************
  SUBROUTINE hpsi(iq,eshift,pinn,pout)
    USE Trivial, ONLY: cmulx, cmuly, cmulz
    USE Levels, ONLY: cdervx,cdervy,cdervz
    INTEGER :: iq
    REAL(db) :: eshift  
    COMPLEX(db),DIMENSION(:,:,:,:) :: pinn,pout
    INTENT(IN) :: iq,eshift
    INTENT(INOUT) :: pinn
    INTENT(OUT) :: pout
    INTEGER :: is,ic
    REAL(db) :: sigis
    COMPLEX(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: pswk,pswk2
    ALLOCATE(pswk(nx,ny,nz,2),pswk2(nx,ny,nz,2))
    ! Step 1: non-derivative parts not involving spin
    DO is=1,2  
       pout(:,:,:,is)=CMPLX(upot(:,:,:,iq)-eshift, &
            -0D0,db)*pinn(:,:,:,is)
    ENDDO
    ! Step 2: the spin-current coupling
    pout(:,:,:,1)=pout(:,:,:,1)  &
         + CMPLX(spot(:,:,:,1,iq),-spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,2)  + spot(:,:,:,3,iq)*pinn(:,:,:,1)
    pout(:,:,:,2)=pout(:,:,:,2) &
         + CMPLX(spot(:,:,:,1,iq),spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,1) - spot(:,:,:,3,iq)*pinn(:,:,:,2)
    ! Step 3: derivative terms in x
    IF(TFFT) THEN
       CALL cdervx(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmulx(der1x,pinn,pswk,0)  
       CALL cmulx(der2x,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,1,iq),0.5D0*aq(:,:,:,1,iq) &
            -sigis*wlspot(:,:,:,2,iq),db)*pswk(:,:,:,is)  &
            -sigis*wlspot(:,:,:,3,iq)*pswk(:,:,:,ic) &
            -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)-wlspot(:,:,:,2,iq))*pinn(:,:,:,1)&
         -0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)+wlspot(:,:,:,2,iq))*pinn(:,:,:,2)&
         +0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervx(pswk2,pswk)  
    ELSE
       CALL cmulx(der1x,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    ! Step 4: derivative terms in y
    IF(TFFT) THEN
       CALL cdervy(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmuly(der1y,pinn,pswk,0)  
       CALL cmuly(der2y,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,2,iq),0.5D0*aq(:,:,:,2,iq) &
            +sigis*wlspot(:,:,:,1,iq),db)*pswk(:,:,:,is) &
            +CMPLX(0.D0,0.5D0*wlspot(:,:,:,3,iq),db)*pswk(:,:,:,ic) &
            -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)+wlspot(:,:,:,1,iq))*pinn(:,:,:,1)&
         +CMPLX(0D0,0.5D0,db)*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)-wlspot(:,:,:,1,iq))*pinn(:,:,:,2)&
         +CMPLX(0D0,0.5D0*wlspot(:,:,:,3,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervy(pswk2,pswk)  
    ELSE
       CALL cmuly(der1y,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    ! Step 5: derivative terms in z
    IF(TFFT) THEN
       CALL cdervz(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmulz(der1z,pinn,pswk,0)  
       CALL cmulz(der2z,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,3,iq),0.5D0*aq(:,:,:,3,iq),db)*pswk(:,:,:,is) &
            +CMPLX(sigis*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)* &
            pswk(:,:,:,ic)-bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,1)&
         +CMPLX(0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,2)&
         +CMPLX(-0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervz(pswk2,pswk)  
    ELSE
       CALL cmulz(der1z,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    !
    DEALLOCATE(pswk,pswk2)
  END SUBROUTINE hpsi
  !***********************************************************************
END Module Meanfield
