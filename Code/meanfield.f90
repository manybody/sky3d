Module Meanfield
  USE Params, ONLY: db,tcoul,wflag,tcoul,tdynamic,todd,toddls,e2,pi
  USE Densities
  USE Forces 
  USE Grids, ONLY: nx,ny,nz,der1x,der2x,der1y,der2y,der1z,der2z,wxyz
  USE Coulomb, ONLY: poisson,wcoul
  IMPLICIT NONE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: upot,bmass,bmass2,dxiq,workden
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) ::  xiq,spot,wlspot,dbmass, &
       laps,ssig,dssigx,dssigy,dssigz, &
       socx,socy,socz,soc,soc2x,soc2y,soc2z,soc0, &
       sfden,dsfx,dsfy,dsfz,divs,a
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_fields
    ALLOCATE(upot(nx,ny,nz,2),bmass(nx,ny,nz,2),bmass2(nx,ny,nz,2), &
         dxiq(nx,ny,nz,2),xiq(nx,ny,nz,3,2),&
         spot(nx,ny,nz,3,2),wlspot(nx,ny,nz,3,2), &
         dbmass(nx,ny,nz,3,2),laps(nx,ny,nz,3,2),  &
         ssig(nx,ny,nz,3,2),dssigx(nx,ny,nz,3,2),dssigy(nx,ny,nz,3,2), &
         dssigz(nx,ny,nz,3,2),socx(nx,ny,nz,3,2), socy(nx,ny,nz,3,2), &
         socz(nx,ny,nz,3,2),soc(nx,ny,nz,3,2),soc0(nx,ny,nz,3,2), &
         soc2x(nx,ny,nz,3,2),soc2y(nx,ny,nz,3,2),soc2z(nx,ny,nz,3,2), &
         sfden(nx,ny,nz,3,2),dsfx(nx,ny,nz,3,2), &
         dsfy(nx,ny,nz,3,2),dsfz(nx,ny,nz,3,2),divs(nx,ny,nz,3,2),&
         a(nx,ny,nz,3,2))
         a=0.0d0
         upot=0.0d0
         bmass=0.0d0 
         bmass2=0.0d0 
         dxiq=0.0d0
         xiq=0.0d0
         wlspot=0.0d0
         dbmass=0.0d0
         laps=0.0d0
         ssig=0.0d0
         dssigx=0.0d0
         dssigy=0.0d0
         dssigz=0.0d0
         socx=0.0d0
         socy=0.0d0
         socz=0.0d0
         soc=0.0d0
         soc0=0.0d0
         soc2x=0.0d0
         soc2y=0.0d0
         soc2z=0.0d0
         sfden=0.0d0
         dsfx=0.0d0
         dsfy=0.0d0
         dsfz=0.0d0
         divs=0.0d0
  END SUBROUTINE alloc_fields

SUBROUTINE skyrme
  USE Params
  USE Densities
  USE Forces
  USE Grids
  USE Coulomb
  USE Trivial
  IMPLICIT NONE
  !
  !***********************************************************************
  !                                                                      *
  !        calculates the most general form of the skyrme force with     *
  !        parameters t0...x3+coulomb+spin-orbit                         *
  !        ayuk=0.0 uses the zero-range form of the skyrme force         *
  !                                                                      *
  !***********************************************************************
  !
  REAL(db),PARAMETER :: epsilon = 1.0d-25  
  LOGICAL,PARAMETER :: tspint = .TRUE.  
  REAL(db) :: spino(2,2),t0c,t0a,ttb,tta,tmb,tma,tdb,tda, &
       t3aa,t3ba,t3a,t3b,t4h,rotspp,rotspn,tsa,tsb, &
       t3sa,t3sb,t3saa,t3sbb,tlsaa,tlsbb,ttka,ttkb,tja,tjb, &
       tjc,tjd,tj2a,tj2b,tlscc,tdivsa,b4q
  INTEGER :: ix,iy,iz,ic,iq,icomp,i
  slate=(3.0d0/pi)**(1.0d0/3.0d0)*e2   !
  wlspot = 0.0d0
  IF(.NOT.todd.AND.toddls) STOP ' todd and toddls need to be the same '
  !
  !***********************************************************************
  !                                                                      *
  !        prepare appropriate combinations of Skyrme parameters         *
  !                                                                      *
  !***********************************************************************
  t0c = f%t0*(1.0d0+0.5d0*f%x0)  
  t0a = t0c-f%t0*(0.5d0+f%x0) 
  tsb = (f%t0*f%x0)/2.0d0  
  tsa = tsb-f%t0/2.0d0   
  ttb =(f%t1+0.5d0*f%x1*f%t1+f%t2+0.5*f%x2*f%t2)/4.0d0  
  tta = ttb-(f%t1+2.0d0*f%x1*f%t1-f%t2-2.0d0*f%x2*f%t2)/8.0d0  
  tmb = ttb  
  tma = tta  
  tdb =-(3.0d0*f%t1+1.5d0*f%x1*f%t1-f%t2-0.5d0*f%x2*f%t2)/8.0d0
  tda = tdb +(3.0d0*f%t1+6.0d0*f%x1*f%t1+f%t2+2.0d0*f%x2*f%t2)/16.0d0
  tlsbb = (f%t2*f%x2-3.0d0*f%t1*f%x1)/16.0d0 
  tlsaa = tlsbb + (3.0d0*f%t1+f%t2)/16.0d0 
  tlscc = (3.0d0*te-toten)/8.0d0
  t3aa = f%t3*f%power/12.0d0 *(1.0d0+0.5d0*f%x3)  
  t3ba = f%t3*f%power/12.0d0 *(0.5d0+f%x3)  
  t3b = f%t3 *(1.0d0+0.5d0*f%x3)/6.0d0  
  t3a = t3b-f%t3 *(0.5d0+f%x3)/6.0d0  
  t3sa = f%t3*f%x3*f%power/24.0d0  
  t3sb = -f%t3*f%power/24.0d0  
  t3saa = f%t3*(f%x3-1.0d0)/12.0d0  
  t3sbb = f%t3*f%x3/12.0d0 
  ttkb = (f%t1*f%x1+f%t2*f%x2)/8.0d0
  ttka = ttkb + (f%t2-f%t1)/8.0d0 
  tdivsa = -3.0d0*(3.0d0*te-toten)/8.0d0 
  tjb = -(f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/4.0d0
  tja = (f%t1-f%t1*f%x1-f%t2-f%t2*f%x2+4.0d0*toten)/4.0d0
  tj2a = -3.0d0*toten/2.0d0
  tj2b = -3.0d0*(te+toten)/4.0d0
  t4h = f%t4/2.0d0  
  b4 = t4h  
  b4q = b4+f%b4p  
  upot=0.0d0
  !***********************************************************************
  !        three-body term                                               *
  !        this step uses dxiq as workspace                              *
  !***********************************************************************
  FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
     dxiq(ix,iy,iz,1) = (rho(ix,iy,iz,1)+rho(ix,iy,iz,2))**f%power  
     upot(ix,iy,iz,1) = dxiq(ix,iy,iz,1)*&
          ( t3aa*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2)) &
          -t3ba*(rho(ix,iy,iz,1)**2+rho(ix,iy,iz,2)**2)/ &
          ((rho(ix,iy,iz,1)+rho(ix,iy,iz,2))+epsilon)    )
     upot(ix,iy,iz,2) = upot(ix,iy,iz,1)  
  END FORALL
  !
  DO iq = 1,2  
     ic = 3-iq  
     upot(:,:,:,iq)=upot(:,:,:,iq)+t3a*dxiq(:,:,:,1)*rho(:,:,:,iq)  
     upot(:,:,:,iq)=upot(:,:,:,iq)+t3b*dxiq(:,:,:,1)*rho(:,:,:,ic)  
  ENDDO
  !***********************************************************************
  !        S^2 parts of three-body term                                  *
  !***********************************************************************
  IF(todd .and. s2on) THEN
    dxiq=0.0d0
    DO icomp = 1,3
      dxiq(:,:,:,1)=dxiq(:,:,:,1) + &
          (sdens(:,:,:,icomp,1)+sdens(:,:,:,icomp,2))**2
      dxiq(:,:,:,2)=dxiq(:,:,:,2) + &
          sdens(:,:,:,icomp,1)**2+sdens(:,:,:,icomp,2)**2
    ENDDO
    DO iq = 1,2
      upot(:,:,:,iq) = upot(:,:,:,iq) + &
          (rho(:,:,:,iq)+rho(:,:,:,3-iq))**f%power* &
          ( t3sa*dxiq(:,:,:,1) + t3sb*dxiq(:,:,:,2))/ &
          ((rho(:,:,:,iq)+rho(:,:,:,3-iq))+epsilon)
    ENDDO    
  ENDIF
  !***********************************************************************
  !        begin calculation of the spin-orbit potential                 *
  !        dbmass used as workspace.                                     *
  !                                                                      *
  !        gradient of spin-orbit current in dbmass                      *
  !***********************************************************************
  DO iq = 1,2
     CALL rmulx(der1x,sodens(:,:,:,1,iq),dbmass(:,:,:,1,iq),0)
     CALL rmuly(der1y,sodens(:,:,:,2,iq),dbmass(:,:,:,2,iq),0)
     CALL rmulz(der1z,sodens(:,:,:,3,iq),dbmass(:,:,:,3,iq),0)
  ENDDO
  !***********************************************************************
  !        add the result to the local hf potential                      *
  !        we sum over all cartesian components to form the dot product  *
  !***********************************************************************
  DO iq = 1,2  
     ic = 3-iq  
     DO icomp = 1,3  
        upot(:,:,:,iq)=upot(:,:,:,iq) &
             -b4q*dbmass(:,:,:,icomp,iq) -b4*dbmass(:,:,:,icomp,ic)
     ENDDO
  ENDDO
  !***********************************************************************
  !        store the gradient of rho in dbmass                           *
  !        these will be used in constructing the bq.(del x sigma) term  *
  !***********************************************************************
  DO iq = 1,2  
     CALL rmulx(der1x,rho(:,:,:,iq),dbmass(:,:,:,1,iq),0)
     CALL rmuly(der1y,rho(:,:,:,iq),dbmass(:,:,:,2,iq),0)
     CALL rmulz(der1z,rho(:,:,:,iq),dbmass(:,:,:,3,iq),0)
  ENDDO
  !***********************************************************************
  !        solve poisson equation for the coulomb potential              *
  !        exchange part is the slater approximation                     *
  !***********************************************************************
  IF(tcoul .eqv. .TRUE.) THEN
     CALL poisson
        bmass(:,:,:,1) = wcoul-slate*rho(:,:,:,2)**(1.0d0/3.0d0)
     upot(:,:,:,2)=upot(:,:,:,2)+bmass(:,:,:,1)
  ENDIF
  !***********************************************************************
  !        calculate the laplacian of rho(store in dxiq)                 *
  !***********************************************************************
  DO iq=1,2
     CALL rmulx(der2x,rho(:,:,:,iq),dxiq(:,:,:,iq),0)  
     CALL rmuly(der2y,rho(:,:,:,iq),dxiq(:,:,:,iq),1)  
     CALL rmulz(der2z,rho(:,:,:,iq),dxiq(:,:,:,iq),1)
  ENDDO
  !***********************************************************************
  !        add all two-body terms                                        *
  !***********************************************************************
  xiq = 0.0d0
  DO iq = 1,2  
     ic = 3-iq  
     !***********************************************************************
     !        t0-dependent part                                             *
     !***********************************************************************
     upot(:,:,:,iq)=upot(:,:,:,iq)+t0a*rho(:,:,:,iq)+t0c*rho(:,:,:,ic)
     !***********************************************************************
     !        t1,t2 and tau-dependent part                                  *
     !***********************************************************************
     upot(:,:,:,iq)=upot(:,:,:,iq)+tta*tau(:,:,:,iq)+ttb*tau(:,:,:,ic)
     !***********************************************************************
     !        two-body laplacian*rho-dependent part                         *
     !***********************************************************************
     upot(:,:,:,iq)=upot(:,:,:,iq)+tda*dxiq(:,:,:,iq)+tdb*dxiq(:,:,:,ic)
     !***********************************************************************
     !        begin storing the effective mass in bmass                     *
     !        initialize bmass with h2m                                     *
     !***********************************************************************
     bmass(:,:,:,iq) = f%h2m(iq)  
     !***********************************************************************
     !        add rho dependence to bmass                                   *
     !***********************************************************************
     bmass(:,:,:,iq)=bmass(:,:,:,iq)+tma*rho(:,:,:,iq)+tmb*rho(:,:,:,ic)
     !***********************************************************************
     !        begin constructing the bq-vector                              *
     !        these will be used in constructing the bq.(del x sigma) term  *
     !        in the action of h on psi                                     *
     !        they are stored in wlspot and used in routine hpsi            *
     !        ****note**** we turned off the non-invariant term(galilean)   *
     !***********************************************************************
     !
     DO icomp = 1,3  
        wlspot(:,:,:,icomp,iq)=wlspot(:,:,:,icomp,iq)+ &
             b4q*dbmass(:,:,:,icomp,iq)+b4*dbmass(:,:,:,icomp,ic)
     ENDDO
     !***********************************************************************
     !        begin constructing the iq-vector                              *
     !        these will be used in constructing the(del.iq+iq.del) term   *
     !        they are stored in xiq                                        *
     !***********************************************************************
     IF(todd) THEN
        DO icomp = 1,3  
           !        xiq(:,:,:,icomp,iq)=xiq(:,:,:,icomp,iq)+ &
           !             tta*current(:,:,:,icomp,iq)+ttb*current(:,:,:,icomp,ic)
           xiq(:,:,:,icomp,iq)=xiq(:,:,:,icomp,iq) - &
                (2.0D0*tta)*current(:,:,:,icomp,iq)-(2.0D0*ttb)*current(:,:,:,icomp,ic)
        ENDDO
     ENDIF
     !
  END DO
  !***********************************************************************
  !        optionally the time-odd parts of the l*s force                *
  !***********************************************************************
  IF(toddls) THEN  
     !                              the rotation of spin density
     !                              on 'spot' as workspace
     DO iq = 1,2  
        CALL rmuly(der1y,sdens(:,:,:,3,iq),spot(:,:,:,1,iq),0)
        CALL rmulz(der1z,sdens(:,:,:,2,iq),spot(:,:,:,1,iq),- 1)
        CALL rmulz(der1z,sdens(:,:,:,1,iq),spot(:,:,:,2,iq),0)
        CALL rmulx(der1x,sdens(:,:,:,3,iq),spot(:,:,:,2,iq),- 1)
        CALL rmulx(der1x,sdens(:,:,:,2,iq),spot(:,:,:,3,iq),0)
        CALL rmuly(der1y,sdens(:,:,:,1,iq),spot(:,:,:,3,iq),- 1)
     ENDDO
     !                                          add to time-odd potential 'xi
     DO iq = 1,2  
        xiq(:,:,:,:,iq)=xiq(:,:,:,:,iq)-b4q*spot(:,:,:,:,iq) &
             -b4*spot(:,:,:,:,3-iq)
     ENDDO
     IF(tspint) THEN  
        spino(1,1) = wxyz*SUM(spot*current)
        spino(2,1) = spino(1,1)+wxyz*&
             SUM(spot(:,:,:,:,1)*current(:,:,:,:,2) &
             +spot(:,:,:,:,2)*current(:,:,:,:,1))
     ENDIF
     !                        the rotation of current
     !                        on 'spot' as workspace
     DO iq = 1,2  
        CALL rmuly(der1y,current(:,:,:,3,iq),spot(:,:,:,1,iq),0)
        CALL rmulz(der1z,current(:,:,:,2,iq),spot(:,:,:,1,iq),-1)
        CALL rmulz(der1z,current(:,:,:,1,iq),spot(:,:,:,2,iq),0)
        CALL rmulx(der1x,current(:,:,:,3,iq),spot(:,:,:,2,iq),-1)
        CALL rmulx(der1x,current(:,:,:,2,iq),spot(:,:,:,3,iq),0)
        CALL rmuly(der1y,current(:,:,:,1,iq),spot(:,:,:,3,iq),-1)
     ENDDO
     IF(tspint) THEN  
        spino(1,2) = wxyz*SUM(spot*sdens)
        spino(2,2) = spino(1,2)+wxyz* &
             SUM(spot(:,:,:,:,1)*sdens(:,:,:,:,2) &
             +spot(:,:,:,:,2)*sdens(:,:,:,:,1))
        IF(wflag.AND.tdynamic)  &
             WRITE(6,'(a,2g13.5)') ' dyn. l*s energies=',  &
             -b4*spino(2,1)-f%b4p*spino(1,1),-b4*spino(2,2)-f%b4p*spino(1,2)
     ENDIF
     !                                           recouple to spin potential
     !?  could one do that with a FORALL ?
     DO icomp = 1,3  
        DO iz = 1,nz  
           DO iy = 1,ny  
              DO ix = 1,nx  
                 rotspp = spot(ix,iy,iz,icomp,1)  
                 rotspn = spot(ix,iy,iz,icomp,2)  
                 spot(ix,iy,iz,icomp,1) =-b4q*rotspp-b4*rotspn  
                 spot(ix,iy,iz,icomp,2) =-b4q*rotspn-b4*rotspp  
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
  ENDIF
  !
  IF (todd) THEN

  !***********************************************************************
  !       time odd S components in t0 term                               *
  !***********************************************************************
  IF (s2on) THEN
    DO iq = 1,2
        spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
           tsa*sdens(:,:,:,:,iq) + tsb*sdens(:,:,:,:,3-iq)
    ENDDO
  ENDIF

  !***********************************************************************
  !       time odd S components in t3 term                               *
  !***********************************************************************
  IF (s2on) THEN
    DO iq = 1,2
      DO icomp = 1,3
         spot(:,:,:,icomp,iq) = spot(:,:,:,icomp,iq) + &
             (rho(:,:,:,1)+rho(:,:,:,2))**f%power* &
             (t3saa*sdens(:,:,:,icomp,iq) + t3sbb*sdens(:,:,:,icomp,3-iq))
      ENDDO
    ENDDO
  ENDIF
  !***********************************************************************
  !        calculate the laplacian of sdens(store in laps)               *
  !***********************************************************************
  DO iq = 1,2
    DO icomp = 1,3
      CALL rmulx(der2x,sdens(:,:,:,icomp,iq),laps(:,:,:,icomp,iq),0)  
      CALL rmuly(der2y,sdens(:,:,:,icomp,iq),laps(:,:,:,icomp,iq),1)  
      CALL rmulz(der2z,sdens(:,:,:,icomp,iq),laps(:,:,:,icomp,iq),1)
    ENDDO
  ENDDO
  !
  IF (lapson) THEN
    DO iq = 1,2
      spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
        (tlsaa -toten/4.0d0)*laps(:,:,:,:,iq)+(tlsbb+tlscc)*laps(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !       time odd components T from S.T                                 *
  !***********************************************************************
  IF (ston) THEN
    DO iq = 1,2
        spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
               (ttka-0.50d0*toten)*tdens(:,:,:,:,iq) + &
               (ttkb-(te+toten)/4.0d0)*tdens(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !       time odd components S from S.T                                 *
  !***********************************************************************
  IF (ston) THEN
    DO iq = 1,2
        ssig(:,:,:,:,iq) = (ttka-0.5d0*toten)*sdens(:,:,:,:,iq) + &
                           (ttkb-(te+toten)/4.0d0)*sdens(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !       time odd components S from S.T div S                           *
  !***********************************************************************
  !
 ! IF (ston) THEN
  !  DO iq  = 1,2
   !   DO icomp = 1,3
    !    CALL rmulx(der1x,ssig(:,:,:,icomp,iq),dssigx(:,:,:,icomp,iq),0)
 !       CALL rmuly(der1y,ssig(:,:,:,icomp,iq),dssigy(:,:,:,icomp,iq),0)
  !      CALL rmulz(der1z,ssig(:,:,:,icomp,iq),dssigz(:,:,:,icomp,iq),0)
  !    ENDDO
  !  ENDDO
  !ENDIF
  !***********************************************************************
  !       time odd components F from S.F                                 *
  !***********************************************************************
  IF (sfon) THEN
    DO iq = 1,2
        spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
                  3.0d0*toten/2.0d0*fdens(:,:,:,:,iq) + &
                  3.0d0*(te+toten)/4.0d0*fdens(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !       time odd components F from S.F                                 *
  !***********************************************************************
  IF (sfon) THEN
    DO iq = 1,2
        sfden(:,:,:,:,iq) = 3.0d0*toten/2.0d0*sdens(:,:,:,:,iq) + &
                  3.0d0*(te+toten)/4.0d0*sdens(:,:,:,:,3-iq)
    ENDDO
  !
    DO iq  = 1,2
      DO icomp = 1,3
        CALL rmulx(der1x,sfden(:,:,:,icomp,iq),dsfx(:,:,:,icomp,iq),0)
        CALL rmuly(der1y,sfden(:,:,:,icomp,iq),dsfy(:,:,:,icomp,iq),0)
        CALL rmulz(der1z,sfden(:,:,:,icomp,iq),dsfz(:,:,:,icomp,iq),0)
      ENDDO
    ENDDO
  ENDIF
  !***********************************************************************
  !        nabla_i(nabla_i S_i) part of divS  (store in laps)            *
  !***********************************************************************
  DO iq = 1,2
    CALL rmulx(der2x,sdens(:,:,:,1,iq),laps(:,:,:,1,iq),0)  
    CALL rmuly(der2y,sdens(:,:,:,2,iq),laps(:,:,:,2,iq),0)  
    CALL rmulz(der2z,sdens(:,:,:,3,iq),laps(:,:,:,3,iq),0)  
  ENDDO
  !
  IF (divson) THEN
    DO iq = 1,2
        spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
          3.0d0*toten/4.0d0*laps(:,:,:,:,iq) + tdivsa*laps(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !        nabla_j(nabla_i S_j) part of divS  (store in laps)            *
  !***********************************************************************
  DO iq = 1,2
    CALL rmulx(der1x,sdens(:,:,:,1,iq)+sdens(:,:,:,3,iq),laps(:,:,:,3,iq),0) 
    CALL rmulz(der1z,laps(:,:,:,3,iq),divs(:,:,:,3,iq),0)  
    CALL rmuly(der1y,sdens(:,:,:,2,iq)+sdens(:,:,:,1,iq),laps(:,:,:,1,iq),0)  
    CALL rmulx(der1x,laps(:,:,:,1,iq),divs(:,:,:,1,iq),0)  
    CALL rmulz(der1z,sdens(:,:,:,3,iq)+sdens(:,:,:,2,iq),laps(:,:,:,2,iq),0)   
    CALL rmuly(der1y,laps(:,:,:,2,iq),divs(:,:,:,2,iq),0)  
  ENDDO
  !
  IF (divson) THEN
    DO iq = 1,2
      spot(:,:,:,:,iq) = spot(:,:,:,:,iq) + &
       3.0d0*toten/4.0d0*divs(:,:,:,:,iq) +tdivsa*divs(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !
  ENDIF
  !***********************************************************************
  !       Weights for J_ij components including tensor                   *
  !***********************************************************************
  !
  IF(jfon) THEN
    soc2x(:,:,:,1,:)=scurrentx(:,:,:,1,:)
    soc2x(:,:,:,2,:)=scurrenty(:,:,:,1,:)
    soc2x(:,:,:,3,:)=scurrentz(:,:,:,1,:)
    soc2y(:,:,:,1,:)=scurrentx(:,:,:,2,:)
    soc2y(:,:,:,2,:)=scurrenty(:,:,:,2,:)
    soc2y(:,:,:,3,:)=scurrentz(:,:,:,2,:)
    soc2z(:,:,:,1,:)=scurrentx(:,:,:,3,:)
    soc2z(:,:,:,2,:)=scurrenty(:,:,:,3,:)
    soc2z(:,:,:,3,:)=scurrentz(:,:,:,3,:)
  ENDIF
    !
  IF(j2on .OR. jfon) THEN
    IF(.NOT. j2on) THEN
      tja=0.0d0
      tjb=0.0d0
    ENDIF	
    !
    DO iq = 1,2
      socx(:,:,:,:,iq)=tja*scurrentx(:,:,:,:,iq)+tj2a*soc2x(:,:,:,:,iq)+ &
                       tjb*scurrentx(:,:,:,:,3-iq)+tj2b*soc2x(:,:,:,:,3-iq)
      socy(:,:,:,:,iq)=tja*scurrenty(:,:,:,:,iq)+tj2a*soc2y(:,:,:,:,iq)+ &
                       tjb*scurrenty(:,:,:,:,3-iq)+tj2b*soc2y(:,:,:,:,3-iq)
      socz(:,:,:,:,iq)=tja*scurrentz(:,:,:,:,iq)+tj2a*soc2z(:,:,:,:,iq)+ &
                       tjb*scurrentz(:,:,:,:,3-iq)+tj2b*soc2z(:,:,:,:,3-iq)
    ENDDO
  ENDIF
  !***********************************************************************
  !       nabla J_ij components including tensor                         *
  !***********************************************************************
  IF(j2on .OR. jfon) THEN
    DO iq = 1,2
      DO icomp = 1,3
        CALL rmulx(der1x,socx(:,:,:,icomp,iq),soc(:,:,:,icomp,iq),0)
        CALL rmuly(der1y,socy(:,:,:,icomp,iq),soc(:,:,:,icomp,iq),1)
        CALL rmulz(der1z,socz(:,:,:,icomp,iq),soc(:,:,:,icomp,iq),1)
      ENDDO
    ENDDO
  ENDIF
  !***********************************************************************
  !       Weights for J_iiJ_jj components for tensor                     *
  !***********************************************************************
  IF(jfon) THEN
    DO iq = 1,2
      soc0(:,:,:,1,iq)=tj2a*(scurrentx(:,:,:,1,iq)+scurrenty(:,:,:,2,iq)+ &
                      scurrentz(:,:,:,3,iq))+tj2b*(scurrentx(:,:,:,1,3-iq)+ &
                      scurrenty(:,:,:,2,3-iq)+scurrentz(:,:,:,3,3-iq))
    ENDDO
  !***********************************************************************
  !       nabla J_ii components including tensor                         *
  !***********************************************************************
    DO iq = 1,2
      CALL rmulx(der1x,soc0(:,:,:,1,iq),soc(:,:,:,1,iq),1)
      CALL rmuly(der1y,soc0(:,:,:,1,iq),soc(:,:,:,2,iq),1)
      CALL rmulz(der1z,soc0(:,:,:,1,iq),soc(:,:,:,3,iq),1)
    ENDDO
  ENDIF
  !***********************************************************************
  !        calculate the del.iq and store in dxiq                        *
  !        dbmass is used as work space                                  *
  !***********************************************************************
  IF(todd) THEN
     DO iq = 1,2  
        CALL rmulx(der1x,xiq(:,:,:,1,iq),dxiq(:,:,:,iq),0)
        CALL rmuly(der1y,xiq(:,:,:,2,iq),dxiq(:,:,:,iq),1)
        CALL rmulz(der1z,xiq(:,:,:,3,iq),dxiq(:,:,:,iq),1)
     ENDDO
  ENDIF
  !***********************************************************************
  !        calculate the gradient of the effective mass and store in     *
  !        dbmass. this is used as part of the kinetic energy term in    *
  !        hpsi                                                          *
  !***********************************************************************
  DO iq = 1,2  
     CALL rmulx(der1x,bmass(:,:,:,iq),dbmass(:,:,:,1,iq),0)
     CALL rmuly(der1y,bmass(:,:,:,iq),dbmass(:,:,:,2,iq),0)
     CALL rmulz(der1z,bmass(:,:,:,iq),dbmass(:,:,:,3,iq),0)
  ENDDO
END SUBROUTINE skyrme




SUBROUTINE hpsi(iq,eshift,pinn,pout)
  USE Params
  USE Grids
  USE Densities
  USE Levels
  USE Trivial, ONLY: cmulx, cmuly, cmulz
  !
  IMPLICIT NONE
  !***********************************************************************
  !        hpsi:  program to form the product of h on the                *
  !        state vector pinn and return the result in pout               *
  !***********************************************************************
  COMPLEX(db),DIMENSION(:,:,:,:) :: pinn,pout 
  REAL(db) :: eshift  
  INTENT(IN) :: eshift
  INTENT(OUT) :: pout
  INTENT(INOUT) :: pinn
  COMPLEX(db),DIMENSION(nx,ny,nz,2) :: pswk1,pswk,pswk2,pswk3,pswk4
  INTEGER :: ix,iy,iz,is,ic,iq
  REAL(db) :: sigis
  INTENT(IN) :: iq
  !***********************************************************************
  !        diagonal part of h times pinn                                 *
  !***********************************************************************
  DO is = 1,2  
     pout(:,:,:,is) = CMPLX(upot(:,:,:,iq)-eshift,&
				0.0d0,db)*pinn(:,:,:,is)
  ENDDO
  !***********************************************************************
  !        action of spin potential                                      *
  !***********************************************************************
  IF(toddls) THEN  
     FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
        pout(ix,iy,iz,1) = pout(ix,iy,iz,1)  &
             + CMPLX(spot(ix,iy,iz,1,iq),-spot(ix,iy,iz,2,iq),db)*pinn(ix,iy,iz,2) &
             + spot(ix,iy,iz,3,iq)*pinn(ix,iy,iz,1)
        pout(ix,iy,iz,2) = pout(ix,iy,iz,2) &
             + CMPLX(spot(ix,iy,iz,1,iq),spot(ix,iy,iz,2,iq),db)*pinn(ix,iy,iz,1) &
             - spot(ix,iy,iz,3,iq)*pinn(ix,iy,iz,2)
     END FORALL
  ENDIF
  !**********************************************************************
  !        J_ij Lesinski way                                            *
  !**********************************************************************
  IF(j2on .OR. jfon) THEN
    DO is = 1,2  
       sigis = 3-2*is  
       pout(:,:,:,is)=pout(:,:,:,is) &
         -0.5d0*CMPLX(sigis*soc(:,:,:,2,iq),soc(:,:,:,1,iq),db)*pinn(:,:,:,3-is) &
         -0.5d0*CMPLX(0.0d0,sigis*soc(:,:,:,3,iq),db)*pinn(:,:,:,is) 
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require d/dx psi                                  *
  !***********************************************************************
  IF(TFFT) THEN
     CALL cdervx(pinn,pswk,d2psout=pswk2)  
     CALL cdervy(pswk,pswk1)  
  ELSE
     CALL cmulx(der1x,pinn,pswk,0)  
     CALL cmulx(der2x,pinn,pswk2,0)  
     CALL cmuly(der1y,pswk,pswk1,0)  
  ENDIF
  DO is = 1,2  
     ic = 3-is  
     sigis = (3-2*is)*0.5d0  
     pout(:,:,:,is)=pout(:,:,:,is) &
          -CMPLX(dbmass(:,:,:,1,iq),0.5d0*xiq(:,:,:,1,iq)-sigis*wlspot(:,:,:,2,iq),db) &
          *pswk(:,:,:,is)  &
          -sigis*wlspot(:,:,:,3,iq)*pswk(:,:,:,ic)
  ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (xiq(:,:,:,1,iq)-wlspot(:,:,:,2,iq))*pinn(:,:,:,1)&
         -0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (xiq(:,:,:,1,iq)+wlspot(:,:,:,2,iq))*pinn(:,:,:,2)&
         +0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervx(pswk2,pswk)  
    ELSE
       CALL cmulx(der1x,pswk2,pswk,0) 
  END IF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
  IF(TFFT) THEN
     CALL cdervx(pinn,pswk,d2psout=pswk2)  
     CALL cdervy(pswk,pswk1)  
  ELSE
     CALL cmulx(der1x,pinn,pswk,0)  
     CALL cmulx(der2x,pinn,pswk2,0)  
     CALL cmuly(der1y,pswk,pswk1,0)  
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.T requiring d/dx psi                       *
  !***********************************************************************
  IF (ston) THEN
    DO is = 1,2
       sigis = 3-2*is
       pswk3(:,:,:,is)=-CMPLX(ssig(:,:,:,1,iq), &
                      -sigis*ssig(:,:,:,2,iq),db)*pswk(:,:,:,3-is) &
                      -sigis*ssig(:,:,:,3,iq)*pswk(:,:,:,is)
    ENDDO
  IF(TFFT) THEN
     CALL cdervx(pswk3,pswk4)  
  ELSE
     CALL cmulx(der1x,pswk3,pswk4,0) 
  ENDIF
     DO is = 1,2
       pout(:,:,:,is)=pout(:,:,:,is)+pswk4(:,:,:,is)
    ENDDO
  ENDIF  
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dx psi                       *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) -CMPLX(dsfx(:,:,:,1,iq)+ &
                       0.5d0*dsfy(:,:,:,2,iq)+0.5d0*dsfz(:,:,:,3,iq), &
                      -0.5d0*sigis*dsfy(:,:,:,1,iq),db)*pswk(:,:,:,3-is) &
                      -0.5d0*sigis*dsfz(:,:,:,1,iq)*pswk(:,:,:,is)
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ij on d/dx psi                                             *
  !**********************************************************************
  IF(j2on .OR. jfon) THEN
    DO is = 1,2  
       sigis = 3-2*is  
       pout(:,:,:,is)=pout(:,:,:,is) &
        -CMPLX(sigis*socx(:,:,:,2,iq),socx(:,:,:,1,iq),db)*pswk(:,:,:,3-is) &
        -CMPLX(0.0d0,sigis*socx(:,:,:,3,iq),db)*pswk(:,:,:,is) 
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ii on d/dx psi                                             *
  !**********************************************************************
  IF(jfon) THEN
    DO is = 1,2  
       pout(:,:,:,is)=pout(:,:,:,is) &
         -CMPLX(0.0d0,soc0(:,:,:,1,iq),db)*pswk(:,:,:,3-is) 
    ENDDO
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dx d/dy psi                  *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) -CMPLX(sfden(:,:,:,2,iq), &
                      -sigis*sfden(:,:,:,1,iq),db)*pswk1(:,:,:,3-is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require(d/dx)**2 psi                              *
  !***********************************************************************
  DO is=1,2
     pout(:,:,:,is)=pout(:,:,:,is) - bmass(:,:,:,iq)*pswk2(:,:,:,is)
  ENDDO
  !***********************************************************************
  !        S.F terms requiring (d/dx)**2 psi                             *
  !***********************************************************************
  IF(sfon) THEN
    DO is=1,2
       pout(:,:,:,is)=pout(:,:,:,is) - sfden(:,:,:,1,iq)*pswk2(:,:,:,3-is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require d/dy psi                                  *
  !***********************************************************************
  IF(TFFT) THEN
     CALL cdervy(pinn,pswk,d2psout=pswk2)  
     CALL cdervz(pswk,pswk1)  
  ELSE
     CALL cmuly(der1y,pinn,pswk,0)  
     CALL cmuly(der2y,pinn,pswk2,0)  
     CALL cmulz(der1z,pswk,pswk1,0)  
  ENDIF
  DO is = 1,2  
     ic = 3-is  
     sigis = 0.5d0*(3-2*is)
     pout(:,:,:,is)=pout(:,:,:,is) &
          -CMPLX(dbmass(:,:,:,2,iq),0.5d0*xiq(:,:,:,2,iq)+sigis*wlspot(:,:,:,1,iq),db) &
          *pswk(:,:,:,is)  &
          +CMPLX(0.d0,0.5d0*wlspot(:,:,:,3,iq),db)*pswk(:,:,:,ic)
  ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (xiq(:,:,:,2,iq)+wlspot(:,:,:,1,iq))*pinn(:,:,:,1)&
         +CMPLX(0D0,0.5D0,db)*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (xiq(:,:,:,2,iq)-wlspot(:,:,:,1,iq))*pinn(:,:,:,2)&
         +CMPLX(0D0,0.5D0*wlspot(:,:,:,3,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervy(pswk2,pswk)  
    ELSE
       CALL cmuly(der1y,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
  IF(TFFT) THEN
     CALL cdervy(pinn,pswk,d2psout=pswk2)  
     CALL cdervz(pswk,pswk1)  
  ELSE
     CALL cmuly(der1y,pinn,pswk,0)  
     CALL cmuly(der2y,pinn,pswk2,0)  
     CALL cmulz(der1z,pswk,pswk1,0)  
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.T requiring d/dy psi                       *
  !***********************************************************************
  IF (ston) THEN
    DO is = 1,2
       sigis = 3-2*is
       pswk3(:,:,:,is)=-CMPLX(ssig(:,:,:,1,iq), &
                      -sigis*ssig(:,:,:,2,iq),db)*pswk(:,:,:,3-is) &
                      -sigis*ssig(:,:,:,3,iq)*pswk(:,:,:,is)
    ENDDO
  IF(TFFT) THEN
     CALL cdervy(pswk3,pswk4)
  ELSE
     CALL cmuly(der1y,pswk3,pswk4,0)  
  ENDIF
     DO is = 1,2
       pout(:,:,:,is)=pout(:,:,:,is)+pswk4(:,:,:,is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dy psi                       *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) -CMPLX(0.5d0*dsfx(:,:,:,2,iq), &
                      -sigis*dsfy(:,:,:,2,iq)-0.5d0*sigis*dsfx(:,:,:,1,iq)-&
                      0.5d0*sigis*dsfz(:,:,:,3,iq),db)*pswk(:,:,:,3-is) &
                      -0.5d0*sigis*dsfz(:,:,:,2,iq)*pswk(:,:,:,is)
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ij on d/dy psi                                             *
  !**********************************************************************
  IF(j2on .OR. jfon) THEN
    DO is = 1,2  
       sigis = 3-2*is  
       pout(:,:,:,is)=pout(:,:,:,is) &
        -CMPLX(sigis*socy(:,:,:,2,iq),socy(:,:,:,1,iq),db)*pswk(:,:,:,3-is) &
        -CMPLX(0.0d0,sigis*socy(:,:,:,3,iq),db)*pswk(:,:,:,is)
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ii on d/dy psi                                             *
  !**********************************************************************
  IF(jfon) THEN
    DO is = 1,2  
       sigis = 3-2*is  
       pout(:,:,:,is)=pout(:,:,:,is) &
         -sigis*soc0(:,:,:,1,iq)*pswk(:,:,:,3-is) 
    ENDDO
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dy d/dz psi                  *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) &
             -sigis*sfden(:,:,:,2,iq)*pswk1(:,:,:,is) &
             +CMPLX(0.0d0,sigis*sfden(:,:,:,3,iq),db)*pswk1(:,:,:,3-is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require(d/dy)**2 psi                              *
  !***********************************************************************
  DO is = 1,2  
     pout(:,:,:,is) = pout(:,:,:,is) - bmass(:,:,:,iq)*pswk2(:,:,:,is)
  ENDDO

  !***********************************************************************
  !        S.F terms requiring (d/dy)**2 psi                             *
  !***********************************************************************
  IF(sfon) THEN
    DO is=1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) &
         +CMPLX(0.0d0,sigis*sfden(:,:,:,2,iq),db)*pswk2(:,:,:,3-is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require d/dz psi                                  *
  !***********************************************************************
  IF(TFFT) THEN
     CALL cdervz(pinn,pswk,d2psout=pswk2)  
     CALL cdervx(pswk,pswk1)  
  ELSE
     CALL cmulz(der1z,pinn,pswk,0)  
     CALL cmulz(der2z,pinn,pswk2,0)  
     CALL cmulx(der1x,pswk,pswk1,0)  
  ENDIF
  DO is = 1,2  
     ic = 3-is  
     sigis = 0.5d0*(3-2*is)  
     pout(:,:,:,is)=pout(:,:,:,is) &
          -CMPLX(dbmass(:,:,:,3,iq),0.5d0*xiq(:,:,:,3,iq),db)*pswk(:,:,:,is) &
          +CMPLX(sigis*wlspot(:,:,:,1,iq),-0.5d0*wlspot(:,:,:,2,iq),db)*pswk(:,:,:,ic)
  ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*xiq(:,:,:,3,iq)*pinn(:,:,:,1)&
         +CMPLX(0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*xiq(:,:,:,3,iq)*pinn(:,:,:,2)&
         +CMPLX(-0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervz(pswk2,pswk)  
    ELSE
       CALL cmulz(der1z,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
  IF(TFFT) THEN
     CALL cdervz(pinn,pswk,d2psout=pswk2)  
     CALL cdervx(pswk,pswk1)  
  ELSE
     CALL cmulz(der1z,pinn,pswk,0)  
     CALL cmulz(der2z,pinn,pswk2,0)  
     CALL cmulx(der1x,pswk,pswk1,0)  
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.T requiring d/dz psi                       *
  !***********************************************************************
  IF(ston) THEN
    DO is = 1,2
       sigis = 3-2*is
       pswk3(:,:,:,is)=-CMPLX(ssig(:,:,:,1,iq), &
                      -sigis*ssig(:,:,:,2,iq),db)*pswk(:,:,:,3-is) &
                      -sigis*ssig(:,:,:,3,iq)*pswk(:,:,:,is)
    ENDDO
  IF(TFFT) THEN
     CALL cdervz(pswk3,pswk4)  
  ELSE
     CALL cmulz(der1z,pswk3,pswk4,0)    
  ENDIF
     DO is = 1,2
       pout(:,:,:,is)=pout(:,:,:,is)+pswk4(:,:,:,is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dz psi                       *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) -CMPLX(0.5d0*dsfx(:,:,:,3,iq), &
                      -0.5d0*sigis*dsfy(:,:,:,3,iq),db)*pswk(:,:,:,3-is) &
                      -(sigis*dsfz(:,:,:,3,iq)+0.5d0*sigis*dsfx(:,:,:,1,iq) &
                      +0.5d0*sigis*dsfy(:,:,:,2,iq))*pswk(:,:,:,is)
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ij on d/dz psi                                             *
  !**********************************************************************
  IF(j2on .OR. jfon) THEN
    DO is = 1,2  
       sigis = 3-2*is  
       pout(:,:,:,is)=pout(:,:,:,is) &
        -CMPLX(sigis*socz(:,:,:,2,iq),socz(:,:,:,1,iq),db)*pswk(:,:,:,3-is) &
        -CMPLX(0.0d0,sigis*socz(:,:,:,3,iq),db)*pswk(:,:,:,is)
    ENDDO
  ENDIF
  !**********************************************************************
  !        J_ii on d/dz psi                                             *
  !**********************************************************************
  IF(jfon) THEN
    DO is = 1,2  
     sigis = 3-2*is  
     pout(:,:,:,is)=pout(:,:,:,is) &
       -CMPLX(0.0d0,sigis*soc0(:,:,:,1,iq),db)*pswk(:,:,:,is) 
    ENDDO
  ENDIF
  !***********************************************************************
  !        S.sigma parts of S.F requiring d/dz d/dx psi                  *
  !***********************************************************************
  IF (sfon) THEN
    DO is = 1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is) &
              -sigis*sfden(:,:,:,1,iq)*pswk1(:,:,:,is) &
              -sfden(:,:,:,3,iq)*pswk1(:,:,:,3-is)
    ENDDO
  ENDIF
  !***********************************************************************
  !        terms which require(d/dz)**2 psi                              *
  !***********************************************************************
  DO is = 1,2  
     pout(:,:,:,is) = pout(:,:,:,is) - bmass(:,:,:,iq)*pswk2(:,:,:,is)
  ENDDO
  !***********************************************************************
  !        S.F terms requiring (d/dz)**2 psi                             *
  !***********************************************************************
  IF(sfon) THEN
    DO is=1,2
       sigis = 3-2*is
       pout(:,:,:,is)=pout(:,:,:,is)-sigis*sfden(:,:,:,3,iq)*pswk2(:,:,:,is)
    ENDDO
  ENDIF
END SUBROUTINE hpsi

END Module Meanfield
