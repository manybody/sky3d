MODULE Twobody
  USE Params
  USE Grids
  USE Densities, ONLY: rho
  USE Moment, ONLY: cmtot
  USE Forces, ONLY: nucleon_mass
  IMPLICIT NONE
  SAVE
  REAL(db) :: roft
  REAL(db) :: rdot
  REAL(db),PRIVATE :: xmin,zmin
  REAL(db),PRIVATE :: slope,slold
  REAL(db),PRIVATE :: bb
  REAL(db) :: centerx(2),centerz(2)
  REAL(db),PRIVATE :: mass(2),charge(2),tke2body(2)
  LOGICAL  :: istwobody
  REAL(db),PARAMETER :: vacuum=0.03D0
  REAL(db) :: xmu,xlf,ecmf,tdotc,tketot,teti,tetc,tetf,tets,xcoul,xcent
  PRIVATE :: getslope, divpoint
CONTAINS
  !
  ! Calculate the two-body quantities from the density at any time,
  ! first checking whther the system can be decomposed
  !
  !*************************************************************** 
  SUBROUTINE twobody_case(xdt) 
    REAL(db),INTENT(IN) :: xdt
    REAL(db),SAVE :: xold(2),zold(2),rold,tetold 
    REAL(db) :: centx(2),centz(2), &
         angle,slopev,xx,zz,diff,rhotot,ratio,vxx(2),vzz(2), & 
         tdotc,epsf,ttt,temp,vol 
    INTEGER :: ix,iy,iz,itcm,ifrag 
    ! *** Step 1
    xold=centerx
    zold=centerz
    rold=SQRT((xold(2)-xold(1))**2 +(zold(2)-zold(1))**2)   
    slold=(zold(2)-zold(1))/(xold(2)-xold(1))   
    tetold=ATAN(slold)   
    ! *** Step 2
    Iteration: DO itcm=1,10 
       ! substep 1
       IF(itcm>1) THEN   
          centx=centerx
          centz=centerz   
          slope=(centz(2)-centz(1))/(centx(2)-centx(1))   
          bb=centz(1)-slope*centx(1)   
       ELSE   ! first iteration
          centx=xold   
          centz=zold   
          CALL getslope 
          bb=cmtot(3)-slope*cmtot(1)   
       ENDIF
       ! substep 2
       istwobody=divpoint() 
       ! substep 3
       angle=ATAN(slope)   
       slopev=dtan(angle+pi/2.0D0)   
       bb=zmin-slopev*xmin   
       ! substep 4
       centerx=0.D0; centerz=0.D0 
       mass=0.0D0   
       charge=0.0D0   
       DO iz=1,nz   
          zz=z(iz)   
          DO ix=1,nx   
             xx=x(ix)   
             diff=zz-slopev*xx-bb   
             DO iy=1,ny   
                vol=wxyz   
                rhotot=rho(ix,iy,iz,1)+rho(ix,iy,iz,2)   
                ifrag=1 ! 1 for left and 2 for right fragment 
                IF(diff<0.0D0) THEN   
                   ifrag=1 
                ELSE 
                   ifrag=2 
                ENDIF
                mass(ifrag)=mass(ifrag)+vol*rhotot   
                charge(ifrag)=charge(ifrag)+vol*rho(ix,iy,iz,2)   
                centerx(ifrag)=centerx(ifrag)+vol*rhotot*xx   
                centerz(ifrag)=centerz(ifrag)+vol*rhotot*zz   
             ENDDO
          ENDDO
       ENDDO
       centerx=centerx/mass   
       centerz=centerz/mass   
       ! substep 5 
       IF(.NOT.istwobody) EXIT 
       IF(MAX(MAXVAL(ABS(centx-centerx)),MAXVAL(ABS(centz-centerz))) & 
            <1.0d-05) EXIT
    ENDDO Iteration
    !  
    ! *** Step 3
    ! 
    xmu=nucleon_mass*mass(1)*mass(2)/(mass(1)+mass(2))   
    ratio=xmu/hbc**2   
    ! 
    vxx=(centx-xold)/xdt   
    vzz=(centz-zold)/xdt   
    tke2body=0.5D0*mass*nucleon_mass*(vxx**2+vzz**2)
    tketot=0.5D0*xmu *((vxx(1)-vxx(2))**2+(vzz(1)-vzz(2))**2)   
    roft=SQRT((centx(2)-centx(1))**2+(centz(2)-centz(1))**2) 
    rdot=(roft-rold)/xdt   
    teti=angle   
    tdotc=(teti-tetold)/xdt   
    ! 
    xlf=tdotc*xmu/hbc*roft**2   
    xcoul=charge(1)*charge(2)*e2/roft   
    xcent=xlf**2/(2.0D0*ratio*roft**2)   
    ecmf=0.5D0*xmu*rdot**2+xcoul+xcent   
    epsf=SQRT(1.0D0+2.0D0*ecmf*xlf**2/ & 
         (ratio*(charge(1)*charge(2)*e2)**2)) 
    ttt=xlf**2/(ratio*charge(1)*charge(2)*e2*roft)   
    IF(ABS(slope)>ABS(slold)) THEN   
       tetc=pi-teti   
       tetf=2.0D0*ACOS(1.0D0/epsf)-teti   
    ELSE   
       tetc=teti   
       teti=pi-tetc   
       tetf=ACOS(1.0D0/epsf)-ACOS(MAX(1.0D0/epsf *(1.0D0+ & 
            ttt),1.D0)) 
    ENDIF
    tets=tetc-tetf   
    IF(ABS(slope) <ABS(slold)) THEN   
       temp=mass(1); mass(1)=mass(2); mass(2)=temp   
       temp=charge(1); charge(1)=charge(2); charge(2)=temp   
       temp=tke2body(1); tke2body(1)=tke2body(2); tke2body(2)=temp   
    ENDIF
  END SUBROUTINE twobody_case
  !*************************************************************** 
  SUBROUTINE getslope
    REAL(db) ::  q2(3,3),xx,yy,zz,vol,denom
    INTEGER :: ix,iy,iz 
    q2=0.D0 
    DO iz=1,nz   
       zz=z(iz)-cmtot(3)   
       DO iy=1,ny   
          yy=y(iy)-cmtot(2)   
          DO ix=1,nx   
             xx=x(ix)-cmtot(1)   
             vol=wxyz*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2)) 
             q2(1,1)=q2(1,1)+(xx*xx+xx*xx-yy*yy-zz*zz)*vol   
             q2(2,2)=q2(2,2)+(yy*yy+yy*yy-xx*xx-zz*zz)*vol   
             q2(3,3)=q2(3,3)+(zz*zz+zz*zz-xx*xx-yy*yy)*vol   
             q2(1,2)=q2(1,2)+3.D0*xx*yy*vol   
             q2(1,3)=q2(1,3)+3.D0*xx*zz*vol   
             q2(2,3)=q2(2,3)+3.D0*yy*zz*vol   
          ENDDO
       ENDDO
    ENDDO
    q2(2,1)=q2(1,2)   
    q2(3,1)=q2(1,3)   
    q2(3,2)=q2(2,3)   
    denom=0.5D0*(q2(1,1)-q2(3,3)+SQRT((q2(1,1) & 
         -q2(3,3))**2+4.D0*q2(3,1)**2)) 
    IF(ABS(denom)<1.D-4) THEN   
       slope=100.D0   
    ELSE   
       slope=q2(1,3)/denom
    ENDIF
  END SUBROUTINE getslope
  !*************************************************************** 
  LOGICAL FUNCTION divpoint() 
    LOGICAL :: online,two,in_vacuum 
    INTEGER :: iyy,iz1,iz2,idz,il,ix,iz,ixl(nx+ny+nz), &
         izl(nx+ny+nz),i,n1,n2,nf 
    REAL(db) :: deltax,deltaz,xx,zz,rhol(nx+ny+nz) 
    ! 
    ! now the calculation of the dividing plane starts by examining the line 
    ! along the axis of largest quadrupole moment 
    ! 
    iyy=ny/2;
    deltax=0.5*dx
    deltaz=0.5*dz
    ! Loop 1
    IF(slope>=0.D0) THEN 
       iz1=1; iz2=nz; idz=1 
    ELSE 
       iz1=nz; iz2=1; idz=-1 
    END IF
    il=0 
    DO ix=1,nx 
       xx=x(ix) 
       DO iz=iz1,iz2,idz 
          zz=z(iz) 
          online=.FALSE. 
          IF(ABS(slope)<=1.D0) THEN 
             online=ABS(zz-slope*xx-bb)<=deltaz 
          ELSE 
             online=ABS(xx-(zz-bb)/slope)<=deltax 
          END IF
          IF(online) THEN 
             il=il+1 
             IF(il>nx+ny+nz) THEN 
                IF(wflag) WRITE(*,*) ' Increase dimensioning in function divpoint' 
                STOP 
             END IF
             ixl(il)=ix; izl(il)=iz 
             rhol(il)=rho(ix,iyy,iz,1)+rho(ix,iyy,iz,2) 
          END IF
       END DO
    END DO
    ! Loop 2
    nf=0; n1=0; n2=0; in_vacuum=.TRUE. 
    DO i=1,il
       IF(rhol(i)>vacuum) THEN 
          IF(in_vacuum) THEN 
             in_vacuum=.FALSE. 
             nf=nf+1 
             IF(nf==2) n2=i 
             IF(nf>2) EXIT 
          END IF
       ELSE 
          IF(.NOT.in_vacuum) THEN 
             IF(nf==1) n1=MAX(1,i-1) 
             in_vacuum=.TRUE. 
          END IF
       END IF
    END DO
    two=nf==2 
    ! final processing
    IF(two) THEN 
       xmin=0.5D0*(x(ixl(n1))+x(ixl(n2))) 
       zmin=0.5D0*(z(izl(n1))+z(izl(n2))) 
    END IF
    divpoint=two 
  END FUNCTION divpoint
  !*************************************************************** 
  SUBROUTINE twobody_print
    REAL(db) :: centx(2),centz(2)
    INTEGER :: i 
    CHARACTER(8),PARAMETER :: name(2)=(/ '  left: ',' right: '/) 
    IF(wflag)  &
         WRITE(*,'(/A,/(4(A12,F12.4)))') & 
         ' Relative motion / Coulomb kinematics:',' red. mass:', & 
         xmu/nucleon_mass,' l/hbar:',xlf,' ecmf(MeV):',ecmf, & 
         ' roft(fm):',roft,' rdot/c:',rdot,' td/c(°):', & 
         tdotc*180.0D0/pi,' trke(MeV):',tketot, & 
         ' teti(°):',teti*180.0D0/pi,' tetf(°):', & 
         tetf*180.0D0/pi,' tetc(°):',tetc*180.0D0/pi, & 
         ' tets(°):',tets*180.0D0/pi,' Vcoul:',xcoul, & 
         ' Vcent:',xcent 
    IF(ABS(slope)>ABS(slold)) THEN   
       centx=centerx; centz=centerz
    ELSE   
       centx(1)=centerx(2); centz(1)=centerz(2)
       centx(2)=centerx(1); centz(2)=centerz(1)
    ENDIF
    ! 
    IF(wflag)  THEN
       WRITE(*,'(/a/a)') ' Collision kinematics:', & 
            '  Side         Mass       Charge     <x>         &
            &<y>        tke'
       DO i=1,2 
          WRITE(*,'(A,2F12.4,1P,3E12.4)') name(i),mass(i),charge(i), & 
               centx(i),centz(i),tke2body(i) 
       ENDDO
    ENDIF
  END SUBROUTINE twobody_print
END MODULE Twobody
