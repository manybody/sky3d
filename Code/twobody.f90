!------------------------------------------------------------------------------
! MODULE: Twobody
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the code to analyze the final (and also, though
!!less useful) initial stage of a heavy-ion reaction. It is applicable
!!only to the dynamic calculations and only if there are essentially two
!!separated fragments. It calculates the fragment masses and charges,
!!their distance, the relative motion kinetic energy, angular momentum,
!!and the scattering angle.
!!
!>
!>@details 
!!The analysis assumes that the two fragments are separated by a region
!!of noticeably lower density, tries to find the line connecting their
!!centers of mass and then divides up space by a plane perpendicular to
!!this line. Of necessity the results are not of high precision, but
!!tend to be useful nevertheless.
!!
!!<b> The user should be critical and always make sure that the real
!!  situation is a two-fragment one before accepting the results of
!!  these routines. It is assumed that the reaction plane is the
!!  (x,z)-plane. The code could be generalized to a three-dimensional
!!  situation, if necessary, but usually the assumption of a fixed
!!  scattering plane should not be a problem. If the initial nuclei have
!!  nonzero internal angular momentum, e.~g., the reaction plane can
!!  rotate and then a more general analysis cannot be avoided.</b>
!------------------------------------------------------------------------------
MODULE Twobody
  USE Params
  USE Grids
  USE Densities, ONLY: rho
  USE Moment, ONLY: cmtot
  USE Forces, ONLY: nucleon_mass
  IMPLICIT NONE
  SAVE
  REAL(db) :: roft                    !<separation distance between fragments in fm.
  REAL(db) :: rdot                    !<relative-motion velocity in units of \f$ c \f$.
  REAL(db),PRIVATE :: xmin            !<x-coordinate of the point in
  !!the (x,z)-plane where minimum density is found between the fragments.
  REAL(db),PRIVATE :: zmin            !<x-coordinate of the point in
  !!the (x,z)-plane where minimum density is found between the fragments.
  REAL(db),PRIVATE :: slope           !<Slope of the line
  !!connecting the fragment centers-of-mass.
  REAL(db),PRIVATE :: slold           !<the value from the previous time step or iteration 
  !!of \c slope where the 2-body analysis was last performed.
  REAL(db),PRIVATE :: bb              !<intercept of the line. The line is thus
  !!given as <tt>z=bb+slope*x</tt>.
  REAL(db) :: centerx(2)              !<x-coordinates of the two fragment centers of mass.
  REAL(db) :: centerz(2)              !<z-coordinates of the two fragment centers of mass.
  REAL(db),PRIVATE :: mass(2)         !<masses of the two fragments.
  REAL(db),PRIVATE :: charge(2)       !<charges of the two fragments
  REAL(db),PRIVATE :: tke2body(2)     !<total kinetic energy in MeV.
  LOGICAL  :: istwobody               !<logical, indicates whether this is a
  !!two-body case (\c TRUE) or not.
  REAL(db),PARAMETER :: vacuum=0.03D0 !<parameter indicating the limiting density below
  !!which vacuum is assumed. This value is not critical, since it is
  !!only used in looking for "empty" regions between the nuclei.
  REAL(db) :: xmu                     !<the reduced mass in MeV.
  REAL(db) :: xlf                     !<angular momentum of relative motion in \f$ \hbar \f$. It
  !!is calculated assuming two point bodies via \f$ \mu R^2\omega \f$.
  REAL(db) :: ecmf                    !<final relative motion energy after extrapolated to
  !!infinite separation.
  REAL(db) :: tdotc                   !<time derivative of scattering angle. Units \f$ c/{\rm fm} \f$.
  REAL(db) :: tketot                  !<relative-motion kinetic energy in MeV.
  REAL(db) :: teti                    !< present scattering angle.
  REAL(db) :: tetc                    !<TODO
  REAL(db) :: tetf                    !<TODO
  REAL(db) :: tets                    !<TODO
  REAL(db) :: xcoul                   !<Coulomb energy, calculated from point charges.
  REAL(db) :: xcent                   !<centrifugal energy, based on angular momentum and distance.
  PRIVATE :: getslope, divpoint
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: twobody_case
!> @brief
!!This routine tries to find a separation of the system into two
!!fragments and determine their properties as well as those of the
!!relative motion. It keeps track of previous analysis results, since
!!the 2-body properties are expected to change slowly and this makes the
!!analysis easier. Also, the motion of the fragments is determined by
!!time-differencing.
!>
!> @details
!!It uses the following steps: 
!!  - <b> Step 1 </b>: the results of the last analysis are partially
!!    saved.  This includes the positions of the fragments, their
!!    distance, the slope of the connecting line, and its angle.
!!  - <b> Step 2 </b>: the fragment division is sought in an iterative
!!    process with (at present) a fixed iteration limit of 10. The reason
!!    for the iterations is that the connecting line, which determines the
!!    dividing plane, which is orthogonal to it at the dividing point,
!!    should link the centers-of-mass, but these can be calculated only
!!    assuming a dividing plane. There is thus a self-consistency problem
!!    which as usual is solved iteratively.  During the iterations the new
!!    value of the center of mass is kept in \c centerx and \c centerz, 
!!    while the previous one is in \c centx and \c centz. Each iteration 
!!    has several steps:
!!     -# The connecting line (\c slope and intercept \c bb) is
!!        determined in one of two ways: if this is the first iteration,
!!        subroutine getslope is called to calculate the slope from the
!!        quadrupole tensor; the intercept is then calculated assuming that
!!        the line passes through the center-of-mass (which it should do in
!!        principle, but it might be off in practice). For the later
!!        iterations the line is calculated directly from the two centers of
!!        mass.
!!     -# Next function \c divpoint is called to find out whether
!!        this is indeed a two-body situation. It also returns the midpoint
!!        between the fragments in \c xmin and \c zmin. This result is
!!        not used immediately, since it might be that shifting the
!!        connecting line could change this situation.
!!     -# The slope of the line is now rotated by \f$ 90^\circ \f$ to get the
!!        line defining the dividing plane in the \f$ (x,z) \f$-plane. This has
!!        \c slopev and \c bb as slope and intercept.
!!     -# Now in a simple loop integrals are done separately over the
!!        two regions, summing up \c charge, \c mass, and \c center
!!        of mass for the two fragments. The assignment to the fragments is
!!        recognized by checking whether the point is above or below the
!!        dividing line (variable \c diff). Note that the \f$ y \f$-component
!!        of the centers of mass is not calculated.
!!     -# Finally the iteration process is stopped if no two-body
!!        situation is found or the center-of-mass vectors have converged.
!!     .
!!  - <b> Step 3 </b>: The two-body analysis is now assumed to be complete
!!    and the physical quantities are evaluated. These are defined above
!!    in the list of module variables.  The following calculation concerns
!!    the final scattering angle \c tets extrapolated to infinity. It
!!    can be understood using the Rutherford trajectories.
!>
!> @param[in] xdt
!> REAL(db), takes the current time step. 
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------  
! DESCRIPTION: getslope
!> @brief
!!This subroutine calculates the slope of the line determined by the
!!eigenvector of largest quadrupole moment in the (x,z)-plane. The
!!first loop sums up the quadrupole tensor, which is dimensioned <tt>
!!q2(3,3)</tt> but of which the index 2 is not actually necessary, since
!!the y-direction is not involved. The definition is kept
!!three-dimensional to reduce confusion and make later generalization
!!easier.
!>
!> @details
!!The two-dimensional eigenvalue problem for \c q2 has the secular
!!equation (remember \f$ q_{31}=q_{13} \f$:
!!\f[ (q_{11}-\lambda)(q_{33}-\lambda)-q_{13}^2=0, \f]
!!with \f$ \lambda \f$ the eigenvalue. For the larger eigenvalue we get
!!\f[ \lambda=\frac1{2}\left(q_{11}+q_{33}+\sqrt{(q_{11}-q_{33})^2
!!-4q_{13}^2}\right), \f] and solving the equation
!!\f$ q_{13}x+(q_{33}-\lambda)z=0 \f$ for z yields
!!\f{eqnarray}{
!!  z&=&\frac{q_{13}}{\lambda-q_{33}}\,x\\
!!  &=& \frac{q_{13}}{\frac1{2}\left(q_{11}-q_{33}+\sqrt{(q_{11}-q_{33})^2
!!        -4q_{13}^2}\right)}.
!!\f}
!!The code calls the denominator \c denom and makes sure no division
!!by zero happens (this could happen for a spherical distribution).  The
!!resulting slope is returned in the module variable \c slope. 
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: divpoint
!> @brief
!!The function divpoint is a helper routine.  It examines the line
!!determined by \c slope and intercept \c bb and finds the point
!!<tt>(xmin,zmin)</tt> which is in the center of the void between the two
!!fragments, returning .true. if this is possible.
!>
!> @details
!!To this end it looks at the behavior of the densities along this line.
!!Since the line has no relation to the numerical grid, this is not
!!trivial.
!!
!!  - <b>First loop:</b> Essentially it looks through the (x,z)-plane to
!!    find points closer than half a grid spacing to the desired line with
!!    equation <tt> z=slope*x+bb</tt> (logical variable \c online).  If the
!!    slope is larger than one, i.~e., if the nuclei are separating
!!    predominantly in the x-direction, we need to take the equation
!!    <tt>x=(z-bb)/slope</tt> instead to get better resolution. To make the
!!    result monotonic along the line, the do loop in z runs backward
!!    for negative slopes. The points found are collected in index vectors
!!    \c ixl, \c izl with accompanying densities <tt> rhol </tt> stored
!!    in arrays of length \c il.
!!
!!  - <b>Second loop:</b> now the number of fragments \c nf is counted by
!!  examining this one-dimensional density curve, looking for
!!  disconnected density humps above \c vacuum density. It is also
!!  recorded where the "vacuum" region starts and ends in variables
!!  \c n1 and \c n2. The logic is as follows:
!!    -# The logical variables \c in_vacuum} keeps track of whether
!!       the search is in a vacuum region at the moment or not. It starts
!!       as \c TRUE.
!!    -# Go to the next point. If its density is above \c vacuum,
!!       and we are in the vacuum, a new fragment is starting and we
!!       increase the number of fragments \c nf by1. If it becomes
!!       bigger than 2, exit, because there are three or more fragments.
!!       If it is now 2, record the starting index for the second fragment
!!       in \c n2.
!!    -# If the density is below \c vacuum and we are not in vacuum,
!!       a fragment is being ended. If it is the first fragment, we record
!!       this index in \c n1.
!!  - <b>Final processing</b> At this point we expect a two-fragment
!!    situation if \c nf=2 and in this case the void region between the
!!    fragments extends from \c n1 to \c n2, which are indices into
!!    arrays \c ix1 and \c iz1 giving the position in the
!!    (x,z)-plane. The code calculates the midpoint between the two
!!    positions and returns \c TRUE in this case.
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: twobody_print
!> @brief
!!This subroutine writes the relevant values determined by the twobody 
!!analysis to the main output.
!---------------------------------------------------------------------------
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
