!------------------------------------------------------------------------------
! MODULE: Twobody
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the routines to distinguish between a one and two-body 
!!system and to perform a two-body analysis in the case of a two-body collision 
!!including scattering angles and momenta.
!------------------------------------------------------------------------------
MODULE Twobody
  USE Params
  USE Grids
  USE Densities
  USE Coulomb
  USE Forces, ONLY: nucleon_mass
  IMPLICIT NONE
  ! Constants
  REAL(db),PARAMETER :: threshold=1.D-5 !< threshold density to recognze separate fragmenst
  ! Data describing initial configuration
  INTEGER :: mass_case !< records the initial mass relation. 0 indicates a symmetric
  !! collision, 1 that fragment #1 is the lighter one, 2 that fragment #2 is heavier.
  !! This is needed to order the fragments the same way in the final state.
  ! Data related to final state
  REAL(db) :: roft !< fragment separation at the current time, zero means connected system
  REAL(db) :: roft_old !< fragment separation at the last time step. Needed to check
  !! whether the motion is outgoing or incoming
  REAL(db) :: velocity(3,0:2) !< c.m. velocity vectors for the total system and the fragments
  ! Data for fragments that are needed globally
  ! index 0 for total, 1:2 for fragments
  !>@name mass,charge, c.m., momentum, and angular momentum of the total system and the fragments
  !>@{
  REAL(db) :: mass(0:2),charge(0:2),cmfrag(3,0:2),momentum(3,0:2),angmom(3,0:2)
  !>@}
  ! Data associated with the Coulomb calculation
  REAL(db),ALLOCATABLE,PRIVATE :: rhotot(:,:,:) !< total density 
  REAL(db),ALLOCATABLE,PRIVATE :: currtot(:,:,:,:) !< total current density
  INTEGER,ALLOCATABLE,PRIVATE :: frag(:,:,:) !< for each mesh pint, tells whether it is in 
  !! fragment 1 or 2, or below the cutoff density, (value 0)
  LOGICAL  :: istwobody !< principal return value. Tells whether the system is separated.
CONTAINS
!--------------------------------------------------------------------------- 
! DESCRIPTION: twobody_analysis
!> @brief
!! For a two-body collision, this subroutine checks whether the system
!! is separated and in that case determines either only the distance of 
!! the fragments or, for the final state, produces a detailed analysis
!>
!> @details
!! The subroutine starts off by allocating arrays \c rhotot and, if the full 
!! analysis is desired, \c currtot, which are then assigned the total density and 
!! current density, respectively. These are set to zero wherever the density
!! is below the \c threshold. Then the following steps are taken:
!! - For the full analysis the total current density and the isospin-dependent
!!   density and current density are also cut off. Then the subroutine 
!!   \c analyze is called for the whole system by setting \c frag such that 
!!   all points have fragment index 0. The coordinate values are corrected
!!   for a shift in the center-of-mass. The total current density is then 
!!   adjusted to remove any overall motion of the system in the grid.
!! - Following this, subroutine \c cut_system is called to calculation
!!   fragment properties in the case of separation. If \c full_analysis is true,
!!   this should be the case, as it is done only for the final state.
!! - For \c full_analysis=false, the analysis is finished, since the fragment
!!   distance has been calculated. Therefore there is a \c RETURN after the arrays
!!   have been deallocated. 
!! - After this the subroutine \c analyze is called separately for the two 
!!   fragments to determine their properties
!! - Based on this information, subroutine \c scattering is called to calculate
!!   relative motion and the overall scattering results.
!>
!> @param[in] full_analysis
!> LOGICAL, if TRUE the full analysis is done, otherwise only the 
!> fragment separation is calculated.
!--------------------------------------------------------------------------- 
  SUBROUTINE twobody_analysis(full_analysis)
    LOGICAL,INTENT(IN) :: full_analysis
    INTEGER :: i,j,k
    !***************************************************************
    ALLOCATE(rhotot(nx,ny,nz),frag(nx,ny,nz))
    rhotot=rho(:,:,:,1)+rho(:,:,:,2)
    ! cut off rhotot to threshold
    WHERE(rhotot<threshold) rhotot=0.D0
    ! calculate properties of total system
    IF(full_analysis) THEN
       ALLOCATE(currtot(nx,ny,nz,3))
       currtot(:,:,:,:)=current(:,:,:,:,1)+current(:,:,:,:,2)
       FORALL(i=1:nx,j=1:ny,k=1:nz,rhotot(i,j,k)<threshold)
          rho(i,j,k,:)=0.D0
          currtot(i,j,k,:)=0.D0
          current(i,j,k,:,:)=0.D0  
       END FORALL
       ! analyze total system
       frag=0
       CALL analyze(0)
       ! correct coordinates to have zero at c.m.
       x=x-cmfrag(1,0)
       y=y-cmfrag(2,0)
       z=z-cmfrag(3,0)
       cmfrag(:,0)=0.D0
       ! correct total momentum to zero
       WRITE(*,'(A,3F10.6)') 'Uncorrected velocity:',momentum(:,0)/mass(0)
       DO i=1,3
          momentum(i,0)=momentum(i,0)/mass(0) ! average velocity
          currtot(:,:,:,i)=currtot(:,:,:,i)-rhotot*momentum(i,0)
       END DO
    END IF
    ! compute dividing plane
    CALL cut_system(full_analysis)
    ! return for reduced analysis
    IF(.NOT.full_analysis) THEN
       DEALLOCATE(rhotot,frag)
       RETURN
    END IF
    ! calculate properties separately for the two fragments
    CALL analyze(1)
    CALL analyze(2)
    ! finally determine relative motion and total scattering effect
    CALL scattering
  END SUBROUTINE twobody_analysis
!--------------------------------------------------------------------------- 
! DESCRIPTION: analyze
!> @brief
!! This subroutine calculates properties of the two separate fragment
!! or of the whole system.
!>
!> @details 
!! It evaluates integrated quantities for a subset of the mesh controlled
!! by its argument \c ifrag in conjunction with the array \c
!! frag. More specifically, \c ifrag=0 is used to calculate on the
!! whole mesh while \c ifrag=1 or 2 select the points belonging to the
!! corresponding fragment. The array \c frag contains the relevant index 
!! for each mesh point.
!! - Most of the calculated quantities have an index \c ifrag to store
!!   these results separately. Thus \c cm(1:3,0:2) is used to store the
!!   total center of mass in \c cm(:,0) and the centers of mass of the
!!   two fragments in \c cm(:,1) or \c cm(:,2), respectively.
!! - The first loop and the following statements calculate the mass number,
!!   charge number, center of mass, and total momentum of the system or
!!   fragment. From the latter two velocities are derived: \c velcorr,
!!   which is in the units used for the \c current density in Sky3D,
!!   namely fm\f$^{-4}\f$, while \c velocity is in units of \f$c\f$.
!! - The next loop calculates the collective orbital angular momentum \c
!!   angmom with respect to the center of mass, subtracting its motion,
!!   and the classical tensor \f$\Theta\f$ of inertia, \c inertia. This allows
!!   calculating also the instantaneous angular velocity vector \c omega
!!   by solving \f$\vec L=\Theta\vec\omega\f$ using the subroutine \c gauss.
!>
!> @param[in] ifrag
!> INTEGER, index of fragment 1 or 2 to be investigated, or 0 for whole system
!--------------------------------------------------------------------------- 
  SUBROUTINE analyze(ifrag)
    INTEGER,INTENT(IN) :: ifrag
    REAL(db) :: velcorr(3,0:2),omega(3),inertia(3,3),curr(3),xx,yy,zz
    INTEGER :: i,j,k
    mass(ifrag)=0.D0
    charge(ifrag)=0.D0
    cmfrag(:,ifrag)=0.D0
    momentum(:,ifrag)=0.D0
    DO k=1,nz
       DO j=1,ny
          DO i=1,nx
             IF(frag(i,j,k)/=ifrag) CYCLE
             mass(ifrag)=mass(ifrag)+rhotot(i,j,k)
             charge(ifrag)=charge(ifrag)+rho(i,j,k,2)
             cmfrag(1,ifrag)=cmfrag(1,ifrag)+rhotot(i,j,k)*x(i)
             cmfrag(2,ifrag)=cmfrag(2,ifrag)+rhotot(i,j,k)*y(j)
             cmfrag(3,ifrag)=cmfrag(3,ifrag)+rhotot(i,j,k)*z(k)
             momentum(:,ifrag)=momentum(:,ifrag)+currtot(i,j,k,:)
          END DO
       END DO
    END DO
    cmfrag(:,ifrag)=cmfrag(:,ifrag)/mass(ifrag)
    ! fragment velocity in units of current density
    velcorr(:,ifrag)=momentum(:,ifrag)/mass(ifrag)
    ! fragment velocity in units of c
    velocity(:,ifrag)=momentum(:,ifrag)*wxyz*hbc/ &
         (nucleon_mass*mass(ifrag))
    ! add volume element to quantities depending on it
    mass(ifrag)=mass(ifrag)*wxyz
    charge(ifrag)=charge(ifrag)*wxyz
    momentum(:,ifrag)=momentum(:,ifrag)*wxyz
    ! Now calculate orbital angular momentum for fragment, relative
    ! to its c.m.
    angmom(:,ifrag)=0.D0
    inertia=0.D0
    DO k=1,nz
       zz=z(k)-cmfrag(3,ifrag)
       DO j=1,ny
          yy=y(j)-cmfrag(2,ifrag)
          DO i=1,nx
             xx=x(i)-cmfrag(1,ifrag)
             IF(frag(i,j,k)/=ifrag) CYCLE
             curr=currtot(i,j,k,:)-rhotot(i,j,k)*velcorr(:,ifrag)
             angmom(1,ifrag)=angmom(1,ifrag)+yy*curr(3)-zz*curr(2)
             angmom(2,ifrag)=angmom(2,ifrag)+zz*curr(1)-xx*curr(3)
             angmom(3,ifrag)=angmom(3,ifrag)+xx*curr(2)-yy*curr(1)
             inertia(1,1)=inertia(1,1)+rhotot(i,j,k)*(yy**2+zz**2)
             inertia(2,2)=inertia(2,2)+rhotot(i,j,k)*(xx**2+zz**2)
             inertia(3,3)=inertia(3,3)+rhotot(i,j,k)*(xx**2+yy**2)
             inertia(1,2)=inertia(1,2)-rhotot(i,j,k)*xx*yy
             inertia(1,3)=inertia(1,3)-rhotot(i,j,k)*xx*zz
             inertia(2,3)=inertia(2,3)-rhotot(i,j,k)*yy*zz
          END DO
       END DO
    END DO
    angmom(:,ifrag)=angmom(:,ifrag)*wxyz
    inertia(2,1)=inertia(1,2)
    inertia(3,1)=inertia(1,3)
    inertia(3,2)=inertia(2,3)
    inertia=inertia*wxyz
    CALL gauss(inertia,angmom,omega)
    WRITE(*,'(A,I2)') '******* Fragment #',ifrag
    WRITE(*,'(" c.m.=",3F10.4," mass=",F10.4," charge=",F10.4)') &
         cmfrag(:,ifrag),mass(ifrag),charge(ifrag)
    Write(*,'(" Velocity/c =",3f10.4)') velocity(:,Ifrag)
    WRITE(*,100) 'Inertia tensor',inertia
    WRITE(*,100) 'Angular momentum',angmom(:,ifrag)
    WRITE(*,100) 'Omega from L',omega
100 FORMAT(A/(3D15.6))
  END SUBROUTINE analyze
!--------------------------------------------------------------------------- 
! DESCRIPTION: gauss
!> @brief
!! This subroutine does a simple direct solution of the 3 by 3 linear
!! system of equations.
!>
!> @details
!! This is an explicit implementation using the determinant algorithm.
!>
!> @param[in] a(3,3)
!> REAL(db), Matrix of the linear system.
!> @param[in] b(3) 
!> REAL(db), right-hand side of the equations.
!> @param[out] c(3)
!> REAL(db), solution vector.
!--------------------------------------------------------------------------- 
  SUBROUTINE gauss(a,b,c)
    REAL(db) :: a(3,3),b(3),c(3),det
    INTENT(IN) :: a,b
    INTENT(OUT) :: c
    det=a(1,3)**2*a(2,2)-2.D0*a(1,2)*a(1,3)*a(2,3)+a(1,1)*a(2,3)**2+ &
         a(1,2)**2*a(3,3)-a(1,1)*a(2,2)*a(3,3)
    IF(ABS(det)<1.D-8) THEN
      c=0.D0
    ELSE
      c(1)=(a(2,3)**2*b(1)-a(2,2)*a(3,3)*b(1)-a(1,3)*a(2,3)*b(2)+ &
           a(1,2)*a(3,3)*b(2)+a(1,3)*a(2,2)*b(3) &
           -a(1,2)*a(2,3)*b(3))/det
      c(2)=(-(a(1,3)*a(2,3)*b(1))+a(1,2)*a(3,3)*b(1)+a(1,3)**2*b(2) &
           -a(1,1)*a(3,3)*b(2)-a(1,2)*a(1,3)*b(3) &
           +a(1,1)*a(2,3)*b(3))/det
      c(3)=(a(1,3)*a(2,2)*b(1)-a(1,2)*a(2,3)*b(1)-a(1,2)*a(1,3)*b(2) &
           +a(1,1)*a(2,3)*b(2)+a(1,2)**2*b(3) &
           -a(1,1)*a(2,2)*b(3))/det
    ENDIF
  END SUBROUTINE gauss
  !*************************************************************** 
!--------------------------------------------------------------------------- 
! DESCRIPTION: scattering
!> @brief
!! This subroutine calculates the relative motion and the various energies 
!! in the final state and, by referring back to the initial state before the 
!! collision, also the net scattering results.
!>
!> @details
!! The subroutine calculates the scattering angle as the sum of three contributions:
!! the ingoing Rutherford angle, the rotation angle during the TDHF calculation, 
!! and the outgoing Rutherford angle.
!! 
!! The Rutherford trajectories are hyperbolas described in polar
!! coordinates as \f[r(\phi)=\frac{k}{\epsilon\cos\phi-1}.\f]
!! with \f$r\f$ and \f$\phi\f$ describing the relative vector of the two-body problem.
!! It is confined to the angular region \f$-\phi_\infty<\phi<\phi_\infty\f$
!! with the asymptotic angle \f$\phi_\infty=\cos^{-1}(\epsilon^{-1})\f$. The
!! angle \f$\phi=0\f$ corresponds to the symmetry axis of the hyperbola.
!! 
!! A problem is that this axis is not directly connected to the
!! coordinates used for the TDHF calculation, which is assumed to have
!! \f$x\f$ and \f$z\f$ as the reaction plane. Since, however, the parameters of
!! the hyperbola \f$k\f$ and \f$\epsilon\f$ can be calculated from the impact
!! parameter and energy of the relative motion, it suffices to know \f$r\f$
!! to calculate the angle corresponding to it in the hyperbola system of
!! coordinates. This is done slightly differently for the incoming and
!! outgoing trajectories. 
!! 
!! The following steps are done:
!! - Examine the initial state before the collision: Since the properties
!! of the fragments may not be available anymore after a restart, we have to 
!! partially restore them from the original input files. <i>This requires 
!! that these files, or symbolic links to them, be present in current directory.</i>
!! From \c for005 it reads the center-of-mass energy \f$E_{\rm cm}\f$ and
!! impact parameter \f$b\f$ as well as the initial positions of the nuclei
!! for the start of the TDHF calculation, \f$\vec r_{0i},\;i=1,2\f$ called
!! \c fcent(:,i) in the code. Then the masses \f$A_i\f$ and charges \f$Z_i\f$
!! are input from the fragment data files.
!! - As the next step, the variable \c mass_case is determined to make sure the
!! mass relation is the same in the ingoing and outgoing case. For symmetric 
!! collisions of course the scattering is not uniquely determined.
!! From these data the following quantities are successively calculated:
!!   -# The initial relative vector for the TDHF calculation 
!!   \f[\vec r(t=0)=\vec r_2(t=0)-\vec r_1(t=0),\f]
!!   -# the initial reduced mass
!!   \f[\mu=m_0\frac{A_1A_2}{A_1+A_2},\f]
!!   -# the dimensionless initial orbital angular momentum of relative motion,
!!   \f[L_{\rm init}=\mu v_{\rm init}b=\mu\sqrt{\frac{2E_{\rm cm}}{\mu}}\,
!!   \frac{b}{\hbar}.\f]
!! From these quantities the parameters of the hyperbola can be
!! calculated:
!!   -# the parameter \f$ k \f$
!!   \f[ k=\frac{L_{\rm init}^2\hbar^2}{e^2Z_1Z_2\mu},\f]
!!   -# the major axis
!!   \f[a=\frac{e^2Z_1Z_2}{2E_{\rm cm}},\f]
!!   -# from which \f$ \epsilon \f$ results:
!!   \f[\epsilon=\frac{\sqrt{a^2+b^2}}{a}.\f]
!! Now the angles can be calculated from the equation of the
!! hyperbola. The asymptotic angle
!! \f[\phi_\infty=\cos^{-1}(1/\epsilon)\f]
!! and the angle corresponding to the initial position,
!! \f[\phi(t=0)=\cos^{-1}\left(\frac{k/|\vec r(t=0)|+1}{\epsilon}\right).\f]
!! In the code both are converted to degrees. The incoming Rutherford
!! scattering angle is just the difference,
!! \f[\phi_{\rm incoming}=\phi_\infty-\phi(t=0).\f]
!! - The rotation of the system during the TDHF calculation is obtained from the
!! initial and final positions of the fragments via the scalar product of the 
!! relative vectors.
!! - The code then calculates the relative motion properties in the final state:
!! reduced mass \c mu, relative velocity \c vrel, kinetic energy \c erel, 
!! and angular momentum of relative motion \c relangmom. 
!! In addition the point-charge Coulomb interaction energy \c ecoul_point,
!! the velocity and kinetic energy in the radial direction \c rdot and \c edot
!! are produced for information.
!! - For the extrapolation of the relative motion energy to infinite separation
!! a more refinde calculation of the Coulomb energy is needed. This is done by 
!! using  subroutine \c poisson three times: once for the total system and then 
!! for each fragment individually by removing the charge density from all space points 
!! associated with the other fragment. This yields \c ecoul(0:2) and the Coulomb
!! interaction energy \c ecoul_int can then be obtained as 
!! <tt>ecoul(0)-ecoul(1)-ecoul(2)</tt>. The asymptotic energy is then the sum of \c erel and
!! \c ecoul_int.
!! - Finally the deflection angle for the outgoing Rutherford trajectory 
!! is obtained. The main difference to
!! the incoming trajectory is that now the quantities determining the
!! hyperbola parameters are taken from the two-body analysis, namely the
!! masses and charges of the fragments, the \f$y\f$-component of the angular
!! momentum of relative motion \c relangmom(1:2), and the final energy
!! of relative motion \c final\_ecm. Otherwise the
!! calculation proceeds identically, though the variable \c b has to
!! be renamed \c b, since the impact parameter may be different for
!! the final state.
!! - The subroutine then prints all the details.
!>
!--------------------------------------------------------------------------- 
  SUBROUTINE scattering
    ! Varaiables needed for the calculation of the initial
    ! Rutherford trajectory
    CHARACTER(64) :: filename(2) ! names of fragment wf files
    REAL(db) :: fcent(3,2)   ! c.m. vectors for the fragments
    REAL(db) :: fboost(3,2)  ! boost vectors for the fragments
    REAL(db) :: b ! impact parameter
    REAL(db) :: ecm  ! relative kinetic energy 
    REAL(db) :: mu ! reduced mass
    REAL(db) :: k,a,epsilon,p,pmax ! hyperbola parameters
    REAL(db) :: time ! needed for the namelist
    REAL(db) :: initial_mass(2) ! masses of the original fragments
    REAL(db) :: initial_charge(2) ! charges of the initial fragments
    REAL(db) :: initial_ecm ! c.m. energy in the asymptotic initial state
    REAL(db) :: initial_L ! initial orbital angular momentum of relative motion
    REAL(db) :: initial_angle ! scattering angle for the Rutherford trajectory
    !! from infinity to the starting position
    REAL(db) :: initial_roft ! fragment distance at start of calculation
    CHARACTER(8) :: forcename
    REAL(db) :: init_relvec(3)  ! relative vector at the beginning of the calculation
    LOGICAL :: fix_boost
    INTEGER :: i,dummy(7),iter
    ! rotation of the interacting system
    REAL(db) :: tdhf_angle ! rotation angle of the system during the TDHF calculation
    ! Variables for the outgoing Rutherford trajectory
    REAL(db) :: rdot ! relative velocity
    REAL(db) :: edot ! relative kinetic energy at the present separation
    REAL(db) :: final_ecm ! asymptotic relative  motion kinetic energy
    REAL(db) :: aa,bb
    REAL(db) :: final_angle ! scattering angle for the outgoing Rutherford trajectory
    REAL(db) :: final_relvec(3) ! relative coordinate at the end of the calculation
    REAL(db) :: ecoul_point,vrel,erel,rhof(nx,ny,nz),ecoul(0:2)
    INTEGER :: ifrag
    REAL(db) :: e_centrifugal ! centrifugal energy
    REAL(db) :: ecoul_int ! Coulomb interaction energy
    REAL(db) :: relangmom(3) ! angular momentum of relative motion
    ! First read the original input file defining the initialization of
    ! the collision
    NAMELIST /fragments/ filename,fcent,fboost,ecm,b,fix_boost
    OPEN(10,FILE='for005')
    READ(10,fragments)
    CLOSE(10)
    initial_ecm=ecm
    DO i=1,2
       OPEN(10,FILE=filename(i),STATUS='old',FORM=&
            'unformatted')
       READ(10) iter,time,forcename,dummy,initial_charge(i),initial_mass(i)
       CLOSE(10)
    END DO
    ! distinguish cases of light vs. heavy fragment
    IF(initial_mass(1)==initial_mass(2)) THEN
       mass_case=0  ! equal masses
    ELSEIF(initial_mass(1)<initial_mass(2)) THEN
       mass_case=1  ! lighter fragment is first
    ELSE
       mass_case=2  ! heavier fragment is first
    END IF
    ! Now compute scattering angle due to incoming Rutherford trajectory
    init_relvec=fcent(:,2)-fcent(:,1)  ! relative vector
    mu=nucleon_mass*initial_mass(1)*initial_mass(2) &
         /(initial_mass(1)+initial_mass(2)) ! reduced mass
    initial_L=mu*SQRT(2.D0*initial_ecm/mu)*b/hbc
    k=(initial_L*hbc)**2/(e2*initial_charge(1)*initial_charge(2)*mu)
    a=e2*initial_charge(1)*initial_charge(2)/(2.D0*initial_ecm)
    epsilon=SQRT(a**2+b**2)/a
    pmax=ACOS(1.D0/epsilon)*180.D0/pi
    initial_roft=SQRT(SUM(init_relvec**2))
    p=ACOS((k/initial_roft+1.D0)/epsilon)*180.D0/pi
    initial_angle=pmax-p
    WRITE(*,*) '*** Data for incoming Rutherford hyperbola ***'
    WRITE(*,'(3(A,f8.3))') 'Hyperbola k=',k,' fm. epsilon=',epsilon, &
         ' Asym. angle=',pmax
    WRITE(*,'(3(A,F8.3))') 'Angle at ',initial_roft,' fm: ',p,' Deflection: ', &
         initial_angle
    ! get initial and final relative vectors in TDHF, calculate angle
    init_relvec=init_relvec/SQRT(SUM(init_relvec**2))
    final_relvec=cmfrag(:,2)-cmfrag(:,1)
    roft=SQRT(SUM(final_relvec**2))
    tdhf_angle=ACOS(SUM(init_relvec*final_relvec)/roft)*180.D0/pi
    WRITE(*,'(A,f8.3)') 'TDHF rotation angle: ',tdhf_angle
    !  
    ! Calculate relative motion quantities-
    ! reduced mass
    mu=nucleon_mass*mass(1)*mass(2)/(mass(1)+mass(2)) ! reduced mass
    ! kinetic energy of relative motion
    vrel=SQRT(SUM((velocity(:,1)-velocity(:,2))**2))
    erel=0.5D0*mu*vrel**2
    ! Coulomb energy in point charge approximation
    ecoul_point=e2*charge(1)*charge(2)/roft
    ! angular momentum of relative motion
    relangmom=0.D0
    DO i=1,2
       relangmom(1)=relangmom(1)+nucleon_mass*mass(i)* &
            (cmfrag(2,i)*velocity(3,i)-cmfrag(3,i)*velocity(2,i))/hbc
       relangmom(2)=relangmom(2)+nucleon_mass*mass(i)* &
            (cmfrag(3,i)*velocity(1,i)-cmfrag(1,i)*velocity(3,i))/hbc
       relangmom(3)=relangmom(3)+nucleon_mass*mass(i)* &
            (cmfrag(1,i)*velocity(2,i)-cmfrag(2,i)*velocity(1,i))/hbc
    END DO
    ! rotational energy assuming fragments as point masses
    e_centrifugal=hbc**2*SUM(relangmom**2)/(2.D0*mu*roft**2)
    ! velocity and kinetic energy of radial motion
    rdot=SUM((velocity(:,2)-velocity(:,1))*(cmfrag(:,2)-cmfrag(:,1)))/roft
    edot=0.5D0*mu*rdot**2
    WRITE(*,*) '*** Data for final configuration ***'
    WRITE(*,'(2(A,F10.5))') 'Fragment distance ',roft,' relative velocity',rdot
    WRITE(*,'(A,3F10.5)') 'Relative angular momentum ',relangmom
    WRITE(*,'(3(A,F10.5))') 'Radial kinetic energy ',edot, &
         ' centrifugal energy',e_centrifugal,' point Coulomb energy ',ecoul_point
    WRITE(*,'(2(A,F10.5))') 'Relative kinetic energy',erel,' edot+e_centrifugal', &
         edot+e_centrifugal
    ! Coulomb energies with more refined calculation
    CALL coulinit
    rhof=rho(:,:,:,2)
    CALL poisson(rhof)
    ecoul(0)=0.5D0*wxyz*SUM(rhof*wcoul)
    DO ifrag=1,2
       WHERE(frag==ifrag)
         rhof=rho(:,:,:,2)
       ELSEWHERE
         rhof=0.D0
       END WHERE
       CALL poisson(rhof)
       ecoul(ifrag)=0.5D0*wxyz*SUM(rhof*wcoul)
    END DO
    ecoul_int=ecoul(0)-SUM(ecoul(1:2))
    WRITE(*,'(3(A,F10.4))') 'Total Coulomb energy: ',ecoul(0), &
         ' Fragment 1: ',ecoul(1),' Fragment 2:',ecoul(2), &
         ' Coulomb interaction energy:',ecoul_int
    ! 
    final_ecm=erel+ecoul_int
    ! now calculate final Rutherford scattering angle
    k=(relangmom(2)*hbc)**2/(e2*charge(1)*charge(2)*mu)
    aa=e2*charge(1)*charge(2)/(2.D0*final_ecm)
    bb=relangmom(2)*hbc/SQRT(2.D0*mu*final_ecm)
    epsilon=SQRT(aa**2+bb**2)/a
    pmax=ACOS(1.D0/epsilon)*180.D0/pi
    roft=SQRT(SUM(final_relvec**2))
    p=ACOS((k/roft+1.D0)/epsilon)*180.D0/pi
    final_angle=pmax-p
    WRITE(*,*) '*** Data for outgoing Rutherford trajectory ***'
    WRITE(*,'(3(A,f8.3))') 'Hyperbola k=',k,' fm. epsilon=',epsilon, &
         ' Asym. angle=',pmax
    WRITE(*,'(3(A,F8.3))') 'Angle at ',roft,' fm: ',p,' Deflection: ', &
         final_angle
    ! produce summary
    WRITE(*,*) '***************** Summary *************************'
    WRITE(*,'(2(A20,2F10.3))') 'Initial masses:',initial_mass, &
         ' Initial charges:',initial_charge
    WRITE(*,'(2(A20,2F10.3))') 'Final masses:',mass(1:2), &
         ' Final charges:',charge(1:2)
    WRITE(*,'(2(A20,F10.3))') 'Initial c.m. energy: ',initial_ecm, &
         ' Final: ',final_ecm
    WRITE(*,'(3(A15,F10.3))') 'Radial: ',edot,' Centrifugal:', &
         e_centrifugal,' Coulomb:',ecoul_int
    WRITE(*,'(2(A20,F10.3))') 'Initial L: ',initial_L,' Final orbital L: ',& 
         relangmom(2)
    WRITE(*,'(A20,2F10.3,A20,F10.3)') 'Fragment L: ',angmom(2,1:2), &
         'Final L_tot: ',relangmom(2)+angmom(2,1)+angmom(2,2)
    WRITE(*,'(A20,F10.3)') 'Scattering angle: ', &
         180.D0-initial_angle-tdhf_angle-final_angle
    WRITE(*,'(3(A15,F10.3))') 'Initial: ',initial_angle, &
         ' TDHF part: ',tdhf_angle,' Final: ',final_angle
  END SUBROUTINE scattering
!--------------------------------------------------------------------------- 
! DESCRIPTION: cut_system
!> @brief
!! This subroutine finds out whether the system can be decomposed into two
!! fragments records the geometry of this decomposition.
!>
!> @details
!! There is an iteration using index \c itcm to make the line connecting
!! the centers of the fragments consistent. It starts with the quadrupole
!! tensor principal axis using function \c getslope. After calculating the
!! fragment centers-of-mass, this is corrected iteratively. Within the loop
!! the following steps are executed:
!! -# If this is not the first iteration, calculate the slope of the 
!! center-connecting line from the centers of mass.
!! -# Call \c divpoint to find whether there is a sufficiently low density
!! between two fragments along that line. If there is, keep track of the 
!! minimum-density location.
!! -# Taking the line from that point, calculate the properties of the 
!! line perpendicular to the inter-center connection: the line dividing space 
!! between the fragments.
!! -# Then assign fragment index \c frag(ix,iy,iz) for each grid point by
!! checking which side of the dividing line they are on. At the same time, 
!! accumulate data for the fragment masses and centers of mass.
!! -# If the centers of mass have not changed significantly, stop the
!! iterations. Otherwise repeat.
!! -# Record the center distance and return if this is not a full analysis.
!! -# Finally, make sure that the assignment of fragment index to heavier
!! or lighter fragment remains the same. This uses array \c swap to interchange 
!! 1 and 2.
!>
!> @param[in] full_analysis
!! LOGICAL, if false, do only what is necessary for calculating the
!! distance
!--------------------------------------------------------------------------- 
  SUBROUTINE cut_system(full_analysis)
    LOGICAL,INTENT(IN) :: full_analysis
    REAL(db) :: cent(3,2),center(3,2),angle,slopev,xx,zz,diff,vol, &
         bb,xmin,zmin,slope
    INTEGER :: ix,iy,iz,itcm,ifrag
    INTEGER,PARAMETER :: swap(0:2)=(/ 0,2,1 /)
    Iteration: DO itcm=1,10 
       ! calculate slope from centers of mass
       IF(itcm>1) THEN   
         cent=center 
         slope=(cent(3,2)-cent(3,1))/(cent(1,2)-cent(1,1))   
         bb=cent(3,1)-slope*cent(1,1)   
       ELSE   ! first iteration
         cent=0.D0
         slope=getslope()
         bb=0.D0
       ENDIF
       ! check whether separation
       istwobody=divpoint(xmin,zmin,bb,slope) 
       ! determine dividing line: slope and intercept
       angle=ATAN(slope)   
       slopev=dtan(angle+pi/2.0D0)   
       bb=zmin-slopev*xmin   
       ! assign fragment cells
       center(:,1:2)=0.D0
       mass(1:2)=0.0D0   
       charge(1:2)=0.0D0
       DO iz=1,nz   
          zz=z(iz)   
          DO ix=1,nx   
             xx=x(ix)   
             diff=zz-slopev*xx-bb   
             DO iy=1,ny   
                vol=wxyz   
                frag(ix,iy,iz)=0
                IF(rhotot(ix,iy,iz)<threshold) CYCLE
                IF(diff<0.0D0) THEN
                  frag(ix,iy,iz)=1
                  ifrag=1 
                ELSE 
                  frag(ix,iy,iz)=2
                  ifrag=2 
                ENDIF
                mass(ifrag)=mass(ifrag)+vol*rhotot(ix,iy,iz)
                center(1,ifrag)=center(1,ifrag)+vol*rhotot(ix,iy,iz)*x(ix)   
                center(2,ifrag)=center(2,ifrag)+vol*rhotot(ix,iy,iz)*y(iy)   
                center(3,ifrag)=center(3,ifrag)+vol*rhotot(ix,iy,iz)*z(iz)
             ENDDO
          ENDDO
       ENDDO
       DO ifrag=1,2
          center(:,ifrag)=center(:,ifrag)/mass(ifrag)
       END DO
       ! stop if this is not a separated case
       IF(.NOT.istwobody) EXIT 
       ! end iterations if the c.m. remains constant
       IF(MAX(MAXVAL(ABS(cent(1,:)-center(1,:))), &
            MAXVAL(ABS(cent(3,:)-center(3,:))))<1.0d-05) EXIT
    ENDDO Iteration
    ! record separation distance
    IF(istwobody) THEN
       roft_old=roft
       roft=SQRT(SUM(center(:,1)-cent(:,2))**2)
       WRITE(*,'(A,F10.2,A,F8.4)') 'Time:',time,' Center distance: ',roft
    ELSE
       roft=0.D0
    END IF
    IF(.NOT.full_analysis) RETURN
    ! correct ordering of fragments to agree with initialization
    SELECT CASE(mass_case)
    CASE(1)
      IF(mass(1)<mass(2)) RETURN
    CASE(2)
      IF(mass(1)>mass(2)) RETURN
    END SELECT
    ! incorrect case: need to exchange indices
    DO iz=1,nz   
       DO iy=1,ny   
          DO ix=1,nx   
             frag(ix,iy,iz)=swap(frag(ix,iy,iz))
          END DO
       END DO
    END DO
  END SUBROUTINE cut_system
!---------------------------------------------------------------------------  
! DESCRIPTION: getslope
!> @brief
!!In this function the slope of the fragment connection line in the 
!! \f$x-z \f$ plane is determined. 
!>
!> @details
!! The slope is determined by the
!! eigenvector of largest quadrupole moment in the \f$(x,z)\f$-plane. The
!! first loop sums up the quadrupole tensor, which is dimensioned <tt>
!!   q2(3,3)</tt> but of which the index 2 is not actually necessary, since
!! the \f$y\f$-direction is not involved. The definition is kept
!! three-dimensional to reduce confusion and make later generalization
!! easier.
!! The two-dimensional eigenvalue problem for \c q2 has the secular
!! equation:
!! \f[(q_{11}-\lambda)(q_{33}-\lambda)-q_{13}^2=0,\f]
!! with \f$\lambda\f$ the eigenvalue. For the larger eigenvalue we get
!! \f[\lambda=\frac1{2}\left(q_{11}+q_{33}+\sqrt{(q_{11}-q_{33})^2
!!     -4q_{13}^2}\right),\f] and solving the equation
!! \f$q_{13}x+(q_{33}-\lambda)z=0\f$ for \f$z\f$ yields
!! \f{eqnarray*}{
!!   z&=&\frac{q_{13}}{\lambda-q_{33}}\,x\\
!!   &=& \frac{q_{13}}{\frac1{2}\left(q_{11}-q_{33}+\sqrt{(q_{11}-q_{33})^2
!!         -4q_{13}^2}\right)}.
!! \f}
!! The code calls the denominator \c denom and makes sure no division
!! by zero happens (this could happen for a spherical distribution).  The
!! resulting slope is returned in the module variable \c slope.
!>
!> @retval newslope REAL(db), calculated slope.
!--------------------------------------------------------------------------- 
  FUNCTION getslope() RESULT(newslope)
    REAL(db) ::  q2(3,3),xx,yy,zz,vol,denom,newslope
    INTEGER :: ix,iy,iz 
    q2=0.D0 
    DO iz=1,nz   
       zz=z(iz)
       DO iy=1,ny   
          yy=y(iy)
          DO ix=1,nx   
             xx=x(ix)
             vol=wxyz*rhotot(ix,iy,iz)
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
      newslope=100.D0   
    ELSE   
      newslope=q2(1,3)/denom
    ENDIF
  END FUNCTION getslope
!--------------------------------------------------------------------------- 
! DESCRIPTION: divpoint
!> @brief
!! The function divpoint examines the line determined by the axis of 
!! largest quadrupole moment to look for a suitable separation
!! point between two fragments. If one is found with sufficiently low 
!! density, it returns \c TRUE and the properties of that point in its 
!! arguments. Otherwise it returns \c FALSE.
!>
!> @details
!!To this end it looks at the behavior of the densities along this line.
!!Since the line has no relation to the numerical grid, this is not
!!trivial.
!!  - <b> First loop: </b> Essentially it looks through the 
!!  \f$(x,z)\f$-plane to   find points closer than half a grid spacing 
!!  to the desired line with equation <tt>z=slope*x+bb</tt> 
!!  (logical variable \c online). If the slope is larger than one, 
!!  i. e., if the nuclei are separating predominantly in the 
!!  \f$x\f$-direction, we need to take the equation
!!  <tt>x=(z-bb)/slope</tt> instead to get better resolution. To make the
!!  result monotonic along the line, the do loop in \f$z\f$ runs backward
!!  for negative slopes. The points found are collected in index vectors
!!  <tt>ixl</tt>, <tt>izl</tt> with accompanying densities <tt>rhol
!! </tt> stored
!!  in arrays of length <tt>il</tt>.
!!  - <b> Second loop: </b> now the number of fragments <tt>nf</tt>
!!  is counted by examining this one-dimensional density curve, looking for
!!  disconnected density humps above <tt>vacuum</tt> density. Where 
!!  the ``vacuum'' region starts and ends is recorded in variables
!!  <tt>n1</tt> and <tt>n2</tt>. The logic is as follows:
!!    -# The logical variable <tt>in\_vacuum</tt> keeps track of whether
!!    the search is in a vacuum region at the moment or not. It starts
!!    as <tt>TRUE</tt>.
!!    -# Go to the next point. If its density is above <tt>vacuum</tt>,
!!    and we are in the vacuum, a new fragment is starting and we
!!    increase the number of fragments <tt>nf</tt> by 1. If it becomes
!!    bigger than 2, exit, because there are three or more fragments.
!!    If it is now 2, record the starting index for the second fragment
!!    in <tt>n2</tt>.
!!    -# If the density is below <tt>vacuum</tt> and we are not in vacuum,
!!    a fragment is being ended. If it is the first fragment, we record
!!    this index in <tt>n1</tt>.
!!  - <b> Final processing: </b> At this point we expect a two-fragment
!!  situation if <tt>nf=2</tt> and in this case the void region between the
!!  fragments extends from <tt>n1</tt> to <tt>n2</tt>, which are indices into
!!  arrays <tt>ix1</tt> and <tt>iz1</tt> giving the position in the
!!  \f$(x,z)\f$-plane. The code calculates the midpoint between the two
!!  positions and returns <tt>TRUE</tt> in this case.
!>
!> @return 
!> LOGICAL, indicates whether there is a fragment division or not.
!> @param[out] xmin 
!> REAL(db), \f$x\f$-coordinate of the dividing point.
!> @param[out] zmin 
!> REAL(db), \f$z\f$-coordinate of the dividing point.
!> @param[in] bb
!> REAL(db), intercept of the line joining the fragments..
!> @param[in] slope
!> REAL(db), slope of the line joining the fragments.
!--------------------------------------------------------------------------- 
  LOGICAL FUNCTION divpoint(xmin,zmin,bb,slope) 
    LOGICAL :: online,two,in_vacuum 
    INTEGER :: iyy,iz1,iz2,idz,il,ix,iz,ixl(nx+ny+nz), &
         izl(nx+ny+nz),i,n1,n2,nf 
    REAL(db) :: deltax,deltaz,xx,zz,rhol(nx+ny+nz) 
    REAL(db),INTENT(OUT) :: xmin,zmin
    REAL(db),INTENT(IN) :: bb,slope
    REAL(db),PARAMETER :: vacuum=0.03D0
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
              WRITE(*,*) ' Increase dimensioning in function divpoint' 
              STOP 
            END IF
            ixl(il)=ix; izl(il)=iz 
            rhol(il)=rhotot(ix,iyy,iz)
          END IF
       END DO
    END DO
    ! Loop 2
    nf=0; in_vacuum=.TRUE. 
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
  !***************************************************
END MODULE Twobody
