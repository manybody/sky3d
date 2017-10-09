!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!Module \c Forces describes the interactions used in the code. The
!!idea is to produce a library of Skyrme forces that can be called up
!!simply by name, but for exploratory purposes a force can also be input
!!using individual parameter values. A force is usually fitted together
!!with a prescription for the pairing and center-of-mass correction, so
!!that these properties are here defined as part of the force.
!!
!!The predefined Skyrme forces are contained in \c forces.data, which
!!contains an array \c pforce of <tt> TYPE(Force) </tt> data. This is a
!!bit unreliable since the Fortran standard restricts the length of
!!statements; it should be replaced by reading from a data file in case
!!this limit causes problems. The present version is, however, still preferred as
!!a data file would have to be replicated in every application directory
!!(or an absolute path would have to be defined in the \c OPEN statement).
!------------------------------------------------------------------------------
MODULE Forces
  USE Params, ONLY: db,scratch,wflag, hbc, e2, pi
  IMPLICIT NONE
  SAVE
  !> This variable contains the parameters for pairing:
  TYPE Pairing
     REAL(db) :: v0prot!<the strength of pairing for protons in MeV.
     REAL(db) :: v0neut!<the strength of pairing for neutrons in MeV.
     REAL(db) :: rho0pr!<the density parameter for the density-dependent delta pairing
  END TYPE Pairing
  !> This variable contains the parameters for the Skyrme force:
  TYPE Force
     CHARACTER(8) :: name         !<the name of the force used to identify it.
     INTEGER :: ex                !<some forces are fitted excluding the Coulomb
     !!exchange term. For <tt> ex=1 </tt> it is included (this is the normal
     !!case), for <tt> ex=0 </tt> not.
     INTEGER :: zpe               !<index for the treatment of the center-of-mass correction
     REAL(db) :: h2m(2)           !<value of \f$ \frac{\hbar^2}{2m} \f$, separately 
     !!for neutrons and protons.
     !>@name Skyrme parameters. 
     !>@{
     REAL(db) :: t0,t1,t2,t3,t4
     !>@}
     !>@name Skyrme exchange parameters.
     !>@{
     REAL(db) :: x0,x1,x2,x3,b4p
     !>@}
     REAL(db) :: power            !< exponent in the nonlinear (originally three-body) term.
     TYPE(Pairing) :: vdi         !< parameter set for the volume-delta pairing case. 
     TYPE(Pairing) :: dddi        !< parameter set for the density-dependent delta pairing case.
  END TYPE Force
  ! include predefined forces
  INCLUDE 'forces.data'
  ! now the structure used in the run itself
  INTEGER :: ipair         !< selects one of several pairing modes.  For
  !!historical reasons the values are 0: no pairing, 5: VDI pairing, and
  !!6: DDDI pairing. In the input the symbolic names are used so these
  !!numerical values are hidden to the user. For details see the input
  !!description and module \c Pairs.   
  LOGICAL :: pair_reg=.FALSE.      !<Enables pairing regularization
  REAL(db):: delta_fit(2)=-1.0d0   !<Allows for fitting to pairing gaps
  REAL(db):: pair_cutoff(2)=-1.0d0 !<Selects a hard pairing cutoff
  REAL(db):: cutoff_factor=0D0  !<The relative "margin" of extra particles defining pairing space.
  REAL(db):: ecut_stab=0D0      ! E_cut for stabilized pairing  !!! PGR
  !!! PGR: ecut_stab to be included in 'TYPE pairing'
  TYPE(Force) :: f         !< this contains parameters for
  !!the Skyrme force actually used in the present calculation, packed
  !into the derived-type \c Force.
  TYPE(Pairing) :: p       !< the pairing parameters used in the
  !!present calculation. This is separate from the force itself: the
  !!force definition usually contains suggestions for the associated
  !!pairing, but this often overridden, e.g, by turning off
  !!pairing.
  ! charge and mass number in static case for pairing
  REAL(db) :: h2ma         !< the average of the two \c h2m values for protons and neutrons.
  REAL(db) :: nucleon_mass !<  the mass of the nucleon (average of
  !!neutron and Proton) in MeV calculated from \c h2ma and \c hbc.
  !
  !>@name these are the coefficients actually used 
  !!for the mean-field and single-particle Hamiltonian calculations
  !!in \c skyrme and \c integ_energy. Note that only \c b4p is also included 
  !!in the Skyrme-force definition; the others are derived from the \c t coefficients.
  !>@{
  REAL(db) :: b0,b0p,b1,b1p,b2,b2p,b3,b3p,b4,b4p,slate
  !>@}
  !>@name derived "C" coefficients
  !>@{
  REAL(db) :: Crho0,Crho1,Crho0D,Crho1D,Cdrho0,Cdrho1,Ctau0,Ctau1,CdJ0,CdJ1
  !>@}
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: read_force
!> @brief
!!The purpose of the subroutine is to read the force and pairing
!!definitions. The \c NAMELIST} \c force contains all the
!!defining values for a Skyrme force (but now as individual variables,
!!not in a derived type) plus a selection of pairing type and strengths
!!(see input description). In addition there is a logical variable
!!\c turnoff_zpe which allows turning off the center-of-mass correction.
!>
!> @details
!!Some quantities are first set negative to see whether required input
!!is missing. Then a predefined Skyrme force with the given name is
!!sought; if it is found, it is simply copied into \c f. If no force
!!of this name is found, a new one is composed from the numbers given in
!!the input.
!!
!!If the input varable \c turnoff_zpe is true, the indicator <tt>
!!f\%zpe </tt> is set to 1, which implies not doing anything about the
!!center-of-mass correction. Actually, in the current version this
!!affects only one statement in the static module.
!!
!!Now the "b" coefficients are
!!calculated straightforwardly.  There is one additional coefficient
!!\c slate used for the Slater approximation to the Coulomb exchange
!!term. This is not a free parameter, but precomputed for convenience in
!!\c forces.f90.  The \f$ b \f$ and \f$ b' \f$ coefficients are used in
!!subroutine \c skyrme (module \c Meanfield) and \c energy
!!(module \c Energies).
!!
!!The variable \c pairing from the namelist then determines the
!!pairing. If it is set to \c 'NONE', no pairing is
!!included. Otherwise the strength parameters are taken from the input
!!or from the predefined force. If this process does not find a
!!reasonable pairing combination, stop with an error message.
!!
!!Finally the routine calculates the values of \c nucleon_mass and
!!\c h2ma. It then prints out a description of the force and pairing
!!parameters.
!--------------------------------------------------------------------------- 
  SUBROUTINE read_force
    CHARACTER(8) :: name,pairing
    INTEGER :: ex,zpe
    REAL(db) :: h2m(2)
    REAL(db) :: t0,t1,t2,t3,t4
    REAL(db) :: x0,x1,x2,x3
    REAL(db) :: power
    REAL(db) :: v0prot,v0neut,rho0pr
    INTEGER :: i
    LOGICAL :: predefined
    LOGICAL :: turnoff_zpe=.FALSE.
    ! read force definition
    NAMELIST /force/ name,pairing, &
         ex,zpe,h2m,t0,t1,t2,t3,t4,x0,x1,x2,x3,b4p,power, &
         ipair,v0prot,v0neut,rho0pr,turnoff_zpe,pair_reg,delta_fit,&
         pair_cutoff,cutoff_factor,ecut_stab
    ! mark force & pairing parameters as undefined
    h2m=-1.0; v0prot=-1.0; v0neut=-1.0; rho0pr=-1.0
    READ(5,force)
    ! seek for force in predefined ones
    predefined=.FALSE.
    DO i=1,nforce
       IF(TRIM(name)==TRIM(pforces(i)%name)) THEN
          predefined=.TRUE.
          f=pforces(i)
       ENDIF
    END DO
    IF(wflag)WRITE(*,*)
    IF(wflag)WRITE(*,*) '***** Force definition *****'
    IF(predefined) THEN
       IF(wflag) WRITE(*,*) 'Using predefined force ',name
    ELSEIF(h2m(1)<0.D0) THEN
       IF(wflag) WRITE(*,*) ' Force not found: ',name,' and none given in input'
       STOP
    ELSE
       IF(wflag) WRITE(*,*) 'Using new force ',name,' defined in input'
       IF(ABS(b4p)<1.0d-9) b4p=f%t4/2.0D0
       f%name=name
       f%ex=ex; f%zpe=zpe
       f%h2m=h2m
       f%t0=t0; f%t1=t1; f%t2=t2; f%t3=t3; f%t4=t4
       f%x0=x0; f%x1=x1;f%x2=x2;f%x3=x3;f%b4p=b4p;
       f%power=power
    ENDIF
    ! turn off zpe if desired
    IF(turnoff_zpe) THEN
       f%zpe=1
       WRITE(*,*) '***** Zero-point-energy correction turned off'
    END IF
    ! calculate "b" and Slater coefficients
    b0=f%t0*(1.0D0+0.5D0*f%x0)  
    b0p=f%t0*(0.5D0+f%x0)  
    b1=(f%t1+0.5D0*f%x1*f%t1+f%t2+0.5*f%x2*f%t2)/4.0D0  
    b1p=(f%t1*(0.5D0+f%x1)-f%t2*(0.5D0+f%x2))/4.0D0  
    b2=(3.0D0*f%t1*(1.D0+0.5D0*f%x1)-f%t2*(1.D0+0.5D0*f%x2))/8.0D0
    b2p=(3.D0*f%t1*(0.5D0+f%x1)+f%t2*(0.5D0+f%x2))/8.D0
    b3=f%t3*(1.D0+0.5D0*f%x3)/4.D0
    b3p=f%t3*(0.5D0+f%x3)/4.D0  
    b4=f%t4/2.D0 
    b4p=f%b4p  
    slate=(3.0D0/pi)**(1.0D0/3.0D0)*e2
    Crho0=0.5d0*b0-0.25d0*b0p
    Crho1=-0.25d0*b0p
    Crho0D=0.3333333333d0*b3-0.1666666666d0*b3p
    Crho1D=-0.1666666666d0*b3p    
    Cdrho0=-0.5d0*b2+0.25d0*b2p
    Cdrho1=0.25*b2p
    Ctau0=b1-0.5d0*b1p
    Ctau1=-0.5d0*b1p
    CdJ0=-b4-0.5d0*b4p
    CdJ1=-0.5d0*b4p
    ! now set up pairing: first case of none
    IF(TRIM(pairing)=='NONE') THEN
       p%v0prot=0.D0; p%v0neut=0.D0; p%rho0pr=0.16D0
       ipair=0
    ELSE
       ! set predefined type of pairing
       IF(TRIM(pairing)=='VDI') ipair=5
       IF(TRIM(pairing)=='DDDI') ipair=6
       !      IF(ipair==6) STOP 'DDDI pairing not implemented in this version'
       IF(ipair==6) p%rho0pr=rho0pr
       ! get predefined pairing if applicable
       IF(predefined) THEN
          IF(ipair==5) p=f%vdi
          IF(ipair==6) p=f%dddi
       ENDIF
       IF(v0prot>0.0d0) p%v0prot=v0prot
       IF(v0neut>0.0d0) p%v0neut=v0neut
       IF(rho0pr>0.0d0) p%rho0pr=rho0pr
       ! stop if this has not yielded meaningful parameters
       IF(p%v0prot<0.D0) STOP 'Pairing not defined properly'
       IF(wflag) WRITE(*,*) 'Using pairing type ',pairing,'=',ipair
    END IF
    ! define average h2m and nucleon mass
    h2ma=0.5D0*(f%h2m(1)+f%h2m(2))
    nucleon_mass=hbc**2/(2.D0*h2ma)
    IF(wflag) THEN
       WRITE(*,"(A)") " Skyrme Potential Parameters"
       WRITE(*,*) 'The force is: ',f%name
       WRITE(*,"(5(A6,F12.5))") "t0",f%t0,"t1",f%t1,"t2",f%t2,"t3",f%t3,"t4",f%t4
       WRITE(*,"(5(A6,F12.5))") "x0",f%x0,"x1",f%x1,"x2",f%x2,"x3",f%x3,"b4p",f%b4p
       WRITE(*,"(A6,F12.5)") "Power",f%power
       WRITE(*,"(A,I2)") " Pairing parameters: Option ipair:",ipair  
       WRITE(*,"(3(A7,F12.5))") "v0prot",p%v0prot,"v0neut",p%v0neut,"rho0pr",p%rho0pr
       WRITE(*,"(A,2F12.5)") "cutoff_factor,ecut_stab=",cutoff_factor,ecut_stab

    ENDIF
  END SUBROUTINE read_force
END MODULE Forces
