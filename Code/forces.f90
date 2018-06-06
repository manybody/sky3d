MODULE Forces
  USE Params, ONLY: db,scratch,wflag, hbc, e2, pi
  IMPLICIT NONE
  SAVE
  ! Record defining pairing
  TYPE Pairing
     REAL(db) :: v0prot,v0neut,rho0pr
  END TYPE Pairing
  ! Record defining Skyrme force
  TYPE Force
     CHARACTER(8) :: name
     INTEGER :: ex,zpe
     REAL(db) :: h2m(2)
     REAL(db) :: t0,t1,t2,t3,t4
     REAL(db) :: x0,x1,x2,x3,b4p
     REAL(db) :: power
     TYPE(Pairing) :: vdi  ! volume-delta  
     TYPE(Pairing) :: dddi ! density-dependent delta
  END TYPE Force
  ! include predefined forces
  INCLUDE 'forces.data'
  ! now the structure used in the run itself
  INTEGER :: ipair
  TYPE(Force) :: f         ! force actually used
  TYPE(Pairing) :: p       !pairing parameters actually used
  ! charge and mass number in static case for pairing
  REAL(db) :: h2ma ! average h2m over p and n, used in some places
  REAL(db) :: nucleon_mass
  ! derived "b" and Slater coefficients
  REAL(db) :: b0,b0p,b1,b1p,b2,b2p,b3,b3p,b4,b4p,slate
  ! derived "C" coefficients
  REAL(db) :: Crho0,Crho1,Crho0D,Crho1D,Cdrho0,Cdrho1,Ctau0,Ctau1,CdJ0,CdJ1
CONTAINS
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
         ipair,v0prot,v0neut,rho0pr,turnoff_zpe
    ! mark force & pairing parameters as undefined
    h2m=-1.0; v0prot=-1.0; v0neut=-1.0; 
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
       f%zpe=-1
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
       p%v0prot=v0prot
       p%v0neut=v0neut
       IF(ipair==6) p%rho0pr=rho0pr
       ! get predefined pairing if applicable
       IF(predefined) THEN
          IF(ipair==5) p=f%vdi
          IF(ipair==6) p=f%dddi
       ENDIF
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
    ENDIF
  END SUBROUTINE read_force
END MODULE Forces
