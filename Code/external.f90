MODULE External
  USE Params, ONLY: db,pi,wflag,extfieldfile,time,scratch
  USE Parallel,ONLY: globalindex
  USE Grids, ONLY: nx,ny,nz,x,y,z,dx,dy,dz,wxyz
  USE Levels, ONLY: nstloc,isospin,charge_number,mass_number
  USE MEANFIELD, ONLY: upot
  IMPLICIT NONE  
  INTEGER,PRIVATE :: isoext=0,ipulse=0
  REAL(db),PRIVATE,ALLOCATABLE,DIMENSION(:,:,:,:) :: extfield
  REAL(db),PRIVATE :: amplq0=0.D0,radext=100.D0,widext=1.D0, &
       tau0,taut,omega=0.D0
  LOGICAL,PRIVATE :: textfield_periodic=.true.
  SAVE
CONTAINS
  !***********************************************************************
  SUBROUTINE getin_external
    NAMELIST/extern/ amplq0,radext,widext,isoext,ipulse,omega,tau0,taut, &
         textfield_periodic
    READ(5,extern)
  END SUBROUTINE getin_external
  SUBROUTINE init_external
    REAL(db) :: facn,facp,facr,xlim,ylim,zlim
    INTEGER :: ix,iy,iz
    CHARACTER(14),PARAMETER :: pulsetype(0:2)=(/ 'Instantaneous ', &
         'Gaussian      ','Cosine squared' /)
    IF(ipulse<0.OR.ipulse>2) STOP &
         ' External field: called with invalid pulse type'
    IF(wflag) THEN
       WRITE(*,*) "***** Parameters of external field *****"  
       WRITE(*,"(a,e12.4)") " Amplitude of axial quad.   =",amplq0  
       WRITE(*,"(2(A,F10.4),A)") " Radial damping: radius ",radext, &
            ' fm,  width ',widext,' fm'
       WRITE(*,"(2(A,I2),2A)") " Isospin of excitation:",isoext, &
            ' Pulse type: ',ipulse,'=',pulsetype(ipulse)
       IF(ipulse>0) &
            WRITE(*,'(A,2F10.4)') ' Peak time and width in fm/c:', &
            tau0,taut
       WRITE(*,"(A,F12.4,A)") " Modulating frequency omega:",omega, &
            ' c/fm'
       WRITE(*,*) 'Periodic (T) or damped (F) external field:',textfield_periodic
    ENDIF
    IF(isoext<0.OR.isoext>1) STOP " INIEXT: called with invalid isospin"
    IF(ipulse==2.AND.tau0<taut) STOP &
         ' External field: tau0<taut is nonsense for cos**2 pulse '
    IF(ipulse==1.AND.tau0<3*taut) STOP &
         ' External field: tau0<3*taut is nonsense for Gaussian pulse '
    IF(isoext==0) THEN  
       facn=1.0D0  
       facp=1.0D0  
    ELSE  
       facn=-1.0D0/(mass_number-charge_number)  
       facp=1.0D0/charge_number  
    ENDIF
    WRITE(*,*) 'EXTERNAL: ',facn,facp
    ALLOCATE(extfield(nx,ny,nz,2))
    xlim=nx*dx
    ylim=ny*dy
    zlim=nz*dz
    DO iz=1,nz  
       DO iy=1,ny  
          DO ix=1,nx  
             IF(textfield_periodic) THEN       ! strictly periodic version
                facr=amplq0 *(2.D0*SIN(z(iz)*PI/zlim)**2 &
                     -SIN(x(ix)*PI/xlim)**2-SIN(y(iy)*PI/ylim)**2) 
             ELSE                              ! damped version
                facr=amplq0 *(2.D0*z(iz)**2-x(ix)**2-y(iy)**2) &
                     /(1.0D0+EXP((SQRT(x(ix)**2+y(iy)**2+z(iz)**2)-radext)/widext))
             END IF
             extfield(ix,iy,iz,1)=facr*facn  
             extfield(ix,iy,iz,2)=facr*facp  
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE init_external
  !***********************************************************************
  SUBROUTINE extfld(time)  
    REAL(db) :: time
    INTENT(IN) :: time
    REAL(db) :: time_factor
    IF(ipulse==1) THEN  
       time_factor=EXP(-((time-tau0)/taut) **2)  
    ELSE  
       IF(time<tau0-taut) THEN  
          time_factor=0.0D0  
       ELSEIF(time>tau0+taut) THEN  
          time_factor=0.0D0  
       ELSE  
          time_factor=COS(0.5D0*pi *(time-tau0)/taut) **2  
       ENDIF
    ENDIF
    IF(omega/=0.0D00) THEN
       time_factor=COS(omega*(time-tau0))*time_factor
    ENDIF
    upot=time_factor*extfield + upot
  END SUBROUTINE extfld
  !***********************************************************************
  SUBROUTINE extboost(noboostflag)
    USE Levels, ONLY: psi
    LOGICAL,INTENT(OUT) :: noboostflag
    INTEGER :: nst,is
    noboostflag=ipulse/=0
    IF(noboostflag.OR.time>0.D0) RETURN
    FORALL(nst=1:nstloc,is=1:2)
       psi(:,:,:,is,nst)=psi(:,:,:,is,nst) &
            *EXP(CMPLX(0.0D0, &
            -extfield(:,:,:,isospin(globalindex(nst))),db))
    END FORALL
  END SUBROUTINE extboost
  !***********************************************************************
  SUBROUTINE print_extfield()
    USE Densities, ONLY: rho
    OPEN(UNIT=scratch,file=extfieldfile,POSITION='APPEND')  
    WRITE(scratch,'(F12.3,1pg15.7)') time,wxyz*SUM(rho*extfield)
    CLOSE(UNIT=scratch)
  END SUBROUTINE print_extfield
END MODULE External
