MODULE External
  USE Params, ONLY: db,pi,wflag,extfieldfile,time,scratch
  USE Parallel,ONLY: globalindex,tabc_av
  USE Grids, ONLY: nx,ny,nz,x,y,z,dx,dy,dz,wxyz
  USE Levels, ONLY: nstloc,isospin,nprot,nneut,charge_number,mass_number
  USE MEANFIELD, ONLY: upot
  USE Densities, ONLY: rho
  USE Energies, ONLY: e_extern
  IMPLICIT NONE  
  INTEGER :: ipulse=0
  INTEGER,PRIVATE :: isoext=0
  REAL(db),PRIVATE,ALLOCATABLE,DIMENSION(:,:,:,:) :: extfield
  REAL(db) :: amplq0=0.D0,amplx=0.D0,amply=0.D0,amplz=0.D0,amplrod=0.0D0,&
                      radext=100.D0,widext=1.D0,tau0,taut,omega=0.D0
  LOGICAL,PRIVATE :: textfield_periodic=.true.
  SAVE
CONTAINS
  !***********************************************************************
  SUBROUTINE getin_external
    NAMELIST/extern/ amplq0,amplx,amply,amplz,radext,widext,isoext,&
                     ipulse,omega,tau0,taut,textfield_periodic,amplrod
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
       WRITE(*,"(a,1pg12.4)") " Amplitude of axial quadrupole =",amplq0  
       WRITE(*,"(a,3(1pg12.4))") " Amplitudes Cartesian dipoles  =",&
               amplx,amply,amplz
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
                     -SIN(x(ix)*PI/xlim)**2-SIN(y(iy)*PI/ylim)**2) &
                    +amplx*SIN(x(ix)*2D0*PI/xlim) &
                    +amply*SIN(y(iy)*2D0*PI/ylim) &
                    +amplz*SIN(z(iz)*2D0*PI/zlim) &
                    +amplrod*x(ix)*SIN(2.0d0*z(iz)*PI/(REAL(nz)*dz))
             ELSE                              ! damped version
               facr=(amplq0 *(2.D0*z(iz)**2-x(ix)**2-y(iy)**2) &
                     +amplx*x(ix)+amply*x(iy)+amplz*x(iz)) &
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
!   variables for computing energy absorbed from external field
    REAL(db),SAVE :: time_factor_old=0D0
    REAL(db),SAVE :: upot_ext_int,upot_ext_int_old=0D0
!
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

    upot_ext_int = wxyz*SUM(extfield*rho)
    e_extern=e_extern+(upot_ext_int+upot_ext_int_old)*&
                      (time_factor-time_factor_old)
    upot_ext_int_old = upot_ext_int
    time_factor_old = time_factor

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
    WRITE(*,*)MAXVAL(extfield)
  END SUBROUTINE extboost
  !***********************************************************************
  SUBROUTINE print_extfield()
    USE Densities, ONLY: rho
    INTEGER :: ix,iz
    REAL(db):: xcm(nz),xcmnorm(nz)
    OPEN(UNIT=scratch,file=extfieldfile,POSITION='APPEND')  
    IF(abs(amplrod)<1D-20) THEN
      WRITE(scratch,'(F12.3,4(1pg15.7))') time,wxyz*SUM(rho*extfield),&
         SUM(rho),SUM(rho**2)
    ELSE
      xcm=0.0d0
      xcmnorm=0.0d0
      DO ix=1,nx; DO iz=1,nz
        xcm(iz)=xcm(iz)+x(ix)*sum(rho(ix,:,iz,:))
        xcmnorm(iz)=xcmnorm(iz)+sum(rho(ix,:,iz,:))
      END DO; END DO
      xcm=xcm/xcmnorm
      WRITE(scratch,*) time,wxyz*SUM(rho*extfield),xcm
    END IF
    CLOSE(UNIT=scratch)
  END SUBROUTINE print_extfield
  !***********************************************************************
  FUNCTION tabc_extfield()
    REAL(db)            :: tabc_extfield
    WRITE(*,*)'----------',wxyz,SUM(rho*extfield),tabc_extfield
    tabc_extfield=tabc_av(wxyz*SUM(rho*extfield))
  END FUNCTION tabc_extfield
END MODULE External
