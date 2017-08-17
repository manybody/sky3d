!------------------------------------------------------------------------------
! MODULE: External
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module allows the coupling of the nucleonic wave functions to an
!!external field. 
!>
!>@details 
!!This can be done either by adding a time-dependent external (i.e., not
!!self-consistent) potential to the single-particle Hamiltonian. 
!!The time-dependence can be of Gaussian form
!!\f[ f(t)=\exp\left(-(t-\tau_0)^2/\Delta\tau^2\right)\cos(\omega(\tau-\tau_0)) \f]
!!or cosine squared form
!!\f[  f(t)=\cos\left(\frac\pi{2}\left(\frac{t-\tau_0}{\Delta\tau}\right)^2\right)
!!   \theta\left(\Delta\tau-|t-\tau_0|\right)\cos(\omega(\tau-\tau_0)) \f]
!!or by
!!giving an initial "boost" 
!!\f[ \psi_k(\vec r,s,t\!=\!0)=
!!\psi_{k,0}(\vec r,s)\,\exp\left(-\I \eta F_q(\vec r)\right) \f]
!!to each wave function. Since this is very
!!easy to modify and will probably have to be adjusted for most
!!applications, the present version just contains a sample for a
!!quadrupole coupling. The logic for different time-dependence
!!assumptions and the isospin are however, fully functional and should
!!be useful in many cases.
!!TODO: What is tau??
!------------------------------------------------------------------------------
MODULE External
  USE Params, ONLY: db,pi,wflag,extfieldfile,time,scratch
  USE Parallel,ONLY: globalindex
  USE Grids, ONLY: nx,ny,nz,x,y,z,dx,dy,dz,wxyz
  USE Levels, ONLY: nstloc,isospin,charge_number,mass_number
  USE MEANFIELD, ONLY: upot
  IMPLICIT NONE  
  INTEGER,PRIVATE :: isoext=0        !< isospin behavior of the external field: \c isoext
                                     !! denotes the same action on protons and neutrons, 
                                     !! <tt> isoext=1 </tt> that with opposing signs.
  INTEGER,PRIVATE :: ipulse=0        !< type of pulse. <tt> ipulse=0</tt> denotes the
                                     !! initial boost configuration, <tt> ipulse=1 </tt> a Gaussian pulse
                                     !! and <tt> ipulse=2 </tt> a cosine squared behavior.
  REAL(db),PRIVATE,ALLOCATABLE,DIMENSION(:,:,:,:) :: extfield !< this is the time-independent field
                                     !! generated according to the parameters \c amplq0, \c textfield_periodic
                                     !! and depending on the latter, possibly \c radext
                                     !! and \c widext. It is used either to calculate the initial boost
                                     !! or is added to the mean field multiplied with the time-dependent factor.
  REAL(db),PRIVATE :: amplq0=0.D0    !< a strength parameter for the perturbation.
                                     !! Its numerical magnitude is usually not important by itself, 
                                     !! but varying it allows studying the effects of different 
                                     !! strengths of the excitation.
  REAL(db),PRIVATE :: radext=100.D0  !< parameter \f$r_0\f$ for the cutoff of the field in fm.
  REAL(db),PRIVATE :: widext=1.D0    !< parameter \f$\Delta_0\f$ for the cutoff of the field in fm.
  REAL(db),PRIVATE :: omega=0.D0     !< parameter \f$ \omega \f$ for time-dependent excitation
  REAL(db),PRIVATE :: taut           !< parameter \f$ \Delta\tau \f$ for time-dependent excitation
  REAL(db),PRIVATE :: tau0           !< parameter \f$ \tau_0 \f$ for time-dependent excitation
  LOGICAL,PRIVATE :: textfield_periodic=.true. !< if this is set to \c true, the
                                     !! external field is made periodic by substituting 
                                     !! \f$ x^2\rightarrow\sin^2\left(\pi x/x_L\right) \f$, 
                                     !! otherwise a damping factor is used \f$ F_q(\vec r)
                                     !! \rightarrow \frac{F_q(\vec r)}{1+e^{(r-r_0)/\Delta r}} \f$.
  SAVE
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: getin_external
!> @brief
!!It reads in all the parameters of the external field.
!--------------------------------------------------------------------------- 
  SUBROUTINE getin_external
    NAMELIST/extern/ amplq0,radext,widext,isoext,ipulse,omega,tau0,taut, &
         textfield_periodic
    READ(5,extern)
  END SUBROUTINE getin_external
!---------------------------------------------------------------------------  
! DESCRIPTION: init_external
!> @brief
!!It does some consistency checks of the external parameters and initializes the external field. 
!>
!> @details
!!The relative strength for the neutron and proton fields is set equal for
!!<tt> isoext=0 </tt> but reduced by the corresponding number of particles in
!!the <tt> isoext=1 </tt>. This avoids a shift of the center of mass in this
!!case.
!!
!!Then the array \c extfield is allocated and the time-independent
!!spatial potential calculated for both isospin cases and in either the
!!periodic or damped versions, depending on the value of 
!!\c textfield_periodic.This is just an illustrative sample field of
!!type \f$ Q_{zz} \f$; in a real calculation there is no need to provide both
!!versions.
!!TODO: This routine is not documented correctly in the "old" pdf manual!!!
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: extfld
!> @brief
!!This is again a very straightforward routine. It calculates the
!!time-dependent prefactor \c time_factor depending on the
!!parameters, and adds the time-independent field \c extfield
!!multiplied by this factor to \c upot. The physical time is here
!!given as an argument, because it needs to be evaluated at both full
!!and half time steps.
!>
!> @param[in] time
!> REAL(db), takes the current time.
!--------------------------------------------------------------------------- 
  SUBROUTINE extfld(time)  
    REAL(db),INTENT(IN) :: time
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
!---------------------------------------------------------------------------  
! DESCRIPTION: extboost
!> @brief
!!This performs the initial boost on all single-particle wave functions. 
!>
!> @details
!!Note that it is always called
!!from \c dynamichf and checks itself whether \c ipulse is zero.
!!This enables \c ipulse to also be made a private variable. Its
!!argument is used to communicate whether it has actually done anything:
!!if it has applied a boost, it sets its argument to \c .FALSE.  so that in
!!\c dynamichf the variable \c text_timedep makes it possible to
!!distinguish this case.
!>
!> @param[out] noboostflag
!> LOGICAL, returns if boost has to be performed.
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: print_extfield
!> @brief
!!This subroutine calculates the expectation of the coupling energy to
!!the external field,
!!\f[ \sum_q\int\,\D^3r\,\rho_q(\vec r)F_q(\vec r) \f]
!!and prints one line containing the present time and this value onto
!!the file \c extfieldfile.
!--------------------------------------------------------------------------- 
  SUBROUTINE print_extfield()
    USE Densities, ONLY: rho
    OPEN(UNIT=scratch,file=extfieldfile,POSITION='APPEND')  
    WRITE(scratch,'(F12.3,1pg15.7)') time,wxyz*SUM(rho*extfield)/amplq0
    CLOSE(UNIT=scratch)
  END SUBROUTINE print_extfield
END MODULE External
