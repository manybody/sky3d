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
!!
!!The operator \f$ F_q \f$ is the excitation operator which can be of any form from monopole 
!! to dipole, quadrupole, octupole and so on. However to avoid the spurious states in the strength
!! the definition in the case of isoscalar dipole (L=1) is different than all the other as
!!\f[ F_{1M} = (r^3 - \frac{5}{3}<r^2> r)Y_{1M}  \f]
!! where \f$ <r^2> \f$ is the average value of \f$r^2\f$ given as input (r2_avg) for this case.
!! For all other cases \f$ F_q \f$ is defined as
!! \f[ F_{LM} = \sqrt{2L+1} r^L Y_{LM} \f]
!! where \f$ Y_{LM} \f$ are the spherical harmonics. Only damped version of boundary conditions is implemented
!! in the present work. 
!------------------------------------------------------------------------------
MODULE External
  USE Params, ONLY: db,pi,wflag,extfieldfile,time,scratch
  USE Parallel,ONLY: globalindex
  USE Grids, ONLY: nx,ny,nz,x,y,z,dx,dy,dz,wxyz
  USE Levels, ONLY: nstloc,isospin,charge_number,mass_number
  USE MEANFIELD, ONLY: upot
  Use Spherical_Harmonics
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
  REAL(db),PRIVATE :: radext=100.D0  !< parameter \f$r_0\f$ for the cutoff of the field in fm.
  REAL(db),PRIVATE :: widext=1.D0    !< parameter \f$\Delta_0\f$ for the cutoff of the field in fm.
  REAL(db),PRIVATE :: omega=0.D0     !< parameter \f$ \omega \f$ for time-dependent excitation
  REAL(db),PRIVATE :: taut           !< parameter \f$ \Delta\tau \f$ for time-dependent excitation
  REAL(db),PRIVATE :: tau0           !< parameter \f$ \tau_0 \f$ for time-dependent excitation
  LOGICAL,PRIVATE :: textfield_periodic=.false. !< if this is set to \c true, the
                                     !! external field is made periodic by substituting 
                                     !! \f$ x^2\rightarrow\sin^2\left(\pi x/x_L\right) \f$, 
                                     !! otherwise a damping factor is used \f$ F_q(\vec r)
                                     !! \rightarrow \frac{F_q(\vec r)}{1+e^{(r-r_0)/\Delta r}} \f$.
   REAL(db) :: ampl_ext=0.D0          !< The strength parameter (\f$ \eta \f$) for the perturbation.
                                     !! Its numerical magnitude is usually not important by itself, 
                                     !! but varying it allows studying the effects of different 
                                     !! strengths of the excitation.
   INTEGER :: L_val=0        !< paramter \f$L_{val}\f$, the orbital angular quantum number(L) used
                                     !! in the definition of spherical harmoninc (\f$Y_{L}^{M}\f$) 
                                     !! which is used while defining the external field
   INTEGER :: M_val=0        !< paramter \f$M_{val}\f$, the projection quantum number(M) used
                                     !! in the definition of spherical harmoninc (\f$Y_{L}^{M}\f$) 
                                     !! which is used while defining the external field
   INTEGER,PRIVATE :: only_P=0       !< If this is set to 1, then the strength is calculated by 
                                     !! only exciting the protons to calculate EM response. It should be used carefully
                                     !! because there is no COM correction implemented for this, however, it works fine in the 
                                     !! cases of big doubly magic nuclei like \f$^{208}Pb\f$.  
   REAL(db) :: r2_avg=0.0d0          !< The value of \f$ <r^2> \f$ used in the case of dipole boost.

   PUBLIC ::L_val,M_val,ampl_ext
  SAVE
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: getin_external
!> @brief
!!It reads in all the parameters of the external field.
!--------------------------------------------------------------------------- 
  SUBROUTINE getin_external
     NAMELIST/extern/ ampl_ext,L_val,M_val,radext,widext,isoext,&
                     ipulse,omega,tau0,taut,textfield_periodic,r2_avg,only_P
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
!! Here only the damped version of the boundary conditions
!! are implemented for the external field \f$F_q\f$. 
!--------------------------------------------------------------------------- 
  SUBROUTINE init_external
    REAL(db) :: facn,facp,facr,xlim,ylim,zlim,dip_f
    INTEGER :: ix,iy,iz
    CHARACTER(14),PARAMETER :: pulsetype(0:2)=(/ 'Instantaneous ', &
         'Gaussian      ','Cosine squared' /)
    IF(ipulse<0.OR.ipulse>2) STOP &
         ' External field: called with invalid pulse type'
    IF(wflag) THEN
       WRITE(*,*) "***** Parameters of external field *****"  
       WRITE(*,"(a,e12.4)") " Amplitude of axial quad.   =",ampl_ext 
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
                                                                           ! damped version
               facr = ampl_ext*SQRT(2.0d0*L_val+1.0d0)
               if (L_val .ge. 2)then
                  facr=facr*Y_lm(L_val,M_val,x(ix),y(iy),z(iz))
                  facr=facr*SQRT(x(ix)**2+y(iy)**2+z(iz)**2)**(L_val)
               else if (L_val .eq. 0)then
                  facr=facr*(0.5d0*SQRT(1.0d0/PI))
                  facr=facr*SQRT(x(ix)**2+y(iy)**2+z(iz)**2)**(2)
               else if (L_val .eq. 1)then
                  ! vol=wxyz*rho(ix,iy,iz,iq) 
                  facr=facr*Y_lm(L_val,M_val,x(ix),y(iy),z(iz))
                  dip_f = (r2_avg**2)*5.0d0/3.0d0
                  ! write(*,*)'__DIPOLE FACTOR__',r2_avg,dip_f
                  facr=facr*(SQRT(x(ix)**2+y(iy)**2+z(iz)**2)**3 - dip_f*SQRT(x(ix)**2+y(iy)**2+z(iz)**2))
               end if
               facr=facr/(1.0D0+EXP((SQRT(x(ix)**2+y(iy)**2+z(iz)**2)-radext)/widext)) !< Damping is done using parameters radext and widext
             if (only_P.eq.1)then
               extfield(ix,iy,iz,1)=0.0d0  
               extfield(ix,iy,iz,2)=facr*facp  
             else
               extfield(ix,iy,iz,1)=facr*facn  
               extfield(ix,iy,iz,2)=facr*facp 
             END IF
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
    WRITE(scratch,'(3x,F12.3,3x,F35.15,2x,F25.18,2I5)') time,wxyz*SUM(rho*extfield)/ampl_ext,ampl_ext,L_val,M_val
    CLOSE(UNIT=scratch)
  END SUBROUTINE print_extfield
END MODULE External
