!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!In this module various moments of the density distribution are
!!calculated. Most of them come in two versions: an isospin-dependent
!!one characterized by a final isospin index of dimension~2, and a
!!summed one distinguished by the ending "\c tot" in its name.
!!Thus, e.g., three variants of center of mass are
!!stored: that of the neutrons <tt> cm(1:3,1) </tt>, of the protons <tt>
!!cm(1:3,2)</tt>, and of the total mass distribution <tt> cmtot(1:3) </tt>.
!>
!>@details 
!!Since the geometrical arrangement in space can be arbitrary with
!!respect to the Cartesian coordinate system - this is certainly true
!!for non-central collision situations - some quantities associated
!!with the axes become meaningless in the general situation. For this
!!reason, the code also calculates the complete Cartesian quadrupole
!!tensor and diagonalizes it to obtain the principal
!!axes of the nucleus. In this frame then we compute the spherical
!!quadrupole moments \f$ Q_{2m} \f$ with their
!!dimensionless counterparts \f$ a_o \f$, \f$ a_2 \f$ and the deformation 
!!parameters \f$ \beta \f$ and \f$ \gamma \f$.

!------------------------------------------------------------------------------
MODULE Moment
  USE Params
  USE Grids, ONLY: nx,ny,nz,x,y,z,wxyz
  USE Spherical_Harmonics

  IMPLICIT NONE
  PRIVATE
  REAL(db) :: pnr(2)    !<the numbers of neutrons <tt>pnr(1)</tt>=\f$ N \f$ and 
  !!protons <tt>pnr(2)</tt>= \f$ Z \f$. These are obtained by a simple integration
  !!of the densities \c rho. 
  REAL(db) :: pnrtot    !<the total particle number <tt> pnrtot</tt>=\f$ A \f$. 
  !!These are obtained by a simple integration of the densities \c rho. 
  REAL(db) :: cm(3,2)   !<the center of mass vectors of
  !!the neutron and proton \f$ \vec R_n \f$, \f$ \vec R_p \f$. Dimension: fm.
  REAL(db) :: cmtot(3)  !<the center of mass vectors of
  !!the total mass distribution \f$ \vec R \f$. Dimension: fm.
  REAL(db) :: pcm(3,2)  !<the integrated momentum vectors, not containing
  !!the nucleon mass. They thus correspond to an integral over the
  !!current density only and have a dimension of velocity: \f$ c \f$.
  REAL(db) :: rms(2)    !< these are the root
  !!mean-square radii of neutron and proton mass distribution.  Dimension: fm.
  REAL(db) :: rmstot    !< these are the root mean-square radii of the total
  !!mass distribution.  Dimension: fm.
  REAL(db) :: q20(2)    !< the \f$ m=0 \f$ components of the spherical quadrupole tensor 
  !!\f$ Q_{20} \f$ in the principal-axes frame for neutrons and protons.
  REAL(db) :: q20tot    !< the \f$ m=0 \f$ components of the spherical quadrupole tensor 
  !!\f$ Q_{20} \f$ in the principal-axes frame for the total mass distribution.
  REAL(db) :: q22(2)    !< the \f$ m=2 \f$ components of the spherical quadrupole tensor 
  !!\f$ Q_{22} \f$ in the principal-axes frame for neutrons and protons.
  REAL(db) :: q22tot    !< the \f$ m=2 \f$ components of the spherical quadrupole tensor 
  !!\f$ Q_{20} \f$ in the principal-axes frame for the total mass distribution.
  REAL(db) :: x2m(3,2)  !<the vectors of second moments of the radii in the three coordinate
  !!directions, i.e.,\f[ \langle r_i^2\rangle=\frac1{A} \int \D^3r\, r_i^2\rho(\vec r). \f]
  !!for neutrons and protons
  !!They are useful to get an idea of the shape of the nucleus in static calculations. 
  !!Dimension: \f${\rm fm}^2 \f$
  REAL(db) :: x2mtot(3) !<the vectors of second moments of the radii in the three coordinate
  !!directions, i.e.,\f[ \langle r_i^2\rangle=\frac1{A} \int \D^3r\, r_i^2\rho(\vec r). \f]
  !!for the total mass distribution
  !!They are useful to get an idea of the shape of the nucleus in static calculations. 
  !!Dimension: \f${\rm fm}^2 \f$
  REAL(db) :: beta20tot !<the quadrupole deformation parameters \f$ a_0 \f$.
  REAL(db) :: beta22tot !<the quadrupole deformation parameters \f$ a_2 \f$.
  REAL(db) :: beta      !<the Bohr-Mottelson deformation parameter \f$ \beta \f$. Calculated
  !!only for the total mass distribution, dimensionless
  REAL(db) :: gamma     !<the Bohr-Mottelson deformation parameter \f$ \gamma \f$. Calculated
  !!only for the total mass distribution, expressed in degrees.
  REAL(db) :: r3(2) !< these are the avg
  !!value of \f$r^3\f$ of neutron and proton mass distribution.  Dimension: fm\f$^3\f$.
  REAL(db) :: r4(2) !< these are the avg
  !!value of \f$r^4\f$ of neutron and proton mass distribution.  Dimension: fm\f$^4\f$.
  REAL(db) :: r2tot !< these are the avg
  !!value of \f$r^2\f$ of total mass distribution.  Dimension: fm\f$^2\f$.
  REAL(db) :: r3tot !< these are the avg
  !!value of \f$r^3\f$ of total mass distribution.  Dimension: fm\f$^3\f$.
  REAL(db) :: r4tot !< these are the avg
  !!value of \f$r^4\f$ of total mass distribution.  Dimension: fm\f$^4\f$.

  REAL(db) :: Mono_tot !< The monopole moment 
  REAL(db) :: Di_tot !< The Dipole moment 
  REAL(db) :: Quad_tot !< The quadrupole moment 
  REAL(db) :: Oct_tot !< The Octupole moment 
  REAL(db) :: HexDec_tot !< The Hexadecapole moment 
  REAL(db) :: DiaTriaConta_tot !< The DiaTriaContapole moment 
  
  

  PUBLIC :: pnr,pnrtot,cm,cmtot,pcm,rmstot,beta,gamma, &
       moments,moment_print,moment_shortprint,Oct_tot, &
       HexDec_tot,DiaTriaConta_tot,Mono_tot,Di_tot,Quad_tot,r3tot,r4tot,r2tot
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: Routinename
!> @brief
!!This is the principal subroutine for calculating the geometric
!!quantities. It consists of two loops, both over isospin and space, and
!!a final analysis Section. Along with that, it also calculates the multipole moments as a function of time.
!>
!> @details
!!In detail:
!!  - The first loop calculates the particle numbers, centers of mass,
!!    and, for the dynamic case only, the momenta divided by nucleon mass.
!!  - The second loop is separate because the center-of mass must be
!!    known to use it as the origin for the vectors. The r.m.s. radii
!!    \c rms and the quantities \c x2m are calculated as well as the
!!    components \c qmat of the Cartesian quadrupole tensor
!!    \f$ Q_{kl}=\int \D^3r \left(3x_kx_l-r^2\delta_{kl}\right) \rho(\vec r) \f$.
!!  - After this, the subroutine \c q2diag is used to determine the
!!    spherical components \f$ Q_{20} \f$ and \f$ Q_{22} \f$ in the 
!!    principal-axes frame also generating some printout in the process.
!!    These are then multiplied with a scale factor to yield the
!!    dimensionless deformation parameters \f$ a_0 \f$ and \f$ a_2 \f$ and 
!!    finally by conversion to polar coordinates the Bohr-Mottelson parameters 
!!    \f$ \beta \f$ and \f$ \gamma \f$.
!!  - The Cartesian and polar deformation parameters are then printed.

!--------------------------------------------------------------------------- 
  SUBROUTINE moments(L_val,M_val,isoext)
    USE Densities, ONLY: rho,current
    INTEGER :: ix,iy,iz,iq,L_val,M_val,isoext
    REAL(db) :: xx(3),x2(3),vol,radius,eta
    REAL(db) :: qmat(3,3,2),qmtot(3,3)
    REAL(db) :: Mono(2),Di(2),Quad(2),Oct(2),HexaDeca(2),DiaTriaConta(2),tmp,facn,facp
    pnr=0.D0
    cm=0.D0
    pcm=0.D0
    DO iq=1,2  
       DO iz=1,nz  
          xx(3)=z(iz)  
          DO iy=1,ny  
             xx(2)=y(iy)  
             DO ix=1,nx  
                xx(1)=x(ix)
                pnr(iq)=pnr(iq)+wxyz*rho(ix,iy,iz,iq)
                cm(:,iq)=cm(:,iq)+wxyz*xx*rho(ix,iy,iz,iq)
                IF(tdynamic) pcm(:,iq)=pcm(:,iq)+wxyz*current(ix,iy,iz,:,iq)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    pnrtot=pnr(1)+pnr(2)  
    cmtot=(cm(:,1)+cm(:,2))/pnrtot  
    DO iq=1,2
       cm(:,iq)=cm(:,iq)/pnr(iq)
    ENDDO
    !***********************************
    rms=0.D0
    qmat=0.D0
    x2m=0.D0
    r3=0.D0
    r4=0.D0

    Mono=0.0D0
    Di = 0.0D0
    Quad = 0.0d0
    Oct = 0.0D0
    HexaDeca = 0.0D0
    DiaTriaConta = 0.0D0

    DO iq=1,2  
       DO iz=1,nz  
          xx(3)=z(iz)-cm(3,iq)  
          x2(3)=xx(3)**2  
          DO iy=1,ny  
             xx(2)=y(iy)-cm(2,iq)  
             x2(2)=xx(2)**2
             DO ix=1,nx  
                xx(1)=x(ix)-cm(1,iq)  
                x2(1)=xx(1)**2  
                vol=wxyz*rho(ix,iy,iz,iq)  
                rms(iq)=vol*SUM(x2)+rms(iq)
                r3(iq)=vol*(SUM(x2)**(1.5d0))+r3(iq)
                r4(iq)=vol*(SUM(x2)**(2))+r4(iq)

                qmat(1,1,iq)=qmat(1,1,iq)+vol*(x2(1)+x2(1)-x2(2)-x2(3))
                qmat(1,2,iq)=qmat(1,2,iq)+3.D0*vol*xx(1)*xx(2)
                qmat(1,3,iq)=qmat(1,3,iq)+3.D0*vol*xx(1)*xx(3)
                qmat(2,2,iq)=qmat(2,2,iq)+vol*(x2(2)+x2(2)-x2(1)-x2(3))
                qmat(2,3,iq)=qmat(2,3,iq)+3.D0*vol*xx(2)*xx(3)
                qmat(3,3,iq)=qmat(3,3,iq)+vol*(x2(3)+x2(3)-x2(1)-x2(2))
                x2m(:,iq)=vol*x2(:)+x2m(:,iq)  

                Mono(iq)=vol*((0.5d0*SQRT(1.0d0/PI))*(SUM(x2)))+Mono(iq)


                IF(isoext==0) THEN  
                eta = vol*(SUM(x2))*5.0d0/3.0d0
                Di(iq)=vol*(Y_lm(1,M_val,x(ix),y(iy),z(iz))*(SQRT(x(ix)**2+y(iy)**2+z(iz)**2)**3- &
                           SQRT(x(ix)**2+y(iy)**2+z(iz)**2)*eta)*SQRT(2*1+1.0d0))+Di(iq)
               !  print*,'inside loop 1'
                else if (isoext==1) THEN
                tmp=(Y_lm(1,M_val,x(ix),y(iy),z(iz))*SQRT(3.0d0))
                tmp = tmp*SQRT(x(ix)**2+y(iy)**2+z(iz)**2)*vol
                Di(iq)=tmp+Di(iq)
               !  print*,Di(iq),tmp,iq,'inside looop 2'
                end if

                Quad(iq)=vol*(Y_lm(2,M_val,xx(1),xx(2),xx(3))*SUM(x2)*SQRT(2*2+1.0d0))+Quad(iq)
                Oct(iq)=vol*(Y_lm(3,M_val,xx(1),xx(2),xx(3))*(SQRT(SUM(x2))**3) &
                           *SQRT(2*3+1.0d0))+Oct(iq)
                HexaDeca(iq)= vol*(Y_lm(4,M_val,xx(1),xx(2),xx(3))*(SQRT(SUM(x2))**4) &
                           *SQRT(2*4+1.0d0))+HexaDeca(iq)
                DiaTriaConta(iq)= vol*(Y_lm(5,M_val,xx(1),xx(2),xx(3))*(SQRT(SUM(x2))**5) &
                           *SQRT(2*5+1.0d0))+DiaTriaConta(iq)
             ENDDO
          ENDDO
       ENDDO
       qmat(2,1,iq)=qmat(1,2,iq)
       qmat(3,1,iq)=qmat(1,3,iq)
       qmat(3,2,iq)=qmat(2,3,iq)
      !  print*,Di(iq),iq
    ENDDO
    

    r2tot=(rms(1)+rms(2))/pnrtot
    rmstot=SQRT((rms(1)+rms(2))/pnrtot)
    rms=SQRT(rms/pnr)  
    x2mtot=(x2m(:,1)+x2m(:,2))/pnrtot  
    r3tot=(r3(1)+r3(2))/pnrtot
    r4tot=(r4(1)+r4(2))/pnrtot

    IF(isoext==0) THEN  
       facn=1.0D0  
       facp=1.0D0  
    ELSE  
     if (L_val .eq. 1)then
       facn=-(pnr(2)/pnrtot)
       facp=(pnr(1))/pnrtot
      else
       facn=-1
       facp=1
       end if

    end if

    Mono_tot = facn*Mono(1)+facp*Mono(2)
    Di_tot = facn*Di(1)+facp*Di(2)

   !  print*,'Dipole',Di_tot,Di(1),Di(2)


    Quad_tot = facn*Quad(1)+facp*Quad(2)
    Oct_tot = facn*Oct(1)+facp*Oct(2)
    HexDec_tot = facn*HexaDeca(1)+facp*HexaDeca(2)
    DiaTriaConta_tot = facn*DiaTriaConta(1)+facp*DiaTriaConta(2)

    DO iq=1,2
       x2m(:,iq)=x2m(:,iq)/pnr(iq)
    ENDDO
    qmtot=qmat(:,:,1)+qmat(:,:,2) 
    IF(printnow.AND.wflag) WRITE(*,'(/A)') 'Cartesian quadrupole tensor,&
         &  principal values, and axes:'
    CALL q2diag(qmat(:,:,1),q20(1),q22(1),'Neutrons ')
    CALL q2diag(qmat(:,:,2),q20(2),q22(2),'Protons  ')
    CALL q2diag(qmtot,q20tot,q22tot,'Total    ')
    radius=r0*pnrtot**(1.D0/3.D0)
    beta20tot=q20tot*(4.0D0*PI/(5.0D0*radius**2*pnrtot))
    beta22tot=q22tot*(4.0D0*PI/(5.0D0*radius**2*pnrtot))
    beta=SQRT(beta20tot**2+2.0*beta22tot**2)
    gamma=ABS(ATAN2(SQRT(2.0)*beta22tot,beta20tot)*180.0D0/PI)
    IF(gamma>120.D0) THEN
       gamma=gamma-120.D0
    ELSEIF(gamma>60.D0) THEN
       gamma=120.D0-gamma
    ENDIF
    IF(printnow.AND.wflag) WRITE(*,'(4(A,F8.4)/)') &
         ' Beta20: ',beta20tot,' Beta22: ',beta22tot,' Beta: ',beta, &
         ' Gamma: ',gamma
  END SUBROUTINE moments
!---------------------------------------------------------------------------  
! DESCRIPTION: moment_shortprint
!> @brief
!!This subroutine simply prints some information into the specialized
!!output files. \c monopolesfile receives the monopole moment, while \c quadrupolesfile
!!receives the quadrupole moment and so on for other multipole files. The physical time
!!starts each line to enable easy time-curve plotting.
!--------------------------------------------------------------------------- 
  SUBROUTINE moment_shortprint()

    OPEN(unit=scratch,file=monopolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,Mono_tot
    CLOSE(unit=scratch)
    
    OPEN(unit=scratch,file=dipolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,Di_tot
    CLOSE(unit=scratch)

    OPEN(unit=scratch,file=quadrupolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,Quad_tot
    CLOSE(unit=scratch)
    
    OPEN(unit=scratch,file=octupolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,Oct_tot
    CLOSE(unit=scratch)

    OPEN(unit=scratch,file=hexadecapolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,HexDec_tot
    CLOSE(unit=scratch)

    OPEN(unit=scratch,file=diatriacontapolesfile,POSITION='APPEND')  
    WRITE(scratch,'(F25.5,2x,F25.15,2x,F15.10)') time,DiaTriaConta_tot
    CLOSE(unit=scratch)
  END SUBROUTINE moment_shortprint
!---------------------------------------------------------------------------  
! DESCRIPTION: moment_print
!> @brief
!!This subroutine prints a somewhat more detailed information. Particle
!!number , r.m.s. radius, \f$ Q_{20} \f$, \c x2m, and center-of-mass are
!!printed for the total distribution and also separately for neutrons
!!and protons. This output goes to the regular output unit.
!--------------------------------------------------------------------------- 
  SUBROUTINE moment_print
    INTEGER :: iq
    CHARACTER(11),PARAMETER :: Name(2)=(/ '  Neutron: ','   Proton: '/)
    Write(*,'(A)') '              Part.Num.   rms-radius   q20         &
         &<x**2>      <y**2>      <z**2>        <x>            &
         &<y>            <z>'    
    WRITE(*,'(a,2f12.4,1p,4e12.4,3e15.7)') '    Total: ',pnrtot,rmstot, &
         q20tot,x2mtot,cmtot
    DO iq=1,2
       WRITE(*,'(a,2f12.4,1p,4e12.4,3e15.7)') name(iq),pnr(iq),rms(iq),q20(iq), &
            x2m(:,iq),cm(:,iq)
    ENDDO
  END SUBROUTINE moment_print
!---------------------------------------------------------------------------  
! DESCRIPTION: q2diag
!> @brief
!!his subroutine diagonalizes the Cartesian quadrupole tensor. To this
!!purpose it calls the \c LAPACK routine \c DSYEV.
!>
!> @details
!!The eigenvalues are obtained as \c q_eig(i) in ascending order of
!!magnitude, and the corresponding eigenvectors as \c q_vec(:,i).
!!Both are printed. Then the corresponding spherical moments are
!!calculated assuming the \f$ z \f$-axis is selected as that of largest
!!quadrupole moment:
!!\f[ Q_{20}=\sqrt{\frac{5}{16\pi}}\,Q_{zz},\qquad
!!Q_{22}=\sqrt{\frac{5}{96\pi}}\,\left(Q_{yy}-Q_{xx}\right). \f] They are
!!returned in \c q20x and \c q22x. 
!>
!> @param[in,out] q_mat
!> REAL(db), array, takes the quadrupole matrix.
!> @param[out] q20x
!> REAL(db), returns the value \f$ Q_{20} \f$.
!> @param[out] q22x
!> REAL(db), returns the value \f$ Q_{22} \f$. 
!> @param[in] title
!> CHARACTER, array, takes a title for printout. 
!--------------------------------------------------------------------------- 
  SUBROUTINE q2diag(q_mat,q20x,q22x,title)
    REAL(db),INTENT(INOUT) :: q_mat(3,3)
    REAL(db),INTENT(OUT) :: q20x,q22x
    CHARACTER(LEN=*),INTENT(IN) :: title
    REAL(db) :: q_eig(3),q_vec(3,3),fv1(20)
    INTEGER :: info,i, j,k

    if(printnow.AND.wflag) write(*,'(3(f12.5,1x))') ((q_mat(j,k),k=1,3),j=1,3)    
    CALL DSYEV('V','U',3,q_mat,3,q_eig,fv1,20,info)
    q_vec=q_mat
    IF(info/=0) STOP 'Quadrupole diagonalization failed'
    IF(printnow.AND.wflag) THEN
       WRITE(*,'(1X,A,3(F10.2,''('',3F8.4,'')''))') &
            title,(REAL(q_eig(i)),REAL(q_vec(:,i)),i=3,1,-1)
    ENDIF
    q20x=SQRT(5.D0/(16.D0*pi))*REAL(q_eig(3))
    q22x=SQRT(5.D0/(96.D0*pi))*(REAL(q_eig(2))-REAL(q_eig(1)))
  END SUBROUTINE q2diag
END MODULE Moment
