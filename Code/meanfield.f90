!------------------------------------------------------------------------------
! MODULE: Meanfield
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module calculates all the ingredients needed for the
!!energy functional and for applying the single-particle Hamiltonian to
!!a wave function.
!>
!>@details
!!The work is done by two subroutines: \c skyrme for the calculation
!!of all the fields, which can be scalar or vector and
!!isospin-dependent. \c hpsi then is the routine applying the
!!single-particle Hamiltonian to one single-particle wave function.
!!
!!Note the division of labor between \c skyrme and
!!\c add_density of module \c Densities: everything that
!!constructs fields - densities and current densities - from the
!!single-particle wave functions is done in \c add_density, which is
!!called in a loop over the states, while \c skyrme does the further
!!manipulations to complete the fields entering the single-particle
!!Hamiltonian by combining the densities and their derivatives. It does
!!not need access to the wave functions.
!------------------------------------------------------------------------------
Module Meanfield
  USE Params, ONLY: db,tcoul
  USE Densities
  USE Forces 
  USE Grids, ONLY: nx,ny,nz,der1x,der2x,der1y,der2y,der1z,der2z,dx,dy,dz,x,y,z
  USE Coulomb, ONLY: poisson,wcoul
  IMPLICIT NONE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: upot   !<this is the local part of the mean field 
  !!\f$ U_q \f$. It is a scalar field with isospin index.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: bmass  !<this is the effective mass \f$ B_q \f$.
  !!It is a scalar, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: divaq  !<this is the divergence of \c aq,
  !!i.e., \f$ \nabla\cdot\vec A_q \f$. Its is a scalar, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: aq     !<This is the vector filed \f$ \vec A_q \f$. 
  !!It is a vector, isospin-dependent field.  
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: spot   !<the field \f$ \vec{S}_q \f$.
  !!It is a vector, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: wlspot !<the field \f$ \vec W_q \f$. 
  !!It is a vector, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: dbmass !<contains the gradient of \c bmass.
  !!It is a vector, isospin-dependent field.
  PRIVATE :: divaq,aq,wlspot,dbmass
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_fields
!> @brief
!!This subroutine has the simple task of allocating all the fields that
!!are local to the module \c Meanfield.
!--------------------------------------------------------------------------- 
  SUBROUTINE alloc_fields
    ALLOCATE(upot(nx,ny,nz,2),bmass(nx,ny,nz,2),divaq(nx,ny,nz,2), &
         aq(nx,ny,nz,3,2),spot(nx,ny,nz,3,2),wlspot(nx,ny,nz,3,2), &
         dbmass(nx,ny,nz,3,2))
    upot=0.D0
    bmass=0.D0
    wlspot=0.D0
    aq=0.D0
    divaq=0.D0
    dbmass=0.D0
  END SUBROUTINE alloc_fields
!---------------------------------------------------------------------------  
! DESCRIPTION: skyrme
!> @brief
!!In this subroutine the various fields are calculated from the
!!densities that were previously generated in module \c Densities.
!>
!> @details
!!The expressions divide up the
!!contributions into an isospin-summed part with \f$ b \f$-coefficients
!!followed by the isospin-dependent one with \f$ b' \f$-coefficients. As it
!!would be a waste of space to store the summed densities and currents,
!!the expressions are divided up more conveniently. If we denote the
!!isospin index \f$ q \f$ by \c iq and the index for the opposite isospin
!!\f$ q' \f$ by \c ic (as \c iq can take the values 1 or 2, it can
!!conveniently be calculated as <tt> 3-iq </tt>), this can be written for
!!example as
!!
!!\f$ b_1\rho-b_1'\rho_q\longrightarrow b_1(\rho_q+\rho_{q'})-b_1'\rho_q \f$
!!\f$ \longrightarrow \f$ <tt>(b1-b1p)*rho(:,:,:,iq)+b1*rho(:,:,:,ic)</tt>
!!This decomposition is used in all applicable cases.
!!
!!For intermediate results the fields \c workden (scalar) and \c workvec 
!!(vector) are used.
!!
!!Now the subroutine proceeds in the following steps:
!!  -# all the parts in \c upot involving an \f$ \alpha \f$-dependent power of the
!!     density are collected. Note that in order to avoid having to
!!     calculate different powers, \f$ \rho^\alpha \f$ is factored out. The
!!     division by the total density uses the small number \c epsilon to
!!     avoid division by zero.
!!  -# the divergence of \f$ \vec J \f$ (\c sodens) is
!!     calculated for both isospins in \c workden and the contributions
!!     are added to \c upot.
!!  -# the Coulomb potential is calculated using
!!     subroutine \c poisson  (see module \c Coulomb). It and the
!!     Slater exchange correction (only if the \c ex parameter in the
!!     force is nonzero) are added to \c upot for protons, <tt> iq=2 </tt>.
!!  -# the Laplacian is applied to the densities and
!!     the result stored in \c workden. Then the remaining terms are constructed.  
!!     Note that the  \c iq -loop is combined with the following steps.
!!  -# the effective mass is calculated.
!!  -# the gradient of the density is calculated and the
!!     spin-orbit vector \f$ \vec W_q \f$ is constructed in \c wlspot.
!!  -# the curl of the spin density vector is calculated
!!     and stored in \c workvec.
!!  -# the vector \f$ \vec A_q \f$ is calculated from the current density 
!!     and the curl of the spin density.
!!  -# the curl of the current density is calculated and stored in \c spot.
!!  -# now the two isospin contributions in \c spot
!!     are combined in the proper way.
!!     This way of handling it avoids the introduction of an additional
!!     work vector for \f$ \nabla\times\vec\jmath_q \f$.
!!  -# the divergence of \f$ \vec A_q \f$ is calculated and stored in \c divaq.
!!  -# finally, the gradient of the effective mass term
!!     \f$ B_q \f$ is calculated and stored in the vector variable \c dbmass.
!!  .
!!This concludes the calculation of all scalar and vector fields needed
!!for the application of the Skyrme force. 
!> @param[in] outpot
!> LOGICAL, if true, only outer potential of type \c outertype is calculated.
!> @param[in] outertype
!> CHARACTER, Determines type of outer potential.
!--------------------------------------------------------------------------- 
  SUBROUTINE skyrme(outpot,outertype)
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    LOGICAL,INTENT(IN)        :: outpot
    CHARACTER(1),INTENT(IN)   :: outertype
    REAL(db),PARAMETER        :: epsilon=1.0d-25  
    REAL(db)                  :: rotspp,rotspn
    REAL(db),ALLOCATABLE      :: workden(:,:,:,:),workvec(:,:,:,:,:)
    REAL(db),ALLOCATABLE,SAVE :: random_k(:,:)
    INTEGER                   :: ix,iy,iz,ic,iq,icomp
    ALLOCATE(workden(nx,ny,nz,2),workvec(nx,ny,nz,3,2))
    !  Step 1: 3-body contribution to upot.
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=(rho(:,:,:,1)+rho(:,:,:,2))**f%power * &
            ((b3*(f%power+2.D0)/3.D0-2.D0*b3p/3.D0)*rho(:,:,:,iq) &
            +b3*(f%power+2.D0)/3.D0*rho(:,:,:,ic) &
            -(b3p*f%power/3.D0)*(rho(:,:,:,1)**2+rho(:,:,:,2)**2)/ &
            (rho(:,:,:,1)+rho(:,:,:,2)+epsilon))
    ENDDO
    ! Step 2: add divergence of spin-orbit current to upot
    DO iq=1,2
       CALL rmulx(der1x,sodens(:,:,:,1,iq),workden(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),workden(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),workden(:,:,:,iq),1)
    ENDDO
    DO iq=1,2  
       ic=3-iq 
       upot(:,:,:,iq)=upot(:,:,:,iq) &
            -(b4+b4p)*workden(:,:,:,iq)-b4*workden(:,:,:,ic)
    ENDDO
    ! Step 3: Coulomb potential
    IF(tcoul) THEN
       CALL poisson
       upot(:,:,:,2)=upot(:,:,:,2)+wcoul
       IF(f%ex/=0) &
            upot(:,:,:,2)=upot(:,:,:,2)-slate*rho(:,:,:,2)**(1.0D0/3.0D0)
    ENDIF
    ! Step 4: remaining terms of upot
    DO iq=1,2
       CALL rmulx(der2x,rho(:,:,:,iq),workden(:,:,:,iq),0)  
       CALL rmuly(der2y,rho(:,:,:,iq),workden(:,:,:,iq),1)  
       CALL rmulz(der2z,rho(:,:,:,iq),workden(:,:,:,iq),1)
    ENDDO
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=upot(:,:,:,iq)+(b0-b0p)*rho(:,:,:,iq)+b0*rho(:,:,:,ic) &
                                ! t1,t2, and tau-dependent part      !
            +(b1-b1p)*tau(:,:,:,iq)+b1*tau(:,:,:,ic) &
                                ! two-body laplacian*rho-dependent part
            -(b2-b2p)*workden(:,:,:,iq)-b2*workden(:,:,:,ic)
       ! Step 5: effective mass
       bmass(:,:,:,iq)=f%h2m(iq)+(b1-b1p)*rho(:,:,:,iq)+b1*rho(:,:,:,ic)
       ! Step 6: calculate grad(rho) and wlspot
       CALL rmulx(der1x,rho(:,:,:,iq),workvec(:,:,:,1,iq),0)
       CALL rmuly(der1y,rho(:,:,:,iq),workvec(:,:,:,2,iq),0)
       CALL rmulz(der1z,rho(:,:,:,iq),workvec(:,:,:,3,iq),0)
    ENDDO
    DO iq=1,2
       ic=3-iq
       wlspot(:,:,:,:,iq)= &
            (b4+b4p)*workvec(:,:,:,:,iq)+b4*workvec(:,:,:,:,ic)
    END DO
    ! Step 7: calculate curl of spin density vector, store in workvec
    DO iq=1,2  
       CALL rmuly(der1y,sdens(:,:,:,3,iq),workvec(:,:,:,1,iq),0)
       CALL rmulz(der1z,sdens(:,:,:,2,iq),workvec(:,:,:,1,iq),-1)
       CALL rmulz(der1z,sdens(:,:,:,1,iq),workvec(:,:,:,2,iq),0)
       CALL rmulx(der1x,sdens(:,:,:,3,iq),workvec(:,:,:,2,iq),-1)
       CALL rmulx(der1x,sdens(:,:,:,2,iq),workvec(:,:,:,3,iq),0)
       CALL rmuly(der1y,sdens(:,:,:,1,iq),workvec(:,:,:,3,iq),-1)
    ENDDO
    ! Step 8: calculate A_q vector
    DO iq=1,2
       ic=3-iq
       aq(:,:,:,:,iq)=-2.0D0*(b1-b1p)*current(:,:,:,:,iq) &
            -2.0D0*b1*current(:,:,:,:,ic) &
            -(b4+b4p)*workvec(:,:,:,:,iq)-b4*workvec(:,:,:,:,ic)
    ENDDO
    ! Step 9: calculate the curl of the current density, stopr in spot
    DO iq=1,2  
       CALL rmuly(der1y,current(:,:,:,3,iq),spot(:,:,:,1,iq),0)
       CALL rmulz(der1z,current(:,:,:,2,iq),spot(:,:,:,1,iq),-1)
       CALL rmulz(der1z,current(:,:,:,1,iq),spot(:,:,:,2,iq),0)
       CALL rmulx(der1x,current(:,:,:,3,iq),spot(:,:,:,2,iq),-1)
       CALL rmulx(der1x,current(:,:,:,2,iq),spot(:,:,:,3,iq),0)
       CALL rmuly(der1y,current(:,:,:,1,iq),spot(:,:,:,3,iq),-1)
    ENDDO
    ! Step 10: combine isospin contributions
    DO icomp=1,3  
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx  
                rotspp=spot(ix,iy,iz,icomp,1)  
                rotspn=spot(ix,iy,iz,icomp,2)  
                spot(ix,iy,iz,icomp,1)=-(b4+b4p)*rotspp-b4*rotspn  
                spot(ix,iy,iz,icomp,2)=-(b4+b4p)*rotspn-b4*rotspp  
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! Step 11: calculate divergence of aq in divaq 
    DO iq=1,2  
       CALL rmulx(der1x,aq(:,:,:,1,iq),divaq(:,:,:,iq),0)
       CALL rmuly(der1y,aq(:,:,:,2,iq),divaq(:,:,:,iq),1)
       CALL rmulz(der1z,aq(:,:,:,3,iq),divaq(:,:,:,iq),1)
    ENDDO
    ! Step 12: calculate the gradient of the effective mass in dbmass
    DO iq=1,2  
       CALL rmulx(der1x,bmass(:,:,:,iq),dbmass(:,:,:,1,iq),0)
       CALL rmuly(der1y,bmass(:,:,:,iq),dbmass(:,:,:,2,iq),0)
       CALL rmulz(der1z,bmass(:,:,:,iq),dbmass(:,:,:,3,iq),0)
    ENDDO
    DEALLOCATE(workden,workvec)
!   external guiding potential
    IF (outpot) THEN
      SELECT CASE(outertype)
      CASE('A')
        IF(.NOT.ALLOCATED(random_k)) THEN
          ALLOCATE(random_k(3,2))
          CALL init_random_seed()
          CALL RANDOM_NUMBER(random_k)
          random_k=(random_k-0.5d0)*2
          IF(wflag)WRITE(*,*) 'Random potential:',random_k
        END IF
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
          upot(ix,iy,iz,iq)=random_k(1,1)*cos(REAL(ix)/nx*2*pi)+random_k(1,2)*sin(REAL(ix)/nx*2*pi)&
                           +random_k(2,1)*cos(REAL(iy)/ny*2*pi)+random_k(2,2)*sin(REAL(iy)/ny*2*pi)&
                           +random_k(3,1)*cos(REAL(iz)/nz*2*pi)+random_k(3,2)*sin(REAL(iz)/nz*2*pi)
        END FORALL
        IF(wflag)WRITE(*,*) 'Random guiding potential'
      CASE('P')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
          upot(ix,iy,iz,iq)=10*(cos(REAL(ix)/nx*2*pi)+cos(REAL(iy)/ny*2*pi)+cos(REAL(iz)/nz*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'P-surface guiding potential'
      CASE('E')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
          upot(ix,iy,iz,iq)=10*abs(cos(REAL(ix)/nx*2*pi)+cos(REAL(iy)/ny*2*pi)+cos(REAL(iz)/nz*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Double P-surface guiding potential'
      CASE('F')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
          upot(ix,iy,iz,iq)=-10*abs(cos(REAL(ix)/nx*2*pi)+cos(REAL(iy)/ny*2*pi)+cos(REAL(iz)/nz*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Reverse double P-surface guiding potential'
      CASE('G')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
              upot(ix,iy,iz,iq)=30*(cos(REAL(ix)/nx*2*pi)*sin(REAL(iy)/ny*2*pi)+&
              cos(REAL(iy)/ny*2*pi)*sin(REAL(iz)/nz*2*pi)+cos(REAL(iz)/nz*2*pi)*sin(REAL(ix)/nx*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Gyroid guiding potential'
      CASE('B')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
              upot(ix,iy,iz,iq)=30*abs(cos(REAL(ix)/nx*2*pi)*sin(REAL(iy)/ny*2*pi)+&
              cos(REAL(iy)/ny*2*pi)*sin(REAL(iz)/nz*2*pi)+cos(REAL(iz)/nz*2*pi)*sin(REAL(ix)/nx*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Double Gyroid guiding potential'
      CASE('C')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
              upot(ix,iy,iz,iq)=-30*abs(cos(REAL(ix)/nx*2*pi)*sin(REAL(iy)/ny*2*pi)+&
              cos(REAL(iy)/ny*2*pi)*sin(REAL(iz)/nz*2*pi)+cos(REAL(iz)/nz*2*pi)*sin(REAL(ix)/nx*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Reverse Double Gyroid guiding potential'
      CASE('D')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
              upot(ix,iy,iz,iq)=30*(cos(REAL(ix)/nx*2.0*pi)*cos(REAL(iy)/ny*2.0*pi)*cos(REAL(iz)/nz*2.0*pi)+&
                                    cos(REAL(ix)/nx*2.0*pi)*sin(REAL(iy)/ny*2.0*pi)*sin(REAL(iz)/nz*2.0*pi)+&
                                    sin(REAL(ix)/nx*2.0*pi)*cos(REAL(iy)/ny*2.0*pi)*sin(REAL(iz)/nz*2.0*pi)+&
                                    sin(REAL(ix)/nx*2.0*pi)*sin(REAL(iy)/ny*2.0*pi)*cos(REAL(iz)/nz*2.0*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Diamond guiding potential'
      CASE('S')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
         upot(ix,iy,iz,iq)=50*(cos(REAL(ix)/nx*2*pi))
        END FORALL
        IF(wflag)WRITE(*,*) 'Slab guiding potential'
      CASE('X')
        DO iq=1,2; DO ix=1,nx; DO iy=1,ny; DO iz=1,nz
!          IF(y(iy)>0) upot(ix,iy,iz,iq)=40*(cos(REAL(ix)/nx*2*pi+pi*REAL(iy)/ny))!+&
            !(1.0d0-abs(REAL(iy-ny/2))/REAL(ny)*2)*REAL(iz)/nz*pi*2))
          upot(ix,iy,iz,iq)=40*(cos(REAL(ix)/nx*2*pi-&
            (1.0d0-abs(REAL(iy-ny/2))/REAL(ny)*2)*REAL(iz)/nz*pi*2))
        END DO; END DO; END DO; END DO
        IF(wflag)WRITE(*,*) 'Parking ramp guiding potential'
      CASE('R')
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
         upot(ix,iy,iz,iq)=exp(-((x(ix))**2+(y(iy))**2)/(dx*nx/4)**2)
        END FORALL  
        upot=upot*(-50)/maxval(upot(:,:,:,:))
        IF(wflag)WRITE(*,*) 'Rod guiding potential'
      CASE('H')
        IF(abs((ny*dy)/(nx*dx)-sqrt(3.0))>0.01) STOP 'Ratio L_y/L_x makes no &
        &sense for hexagonal'
        FORALL(iq=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
         upot(ix,iy,iz,iq)=exp(-((x(ix)-dx*nx/4)**2+(y(iy)-dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+dx*nx/4)**2+(y(iy)+dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-5*dx*nx/4)**2+(y(iy)-dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+5*dx*nx/4)**2+(y(iy)+dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-dx*nx/4)**2+(y(iy)-5*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+dx*nx/4)**2+(y(iy)+5*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-5*dx*nx/4)**2+(y(iy)-5*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+5*dx*nx/4)**2+(y(iy)+5*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+3*dx*nx/4)**2+(y(iy)-dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-3*dx*nx/4)**2+(y(iy)+dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-dx*nx/4)**2+(y(iy)+3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+dx*nx/4)**2+(y(iy)-3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+3*dx*nx/4)**2+(y(iy)+3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-3*dx*nx/4)**2+(y(iy)-3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-5*dx*nx/4)**2+(y(iy)+3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+5*dx*nx/4)**2+(y(iy)-3*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)+3*dx*nx/4)**2+(y(iy)-5*dy*ny/4)**2)/(dx*nx/4)**2)&
                          +exp(-((x(ix)-3*dx*nx/4)**2+(y(iy)+5*dy*ny/4)**2)/(dx*nx/4)**2)
        END FORALL  
        upot=upot*(-50)/maxval(upot(:,:,:,:))
        IF(wflag)WRITE(*,*) 'Hexagonal rod guiding potential'
      END SELECT
    END IF
  END SUBROUTINE skyrme
!---------------------------------------------------------------------------  
! DESCRIPTION: hpsi
!> @brief
!!This subroutine applies the single-particle Hamiltonian to a
!!single-particle wave function \c pinn to produce an output wave
!!function \c pout. The argument \c iq indicates the isospin for
!!the wave function and \c eshift is an energy shift which is zero in
!!the dynamic calculation but crucial to the static algorithm (see 
!!\c grstep in module \c Static).
!>
!> @details
!!For an understanding of this module the role of the following local
!!variables is crucial.  They are
!!  - \c is: this is used in the loops over spin to indicate the
!!    spin component: <tt> is=1 </tt> for spin up and <tt> is=2 </tt> for spin down.
!!  - \c ic: denotes the index for the opposite spin; it is
!!    calculated as <tt> ic=3-is </tt>. Note the similar handling of
!!    the two isospin projections using \c iq and \c icomp in
!!    subroutine \c skyrme.
!!  - \c sigis: this variable denotes the sign of the spin
!!    projection. It is calculated as <tt> sigis=3-2*is </tt> and thus is \c +1
!!    for spin up (<tt> is=1 </tt>) and \c - for spin down (<tt> is=2 </tt>).
!!
!!The general structure of the subroutine is as follows: first the part
!!of the Hamiltonian not involving derivatives is applied, followed by the
!!terms involving derivatives in order \f$ x \f$, \f$ y \f$, $z$.
!!Since the structure of the Hamiltonian involves only first or second
!!derivatives in one spatial direction in each term, the derivatives can
!!be calculated for one direction and then the working space can be
!!reused for the next one.
!!
!!The expressions for the different spatial derivatives are quite
!!analogous, so that only the $x$-direction will be discussed at length
!!below.
!!
!!The expressions is repeated here:
!!\f[
!!  \hat h=U_q(\vec r)-\nabla\cdot\left[B_q(\vec r)\nabla\right]
!!  +\I\vec W_q\cdot(\vec\sigma\times\nabla)
!!  +\vec S_q\cdot\vec\sigma
!!  -\frac{\I}{2} \left[(\nabla\cdot\vec A_q)+2\vec A_q\cdot\nabla\right].
!!\f]
!! 
!!  -# the non-derivative parts not involving spin. These
!!     arise from \f$ U_q \f$ and \f$ -\tfrac{\I}{2}\,\nabla\cdot\vec A_q \f$, which
!!     are combined into a complex expression. The energy shift \c eshift is also included.
!!  -# the spin current coupling is constructed by simply
!!     using the explicit definition of the Pauli matrices and multiplying
!!     the resulting matrix onto the spinor wave function.
!!  -# the first and second derivative in the \f$ x \f$-direction are evaluated 
!!     and stored in the arrays \c pswk and \c pswk2. The last term in the Hamiltonian 
!!     gives rise to the two contributions
!!     \f[ -\frac{\partial B_q}{\partial x}\frac{\partial}{\partial x}-B_q 
!!     \frac{\partial^2}{\partial x^2}, \f] 
!!     of which the second is evaluated
!!     straightforwardly, while the first one is combined with the spin-orbit
!!     contribution. The part of \f$ \I\vec W_q\cdot(\vec\sigma\times\nabla) \f$
!!     that contains an \f$ x \f$-derivative is
!!     \f[ (\I W_y\sigma_z-\I W_z\sigma_y)\frac{\partial}{\partial x}=
!!     \begin{pmatrix} \I W_y&-W_z\\ W_z&-\I W_y
!!     \end{pmatrix}\frac{\partial}{\partial x} \f]
!!     This is programmed employing the variable \c sigis to account for
!!     the different signs in the rows of the matrix.
!!  -# for the derivatives in the $y$-direction the
!!     procedure is similar; the spin-orbit part is now
!!     \f[ (\I W_z\sigma_x-\I W_x\sigma_z)\frac{\partial}{\partial y}=
!!     \begin{pmatrix} -\I W_x&\I W_z\\ \I W_z&\I W_x
!!     \end{pmatrix}\frac{\partial}{\partial y} \f]
!!  -# for the derivatives in the \f$ z \f$-direction the
!!     procedure is again similar; the spin-orbit part is now
!!     \f[ (\I W_x\sigma_y-\I W_y\sigma_x)\frac{\partial}{\partial z}=
!!     \begin{pmatrix} 0&W_x-\I W_y\\ -W_x-\I W_y & 0
!!     \end{pmatrix}\frac{\partial}{\partial z} \f]
!>
!> @param[in] iq
!> INTEGER, takes the isospin.
!> @param[in] eshift
!> REAL(db), takes the energy shift.
!> @param[in,out] pinn
!> COMPLEX(db), array, takes the wave function the Hamiltonian is supposed to be applied to.
!> @param[out] pout
!> COMPLEX(db), array, returns the Hamiltonian applied to the wave function.
!--------------------------------------------------------------------------- 
  SUBROUTINE hpsi(iq,eshift,pinn,pout)
    USE Trivial, ONLY: cmulx, cmuly, cmulz
    USE Levels, ONLY: cdervx,cdervy,cdervz
    INTEGER :: iq
    REAL(db) :: eshift  
    COMPLEX(db),DIMENSION(:,:,:,:) :: pinn,pout
    INTENT(IN) :: iq,eshift
    INTENT(INOUT) :: pinn
    INTENT(OUT) :: pout
    INTEGER :: is,ic
    REAL(db) :: sigis
    COMPLEX(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: pswk,pswk2
    ALLOCATE(pswk(nx,ny,nz,2),pswk2(nx,ny,nz,2))
    ! Step 1: non-derivative parts not involving spin
    DO is=1,2  
       pout(:,:,:,is)=CMPLX(upot(:,:,:,iq)-eshift, &
            -0D0,db)*pinn(:,:,:,is)
    ENDDO
    ! Step 2: the spin-current coupling
    pout(:,:,:,1)=pout(:,:,:,1)  &
         + CMPLX(spot(:,:,:,1,iq),-spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,2)  + spot(:,:,:,3,iq)*pinn(:,:,:,1)
    pout(:,:,:,2)=pout(:,:,:,2) &
         + CMPLX(spot(:,:,:,1,iq),spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,1) - spot(:,:,:,3,iq)*pinn(:,:,:,2)
    ! Step 3: derivative terms in x
    IF(TFFT) THEN
       CALL cdervx(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmulx(der1x,pinn,pswk,0)  
       CALL cmulx(der2x,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,1,iq),0.5D0*aq(:,:,:,1,iq) &
            -sigis*wlspot(:,:,:,2,iq),db)*pswk(:,:,:,is)  &
            -sigis*wlspot(:,:,:,3,iq)*pswk(:,:,:,ic) &
            -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)-wlspot(:,:,:,2,iq))*pinn(:,:,:,1)&
         -0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)+wlspot(:,:,:,2,iq))*pinn(:,:,:,2)&
         +0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervx(pswk2,pswk)  
    ELSE
       CALL cmulx(der1x,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    ! Step 4: derivative terms in y
    IF(TFFT) THEN
       CALL cdervy(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmuly(der1y,pinn,pswk,0)  
       CALL cmuly(der2y,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,2,iq),0.5D0*aq(:,:,:,2,iq) &
            +sigis*wlspot(:,:,:,1,iq),db)*pswk(:,:,:,is) &
            +CMPLX(0.D0,0.5D0*wlspot(:,:,:,3,iq),db)*pswk(:,:,:,ic) &
            -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)+wlspot(:,:,:,1,iq))*pinn(:,:,:,1)&
         +CMPLX(0D0,0.5D0,db)*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)-wlspot(:,:,:,1,iq))*pinn(:,:,:,2)&
         +CMPLX(0D0,0.5D0*wlspot(:,:,:,3,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervy(pswk2,pswk)  
    ELSE
       CALL cmuly(der1y,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    ! Step 5: derivative terms in z
    IF(TFFT) THEN
       CALL cdervz(pinn,pswk,d2psout=pswk2)  
    ELSE
       CALL cmulz(der1z,pinn,pswk,0)  
       CALL cmulz(der2z,pinn,pswk2,0)  
    ENDIF
    DO is=1,2  
       ic=3-is  
       sigis=(3-2*is)*0.5D0  
       pout(:,:,:,is)=pout(:,:,:,is) &
            -CMPLX(dbmass(:,:,:,3,iq),0.5D0*aq(:,:,:,3,iq),db)*pswk(:,:,:,is) &
            +CMPLX(sigis*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)* &
            pswk(:,:,:,ic)-bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,1)&
         +CMPLX(0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*&
            pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,2)&
         +CMPLX(-0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*&
            pinn(:,:,:,1)
    IF(TFFT) THEN
       CALL cdervz(pswk2,pswk)  
    ELSE
       CALL cmulz(der1z,pswk2,pswk,0)  
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    !
    DEALLOCATE(pswk,pswk2)
  END SUBROUTINE hpsi
!---------------------------------------------------------------------------  
! DESCRIPTION: init_random_seed
!> @brief
!!Initializes random number generator for random potential.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
  END SUBROUTINE
END Module Meanfield
