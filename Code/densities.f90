!------------------------------------------------------------------------------
! MODULE: Densities
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module has two purposes: it defines and allocates the densities
!!and currents making up the mean field,
!!and also contains the subroutine \c add_density which accumulates
!!the basic densities over the single-particle wave functions.
!!Subroutine \c skyrme in module \c Meanfield then uses these
!!densities to build up the components of the single-particle Hamiltonian.
!------------------------------------------------------------------------------
MODULE Densities
  USE Params, ONLY: db,tfft
  USE Grids, ONLY: nx,ny,nz,der1x,der1y,der1z
  USE Levels, ONLY: cdervx,cdervy,cdervz
  USE Trivial, ONLY: cmulx,cmuly,cmulz
  IMPLICIT NONE
  SAVE
  !>@name Scalar densities: 
  !>These are dimensioned <tt>(nx,ny,nz,2)</tt>,
  !>where the last index is 1 for neutrons and 2 for protons.
  !>@{
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: rho
  !< the density, separately for each isospin (in
  !!\f${\rm fm}^{-3} \f$).  The definition is:
  !!\f[ \rho_q(\vec r)=\sum_{k\in q}w_k^2\sum_s|\phi_k(\vec r,s)|^2,\qquad q=n,p \f]
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: tau
  !<the kinetic energy density, separately
  !!for each isospin. It is defined as the sum of the spin
  !!contributions and all particles of the given isospin
  !!\f[ \tau_q(\vec r)=\sum_{k\in q}w_k^2\sum_s|\nabla\phi_k(\vec r,s)|^2,
  !!\qquad q=n,p\f] Note that it does not include the factor
  !!\f$ \hbar^2/2m \f$. Units: \f$ {\rm fm}^{-5}$.
  !>@}
  !
  !>@name Vector densities: 
  !>These are dimensioned <tt>(nx,ny,nz,3,2)</tt>, where the last index is 1 for neutrons and 2 for
  !>protons, and the next-to-last stands for the Cartesian direction.
  !>@{
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: current
  !<this is the total probability current
  !!density, defined in the familiar way as
  !!\f[ \vec\jmath_q(\vec r)
  !!= \frac{1}{2\I}\sum_{\alpha\in q}w_\alpha^2
  !!\sum_s\left(\psi_\alpha^*(\vec r,s)\nabla\psi_\alpha(\vec
  !!r,s)-\psi_\alpha(\vec r,s) \nabla\psi_\alpha^*(\vec r,s)\right). \f]
  !!Note that the factor \f$ \frac{\hbar}{m} \f$ is not included. Its units
  !!are therefore \f$ {\rm fm}^{-4} \f$.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: sdens
  !<the spin density. It is defined as
  !!\f[ \vec\sigma_q(\vec r)=\sum_{\alpha\in q} w_\alpha^2 \sum_{ss'}\psi_\alpha^*(\vec
  !!r,s)\,\sigma_{ss'} \,\psi_\alpha(\vec r,s').\f] 
  !!Note that it does not include the factor \f$ \hbar/2 \f$. Units: \f$ {\rm fm}^{-3} \f$.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: sodens
  !<the spin-orbit density, defined as
  !!\f[ \vec J_q(\vec r)=\frac{1}\I\sum_{\alpha\in q}w_\alpha^2
  !!\sum_{ss'}\left(\psi_\alpha^*(\vec r,s)\nabla\times\sigma_{ss'}
  !!\psi_\alpha(\vec r,s')\right).\f]
  !!Its units are also \f$ {\rm fm}^{-4} \f$.
  !>@}
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_densities
!> @brief
!!This is simply a short routine to allocate all the arrays defined in this module.
!--------------------------------------------------------------------------- 
  SUBROUTINE alloc_densities
    ALLOCATE(rho(nx,ny,nz,2),tau(nx,ny,nz,2),current(nx,ny,nz,3,2), &
         sdens(nx,ny,nz,3,2),sodens(nx,ny,nz,3,2))
  END SUBROUTINE alloc_densities
!---------------------------------------------------------------------------  
! DESCRIPTION: add_density
!> @brief
!!This subroutine is given a single-particle wave function \c psin
!!with its isospin index \c iq and occupation \f$ w_\alpha^2= \f$ \c weight and adds
!!its contribution to the density arrays.
!>
!> @details
!!The reason for not including the loop over states in the subroutine is
!!that in the dynamic code, the contribution of a new single-particle
!!wave function (calculated by tstep) to the densities is added
!!without saving that wave function, eliminating the requirement for a
!!second huge wave-function array.
!!
!!It may seem strange that \c add_density has the densities
!!themselves as parameters, which are readily available in the module.
!!The reason for this is \c OPENMP parallelization. The loop over
!!wave functions is done in parallel under \c OPENMP. Since any of
!!the parallel tasks must add up the contributions of its assigned wave
!!functions, each task must have a copy of the densities to work on;
!!otherwise they would try to update the same density at the same time.
!!The separate copies are then combined using the \c OPENMP <tt>REDUCE(+)</tt> directive.
!!
!!The local copies of the densities passed as arrays are denoted with
!!the prefixed letter "l" for \a local; they are \c lrho, \c ltau, 
!!\c lcurrent, \c lsdens, and \c lsodens.
!!
!!If the weight is zero, there is nothing to do and the subroutine
!!returns immediately. Otherwise, the contributions not involving
!!derivatives are first computed and added to the affected densities,
!!i.~e., number and spin density.
!!
!!After this the derivative terms are evaluated by computing each
!!Cartesian direction separately. In all three cases the derivative is
!!evaluated first and put into \c ps1, after which the contributions
!!are added straightforwardly. They involve the wave function itself,
!!the derivative, and for the spin-orbit density also a Pauli matrix, so
!!that different spin projections have to be combined properly.
!!
!!The complex products always in the end evaluate to something real and
!!the expressions are simplified to take this into account. For example,
!!the following transformation is done: 
!!\f{eqnarray*}{
!!\frac{1}{2\I}(\psi^*\nabla\psi-\psi\nabla\psi^*)&=&
!!\frac{1}{2\I}\left(\psi^*\nabla\psi-(\psi^*\nabla\psi)^*\right)\\
!!&=&\frac{1}{2\I}\left(2\I\Im(\psi^*\nabla\psi)\right)\\
!!&\rightarrow&{\tt AIMAG(CONJG(psin)*psi1)}
!!\f}
!!and similarly for the other expressions.
!!
!!The efficiency of this relies on the FORTRAN compiler recognizing that
!!only the imaginary part of the complex product is needed and not
!!computing the real part at all. This seems to be the case with all
!!present compilers.
!>
!> @param[in] iq
!> INTEGER, takes the isospin of the wave function.
!> @param[in] weight
!> REAL(db), takes the BCS weight of the wave function.
!> @param[in, out] psin
!> CMPLEX(db), array, takes wave function .
!> @param[in, out] lrho
!> REAL(db), array, takes and adds the density.
!> @param[in, out] ltau
!> REAL(db), array, takes and adds the kinetic density.
!> @param[in, out] lcurrent
!> REAL(db), array, takes and adds the current density.
!> @param[in, out] lsdens
!> REAL(db), array, takes and adds the spin density.
!> @param[in, out] lsodens
!> REAL(db), array, takes and adds the spin-orbit density.
!--------------------------------------------------------------------------- 
  SUBROUTINE add_density(iq,weight,psin,lrho,ltau,lcurrent,lsdens,lsodens)  
    COMPLEX(db),INTENT(INOUT) :: psin(nx,ny,nz,2)
    REAL(db),DIMENSION(:,:,:,:),INTENT(INOUT) :: lrho,ltau
    REAL(db),DIMENSION(:,:,:,:,:),INTENT(INOUT) :: lcurrent,lsdens,lsodens
    INTEGER,INTENT(IN) :: iq
    REAL(db),INTENT(IN) :: weight
    COMPLEX(db),ALLOCATABLE :: ps1(:,:,:,:)  
    INTEGER :: ix,iy,iz
    ALLOCATE(ps1(nx,ny,nz,2))
    IF(weight<=0.D0) RETURN
    !***********************************************************************
    ! non-derivative terms
    !***********************************************************************
    lrho(:,:,:,iq)=lrho(:,:,:,iq)+weight* &
         (psin(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         psin(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(iz=1:nz,iy=1:ny,ix=1:nx)
       lsdens(ix,iy,iz,1,iq)=lsdens(ix,iy,iz,1,iq)+2.D0*weight &
            *REAL(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,2))
       lsdens(ix,iy,iz,2,iq)=lsdens(ix,iy,iz,2,iq)+2.D0*weight &
            *AIMAG(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,2))
       lsdens(ix,iy,iz,3,iq)=lsdens(ix,iy,iz,3,iq)+weight &
            *(REAL(CONJG(psin(ix,iy,iz,1))*psin(ix,iy,iz,1)) &
            -REAL(CONJG(psin(ix,iy,iz,2))*psin(ix,iy,iz,2)))
    END FORALL
    !***********************************************************************
    ! x-derivatives
    !***********************************************************************
    IF(TFFT) THEN
       CALL cdervx(psin,ps1)  
    ELSE
       CALL cmulx(der1x,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,1,iq)=lcurrent(:,:,:,1,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
       lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)-weight &
            *(  AIMAG(ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1))) &
            - AIMAG(ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2))) )
       lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
            *(  REAL(psin(ix,iy,iz,1)*CONJG(ps1(ix,iy,iz,2))) &
            - REAL(psin(ix,iy,iz,2)*CONJG(ps1(ix,iy,iz,1))) )
    END FORALL
    !***********************************************************************
    ! y-derivatives
    !***********************************************************************
    IF(TFFT) THEN
       CALL cdervy(psin,ps1)  
    ELSE
       CALL cmuly(der1y,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,2,iq)=lcurrent(:,:,:,2,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
       lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
            *AIMAG( ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,1)) &
            -ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,2)) )
       lsodens(ix,iy,iz,3,iq)=lsodens(ix,iy,iz,3,iq)-weight &
            *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
            +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
    END FORALL
    !***********************************************************************
    ! z-derivatives
    !***********************************************************************
    IF(TFFT) THEN
       CALL cdervz(psin,ps1)  
    ELSE
       CALL cmulz(der1z,psin,ps1,0)  
    ENDIF
    ltau(:,:,:,iq)=ltau(:,:,:,iq)+weight* &
         (ps1(:,:,:,1)*CONJG(ps1(:,:,:,1))+ps1(:,:,:,2)*CONJG(ps1(:,:,:,2)))
    lcurrent(:,:,:,3,iq)=lcurrent(:,:,:,3,iq)+weight* &
         AIMAG(ps1(:,:,:,1)*CONJG(psin(:,:,:,1))+ &
         ps1(:,:,:,2)*CONJG(psin(:,:,:,2)))
    FORALL(ix=1:Nx,iy=1:ny,iz=1:nz)
       lsodens(ix,iy,iz,1,iq)=lsodens(ix,iy,iz,1,iq)+weight &
            *REAL( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
            -ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
       lsodens(ix,iy,iz,2,iq)=lsodens(ix,iy,iz,2,iq)+weight &
            *AIMAG( ps1(ix,iy,iz,2)*CONJG(psin(ix,iy,iz,1)) &
            +ps1(ix,iy,iz,1)*CONJG(psin(ix,iy,iz,2)) )
    END FORALL
    DEALLOCATE(ps1)
  END SUBROUTINE add_density
END MODULE Densities
