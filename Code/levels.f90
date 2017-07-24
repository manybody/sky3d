!------------------------------------------------------------------------------
! MODULE: Levels
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module is concerned with the wave function data: definition of
!!the pertinent arrays, allocating their storage and simple operations on them.
!------------------------------------------------------------------------------
MODULE Levels
  USE Params, ONLY: db,pi
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz,bangx,bangy,bangz
  USE Fourier
  IMPLICIT NONE
  SAVE
  INTEGER :: nstmax                        !<is the total number of wave functions
  !!present in the calculation. For the MPI version only \c nstloc
  !!wave functions are present on each node. Note that for all other
  !!wave-function related arrays, such as single-particle energies, the
  !!full set is stored on each node.
  INTEGER :: nstloc                         !<the number of wave functions stored 
  !!on the present node. In principle this should be defined in module 
  !!\c Parallel, but this would lead to a circular module dependence.
  INTEGER :: nneut                          !<the physical numbers of
  !!neutrons. These may be smaller than the number of single-particle states.
  INTEGER :: nprot                          !<the physical numbers of
  !!protons. These may be smaller than the number of single-particle states.
  !
  !>@name neutron and proton start and ending indices
  !!The neutron states are
  !!numbered <tt> npmin(1)</tt> through <tt> npsi(1)</tt> and the proton states
  !>run from <tt> npmin(2)</tt> through <tt> npsi(2)</tt>. Protons follow
  !>neutrons, so <tt> npmin(1)=1</tt> and <tt> npmin(2)=npsi(1)+1</tt>. 
  !><em>Note that for each particle type the number of states can be
  !>larger than the particle number, as states may be fractionally
  !>occupied or even empty.</em>  If initialization is not from fragments,
  !><tt> npsi(2)</tt> <em> as an input value</em> refers to the total number of
  !>proton states, it is later updated (in <tt> init.f90</tt>) to its normal
  !>meaning as the final index for proton states, which coincides with
  !>the total number of states, <tt> npsi(2)=nstmax</tt>.
  !>@{
  INTEGER :: npmin(2)                     
  INTEGER :: npsi(2)                        
  !>@}
  REAL(db) :: charge_number,&               !<the physical charge numbers.
              mass_number                   !<the physical mass numbers.
  COMPLEX(db),POINTER,CONTIGUOUS :: psi(:,:,:,:,:)     !<this is the main array for the wave functions.
  !!It has an additional last index counting the states. In the
  !!sequential case it runs from \f$ 1\ldots{\tt nstmax} \f$, in the MPI
  !!version each node has only \c nstloc wave functions.
  COMPLEX(db),POINTER,CONTIGUOUS :: hampsi(:,:,:,:,:)     !<this is the array for the Hamiltonian 
  !!applied to the wave functions.
  !!It has an additional last index counting the states. In the
  !!sequential case it runs from \f$ 1\ldots{\tt nstmax} \f$, in the MPI
  !!version each node has only \c nstloc wave functions.
  REAL(db),ALLOCATABLE :: sp_orbital(:,:)                !<dimensioned <tt> (3,nstmax)</tt>:
  !!expectation values of three components of single-particle orbital
  !!angular momentum, in units of \f$ \hbar \f$.
  REAL(db),ALLOCATABLE :: sp_spin(:,:)                   !<dimensioned <tt> (3,nstmax)</tt>:
  !!expectation values of three components of single-particle spin, in
  !!units of \f$ \hbar \f$.
  INTEGER, ALLOCATABLE :: isospin(:)                     !<keeps track of isospin of particle,
  !!1=neutron, 2=proton.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_energy       !<single-particle energy in MeV.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_efluct1      !<single-particle energy fluctuation in MeV calculated as
  !!\f[ \sqrt{\langle\psi|\hat h^2|\psi\rangle-\langle\psi|\hat h|\psi\rangle^2}. \f] 
  !!Used only as informational printout in the static part.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_kinetic      !<single-particle kinetic energy in units of MeV.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_norm         !<norm of single-particle wave function; should be unity normally.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_efluct2      !<single-particle energy fluctuation in MeV calculated as
  !!\f[ \sqrt{\langle\hat h\psi|\hat h\psi\rangle-\langle\psi|\hat h|\psi\rangle^2}. \f] 
  !!Used only as informational printout in the static part.
  REAL(db), ALLOCATABLE, DIMENSION(:) :: sp_parity       !<single-particle parity w.r.t.
  !!three-dimensional reflection at the origin; calculated as
  !!\f[ \sum_s\int\D^3r\, \psi^*(\vec r,s)\psi_s(-\vec r,s). \f]
  REAL(db), ALLOCATABLE, DIMENSION(:) :: wocc            !<occupation probability of single-particle
  !!state, may be fractional because of pairing. In the equations this
  !!is usually denoted as \f$ w_k^2 \f$, the square added because of the pairing notation.
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_levels
!> @brief
!!This subroutine allocates all the arrays associated with
!!single-particle wave functions. 
!>
!> @details
!!Note that while most have dimension
!!\c nstmax, \c psi itself is dimensioned for the number \c nstloc of wave 
!!functions on one specific processor. It also records the isospin value.
!--------------------------------------------------------------------------- 
  SUBROUTINE alloc_levels
    ALLOCATE(psi(nx,ny,nz,2,nstloc),hampsi(nx,ny,nz,2,nstloc), &
         sp_energy(nstmax),sp_efluct1(nstmax),sp_kinetic(nstmax),& 
         sp_norm(nstmax),sp_efluct2(nstmax),sp_parity(nstmax), &
         sp_orbital(3,nstmax),sp_spin(3,nstmax),wocc(nstmax), &
         isospin(nstmax))
    isospin(1:npsi(1))=1
    isospin(npmin(2):npsi(2))=2
  END SUBROUTINE alloc_levels
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervx
!> @brief
!!This routine calculates
!!derivatives of wave functions using the FFT method, in the \f$ x \f$-direction. 
!>
!> @details
!!The last argument can be omitted
!!and no second derivative is calculated in this case. The derivatives
!!add the proper dimension of \f$ {\rm fm}^{-1} \f$ and \f$ {\rm fm}^{-2} \f$,
!!respectively. Note the dependence of the \f$ k \f$ -value on index.
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!--------------------------------------------------------------------------- 
  SUBROUTINE cdervx0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: ix
    kfac=(PI+PI)/(dx*nx)
    CALL dfftw_execute_dft(xforward,psin,d1psout)
    IF(PRESENT(d2psout)) THEN
       DO ix=1,nx/2
          d2psout(ix,:,:,:)=-((ix-1)*kfac)**2*d1psout(ix,:,:,:)/REAL(nx)
          d2psout(nx-ix+1,:,:,:)=-(ix*kfac)**2*d1psout(nx-ix+1,:,:,:)/REAL(nx)
       ENDDO
       CALL dfftw_execute_dft(xbackward,d2psout,d2psout)
    ENDIF
    d1psout(1,:,:,:)=(0.D0,0.D0)
    DO ix=2,nx/2
       d1psout(ix,:,:,:)=CMPLX(0.0D0,((ix-1)*kfac),db)*d1psout(ix,:,:,:)/REAL(nx)
    ENDDO
    DO ix=1,nx/2-1
       d1psout(nx-ix+1,:,:,:)=CMPLX(0.0D0,-(ix*kfac),db)*d1psout(nx-ix+1,:,:,:) &
            /REAL(nx)
    ENDDO
    d1psout(nx/2+1,:,:,:)=(0.D0,0.D0)
    CALL dfftw_execute_dft(xbackward,d1psout,d1psout)
  END SUBROUTINE cdervx0
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervy
!> @brief
!!This routine is a wrapper function for cdervx0. 
!>
!> @details
!!It calculates the corrections for derivatives with TABC
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!---------------------------------------------------------------------------   
  SUBROUTINE cdervx(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangx)>0.00001) THEN
       DO i=1,nx
          p_temp(i,:,:,:)=psin(i,:,:,:)*exp(CMPLX(0.0d0,-bangx/nx*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervx0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervx0(p_temp,d1psout)
    END IF
    IF(abs(bangx)>0.00001) THEN
       DO i=1,nx
          d1psout(i,:,:,:)=d1psout(i,:,:,:)*exp(CMPLX(0.0d0,bangx/nx*i,db))&
               +CMPLX(0.0d0,bangx/(nx*dx),db)*psin(i,:,:,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,nx
             d2psout(i,:,:,:)=d2psout(i,:,:,:)*exp(CMPLX(0.0d0,bangx/nx*i,db))&
                  +(bangx/(nx*dx))**2*psin(i,:,:,:)&
                  +CMPLX(0.0d0,2.0d0*bangx/(nx*dx),db)*d1psout(i,:,:,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervx
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervy
!> @brief
!!This routine calculates
!!derivatives of wave functions using the FFT method, in the \f$ y \f$-direction. 
!>
!> @details
!!The last argument can be omitted
!!and no second derivative is calculated in this case. The derivatives
!!add the proper dimension of \f$ {\rm fm}^{-1} \f$ and \f$ {\rm fm}^{-2} \f$,
!!respectively. Note the dependence of the \f$ k \f$ -value on index.
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!--------------------------------------------------------------------------- 
  SUBROUTINE cdervy0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: iy,is,k
    kfac=(PI+PI)/(dy*ny)
!guru    CALL dfftw_execute_dft(yforward,psin,d1psout)
    DO is=1,2
       DO k=1,nz
          CALL dfftw_execute_dft(yforward,psin(:,:,k,is),d1psout(:,:,k,is))
       END DO
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iy=1,ny/2
          d2psout(:,iy,:,:)=-((iy-1)*kfac)**2*d1psout(:,iy,:,:)/REAL(ny)
          d2psout(:,ny-iy+1,:,:)=-(iy*kfac)**2*d1psout(:,ny-iy+1,:,:)/REAL(ny)
       ENDDO
!guru       CALL dfftw_execute_dft(ybackward,d2psout,d2psout)
       DO is=1,2
          DO k=1,nz
             CALL dfftw_execute_dft(ybackward,d2psout(:,:,k,is),d2psout(:,:,k,is))
          END DO
       END DO
    ENDIF
    d1psout(:,1,:,:)=(0.D0,0.D0)
    DO iy=2,ny/2
       d1psout(:,iy,:,:)=CMPLX(0.0D0,((iy-1)*kfac),db)*d1psout(:,iy,:,:)/REAL(ny)
    ENDDO
    DO iy=1,ny/2-1
       d1psout(:,ny-iy+1,:,:)=CMPLX(0.0D0,-(iy*kfac),db)*d1psout(:,ny-iy+1,:,:) &
            /REAL(ny)
    ENDDO
    d1psout(:,ny/2+1,:,:)=(0.D0,0.D0)
!guru    CALL dfftw_execute_dft(ybackward,d1psout,d1psout)
    DO is=1,2
       DO k=1,nz
          CALL dfftw_execute_dft(ybackward,d1psout(:,:,k,is),d1psout(:,:,k,is))
       END DO
    END DO
  END SUBROUTINE cdervy0
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervy
!> @brief
!!This routine is a wrapper function for cdervy0. 
!>
!> @details
!!It calculates the corrections for derivatives with TABC
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!--------------------------------------------------------------------------- 
  SUBROUTINE cdervy(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangy)>0.00001) THEN
       DO i=1,ny
          p_temp(:,i,:,:)=psin(:,i,:,:)*exp(CMPLX(0.0d0,-bangy/ny*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervy0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervy0(p_temp,d1psout)
    END IF
    IF(abs(bangy)>0.00001) THEN
       DO i=1,ny
          d1psout(:,i,:,:)=d1psout(:,i,:,:)*exp(CMPLX(0.0d0,bangy/ny*i,db))&
               +CMPLX(0.0d0,bangy/(ny*dy),db)*psin(:,i,:,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,ny
             d2psout(:,i,:,:)=d2psout(:,i,:,:)*exp(CMPLX(0.0d0,bangy/ny*i,db))&
                  +(bangy/(ny*dy))**2*psin(:,i,:,:)&
                  +CMPLX(0.0d0,2.0d0*bangy/(ny*dy),db)*d1psout(:,i,:,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervy
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervz
!> @brief
!!This routine calculates
!!derivatives of wave functions using the FFT method, in the \f$ z \f$-direction. 
!>
!> @details
!!The last argument can be omitted
!!and no second derivative is calculated in this case. The derivatives
!!add the proper dimension of \f$ {\rm fm}^{-1} \f$ and \f$ {\rm fm}^{-2} \f$,
!!respectively. Note the dependence of the \f$ k \f$ -value on index.
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!--------------------------------------------------------------------------- 
  SUBROUTINE cdervz0(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:)  
    REAL(db) :: kfac
    INTEGER :: iz,is
    kfac=(PI+PI)/(dz*nz)
!guru    CALL dfftw_execute_dft(zforward,psin,d1psout)
    DO is=1,2
       CALL dfftw_execute_dft(zforward,psin(:,:,:,is),d1psout(:,:,:,is))
    END DO
    IF(PRESENT(d2psout)) THEN
       DO iz=1,nz/2
          d2psout(:,:,iz,:)=-((iz-1)*kfac)**2*d1psout(:,:,iz,:)/REAL(nz)
          d2psout(:,:,nz-iz+1,:)=-(iz*kfac)**2*d1psout(:,:,nz-iz+1,:)/REAL(nz)
       ENDDO
!guru       CALL dfftw_execute_dft(zbackward,d2psout,d2psout)
       DO is=1,2
          CALL dfftw_execute_dft(zbackward,d2psout(:,:,:,is),d2psout(:,:,:,is))
       END DO
    ENDIF
    d1psout(:,:,1,:)=(0.D0,0.D0)
    DO iz=2,nz/2
       d1psout(:,:,iz,:)=CMPLX(0.0D0,((iz-1)*kfac),db)*d1psout(:,:,iz,:)/REAL(nz)
    ENDDO
    DO iz=1,nz/2-1
       d1psout(:,:,nz-iz+1,:)=CMPLX(0.0D0,-(iz*kfac),db)*d1psout(:,:,nz-iz+1,:) &
            /REAL(nz)
    ENDDO
    d1psout(:,:,nz/2+1,:)=(0.D0,0.D0)
!guru    CALL dfftw_execute_dft(zbackward,d1psout,d1psout)
    DO is=1,2
       CALL dfftw_execute_dft(zbackward,d1psout(:,:,:,is),d1psout(:,:,:,is))
    END DO
  END SUBROUTINE cdervz0
!---------------------------------------------------------------------------  
! DESCRIPTION: cdervy
!> @brief
!!This routine is a wrapper function for cdervz0. 
!>
!> @details
!!It calculates the corrections for derivatives with TABC
!>
!> @param[in] psin
!> COMPLEX(db), array, takes the wave function to be differantiated.
!> @param[out] d1psout
!> COMPLEX(db), array, returns the first derivative.
!> @param[out] d2psout
!> COMPLEX(db), array, OPTIONAL, returns the first derivative.
!--------------------------------------------------------------------------- 
  SUBROUTINE cdervz(psin,d1psout,d2psout)  
    COMPLEX(db), INTENT(IN) ::  psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT):: d1psout(:,:,:,:)
    COMPLEX(db), INTENT(OUT), OPTIONAL :: d2psout(:,:,:,:) 
    COMPLEX(db),ALLOCATABLE :: p_temp(:,:,:,:)
    INTEGER :: shapep(4),i
    shapep=shape(psin)
    ALLOCATE(p_temp(shapep(1),shapep(2),shapep(3),shapep(4))) 
    IF(abs(bangz)>0.00001) THEN
       DO i=1,nz
          p_temp(:,:,i,:)=psin(:,:,i,:)*exp(CMPLX(0.0d0,-bangz/nz*i,db))
       END DO
    ELSE
       p_temp=psin
    END IF
    IF(PRESENT(d2psout)) THEN
       CALL cdervz0(p_temp,d1psout,d2psout)
    ELSE
       CALL cdervz0(p_temp,d1psout)
    END IF
    IF(abs(bangz)>0.00001) THEN
       DO i=1,nz
          d1psout(:,:,i,:)=d1psout(:,:,i,:)*exp(CMPLX(0.0d0,bangz/nz*i,db))&
               +CMPLX(0.0d0,bangz/(nz*dz),db)*psin(:,:,i,:)
       END DO
       IF(PRESENT(d2psout)) THEN
          DO i=1,nz
             d2psout(:,:,i,:)=d2psout(:,:,i,:)*exp(CMPLX(0.0d0,bangz/nz*i,db))&
                  +(bangz/(nz*dz))**2*psin(:,:,i,:)&
                  +CMPLX(0.0d0,2.0d0*bangz/(nz*dz),db)*d1psout(:,:,i,:)
          END DO
       END IF
    END IF
    DEALLOCATE(p_temp)
  END SUBROUTINE cdervz
!---------------------------------------------------------------------------  
! DESCRIPTION: laplace
!> @brief
!!depending on the presence of the third argument \c e0inv, it can
!!calculate two things using FFT methods:
!!
!!  - if \c e0inv is not present, it calculates the Laplacian
!!    \f[\frac{\partial^2\psi}{\partial x^2}+\frac{\partial^2\psi}
!!    {\partial y^2}+\frac{\partial^2\psi} {\partial z^2}. \f]
!!  - if \c e0inv is present and positive, it calculates
!!    \f[ \frac1{E_{0\rm inv}+\hat t}\psi, \f]
!!    with the kinetic-energy operator
!!    \f[ \hat t=-\frac{\hbar^2}{2m}\left(\frac{\partial^2}
!!    {\partial x^2}+\frac{\partial^2} {\partial
!!    y^2}+\frac{\partial^2}{\partial z^2}\right), \f]
!!    which is of course expressed through \f$ \vec k \f$ in momentum space.
!!  .
!!The methods used in this routine are similar to those for the
!!derivatives \c cdervx etc.
!> @param[in] psin
!> COMPLEX(db), array, takes wave function on which operation is performed.
!> @param[out] psout
!> COMPLEX(db), array, returns the resulting wave function.
!> @param[in] e0inv
!> REAL(db), OPTIONAL, takes the damping factor.
!--------------------------------------------------------------------------- 
  SUBROUTINE laplace(psin,psout,e0inv)  
    USE Forces, ONLY: h2ma
    USE Grids, ONLY: dx,dy,dz
    COMPLEX(db), INTENT(IN)   :: psin(:,:,:,:)
    COMPLEX(db), INTENT(OUT)  :: psout(:,:,:,:)
    REAL(db), INTENT(IN), OPTIONAL :: e0inv
    REAL(db) :: kfacx, kfacy, kfacz
    REAL(db) :: k2facx(nx),k2facy(ny),k2facz(nz)
    INTEGER :: ix, iy, iz, is
    IF(.NOT.PRESENT(e0inv).AND.(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001))&
      STOP 'Laplace does not work with Bloch boundaries'
    kfacz=(PI+PI)/(dz*nz)
    kfacy=(PI+PI)/(dy*ny)
    kfacx=(PI+PI)/(dx*nx)
    DO iz=1,nz/2
       k2facz(iz)=-((iz-1)*kfacz)**2
       k2facz(nz-iz+1)=-(iz*kfacz)**2
    ENDDO
    DO iy=1,ny/2
       k2facy(iy)=-((iy-1)*kfacy)**2
       k2facy(ny-iy+1)=-(iy*kfacy)**2
    ENDDO
    DO ix=1,nx/2
       k2facx(ix)=-((ix-1)*kfacx)**2
       k2facx(nx-ix+1)=-(ix*kfacx)**2
    ENDDO
    IF(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          psout(ix,iy,iz,:)=psin(ix,iy,iz,:)*exp(CMPLX(0.0d0,-bangx/nx*ix,db))*exp(CMPLX(0.0d0,-bangy/ny*iy,db))*&
               exp(CMPLX(0.0d0,-bangz/nz*iz,db))
       END FORALL
    ELSE
       psout=psin
    END IF
    DO is=1,2
       CALL dfftw_execute_dft(pforward,psout(:,:,:,is),psout(:,:,:,is))
    END DO
    IF(PRESENT(e0inv)) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz,is=1:2)
          psout(ix,iy,iz,is)=psout(ix,iy,iz,is)  &
               /(e0inv-h2ma*(k2facx(ix)+k2facy(iy)+k2facz(iz))) &
               / DBLE(nx*ny*nz)
       END FORALL
    ELSE
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz,is=1:2)
          psout(ix,iy,iz,is)=psout(ix,iy,iz,is)  &
               *(k2facx(ix)+k2facy(iy)+k2facz(iz)) &
               /DBLE(nx*ny*nz)
       END FORALL
    ENDIF
    DO is=1,2
       CALL dfftw_execute_dft(pbackward,psout(:,:,:,is),psout(:,:,:,is))
    END DO
    IF(bangx>0.000001 .OR. bangy>0.000001 .OR. bangz>0.000001) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          psout(ix,iy,iz,:)=psout(ix,iy,iz,:)*exp(CMPLX(0.0d0,bangx/nx*ix,db))*exp(CMPLX(0.0d0,bangy/ny*iy,db))*&
               exp(CMPLX(0.0d0,bangz/nz*iz,db))
       END FORALL
    END IF
  END SUBROUTINE laplace
END MODULE Levels
