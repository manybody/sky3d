!------------------------------------------------------------------------------
! MODULE: Static
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the routines needed for the static iterations. All the 
!!logic needed for this case is concentrated here; all other modules except 
!!for \c External are equally used in the dynamic calculation.
!------------------------------------------------------------------------------
MODULE Static
  USE Params
  USE Densities
  USE Meanfield, ONLY: skyrme, hpsi, upot, bmass
  USE Levels
  USE Grids
  USE Moment
  USE Energies
  USE Inout, ONLY: write_wavefunctions, write_densities, plot_density, &
       sp_properties,start_protocol
  USE Pairs, ONLY: pair,epair,avdelt
  IMPLICIT NONE
  LOGICAL  :: tdiag=.FALSE.    !< if \c true, there is a diagonalization of
  !!the Hamiltonian during the later (after the 20th) static iterations.
  !!The 20 is hard coded in \c static.f90.
  LOGICAL  :: tlarge=.FALSE.   !< if \c true, during the diagonalization
  !!the new wave functions are temporarily written on disk to avoid
  !doubling the memory requirements.
  LOGICAL  :: tvaryx_0=.FALSE. !< it <tt>.TRUE.</tt> the parameter \f$ x_0 \f$
  !!is changed in every iteration in order to achieve faster convergence.
  INTEGER  :: maxiter          !< maximum number of iterations allowed.
  REAL(db) :: radinx           !< radius parameters (x-direction) used in the 
  !!harmonic-oscillator initialization (see subroutine \c harmosc).
  REAL(db) :: radiny           !< radius parameters (y-direction) used in the 
  !!harmonic-oscillator initialization (see subroutine \c harmosc).
  REAL(db) :: radinz           !< radius parameters (z-direction) used in the 
  !!harmonic-oscillator initialization (see subroutine \c harmosc).
  REAL(db) :: serr             !< convergence criterion. Iterations are stopped
  !!if the fluctuation in single-particle energies falls below this
  !!value (see near the end of subroutines {\c statichf).
  REAL(db) :: delesum          !< used to sum up the changes in
  !!single-particle energies during one iteration; it is calculated in
  !!\c statichf but printed in \c sinfo so that it is a module variable.
  REAL(db) :: x0dmp=0.2D0      !< corresponds to the parameter \f$ x_0 \f$ appearing 
  !!in the damped gradient iteration.
  REAL(db) :: e0dmp=100.D0     !< corresponds to the parameter \f$ E_0 \f$ appearing 
  !!in the damped gradient iteration.
  REAL(db) :: x0dmpmin=0.2d0   !< corresponds to minimal value of \f$ x_0 \f$ appearing 
  !!in the damped gradient iteration, if tvaryx_0 is set to <tt>.TRUE.</tt>.
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: getin_static
!> @brief
!!This routine reads the input for the static calculation using namelist
!!\c static. 
!>
!> @details
!!This includes the module variables of this module, but
!!also the numbers of particles, which are input quantities only in the
!!static mode if initialization is not done via fragment files. In the
!!case of user initialization they are also needed to correctly allocate
!!the fields; only the wave functions are then calculated in 
!!\c user_init. The values given in the input are overwritten by
!!fragment data otherwise.
!!
!!Thus only if <tt>nof<=0</tt> the input numbers are used. If \c npsi is
!!not given in the input, the values of \c nneut and \c nprot are
!!used for the number of wave functions, except for the pairing case,
!!when they are computed from a formula.
!!
!!The variables \c charge_number and \c mass_number are also set for this case.
!--------------------------------------------------------------------------- 
  SUBROUTINE getin_static
    NAMELIST/static/ tdiag,tlarge,maxiter, &
         radinx,radiny,radinz,serr,x0dmp,e0dmp,nneut,nprot,npsi,tvaryx_0
    npsi=0
    READ(5,static)
    IF(nof<=0) THEN
       IF(npsi(1)==0) THEN  
          IF(ipair==0.OR.nof<0) THEN  
             npsi(1)=nneut  
          ELSE  
             npsi(1)=NINT(nneut+1.65*FLOAT(nneut)**0.666667D0)  
             IF(MOD(npsi(1),2)/=0) npsi(1)=npsi(1)+1
          ENDIF
       ENDIF
       IF(npsi(2)==0) THEN  
          IF(ipair==0.OR.nof<0) THEN  
             npsi(2)=nprot  
          ELSE  
             npsi(2)=NINT(nprot+1.65*FLOAT(nprot)**0.666667D0)  
             IF(MOD(npsi(2),2)/=0) npsi(2)=npsi(2)+1
          ENDIF
       ENDIF
       IF(nneut>npsi(1).OR.nprot>npsi(2)) & 
            STOP 'Particle & state numbers in conflict'
       nstmax=npsi(1)+npsi(2)
       charge_number=nprot  
       mass_number=nneut+nprot  
       x0dmpmin=x0dmp
       WRITE(*,*) "x0dmpmin=", x0dmpmin
    END IF
  END SUBROUTINE getin_static
!---------------------------------------------------------------------------  
! DESCRIPTION: init_static
!> @brief
!!This subroutine essentially just prints the static input and then
!!initializes the header files with their header lines. 
!>
!> @details
!!This is not included in \c getin_static, because the particle 
!!and state numbers may have been changed by fragment input.
!!
!!In addition, the damping matrices are constructed by calling 
!!\c setup_damping. Finally for some Skyrme forces the effective mass
!!is changed to account for the center-of-mass correction. This should
!!only be used if necessary and not in dynamic calculations.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_static
    IF(wflag) THEN
       WRITE(*,*) '***** Parameters for static calculation *****'
       WRITE(*,"(3(A,I4))") ' Neutrons:',nneut,' in ',npsi(1),' levels'
       WRITE(*,"(3(A,I4))") ' Protons :',nprot,' in ',npsi(2)-npsi(1),' levels'
       WRITE(*,"(A,I6)") " Maximum number of iterations: ",maxiter  
       WRITE(*,"(2(A,G12.5))") ' Damping coefficient:',x0dmp, &
            " Damping energy scale: ",e0dmp
       WRITE(*,"(A,1PE12.4)") " Convergence limit: ",serr  
       ! initialize *.res files
       CALL start_protocol(converfile, &
            '# Iter   Energy  d_Energy    h**2        h*h        rms      r^2     r^3     r^4    &
            &beta2   gamma      Avdelta(N)   Avdelta(P)')
       CALL start_protocol(dipolesfile, &
            '# Iter    c.m. x-y-z                                  Isovector&
            &dipoles x-y-z')
       CALL start_protocol(spinfile, &
            '# Iter      Lx        Ly        Lz        Sx        Sy        &
            &Sz        Jx        Jy        Jz')
    ENDIF
    ! calculate damping matrices
    IF(e0dmp>0.0D0) CALL setup_damping(e0dmp)
    ! c.m. fixing term
    IF(f%zpe==0) THEN
       f%h2m=f%h2m*(mass_number-1.0D0)/mass_number
       WRITE(*,*) '***** Nucleon mass modified for z.p.e. correction'
    END IF
  END SUBROUTINE init_static
!---------------------------------------------------------------------------  
! DESCRIPTION: statichf
!> @brief
!!This is the principal routine for the static iterations. It applies
!!the gradient step repeatedly until convergence is achieved or the
!!maximum number of iterations is reached. 
!>
!> @details
!!The following local variables
!!are worth defining:
!!  - <b>\c sumflu</b> keeps track of the fluctuations 
!!    in single-particle energies summed over the states. It is used as
!!    the convergence criterion.
!!  - <b>\c addnew, \c addco</b> are simple factors with standard
!!    values 0.2 and 0.8 used for relaxation (see Step 8).  They are
!!    defined as parameter variables so they can be changed easily if
!!    desired.
!!  - <b>\c denerg</b> is used as an argument to \c grstep to contain
!!    the relative change in energy of the single-particle state.
!!  .
!!  - <b> Step 1: Initialization</b>: if this is a restart, the number
!!    of the initial iteration is set to the value of <tt> iter+1</tt>
!!    obtained from \c wffile. In this case the single-particle
!!    quantities do not have to be set to zero and orthogonalization is
!!    not necessary.  If this is not a restart, the initialization is done
!!    as zeroth iteration and the first iteration number for the loop is
!!    set to one.  some variables are initialized to zero and the matrix
!!    for the diagonalization \c hmatr is allocated. Then \c schmid
!!    is called for initial orthogonalization.
!!  - <b> Step2: calculating densities and mean fields</b>: the
!!    densities are reset to zero and then in a loop over states the
!!    contributions of the single-particle states are added up. The
!!    subroutine \c skyrme is called to compute the mean-field
!!    components.
!!  - <b> Step 3: initial gradient step</b>: in a loop over the
!!    single-particle wave functions the gradient step is applied - 
!!    see subroutine \c grstep. The
!!    sum of relative changes in single-particle energies and fluctuations
!!    are accumulated in \c delesum and \c sumflu. This loop is
!!    followed by the pairing calculation (which needs the single-particle
!!    energies calculated in the gradient step) and renewed
!!    orthogonalization. Then the detailed single-particle properties and
!!    energies are calculated using \c sp_properties} and \c sinfo.
!!  - <b> Step 4: start iteration</b>: this is
!!    the principal loop for the static calculation. The iteration number
!!    is printed.
!!  - <b> Step 5: gradient step</b>: this is identical to the initial
!!    gradient step in "Step 3".
!!  - <b> Step 6: diagonalization</b>: after 20 iterations and only if
!!    the switch \c tdiag is true, \c diagstep is called to
!!    diagonalize the single-particle Hamiltonian.
!!  - <b> Step 7: pairing and orthogonalization</b>: these are called for
!!    the new wave functions.
!!  - <b> Step 8: calculate densities and fields with relaxation</b>: the
!!    old density \c rho and kinetic energy density \c tau are saved
!!    in \c upot and \c bmass, which are here used purely as work
!!    arrays. Then the new densities are accumulated from the wave
!!    functions and for \c rho and \c tau they are mixed with the
!!    old densities in a ratio given by \c addnew and \c addco, if
!!    this is turned on by \c taddnew. \c skyrme is called to
!!    calculate the new fields. Then the detailed single-particle
!!    properties and energies are calculated and printed using 
!!    \c sp_properties} and \c sinfo.
!!  - <b> Step 9: finalizing the loop</b>: convergence is checked by
!!    comparing \c sumflu per particle to \c serr, if it is smaller,
!!    the job terminates after writing the final wave functions.
!!    Otherwise, the wave functions are written if indicated by
!!    \c mrest.
!!  . 
!!and the loop continues.
!--------------------------------------------------------------------------- 
  SUBROUTINE statichf
    LOGICAL, PARAMETER :: taddnew=.TRUE. ! mix old and new densities
    INTEGER :: iq,nst,firstiter
    REAL(db) :: sumflu,denerg
    REAL(db),PARAMETER :: addnew=0.2D0,addco=1.0D0-addnew      
    ! Step 1: initialization
    IF(tdiag) ALLOCATE(hmatr(nstmax,nstmax))
    IF(trestart) THEN
       firstiter=iter+1
    ELSE
       iter=0
       firstiter=1
       sp_energy=0.0D0  
       sp_efluct1=0.0D0  
       sp_efluct2=0.D0
       sp_norm=0.0D0  
       sumflu=0.D0
       CALL schmid
    END IF
    ! Step 2: calculate densities and mean field
    rho=0.0D0
    tau=0.0D0
    current=0.0D0
    sdens=0.0D0
    sodens=0.0D0
    DO nst=1,nstmax
       CALL add_density(isospin(nst),wocc(nst),psi(:,:,:,:,nst), &
            rho,tau,current,sdens,sodens)  
    ENDDO
    CALL skyrme
    ! Step 3: initial gradient step
    delesum=0.0D0  
    sumflu=0.0D0  
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
    !$OMP    SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
    DO nst=1,nstmax
       CALL grstep(nst,isospin(nst),sp_energy(nst),denerg, &
            psi(:,:,:,:,nst))
       sumflu=sumflu+wocc(nst)*sp_efluct1(nst)  
       delesum=delesum+wocc(nst)*denerg  
    ENDDO
    !$OMP END PARALLEL DO
    ! pairing and orthogonalization
    IF(ipair/=0) CALL pair
    CALL schmid
    ! produce and print detailed information
    CALL sp_properties
    CALL sinfo
    !set x0dmp to 3* its value to get faster convergence
    IF(tvaryx_0) x0dmp=3.0d0*x0dmp
    ! step 4: start static iteration loop
    Iteration: DO iter=firstiter,maxiter  
       WRITE(*,'(a,i6,a,F12.4)') ' Static Iteration No.',iter,'  x0dmp= ',x0dmp
       ! Step 5: gradient step
       delesum=0.0D0  
       sumflu=0.0D0
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
       !$OMP    SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
       DO nst=1,nstmax
          CALL grstep(nst,isospin(nst),sp_energy(nst),denerg, &
               psi(:,:,:,:,nst))
          sumflu=sumflu+wocc(nst)*sp_efluct1(nst)  
          delesum=delesum+wocc(nst)*denerg  
       ENDDO
       !$OMP END PARALLEL DO
       ! Step 6: diagonalize if desired
       IF(tdiag.AND.iter>20) THEN
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iq) SCHEDULE(STATIC) 
          DO iq=1,2
             CALL diagstep(iq,npsi(iq)-npmin(iq)+1)
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
       ! Step 7: do pairing and orthogonalization
       IF(ipair/=0) CALL pair
       CALL schmid
       ! Step 8: get new densities and fields with relaxation
       IF(taddnew) THEN
          upot=rho
          bmass=tau
       ENDIF
       rho=0.0D0
       tau=0.0D0
       current=0.0D0
       sdens=0.0D0
       sodens=0.0D0
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst) SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:rho, tau, current, sdens, sodens)
       DO nst=1,nstmax
          CALL add_density(isospin(nst),wocc(nst),psi(:,:,:,:,nst), &
               rho,tau,current,sdens,sodens)  
       ENDDO
       !$OMP END PARALLEL DO
       IF(taddnew) THEN
          rho=addnew*rho+addco*upot
          tau=addnew*tau+addco*bmass
       ENDIF
       CALL skyrme
       ! calculate and print information
       IF(mprint>0.AND.MOD(iter,mprint)==0) THEN
          CALL sp_properties
          CALL sinfo
       ELSEIF(tvaryx_0) THEN
          CALL sp_properties
          CALL sum_energy
       ENDIF
       ! Step 9: check for convergence, saving wave functions
       IF(sumflu/nstmax<serr.AND.iter>1) THEN
          CALL write_wavefunctions
          EXIT Iteration  
       END IF
       IF(MOD(iter,mrest)==0) THEN  
          CALL write_wavefunctions
       ENDIF
       ! Step 10: update step size for the next iteration
       IF(tvaryx_0) THEN
          IF(ehf<ehfprev .OR. efluct1<(efluct1prev*(1.0d0-1.0d-5)) &
               .OR. efluct2<(efluct2prev*(1.0d0-1.0d-5))) THEN
             x0dmp=x0dmp*1.005
          ELSE
             x0dmp=x0dmp*0.8
          END IF
          IF(x0dmp<x0dmpmin) x0dmp=x0dmpmin
          efluct1prev=efluct1
          efluct2prev=efluct2
          ehfprev=ehf
       END IF
    END DO Iteration
    IF(tdiag) DEALLOCATE(hmatr)
  END SUBROUTINE statichf
!---------------------------------------------------------------------------  
! DESCRIPTION: grstep
!> @brief
!!This subroutine applies the damped gradient step to a wave
!!function \c psi. Its index is \c nst and isospin \c iq - 
!!these data are needed for the construction of the
!!Hamiltonian matrix \c hmatr. The argument \c spe
!!\f$ \rightarrow\epsilon \f$ represents the single-particle energy which is
!!used a an energy shift in the calculation.
!>
!> @details
!!The work is done in the following steps:
!!  -# apply \f$ \hat h-\epsilon \f$ to the wave function \c psin to obtain
!!     \c ps1.
!!  -# calculate the diagonal matrix element of \f$ \hat h \f$
!!     \f[ {\tt xnormb}=\langle {\tt psin}|{\tt ps1}\rangle=\langle {\tt
!!     psin}|\hat h|{\tt psin}\rangle\f] 
!!     and the squared norm of \c psin (in \c xnorm). The expression 
!!     <tt>xnormb/xnorm</tt> thus corresponds
!!     to the expectation value \f$ \|\hat h\| \f$ in \c psin.
!!
!!     Then the matrix elements of \f$ \hat h \f$ with all other states of the same
!!     isospin are calculated and inserted into \c hmatr; in the diagonal
!!     matrix elements the energy shift \f$ \epsilon \f$ is added back.
!!  -# for the calculation of the fluctuation in the
!!     single-particle energy, the quantity
!!     \f[ {\tt exph2}=\langle {\tt psin}|\hat h^2|{\tt psin}\rangle \f]
!!     is used to compute
!!     \f[ {\tt sp\_efluct1}=\sqrt{\|\hat h^2\|-\|\hat h\|^2} \f]
!!     and the squared norm
!!     \f[ {\tt varh2}=\|\hat h\,{\tt psin}\|^2 \f]
!!     similarly for the second fluctuation measure
!!     \f[ {\tt sp\_efluct2}=\sqrt{\|\hat h\,{\tt psin}\|^2/\|{\tt psin}\|^2-\|\hat
!!     h\|^2}. \f]
!!     They are calculated only for time steps with output turned on.
!!  -# now the damping is performed. We first compute
!!     \f[ |{\tt ps1}\rangle-{\tt xnormb}\,|{\tt psin}\rangle=\left(\hat h
!!     -\langle {\tt psin}|\hat h|{\tt psin}\rangle\right)\,|{\tt
!!     psin}\rangle \f]
!!     replacing \c ps1, on which then the real damping
!!     operator acts. If \c FFT is being used for the derivatives, we use
!!     the routine \c laplace from module \c levels to compute
!!     \f[ \frac{x_{0\rm dmp}}{E_{0\rm inv}+\hat t}\,|{\tt ps1}\rangle. \f]
!!     If derivatives are to done by matrices, the damping matrices 
!!     \c cdmpx etc. are used (see subroutine \c setdmc in module \c Grids.  
!!     The factors \c x0dmp and \c e0dmp have to be
!!     manipulated a bit in this case.
!!     Finally we subtract this result multiplied by \c x0act from the
!!     original wave function to get the damped one.
!!
!!  -# the single-particle energy is calculated from its
!!     new expectation value with the energy shift restore, and the
!!     comparison with the initial value yields the relative change 
!!     \c denerg, which is passed back to the caller.
!>
!> @param[in] nst
!> INTEGER, takes the index of the wave function.
!> @param[in] iq
!> INTEGER, takes the isospin.
!> @param[in,out] spe
!> REAL(db), takes and returns single-particle energy. 
!> @param[out] denerg
!> REAL(db), returns the difference in energy of the particle before and 
!! after the iteration.
!> @param[in,out] psin
!> REAL(db), array, wave function to be iterated. 
!--------------------------------------------------------------------------- 
  SUBROUTINE grstep(nst,iq,spe,denerg,psin)
    USE Trivial, ONLY: cmulx,cmuly,cmulz,rpsnorm,overlap
    !***********************************************************************
    !                                                               
    !     grstep=one damped gradient iteration step for a given  
    !                 wave function psin with isospin iq.          
    !        psi=o[ psi - x0*damp*[(h-spe)psi] ]                
    !                                                              
    !***********************************************************************
    INTEGER,INTENT(IN) :: nst,iq
    REAL(db) :: spe,denerg
    COMPLEX(db) :: psin(:,:,:,:)  
    INTENT(OUT) :: denerg
    INTENT(INOUT) :: spe,psin
    REAL(db) :: x0act,esf,enrold,xnorm,xnormb,exph2,varh2
    COMPLEX(db) :: ps1(nx,ny,nz,2),ps2(nx,ny,nz,2)
    INTEGER :: nst2
    ! Step 1:(h-esf) on psin yields ps1.                        *
    esf=spe 
    CALL hpsi(iq,esf,psin,ps1)
    ! Step 2: calculate matrix elements
    xnorm=rpsnorm(psin)
    xnormb=overlap(psin,ps1)
    DO nst2=1,nstmax
       IF(tdiag.AND.isospin(nst2)==isospin(nst))   &
            hmatr(nst2,nst)=overlap(psi(:,:,:,:,nst2),ps1)
    ENDDO
    IF(tdiag) hmatr(nst,nst)=hmatr(nst,nst)+spe
    ! Step 3: calculate fluctuation, i.e. <h*h> and |h|**2
    IF((mprint>0.AND.MOD(iter,mprint)==0).OR.tvaryx_0) THEN  
       CALL hpsi(iq,esf,ps1,ps2)
       exph2=overlap(psin,ps2)
       varh2=rpsnorm(ps1)
       sp_efluct1(nst)=SQRT(ABS(exph2/xnorm-(xnormb/xnorm)**2))  
       sp_efluct2(nst)=SQRT(ABS(varh2/xnorm-(xnormb/xnorm)**2))  
    ENDIF
    ! Step 4: the damping step
    IF(e0dmp>0.0D0) THEN  
       ps1=ps1 - xnormb*psin
       x0act=x0dmp
       IF(TFFT) THEN
          CALL laplace(ps1,ps2,e0inv=e0dmp)  
       ELSE
          CALL cmulz(cdmpz,ps1,ps2,0)  
          CALL cmuly(cdmpy,ps2,ps1,0)  
          CALL cmulx(cdmpx,ps1,ps2,0)  
          x0act=x0act/e0dmp
       ENDIF
       psin=psin - x0act*ps2
    ELSE  
       psin=(1.0+x0dmp*xnormb)*psin-x0dmp*ps1
    ENDIF
    ! Step 5: energy convergence criterion
    enrold=spe  
    spe=xnormb+esf  
    denerg=(enrold-spe)/ABS(spe)  
  END SUBROUTINE grstep
!---------------------------------------------------------------------------  
! DESCRIPTION: diagstep
!> @brief
!!This subroutine performs a diagonalization of the single-particle
!!Hamiltonian using the LAPACK routine \c ZHEEVD. This is of course
!!done separately for protons and neutrons indexed by \c iq. The
!!matrix \c hmatr is produced in \c grstep.
!>
!> @details
!!We do not give excessive detail here but summarize the main
!!points, which should be easy to analyze. The steps are:
!!
!!  - <b>Step 1</b>: the matrix is copied into
!!    array \c unitary. The \c ZHEEVD is called, which leaves as
!!    main results the vector of eigenvalues \c eigen and the unitary
!!    matrix \c unitary describing the transformation from the
!!    original states to the diagonalized ones. The latter matrix is
!!    stored in lower-diagonal form.
!!  - <b>Step 2: transform states</b>: the matrix \c unitary is used
!!    to form the appropriate linear combinations of the original
!!    single-particle states. If \c tlarge is true, each state if
!!    formed independently and written onto a scratch file, otherwise an
!!    intermediate array \c ps1, dimensioned to hold either only
!!    protons or only neutrons (number of states through the argument 
!!    \c nlin), is allocated to receive the new states, which are then
!!    copied back into \c psi.
!>
!> @param[in] iq
!> INTEGER, takes tisospin.
!> @param[in] nlin
!> INTEGER, takes the number of particles with isospin \c iq.
!--------------------------------------------------------------------------- 
  SUBROUTINE diagstep(iq,nlin)
    !***********************************************************************
    !                                                                      *
    !     diagstep= diagonalize Hamiltonian matrix of active shells      *
    !                                                                      *
    !***********************************************************************
    INTEGER,INTENT(IN) :: iq,nlin
    INTEGER :: nst,nst2,noffset
    INTEGER :: infoconv
    REAL(db) :: eigen(nstmax)
    COMPLEX(db) :: unitary(nlin,nlin)
    COMPLEX(db), ALLOCATABLE :: psiw(:,:,:,:)          ! work space
    COMPLEX(db), ALLOCATABLE :: ps1(:,:,:,:,:)
    COMPLEX(db) :: cwork(2*nlin*nlin)
    REAL(db)    :: rwork(2*nlin*nlin+5*nlin+1)
    INTEGER     :: iwork(5*nlin+3)
    INTERFACE
       SUBROUTINE zheevd( jobz, uplo, n, a, lda, w, work, &
            lwork, rwork, lrwork, iwork, liwork, info )
         USE Params, ONLY: db
         CHARACTER(1) :: jobz, uplo
         INTEGER :: info, ldab, liwork, lrwork, lwork, n, iwork(*)
         DOUBLE PRECISION ::  rwork( * ), w( * )
         COMPLEX(8) :: a( lda, * ), work( * )
         INTENT(IN) :: jobz,uplo,n,lda,lwork,lrwork,liwork
         INTENT(INOUT) :: a
         INTENT(OUT) :: w,work,rwork,iwork,info
       END SUBROUTINE zheevd
    END INTERFACE
    ! Step 1: copy matrix, then diagonalize
    noffset=npmin(iq)-1
    unitary=hmatr(npmin(iq):npsi(iq),npmin(iq):npsi(iq))
    CALL ZHEEVD('V','L',nlin,unitary,nlin,eigen, &
         cwork,nlin*nlin*2,rwork,2*nlin*nlin+5*nlin+1, &
         iwork,5*nlin+3,infoconv)
    !  Step 2: transform states, replace original ones
    IF(tlarge) THEN
       OPEN(scratch,status='scratch',form='unformatted')
       ALLOCATE(psiw(nx,ny,nz,2))
       noffset=npmin(iq)-1
       DO nst=npmin(iq),npsi(iq)
          psiw=CMPLX(0.0D0,0.0D0,db)
          DO nst2=npmin(iq),npsi(iq)
             psiw(:,:,:,:)=psiw(:,:,:,:) &
                  + unitary(nst2-noffset,nst-noffset)*psi(:,:,:,:,nst2)
          ENDDO
          WRITE(scratch) psiw
       ENDDO
       REWIND(scratch)
       DO nst=npmin(iq),npsi(iq)
          READ(scratch) psi(:,:,:,:,nst)
       ENDDO
       DEALLOCATE(psiw)
       CLOSE(scratch)
    ELSE
       ALLOCATE(ps1(nx,ny,nz,2,nlin))
       noffset=npmin(iq)-1
       ps1=0.0D0
       DO nst=npmin(iq),npsi(iq)
          DO nst2=npmin(iq),npsi(iq)
             ps1(:,:,:,:,nst-noffset)=ps1(:,:,:,:,nst-noffset) &
                  + unitary(nst2-noffset,nst-noffset)*psi(:,:,:,:,nst2)
          ENDDO
       ENDDO
       psi(:,:,:,:,npmin(iq):npsi(iq))=ps1
       DEALLOCATE(ps1)
    ENDIF
  END SUBROUTINE diagstep
!---------------------------------------------------------------------------  
! DESCRIPTION: sinfo
!> @brief
!!This subroutine computes the data that are not needed for the
!!calculation itself, but only for informative output and writes them
!!onto the appropriate output files.
!>
!> @details
!!First \c moments, \c integ_energy, and
!!\c sum_energy are called to compute the relevant physical
!!quantities. Output lines are added to \c converfile,
!!\c dipolesfile, \c momentafile, and
!!\c spinfile. Information on the current energy \c ehf, the
!!total kinetic energy \c tke, the relative change in energy over
!!the last iteration \c delesum, the fluctuations in
!!single-particle energy, and the energy corrections are printed on
!!standard output.
!!
!!At intervals of \c mplot iterations the density printer plot is
!!produced and the <tt> *.tdd</tt> file containing the densities written.
!!
!!Finally, on standard output more detailed output is given: an overview
!!of the different contributions to the energy, a list of
!!single-particle state properties, and a listing of various moments
!!using \c moment_print.
!--------------------------------------------------------------------------- 
  SUBROUTINE sinfo
    INTEGER :: il
    CHARACTER(*),PARAMETER :: &
         header='  #  Par   v**2   var_h1   var_h2    Norm     Ekin    Energy &
         &    Lx      Ly      Lz     Sx     Sy     Sz  '   
    ! calculate static observables for printout                       *
    CALL moments(0,0,0)
    CALL integ_energy
    CALL sum_energy
    ! add information to summary files
    OPEN(unit=scratch,file=converfile,POSITION='APPEND')  
   WRITE(scratch,'(2x,i5,f9.2,3(1pg11.3),6(0pf15.8),f15.8,f15.8)') &
         iter,ehf,delesum/pnrtot,efluct1,efluct2,rmstot,r2tot,r3tot,r4tot,beta,gamma,avdelt(1),avdelt(2)
    CLOSE(scratch)
    OPEN(unit=scratch,file=dipolesfile, POSITION='APPEND')  
    WRITE(scratch,'(1x,i5,6E14.4)') iter,cmtot,cm(:,2)-cm(:,1)
    CLOSE(unit=scratch)
    OPEN(unit=scratch,file=spinfile, POSITION='APPEND')  
    WRITE(scratch,'(1x,i5,9F10.4)') iter,orbital,spin,total_angmom 
    CLOSE(unit=scratch)
    WRITE(*,'(/,A,I7,A/2(A,F12.4),A/(3(A,E12.5),A),3(A,E12.5))') &
         ' ***** Iteration ',iter,' *****',' Total energy: ',ehf,' MeV  Total kinetic energy: ', &
         tke,' MeV',' de/e:      ',delesum,'      h**2  fluct.:    ',efluct1, &
         ' MeV, h*hfluc.:    ',efluct2,' MeV', &
         ' MeV. Rearrangement E: ',e3corr,' MeV. Coul.Rearr.: ', &
         ecorc,' MeV x0dmp: ',x0dmp
    ! detail printout
    WRITE(*,'(/A)') ' Energies integrated from density functional:'
    WRITE(*,'(4(A,1PE14.6),A/26X,3(A,1PE14.6),A)') &
         ' Total:',ehfint,' MeV. t0 part:',ehf0,' MeV. t1 part:',ehf1, &
         ' MeV. t2 part:',ehf2,' MeV',' t3 part:',ehf3,' MeV. t4 part:',ehfls, &
         ' MeV. Coulomb:',ehfc,' MeV.'
    IF(ipair/=0) WRITE(*,'(2(A,1PE14.6))') ' Pairing energy neutrons: ', &
         epair(1),' protons: ',epair(2)
    ! output densities
    IF(mplot/=0) THEN  
       IF(MOD(iter,mplot)==0) THEN
          CALL plot_density
          CALL write_densities
       ENDIF
    ENDIF
    IF(.NOT.wflag) RETURN
    ! print details of s.p. levels
    WRITE(*,'(A)') ' Neutron Single Particle States:',header
    DO il=1,nstmax
       IF(il==npmin(2)) THEN
          WRITE(*,'(A)') ' Proton Single Particle States:',header  
       END IF
       WRITE(*,'(1X,I3,F4.0,F8.5,2F9.5,F9.6,F8.3,F10.3,3F8.3,3F7.3)') &
            il,sp_parity(il),wocc(il),sp_efluct1(il),sp_efluct2(il), &
            sp_norm(il),sp_kinetic(il),sp_energy(il), &
            sp_orbital(:,il),sp_spin(:,il)
    ENDDO
    CALL moment_print
  END SUBROUTINE sinfo
!---------------------------------------------------------------------------  
! DESCRIPTION: harmosc
!> @brief
!!This subroutine initializes all wave function with harmonic oscillator states. 
!>
!> @details
!!As a first step it calculates the first state, which is a Gaussian. In the 
!!second step the other wave functions are obtained by multiplying the Gaussian 
!!with polynomials. The subroutine does not return normalized states. In combination 
!!with \c schmid the normalized lowest oscillator states are obtained. 
!--------------------------------------------------------------------------- 
  SUBROUTINE harmosc
    USE Trivial, ONLY: rpsnorm
    REAL(db) :: xx,yy,zz,xx2,zz2,y2,anorm,temp
    INTEGER  :: nst,iq,is,ix,iy,iz,nps,i,j,k,ka,nshell(3,nstmax)
    wocc=0.D0
    wocc(1:nneut)=1.D0
    wocc(npmin(2):npmin(2)+nprot-1)=1.D0
    !*************************************************************************
    ! Lowest state: Gaussian
    !*************************************************************************
    nst=0  
    DO iq=1,2  
       IF(iq==1) THEN
          nps=npsi(1)
       ELSE
          nps=npsi(2)
       ENDIF
       ka_loop: DO ka=0,nps  
          DO k=0,nps  
             DO j=0,nps  
                DO i=0,nps  
                   IF(ka==i+j+k) THEN  
                      DO is=1,2
                         nst=nst+1  
                         IF(nst>nps) EXIT ka_loop
                         nshell(1,nst)=i  
                         nshell(2,nst)=j  
                         nshell(3,nst)=k  
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO ka_loop
       nst=nst-1
    END DO
    DO iq=1,2
       nst=npmin(iq)
       DO iz=1,nz  
          zz2=(z(iz)/radinz)**2
          DO ix=1,nx  
             xx2=(x(ix)/radinx)**2  
             DO iy=1,ny  
                y2=(y(iy)/radiny)**2  
                temp=xx2+y2+zz2
                psi(ix,iy,iz,1,nst)=EXP(-(temp))  
             ENDDO
          ENDDO
       ENDDO
       anorm=rpsnorm(psi(:,:,:,:,nst))
       psi(:,:,:,:,nst)=psi(:,:,:,:,nst)/SQRT(anorm)
       !*************************************************************************
       ! Higher states: lowest * polynomial
       !*************************************************************************
       DO nst=npmin(iq)+1,npsi(iq)  
          is=MOD(nst-npmin(iq),2)+1  
          DO iz=1,nz  
             IF(nshell(3,nst)/=0) THEN  
                zz=z(iz)**nshell(3,nst)  
             ELSE  
                zz=1.0D0  
             ENDIF
             DO iy=1,ny  
                IF(nshell(2,nst)/=0) THEN  
                   yy=y(iy)**nshell(2,nst)  
                ELSE  
                   yy=1.0D0  
                ENDIF
                DO ix=1,nx  
                   IF(nshell(1,nst)/=0) THEN  
                      xx=x(ix)**nshell(1,nst)  
                   ELSE  
                      xx=1.0D0  
                   ENDIF
                   psi(ix,iy,iz,is,nst)=psi(ix,iy,iz,1,npmin(iq))*xx*yy*zz
                ENDDO
             ENDDO
          ENDDO
       END DO
    END DO
    WRITE(*,*) "Harmonic oscillators widths (x-y-z):"
    WRITE(*,"(3F12.4)") radinx,radiny,radinz
    WRITE(*,*) '***** Harmonic oscillator initialization complete *****'
  END SUBROUTINE harmosc
END MODULE Static
