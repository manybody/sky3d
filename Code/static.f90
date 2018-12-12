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
  USE Parallel
  USE Constraint, ONLY: tconstraint,before_constraint,init_constraint,tune_constraint,add_constraint
  USE Temperature, ONLY: kbT,temp_dist
  USE Inout, ONLY: write_wavefunctions, write_densities, plot_density, &
       sp_properties,start_protocol
  USE Pairs, ONLY: pair,epair,avdelt,avdeltv2,avg,eferm,eferm_cutoff,partnum_cutoff,pairwg,deltaf
  USE Formfactor, ONLY: radius_print
  IMPLICIT NONE
  LOGICAL  :: tdiag=.FALSE.    !< if \c true, there is a diagonalization of
  !!the Hamiltonian during the later (after the 20th) static iterations.
  !!The 20 is hard coded in \c static.f90.
  LOGICAL  :: tlarge=.FALSE.   !< if \c true, during the diagonalization
  !!the new wave functions are temporarily written on disk to avoid
  !doubling the memory requirements.
  LOGICAL  :: tvaryx_0=.FALSE. !< it <tt>.TRUE.</tt> the parameter \f$ x_0 \f$
  !!is changed in every iteration in order to achieve faster convergence.
  LOGICAL :: ttime=.FALSE.     !< If <tt>.TRUE.</tt> and in MPI mode, outputs 
  !!timings for routines.
  INTEGER  :: maxiter          !< maximum number of iterations allowed.
  INTEGER  :: outerpot=0       !< dtermines how many iterations will be 
  !!performed with outer potential.
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
  CHARACTER(1) :: outertype='N'!<determines type of outer potential.
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
         radinx,radiny,radinz,serr,x0dmp,e0dmp,nneut,nprot,npsi,tvaryx_0,&
         outerpot,outertype,ttime,kbT
    npsi=0
    READ(5,static)

    IF(.NOT.tdiag) STOP "static code runs only with TDIAG=.TRUE."

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
    END IF
    CALL init_constraint()
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
  !***********************************************************************
  !begins protocols, inits damping and calculates zpe correction
  !***********************************************************************
    IF(wflag) THEN
       WRITE(*,*)
       WRITE(*,*) '***** Parameters for static calculation *****'
       WRITE(*,"(3(A,I4))") ' Neutrons:',nneut,' in ',npsi(1),' levels'
       WRITE(*,"(3(A,I4))") ' Protons :',nprot,' in ',npsi(2)-npsi(1),' levels'
       WRITE(*,"(A,I6)") " Maximum number of iterations: ",maxiter
       IF(tvaryx_0) THEN
       WRITE(*,"(2(A,G12.5))") ' Min. damping coefficient:',x0dmpmin, &
            " Damping energy scale: ",e0dmp
       ELSE
       WRITE(*,"(2(A,G12.5))") ' Damping coefficient:',x0dmp, &
            " Damping energy scale: ",e0dmp
       END IF
       WRITE(*,"(A,1PE12.4)") " Convergence limit: ",serr  
       ! initialize *.res files
       CALL start_protocol(converfile, &
            '# Iter   Energy  d_Energy    h**2        h*h        rms    &
            &beta2  gamma      x_0')
       CALL start_protocol(dipolesfile, &
            '# Iter    c.m. x-y-z                                  Isovector&
            &dipoles x-y-z')
       CALL start_protocol(spinfile, &
            '# Iter      Lx        Ly        Lz        Sx        Sy        &
            &Sz        Jx        Jy        Jz')
       CALL start_protocol(energiesfile, &
            '#Iter  N(n)      N(p)      E(sum)         E(integ)       &
            &Ekin           E_Coul         ehfCrho0       ehfCrho1       &
            &ehfCdrho0      ehfCdrho1      & 
            &ehfCtau0       ehfCtau1       ehfCdJ0        ehfCdJ1        Ecm_corr       -TS')
       IF(tabc_nprocs>1.AND.tabc_myid==0) CALL start_protocol(tabcfile, &
            '# Iter   Energy         E_kin          E_Coul         E_Skyrme ')
    ENDIF
    ! calculate damping matrices
    IF(e0dmp>0.0D0) CALL setup_damping(e0dmp)
    ! c.m. fixing term
    IF(f%zpe==0) THEN
       f%h2m=f%h2m*(mass_number-1.0D0)/mass_number
       WRITE(*,*) '***** Nucleon mass modified for z.p.e. correction'
    END IF
    IF(.NOT.ALLOCATED(pairwg)) THEN
      ALLOCATE(pairwg(nstmax))
      pairwg=1D0                        ! default if no soft cutoff is invoked
    END IF
    IF(.NOT.ALLOCATED(deltaf)) ALLOCATE(deltaf(nstmax))
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
    USE Linalg, ONLY: init_linalg
!    USE Constraint                   
    LOGICAL, PARAMETER   :: taddnew=.TRUE. ! mix old and new densities
    INTEGER              :: iq,nst,firstiter,number_threads
    REAL(db)             :: sumflu,denerg
    REAL(db) , PARAMETER :: addnew=0.2D0,addco=1.0D0-addnew  
    INTEGER,  EXTERNAL   :: omp_get_num_threads 
    !***********************************************************************
    !
    !performs static iterations
    !
    !***********************************************************************
    !  
    !***********************************************************************
    ! Step 1: initialization
    !*********************************************************************** 
    number_threads=1
    !$OMP PARALLEL
    !$ number_threads=OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    IF(wflag)WRITE(*,*)'number of threads= ',number_threads
    IF(wflag)WRITE(*,*)
    CALL init_linalg(0)
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
       IF(wflag)WRITE(*,'(A29)',advance="no") 'Initial orthogonalization... '
       IF(tmpi)THEN
          CALL diagstep(my_iso,.FALSE.)
       ELSE
         DO iq=1,2
           CALL diagstep(iq,.FALSE.)
          END DO
       END IF 
       IF(wflag)WRITE(*,*)'DONE'
    END IF
    !****************************************************  
    ! Step 2: calculate densities and mean field
    !****************************************************  
    rho=0.0D0
    tau=0.0D0
    current=0.0D0
    sdens=0.0D0
    sodens=0.0D0
    IF(wflag)WRITE(*,'(A25)',advance="no")'Initial add_density... '
    DO nst=1,nstloc
       CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)),&
                        psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
    ENDDO
    IF(tmpi) CALL collect_densities!sum densities over all nodes 
    IF(wflag)WRITE(*,*) 'DONE'
!
    IF(wflag)WRITE(*,'(A25)',advance="no")'Initial skyrme... '
    CALL skyrme(iter<=outerpot,outertype)
    IF(wflag)WRITE(*,*) 'DONE'
    !****************************************************  
    ! Step 3: initial gradient step
    !****************************************************  
    delesum=0.0D0  
    sumflu=0.0D0  
    IF(wflag)WRITE(*,'(A25)',advance="no")'Initial grstep... '
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
    !$OMP SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
    DO nst=1,nstloc
      CALL grstep(globalindex(nst),isospin(globalindex(nst)),sp_energy(globalindex(nst)),denerg,psi(:,:,:,:,nst))
      sumflu=sumflu+wocc(globalindex(nst))*sp_efluct1(globalindex(nst))  
      delesum=delesum+wocc(globalindex(nst))*denerg  
    ENDDO
    !$OMP END PARALLEL DO
    IF(tmpi) CALL collect_energies(delesum,sumflu)!collect fluctuations and change in energy 
    IF(wflag)WRITE(*,*) 'DONE'
    ! pairing and orthogonalization
    IF(kbT>0.0d0) CALL temp_dist
    IF(ipair/=0) CALL pair
    IF(kbT>0.0d0.AND.ipair/=0) STOP 'Pairing and finite temperature is not implemented'
    IF(wflag)WRITE(*,'(A25)',advance="no") 'Initial ortho2... '
    IF (my_iso>0) THEN
      CALL diagstep(my_iso,.FALSE.)
    ELSE
      DO iq=1,2
        CALL diagstep(iq,.FALSE.)    
      END DO 
    END IF
    IF(wflag)WRITE(*,*) 'DONE'
    ! produce and print detailed information
    CALL sp_properties
    IF(tmpi) THEN
      DO nst=1,nstmax
        IF(node(nst)/=mpi_myproc) sp_energy(nst)=0.0d0
      END DO
    CALL collect_sp_properties!collect single particle properties
    END IF
    CALL sinfo(wflag)
    !set x0dmp to 3* its value to get faster convergence
    IF(tvaryx_0) x0dmp=3.0d0*x0dmp
    !****************************************************  
    ! step 4: start static iteration loop
    !****************************************************  
    Iteration: DO iter=firstiter,maxiter
       IF(tmpi) CALL mpi_start_timer(1)
!      compute expectation value of constraint before step
       IF(tconstraint) CALL before_constraint(rho) 
       IF(wflag)WRITE(*,'(a,i6)') ' Static Iteration No.',iter
       !****************************************************  
       ! Step 5: gradient step
       !****************************************************  
       delesum=0.0D0  
       sumflu=0.0D0
       IF(ttime.AND.tmpi) CALL mpi_start_timer_iq(2)
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
       !$OMP SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
       DO nst=1,nstloc
         CALL grstep(globalindex(nst),isospin(globalindex(nst)),sp_energy(globalindex(nst)),denerg,psi(:,:,:,:,nst))
         sumflu=sumflu+wocc(globalindex(nst))*sp_efluct1(globalindex(nst))  
         delesum=delesum+wocc(globalindex(nst))*denerg  
       ENDDO
       !$OMP END PARALLEL DO
       IF(ttime.AND.tmpi) CALL mpi_stop_timer_iq(2,'grstep: ')
       IF(tmpi) CALL collect_energies(delesum,sumflu)!collect fluctuation and change in energy
       !****************************************************
       ! Step 6: diagonalize and orthonormalize
       !****************************************************
       sp_norm=0.0d0
       IF (my_iso>0) THEN
         CALL diagstep(my_iso,tdiag)
       ELSE
         DO iq=1,2
           CALL diagstep(iq,tdiag)    
         END DO 
       END IF
       !****************************************************
       ! Step 7: do pairing
       !****************************************************
       IF(kbT>0.0d0) CALL temp_dist
       IF(ipair/=0) CALL pair
       IF(kbT>0.0d0.AND.ipair/=0) STOP 'Pairing and finite temperature is not implemented'
       !****************************************************
       ! Step 8: get new densities and fields with relaxation
       !****************************************************
       IF(taddnew) THEN
          upot=rho
          bmass=tau
       ENDIF
       rho=0.0D0
       tau=0.0D0
       current=0.0D0
       sdens=0.0D0
       sodens=0.0D0
       IF(ttime.AND.tmpi) CALL mpi_start_timer_iq(2)
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst) SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:rho, tau, current, sdens, sodens)
       DO nst=1,nstloc
         CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)),&
                          psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
       ENDDO
       !$OMP END PARALLEL DO
       IF(tmpi) CALL collect_densities!collect densities from all nodes
       IF(ttime.AND.tmpi) CALL mpi_stop_timer_iq(2,'add density: ')
       !****************************************************
       ! Step 8a: optional constraint step
       !****************************************************      
       IF(tconstraint) THEN
         IF(ttime.AND.tmpi) CALL mpi_start_timer(2)
         CALL tune_constraint(e0dmp,x0dmp) 
         IF(tmpi)THEN
           CALL diagstep(my_iso,.FALSE.)
         ELSE
           DO iq=1,2
             CALL diagstep(iq,.FALSE.)
           END DO
         END IF 
         rho=0.0D0
         tau=0.0D0
         current=0.0D0
         sdens=0.0D0
         sodens=0.0D0
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst) SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:rho, tau, current, sdens, sodens)
       DO nst=1,nstloc
         CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)),&
                          psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
       ENDDO
       !$OMP END PARALLEL DO
       IF(tmpi) CALL collect_densities!collect densities from all nodes
       IF(ttime.AND.tmpi) CALL mpi_stop_timer(2,'constraint: ')
       END IF
       IF(taddnew) THEN
          rho=addnew*rho+addco*upot
          tau=addnew*tau+addco*bmass
       ENDIF  
       !****************************************************
       ! Step 8b: construct potentials
       !****************************************************  
       IF(ttime.AND.tmpi) CALL mpi_start_timer(2)
       CALL skyrme(iter<=outerpot,outertype)
       IF(tconstraint) CALL add_constraint(upot) 
       IF(ttime.AND.tmpi) CALL mpi_stop_timer_iq(2,'skyrme: ')
       !****************************************************
       ! calculate and print information
       !**************************************************** 
       IF(ttime.AND.tmpi)CALL mpi_start_timer_iq(2)
       CALL sp_properties
       IF(tmpi) THEN
         DO nst=1,nstmax
           IF(node(nst)/=mpi_myproc) sp_energy(nst)=0.0d0
         END DO
       IF(ttime.AND.tmpi) CALL mpi_stop_timer_iq(2,'sp properties: ')
       CALL collect_sp_properties!collect single particle properties
       END IF
       CALL sinfo(mprint>0.AND.MOD(iter,mprint)==0.AND.wflag)
       !****************************************************
       ! Step 9: check for convergence, saving wave functions
       !****************************************************
       IF(sumflu/nstmax<serr.AND.iter>1.AND..NOT.ttabc) THEN
          CALL write_wavefunctions
          EXIT Iteration  
       END IF
       IF(MOD(iter,mrest)==0) THEN  
          CALL write_wavefunctions
       ENDIF
       !***********************************************************************
       ! Step 10: calculate new step size
       !***********************************************************************
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
    !***********************************************************************
    ! Step 1:(h-esf) on psin yields ps1.
    !***********************************************************************
    esf=spe
    CALL hpsi(iq,esf,psin,ps1)
    !***********************************************************************
    ! Step 2: store ps1 in hampsi
    !***********************************************************************
    xnorm=rpsnorm(psin)
    xnormb=REAL(overlap(psin,ps1))
    hampsi(:,:,:,:,localindex(nst))=ps1!store h|psi> in hampsi
    !***********************************************************************
    ! Step 3: calculate fluctuation, i.e. <h*h> and |h|**2
    !***********************************************************************
     CALL hpsi(iq,esf,ps1,ps2)
     exph2=REAL(overlap(psin,ps2))
     varh2=rpsnorm(ps1)
     sp_efluct1(nst)=SQRT(ABS(exph2/xnorm-(xnormb/xnorm)**2))  
     sp_efluct2(nst)=SQRT(ABS(varh2/xnorm-(xnormb/xnorm)**2))
    !***********************************************************************
    ! Step 4: the damping step
    !***********************************************************************
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
    !***********************************************************************
    ! Step 5: energy convergence criterion
    !***********************************************************************
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
!> @param[in] diagonalize
!> INTEGER, If <tt>.TRUE.</tt> diagonalization is performed.
!--------------------------------------------------------------------------- 
  SUBROUTINE diagstep(iq,diagonalize)
    USE Trivial, ONLY: overlap,rpsnorm
    USE Linalg,  ONLY: eigenvecs,loewdin,comb_orthodiag,recombine,calc_matrix,wf_1dto2d,wf_2dto1d,&
                       matrix_split,matrix_gather,init_linalg,nlin
    
    INTEGER,INTENT(IN)             :: iq
    LOGICAL,INTENT(IN)             :: diagonalize
    INTEGER                        :: nst,nst2,noffset,ix,iy
    COMPLEX(db),ALLOCATABLE,TARGET :: unitary(:,:),hmatr_lin(:,:),unitary_h(:,:), rhomatr_lin(:,:),&
                                      rhomatr_lin_eigen(:,:), unitary_rho(:,:)
    COMPLEX(db),POINTER            :: unitary_d(:,:),hmatr_lin_d(:,:),unitary_h_d(:,:),rhomatr_lin_d(:,:),&
                                      unitary_rho_d(:,:),psi_2d(:,:),hampsi_2d(:,:)
    INTEGER,EXTERNAL               :: omp_get_num_threads,omp_get_thread_num
    !***********************************************************************
    ! Step 1: Copy |psi> and h|psi> to 2d storage mode
    !***********************************************************************
    IF(.NOT.tmpi) CALL init_linalg(iq)
    noffset=npmin(iq)-1
    ALLOCATE(unitary(nstloc_x,nstloc_y),           hmatr_lin(nstloc_x,nstloc_y),&
             unitary_h(nstloc_x,nstloc_y),         rhomatr_lin(nstloc_x,nstloc_y),&
             rhomatr_lin_eigen(nstloc_x,nstloc_y), unitary_rho(nstloc_x,nstloc_y))
    unitary_h=0.0d0
    hmatr_lin=0.0d0
    rhomatr_lin=0.0d0
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(tmpi) THEN
      ALLOCATE(psi_2d(psiloc_x,psiloc_y),hampsi_2d(psiloc_x,psiloc_y))
      CALL wf_1dto2d(psi(:,:,:,:,:),psi_2d)
      IF(diagonalize) CALL wf_1dto2d(hampsi(:,:,:,:,:),hampsi_2d)
    ELSE
      psi_2d(1:nx*ny*nz*2,1:nlin) => psi(:,:,:,:,npmin(iq):npsi(iq))
      hampsi_2d(1:nx*ny*nz*2,1:nlin) => hampsi(:,:,:,:,npmin(iq):npsi(iq))
    END IF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'1d to 2d: ')
    !***********************************************************************
    ! Step 2: Calculate lower tringular of h-matrix and overlaps.
    !***********************************************************************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    CALL calc_matrix(psi_2d,psi_2d,rhomatr_lin)
    IF(diagonalize)&
    CALL calc_matrix(psi_2d,hampsi_2d,hmatr_lin)
    unitary_h=0.0d0                                                                                        
    DO nst=1,nstloc_x
      DO nst2=1,nstloc_y
        ix=globalindex_x(nst)
        iy=globalindex_y(nst2)
        IF(ix==iy) THEN
          sp_norm(ix)=REAL(rhomatr_lin(nst,nst2))
          IF(diagonalize) THEN
            hmatr_lin(nst,nst2)=sp_energy(ix)!account for hampsi=(h-spe)|psi>
          ELSE
            unitary_h(nst,nst2)=CMPLX(1.0d0,0.0d0)
          END IF
        END IF
      ENDDO    !for nst2
    ENDDO    !for nst
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'CalcMatrix: ')
   !**************switch to diag contexts**************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(tmpi) THEN
      ALLOCATE(unitary_d(nstloc_diag_x,nstloc_diag_y),     hmatr_lin_d(nstloc_diag_x,nstloc_diag_y),&
               unitary_h_d(nstloc_diag_x,nstloc_diag_y),   rhomatr_lin_d(nstloc_diag_x,nstloc_diag_y),&
               unitary_rho_d(nstloc_diag_x,nstloc_diag_y))
      CALL matrix_split(rhomatr_lin,hmatr_lin,rhomatr_lin_d,hmatr_lin_d)
    ELSE
      unitary_d(1:nlin,1:nlin)      => unitary(:,:)
      unitary_rho_d(1:nlin,1:nlin)  => unitary_rho(:,:)
      unitary_h_d(1:nlin,1:nlin)    => unitary_h(:,:)
      hmatr_lin_d(1:nlin,1:nlin)    => hmatr_lin(:,:)
      rhomatr_lin_d(1:nlin,1:nlin)  => rhomatr_lin(:,:)
    END IF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'Diag Ortho Comm: ')
    !***********************************************************************
    ! Step 3: Calculate eigenvectors of h if wanted
    !***********************************************************************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(diagonalize.AND. (my_diag==1.OR.my_diag==3.OR.my_diag==0)) THEN
      CALL eigenvecs(hmatr_lin_d,unitary_h_d)
    END IF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'Diag matrix: ')
    !***********************************************************************
    ! Step 4: Calculate matrix for Loewdin
    !***********************************************************************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(my_diag==2.OR.my_diag==4.OR.my_diag==0)THEN
      CALL loewdin(rhomatr_lin_d,unitary_rho_d)
    ENDIF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'Ortho matrix: ')
    !*********************************switch back to 2d contxts
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(tmpi) CALL matrix_gather(unitary_rho,unitary_h,unitary_rho_d,unitary_h_d)
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'Diag ortho rev comm: ')
    !***********************************************************************
    ! Step 5: Combine h and diagonalization matrix and transpose them
    !***********************************************************************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(diagonalize) THEN
      CALL comb_orthodiag(unitary_h,unitary_rho,unitary)
    ELSE
      unitary=unitary_rho
    END IF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'Combine: ')
    !***********************************************************************
    ! Step 6: Recombine |psi> and write them into 1d storage mode
    !***********************************************************************
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    CALL recombine(psi_2d,unitary,hampsi_2d)
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'recombine: ')
    IF(tmpi.AND.ttime) CALL mpi_start_timer_iq(2)
    IF(tmpi)THEN
      CALL wf_2dto1d(hampsi_2d,psi(:,:,:,:,:))
    ELSE
      psi_2d=hampsi_2d
    END IF
    IF(tmpi.AND.ttime) CALL mpi_stop_timer_iq(2,'2d to 1d: ')
    DEALLOCATE(unitary,hmatr_lin,unitary_h,rhomatr_lin,&
               rhomatr_lin_eigen,unitary_rho)
    IF(tmpi) THEN
      DEALLOCATE(unitary_d,hmatr_lin_d,unitary_h_d,rhomatr_lin_d,&
                 unitary_rho_d,hampsi_2d,psi_2d)
    ELSE
      NULLIFY(unitary_d,hmatr_lin_d,unitary_h_d,rhomatr_lin_d,&
              unitary_rho_d,hampsi_2d,psi_2d)
    END IF
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
  SUBROUTINE sinfo(printing)
    INTEGER :: il,iq
    LOGICAL :: printing
    REAL(db):: tabc_energy, tabc_ekin, tabc_ecoul, tabc_eskyrme,entro
    CHARACTER(*),PARAMETER :: &
         header='  #  Par   v**2   var_h1   var_h2    Norm     Ekin    Energy &
         &    Lx      Ly      Lz     Sx     Sy     Sz    pairwg'   
    ! calculate static observables for printout                       *
    CALL moments
    CALL integ_energy
    CALL sum_energy
    IF(f%zpe==1 .AND. printing) THEN
       CALL cm_correction()
    END IF
    ! add information to summary files
    IF(printing) THEN
       IF(tabc_nprocs>1) THEN
          tabc_energy=tabc_av(ehf)
          tabc_ekin=tabc_av(tke)
          tabc_ecoul=tabc_av(ehfc)
          tabc_eskyrme=tabc_energy-tabc_ekin-tabc_ecoul
          IF(tabc_myid==0) THEN
             OPEN(unit=scratch,file=tabcfile,POSITION='APPEND')  
               WRITE(scratch,'(1x,i5,4F15.7)') &
                    iter, tabc_energy, tabc_ekin, tabc_ecoul, tabc_eskyrme
             CLOSE(unit=scratch)
          END IF
       END IF
       OPEN(unit=scratch,file=energiesfile,POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,1x,2(F9.3,1x),14(F14.5,1x))') &
            iter,pnr,ehf,ehfint,tke,ehfc,ehfCrho0,ehfCrho1,ehfCdrho0,ehfCdrho1,ehfCtau0,&
            ehfCtau1,ehfCdJ0,ehfCdJ1,ecmcorr,-entropy()*kbT
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=converfile,POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,f9.2,3(1pg11.3),2(0pf8.3),f6.1,f10.7)') &
            iter,ehf-ecmcorr,delesum/pnrtot,efluct1,efluct2,rmstot,beta,gamma,x0dmp
       CLOSE(scratch)
       OPEN(unit=scratch,file=dipolesfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,6E14.4)') iter,cmtot,cm(:,2)-cm(:,1)
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=spinfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,9F10.4)') iter,orbital,spin,total_angmom 
       CLOSE(unit=scratch)
       entro=entropy()
       WRITE(*,'(/,A,I7,A/A,2(F12.4,A)/2(A,F12.4),A/(3(A,E12.5),A))') &
            ' ***** Iteration ',iter,' *************************************************&
            &***********************************',&
            ' Free energy: ',ehf-ecmcorr-entro*kbT,' MeV Entropy: ',entro,'.',&
            ' Total energy: ',ehf-ecmcorr,&
            ' MeV  Total kinetic energy: ', tke,' MeV',&
            ' de/e:      ',delesum,'      h**2  fluct.:    ',efluct1,&
            ' MeV, h*hfluc.:    ',efluct2,' MeV', &
            ' MeV. Rearrangement E: ',e3corr,' MeV. Coul.Rearr.: ', &
             ecorc,' MeV   c.m.correction:',ecmcorr,' MeV'
       ! detail printout
       WRITE(*,'(/A)') ' Energies integrated from density functional:********************&
                  &********************************************'
       WRITE(*,'(4(A,1PE14.6),A/26X,3(A,1PE14.6),A)') &
            ' Total:',ehfint-ecmcorr,' MeV. t0 part:',ehf0,' MeV. t1 part:',ehf1, &
            ' MeV. t2 part:',ehf2,' MeV.',' t3 part:',ehf3,' MeV. t4 part:',ehfls, &
            ' MeV. Coulomb:',ehfc,' MeV.'
       WRITE(*,*)'                          *********************************************&
                  &**************************************'
       WRITE(*,'(26X,3(A,1PE14.6),A)')' Crho0:  ',ehfCrho0,' MeV. Crho1:  ',ehfCrho1,' MeV.'
       WRITE(*,'(26X,3(A,1PE14.6),A)')' Cdrho0: ',ehfCdrho0,' MeV. Cdrho1: ',ehfCdrho1,' MeV.'
       WRITE(*,'(26X,3(A,1PE14.6),A)')' Ctau0:  ',ehfCtau0,' MeV. Ctau1:  ',ehfCtau1,' MeV.'
       WRITE(*,'(26X,3(A,1PE14.6),A)')' CdJ0:   ',ehfCdJ0,' MeV. CdJ1:   ',ehfCdJ1,' MeV.'
       WRITE(*,*)'**********************************************************************&
                  &**************************************'
       IF(ipair/=0) THEN
         WRITE(*,'(a)') '          e_ferm      e_pair     <uv delta>   <v2 delta>   aver_force '
         DO il=1,2  
            WRITE(*,'(a,i2,a,5(1pg12.4))') 'iq=',il,': ',eferm(il) , &
               epair(il) ,avdelt(il),avdeltv2(il),avg(il)
         ENDDO
         IF(cutoff_factor>0D0) THEN
           WRITE(*,'(/7x,a)') '  e_ferm_cut    pairing-band   pairing space '
           DO iq=1,2  
             WRITE(*,'(a,i2,a,5(1pg12.4))') 'iq=',iq,': ',eferm_cutoff(iq) , &
               eferm_cutoff(iq)-eferm(iq),partnum_cutoff(iq)
           ENDDO
         END IF
         WRITE(*,*)'**********************************************************************&
                  &**************************************'
       END IF
       CALL flush(6)
       ! output densities
       IF(mplot/=0) THEN  
          IF(MOD(iter,mplot)==0) THEN
             !CALL plot_density
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
          IF(cutoff_factor>0D0) THEN
            WRITE(*,'(1X,I3,F4.0,F8.5,2F9.5,F9.6,F8.3,F10.3,3F8.3,5F7.3)') &
               il,sp_parity(il),wocc(il),sp_efluct1(il),sp_efluct2(il), &
               sp_norm(il),sp_kinetic(il),sp_energy(il), &
               sp_orbital(:,il),sp_spin(:,il),pairwg(il),deltaf(il)
          ELSE
            WRITE(*,'(1X,I3,F4.0,F8.5,2F9.5,F9.6,F8.3,F10.3,3F8.3,5F7.3)') &
               il,sp_parity(il),wocc(il),sp_efluct1(il),sp_efluct2(il), &
               sp_norm(il),sp_kinetic(il),sp_energy(il), &
               sp_orbital(:,il),sp_spin(:,il),pairwg(il),deltaf(il)
          END IF
       ENDDO
       CALL moment_print
       CALL radius_print
       CALL flush(6)
    END IF
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
    USE Trivial, ONLY: rpsnorm,overlap
    REAL(db) :: xx,yy,zz,xx2,zz2,y2,anorm,temp
    INTEGER  :: nst,iq,is,ix,iy,iz,nps,i,j,k,ka,nshell(3,nstmax)
    COMPLEX(db) :: psitemp(nx,ny,nz,2)
    IF(wflag)WRITE(*,*) "Harmonic oscillators widths (x-y-z):"
    IF(wflag)WRITE(*,"(3F12.4)") radinx,radiny,radinz
    psitemp=0.0d0
    psi=0.0d0
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
        DO k=0,ka  
          DO j=0,ka  
            DO i=0,ka  
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
            IF(node(nst)==mpi_myproc) psi(ix,iy,iz,1,localindex(nst))=EXP(-(temp))  
            psitemp(ix,iy,iz,1)=EXP(-(temp))  
          ENDDO
        ENDDO
      ENDDO
      IF(node(nst)==mpi_myproc) THEN
        anorm=rpsnorm(psi(:,:,:,:,localindex(nst)))
        psi(:,:,:,:,localindex(nst))=psi(:,:,:,:,localindex(nst))/SQRT(anorm)
      END IF
        anorm=rpsnorm(psitemp(:,:,:,:))
        psitemp(:,:,:,:)=psitemp(:,:,:,:)/SQRT(anorm)
      !*************************************************************************
      ! Higher states: lowest * polynomial
      !*************************************************************************
      DO nst=npmin(iq)+1,npsi(iq)
       IF (node(nst)==mpi_myproc) THEN
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
                psi(ix,iy,iz,is,localindex(nst))=psitemp(ix,iy,iz,1)*xx*yy*zz
            ENDDO
          ENDDO
        ENDDO
       psi(:,:,:,:,localindex(nst))=psi(:,:,:,:,localindex(nst))/sqrt(rpsnorm(psi(:,:,:,:,localindex(nst))))
       END IF
      END DO
    END DO
    IF(wflag)WRITE(*,*) '***** Harmonic oscillator initialization complete *****'
  END SUBROUTINE harmosc
!---------------------------------------------------------------------------  
! DESCRIPTION: planewaves
!> @brief
!!This subroutine initializes all wave function with plane wave states. 
!--------------------------------------------------------------------------- 
  SUBROUTINE planewaves
  IMPLICIT NONE
  INTEGER,PARAMETER :: npwmax=10
  INTEGER  :: nst, iq
  INTEGER :: i,j,l,ii,jj,kk
  INTEGER :: kf(3,npsi(2)),temp_k(3)
  LOGICAL :: check
  INTEGER :: ki(3,8*npwmax**3),ki_t(3,npwmax**3)
  REAL(db) :: temp_e,temp_energies(8*npwmax**3)
  WRITE(*,*)
  WRITE(*,*)'*****init plane waves:*****'
  psi=(0.d0,0.d0)
  wocc=0.D0
  wocc(1:nneut)=1.D0
  wocc(npmin(2):npmin(2)+nprot-1)=1.D0
  !***********************************************************************
  !                           calculate all k                            *
  !***********************************************************************
  j=0
  ii=0
  jj=0
  kk=0
  DO i=1,npwmax**3
    IF (ii==7) THEN
      ii=0
      jj=jj+1
    END IF
    IF (jj==7) THEN
      jj=0
      kk=kk+1
    END IF
    ki(1,i)=ii
    ki(2,i)=jj
    ki(3,i)=kk
    ii=ii+1
  END DO
  ki_t(:,1:npwmax**3)=ki(:,1:npwmax**3)
  l=1
  DO i=1,npwmax**3
    DO j=1,8
      ki(:,l)=ki_t(:,i)
      IF(j==2.OR.j==4.OR.j==6.OR.j==8) THEN 
        ki(1,l)=-ki_t(1,i)
      END IF
      IF(j==3.OR.j==4.OR.j==7.OR.j==8) THEN
        ki(2,l)=-ki_t(2,i)
      END IF
      IF(j==5.OR.j==6.OR.j==7.OR.j==8) THEN
        ki(3,l)=-ki_t(3,i)
      END IF
      temp_energies(l)=epw(ki(1,l),ki(2,l),ki(3,l))
      l=l+1
    END DO
  END DO
  !insertion_sort
  DO i=2,8*npwmax**3
    temp_e=temp_energies(i)
    temp_k(:)=ki(:,i)
    j=i
    DO WHILE (j>1 .AND. temp_energies(j-1)>temp_e)
      temp_energies(j)=temp_energies(j-1)
      ki(:,j)=ki(:,j-1)
      j=j-1
    END DO
    temp_energies(j)=temp_e
    ki(:,j)=temp_k(:)
  END DO
  nst = 1
  DO iq = 1,2  
    i=1 !counts nobs/2 (spin)
    j=1 !counts "-"-signs
    l=1 !counts ki
    DO WHILE(nst<=npsi(iq))
      IF(l>=8*npwmax**3)STOP 'PROBLEM'
      kf(:,i)=ki(:,l)
      CALL check_kf(kf,i,check)
      IF(check) THEN 
        CALL pw(nst,kf(1,i),kf(2,i),kf(3,i),1)
        nst=nst+1
        CALL pw(nst,kf(1,i),kf(2,i),kf(3,i),-1)
        nst=nst+1
        i=i+1
      END IF
      l=l+1
    END DO
  END DO
  END SUBROUTINE planewaves
!---------------------------------------------------------------------------  
! DESCRIPTION: epw
!> @brief
!!This function calculates the energy of a plane wave state with wave vectors kx,ky,kz. 
!> @param[in] kx
!> REAL(db), takes wave number in x-direction.
!> @param[in] ky
!> REAL(db), takes wave number in y-direction.
!> @param[in] kz
!> REAL(db), takes wave number in z-direction.
!--------------------------------------------------------------------------- 
REAL(db) FUNCTION epw(kx,ky,kz) RESULT(e)
  USE FORCES, ONLY: nucleon_mass
  INTEGER,INTENT(IN) :: kx,ky,kz
  REAL(db) :: dx,dy,dz
  dx=x(2)-x(1)
  dy=y(2)-y(1)
  dz=z(2)-z(1)
  e=(hbc**2)/(2*nucleon_mass)*(((2*pi*kx+bangx)/(nx*dx))**2&
  +((2*pi*ky+bangy)/(ny*dy))**2+((2*pi*kz+bangz)/(nz*dz))**2)
END FUNCTION
!---------------------------------------------------------------------------  
! DESCRIPTION: pw
!> @brief
!!This subroutine initializes state \c n with a plane wave with wave numbers 
!!\c kx, \c ky ,\c kz and spin \c s. 
!> @param[in] n
!> REAL(db), takes the state index n.
!> @param[in] kx
!> REAL(db), takes wave number in x-direction.
!> @param[in] ky
!> REAL(db), takes wave number in y-direction.
!> @param[in] kz
!> REAL(db), takes wave number in z-direction.
!> @param[in] s
!> REAL(db), takes spin (s>0 for up, s<0 for down).
!--------------------------------------------------------------------------- 
SUBROUTINE pw(n,kx,ky,kz,s)
  USE Trivial, ONLY : rpsnorm
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: n,kx,ky,kz,s
  INTEGER :: ix,iy,iz,nst
  REAL(db) :: facx,facy,facz,norm
  COMPLEX(db) :: fy,fz
  IF(mpi_myproc/=node(n)) RETURN
  nst=localindex(n)
  DO iz = 1,nz  
     facz=REAL(iz-1)*((2.D0*pi*REAL(kz)+bangz)/FLOAT(nz))
     fz=CMPLX(COS(facz),SIN(facz),db)
     DO iy=1,ny
        facy=REAL(iy-1)*((2.D0*pi*REAL(ky)+bangy)/FLOAT(ny))
        fy=CMPLX(COS(facy),SIN(facy),db)
        DO ix=1,nx
           facx=REAL(ix-1)*((2.D0*pi*REAL(kx)+bangx)/FLOAT(nx))
           IF(s>0) THEN
              psi(ix,iy,iz,1,nst)=fz*fy*CMPLX(COS(facx),SIN(facx),db)
              psi(ix,iy,iz,2,nst)=0.D0
           ELSE
              psi(ix,iy,iz,2,nst)=fz*fy*CMPLX(COS(facx),SIN(facx),db)
              psi(ix,iy,iz,1,nst)=0.D0
           END IF
        ENDDO
     ENDDO
  ENDDO
  norm=SQRT(rpsnorm(psi(:,:,:,:,nst)))
  psi(:,:,:,:,nst)=psi(:,:,:,:,nst)/norm
  WRITE(*,'(A14,3I2,A7,I2,A12,I4,A8,F9.5,A10,F6.3)')'state with k=(',kx,ky,kz,&
  '), spin= ',s,' at position',globalindex(nst),' energy ',epw(kx,ky,kz),' and wocc=',wocc(nst)
END SUBROUTINE pw
!---------------------------------------------------------------------------  
! DESCRIPTION: check_kf
!> @brief
!!This subroutine checks if the plane wave state \c i is already occupied
!> @param[in] k
!> REAL(db), takes wave vectors k.
!> @param[in] i
!> INTEGER, state which is searched for.
!> @param[out] check
!> LOGICAL, if true, state i has not been occupied, yet.
!--------------------------------------------------------------------------- 
SUBROUTINE check_kf(k,i,check)
  INTEGER,INTENT(IN) :: k(3,npsi(2)),i
  LOGICAL,INTENT(OUT) :: check
  INTEGER :: j
  check=.TRUE.
  IF(i==1) RETURN
  DO j=1,i-1
    IF(k(1,j)==k(1,i).AND.k(2,j)==k(2,i).AND.k(3,j)==k(3,i)) check=.FALSE.
  END DO
END SUBROUTINE
END MODULE Static
