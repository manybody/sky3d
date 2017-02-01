!------------------------------------------------------------------------------
! MODULE: Dynamic
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the routines needed for time propagation of the
!!system. All the logic needed for this
!!case is concentrated here; all other modules except for \c External
!!are equally used in the static calculation.
!------------------------------------------------------------------------------
MODULE DYNAMIC
  USE Params
  USE Grids, ONLY: nx,ny,nz,wxyz
  USE Densities
  USE Levels
  USE Energies
  USE Moment
  USE Twobody, ONLY: twobody_case,istwobody,roft,rdot, twobody_print
  USE Parallel
  USE Meanfield, ONLY: skyrme, hpsi, spot
  USE Trivial, ONLY: overlap
  USE Inout, ONLY: write_wavefunctions,write_densities, plot_density, &
       sp_properties,start_protocol
  USE External
  IMPLICIT NONE
  SAVE
  INTEGER            :: nt                !< the number of the final time step to be
  !! calculated. In case of a restart this is smaller than the total number of time steps.
  REAL(db)           :: dt                !< the physical time increment in units of fm/c.
  INTEGER            :: mxpact=6          !< the number of terms to be taken in the
  !! expansion of the potential
  INTEGER            :: mrescm=0          !< frequency of c.m. motion correction
  REAL(db)           :: rsep              !<  the final separation distance. The calculation
  !! is stopped if there has been a reseparation into two fragments and their distance exceeds \c rsep. 
  LOGICAL            :: texternal=.FALSE. !< this logical variable indicates that an
  !! external field is present. See module \c External.
  LOGICAL            :: text_timedep      !< this logical variable indicates that the
  !! external field is time-dependent and does not describe an instantaneous boost.
  REAL(db),PARAMETER :: esf=0.0D0         !< this is the energy shift for the call to \c hpsi. 
  !! Since it is not used in the dynamics part of the code, it
  !! is here set to the constant value of zero.
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: getin_dynamic
!> @brief
!!This is a relatively simple routine that reads the input for namelist
!!\c dynamic and prints it on standard output. If \c texternal is
!!true, it also calls \c getin_external to read the external field
!!parameters.
!--------------------------------------------------------------------------- 
  SUBROUTINE getin_dynamic
    NAMELIST /dynamic/ nt,dt,mxpact,mrescm,rsep,texternal
    READ(5,dynamic)  
    IF(wflag) THEN
       WRITE(*,*) '***** Parameters for the dynamic calculation *****'
       WRITE(*,"(A,I8,A,F10.5,A)") " Number of time steps:",nt, &
            ' Time step size: ',dt,' fm/c'
       WRITE(*,'(A,F7.2,A)') ' The calculation stops at ',rsep, &
            ' fm fragment separation'
       WRITE(*,'(A,I3)') ' Power limit in operator expansion:',mxpact
    ENDIF
    IF(texternal) CALL getin_external
  END SUBROUTINE getin_dynamic
!---------------------------------------------------------------------------  
! DESCRIPTION: dynamichf
!> @brief
!!This subroutine performs the main time-integration algorithm starting
!!at time step 0 and then iterating the desired number of steps.
!>
!> @details
!!Its building blocks are:
!!  - <b> Step 1: preparation phase:</b> this phase consists of several substeps. 
!!     -# If this is not a restart, the time and iteration number are
!!        zeroed (for a restart only the physical time is taken from the
!!        \c wffile). The wave functions are saved in the
!!        \c wffile, to save setup time in case the calculation has to
!!        be restarted from this initial point.
!!     -# The instantaneous external boost is applied using 
!!        \c extboost. If this subroutine has applied a boost, it sets 
!!        \c text_timedep to \c .FALSE. so no further calls to
!!        external-field routines are made (except for \c print_extfield).
!!     -# The protocol files <tt> *.res </tt> are initialized with
!!        their header lines.  
!!     -# The densities and current are calculated in a loop over wave
!!        function by calling \c add_densities. 
!!        They are first set to zero and then accumulated in a
!!        loop over the set on the local node, followed by collecting them
!!        over all nodes.
!!     -# The mean field and (if \c texternal is true) the
!!        external field are calculated for time zero using routines 
!!        \c skyrme and \c extfld.
!!     -# Then \c tinfo is called to calculate and print the
!!        single-particle quantities and the total energies at time 0 or
!!        iteration 0.
!!     -# Finally preparations are made for the time-stepping loop: the
!!        starting index is set to \c iter+1: this is either after the end of a
!!        previous job in the case of a restart, or just one for a new
!!        calculation. The physical time is either 0 or the time taken from
!!        a restart file.
!!
!!  - <b> Step 2: predictor time step:</b> the loop over the iteration
!!    index \c iter is started. Then the densities and mean-field
!!    components are estimated, i. e., in effect the Hamiltonian 
!!    \f$ \hat h(t+\tfrac1{2}\Delta t) \f$ required by the numerical method. This is done
!!    by evolving the wave functions for a full \c dt using the old
!!    Hamiltonian and averaging the densities between old and new ones to
!!    obtain the mid-time Hamiltonian for propagation. In detail the
!!    procedure is as follows:
!!     -# The densities are not set to zero in order to allow  adding the
!!        contributions of the wave functions at the end of the time step.
!!     -# In the MPI version, the densities are divided by the number of
!!        nodes. In this way, adding up contributions from all nodes, the
!!        densities from the beginning of the time step will be included
!!        correctly.
!!     -# Subroutine \c tstep is used to propagate the wave functions
!!        to the end of the time step. Note that truncation in the
!!        exponential happens at \c mxpact/2, since the accuracy need not
!!        be as high as in the full step. For each wave function its
!!        contribution is added to the densities. The wave functions
!!        themselves do not need to be saved as they are not used for
!!        anything else.
!!     -# After the loop, contributions from all the nodes are added up
!!        for the MPI case using subroutine \c collect_densities.
!!     -# The densities are multiplied by one half to form the average
!!        of the values at \f$ t \f$ and \f$ t+\Delta t \f$.
!!     -# These average densities are then used to calculate the mean
!!        field at half time using subroutine \c skyrme}; also the
!!        external field is obtained for the half time using \c extfld.
!!
!!  - <b> Step 3: full time step:</b> Now that the single-particle
!!    Hamiltonian has been estimated for the middle of the time step, the
!!    propagation can be carried out to the end of the time step. This is
!!    quite analogous to the half step with only three crucial
!!    differences:
!!     -# The densities are reset to zero before the wave function loop,
!!        so the densities summed up are the purely the densities at the end
!!        of the time step,
!!     -# the series expansion in \c tstep now uses the full 
!!        \c mxpact terms, and
!!     -# the new wave functions are copied back into \c psi to be
!!        available for the next time step.
!!
!!  - <b> Step 4: Center-of-mass correction: </b>
!!    If a center-of-mass correction is desired by the user by setting
!!    <tt> mrescm /=0</tt>, 
!!    subroutine \c resetcm is called every \c mrescm'th time step to
!!    reset the center-of-mass velocity to zero.  
!!
!!  - <b> Step 5: generating some output </b> At this point the time is
!!    advanced by \c dt because the physical time is now the end of
!!    the time step, and this must be printed out correctly by the
!!    following output routines. \c tinfo is called to calculate
!!    single-particle properties, total energies, and so on.
!!
!!  - <b> Step 6: finishing up the time step:</b> \c tinfo is called
!!  to output the calculated data, then \c skyrme and \c extfld
!!  calculate the mean field and the external field, respectively, for the
!!  end of the time step, after which the wave functions are written
!!  onto \c wffile depending on \c mrest.
!!
!!This ends the time loop and subroutine \c dynamichf itself.
!--------------------------------------------------------------------------- 
  SUBROUTINE dynamichf
    INTEGER :: nst,istart
    COMPLEX(db),ALLOCATABLE :: ps4(:,:,:,:)
    ALLOCATE(ps4(nx,ny,nz,2))
    ! Step 1: Preparation phase
    IF(.NOT.trestart) THEN
       iter=0
       time=0.0D0  
       ! save wave functions
       CALL write_wavefunctions
       IF(wflag) WRITE(*,*) 'Wrote wave function file after initialization'
    END IF
    ! external boost
    IF(texternal) THEN
       CALL init_external
       CALL extboost(text_timedep)
    END IF
    ! Create protocol files
    IF(wflag) THEN
       ! Initialize *.res files
       CALL start_protocol(energiesfile, &
            '#    Time    N(n)    N(p)       E(sum)        E(integ)      Ekin &
            &      Ecoll(n)     Ecoll(p)')
       CALL start_protocol(monopolesfile, &
            '#    Time      rms_n     rms_p   rms_tot   rms_n-rms_p')
       CALL start_protocol(dipolesfile, &
            '# Iter    c.m. x-y-z                                  Isovector&
            &dipoles x-y-z')
       CALL start_protocol(quadrupolesfile, &
            '#     Time     Q(n)          Q(p)         Q(n+p)        x²(n)         &
            &y²(n)         z²(n)         x²(p)         y²(p)         z²(p)')
       CALL start_protocol(spinfile, &
            '# Iter      Lx        Ly        Lz        Sx        Sy        &
            &Sz        Jx        Jy        Jz')
       CALL start_protocol(momentafile, &
            '#     Time      Px            Py            Pz')
       IF(texternal) CALL start_protocol(extfieldfile, &
            '#   time       average_extfield')

    END IF
    ! calculate densities and currents
    rho=0.0D0
    tau=0.0D0
    current=0.0D0
    sdens=0.0D0
    sodens=0.0D0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst) SCHEDULE(STATIC) &
    !$OMP REDUCTION(+:rho,tau,current,sdens,sodens)
    DO nst=1,nstloc
       CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)), &
            psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
    ENDDO
    !$OMP END PARALLEL DO
    IF(tmpi) CALL collect_densities
    ! calculate mean fields and external fields
    CALL skyrme
    IF(text_timedep) CALL extfld(0.D0)
    CALL tinfo
    !***********************************************************************
    istart=iter+1
    ! Step 2: start loop and do half-time step
    Timestepping:  DO iter=istart,nt  
       IF(wflag) WRITE(*,'(/A,I6,A,F8.2,A)') ' Starting time step #',iter, &
            ' at time=',time,' fm/c'
       ! correction for parallel version
       IF(tmpi) THEN
          rho=rho/mpi_nprocs
          tau=tau/mpi_nprocs
          current=current/mpi_nprocs
          sodens=sodens/mpi_nprocs
          sdens=sdens/mpi_nprocs
       ENDIF
       ! propagate to end of time step and add to densities
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,ps4) SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:rho,tau,current,sdens,sodens)
       DO nst=1,nstloc
          ps4=psi(:,:,:,:,nst) 
          CALL tstep(isospin(globalindex(nst)),mxpact/2,ps4)
          CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)), &
               ps4,rho,tau,current,sdens,sodens)  
       ENDDO
       !$OMP END PARALLEL DO
       IF(tmpi) CALL collect_densities
       ! average over time step
       rho=0.5D0*rho
       tau=0.5D0*tau
       current=0.5D0*current
       sodens=0.5D0*sodens
       sdens=0.5D0*sdens
       ! compute mean field and add external field
       CALL skyrme  
       IF(text_timedep) CALL extfld(time+dt/2.0D0)
       ! Step 3: full time step
       ! reset densities
       rho=0.0D0
       tau=0.0D0
       current=0.0D0
       sdens=0.0D0
       sodens=0.0D0
       ! propagate to end of step, accumulate densities
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,ps4) SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:rho,tau,current,sdens,sodens)
       DO nst=1,nstloc
          ps4=psi(:,:,:,:,nst) 
          CALL tstep(isospin(globalindex(nst)),mxpact,ps4)
          CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)), &
               ps4,rho,tau,current,sdens,sodens)  
          psi(:,:,:,:,nst)=ps4
       ENDDO
       !$OMP END PARALLEL DO
       ! sum up over nodes
       IF(tmpi) CALL collect_densities
       ! Step 4: eliminate center-of-mass motion if desired
       IF(mrescm/=0) THEN  
          IF(MOD(iter,mrescm)==0) THEN  
             CALL resetcm
             !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst) SCHEDULE(STATIC) &
             !$OMP REDUCTION(+:rho,tau,current,sdens,sodens)
             DO nst=1,nstloc
                CALL add_density(isospin(globalindex(nst)), &
                     wocc(globalindex(nst)), &
                     psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
             ENDDO
             !$OMP END PARALLEL DO
             IF(tmpi) CALL collect_densities
          ENDIF
       ENDIF
       ! Step 5: generating some output
       time=time+dt
       CALL tinfo
       ! Step 6: finishing up
       ! compute densities, currents, potentials etc.                  *
       CALL skyrme  
       IF(text_timedep) CALL extfld(time+dt)
       IF(MOD(iter,mrest)==0) THEN  
          CALL write_wavefunctions
          IF(wflag) WRITE(*,*) ' Wrote restart file at end of  iter=',iter
       ENDIF
    END DO Timestepping
    DEALLOCATE(ps4)
  END SUBROUTINE dynamichf
!---------------------------------------------------------------------------  
! DESCRIPTION: tstep
!> @brief
!!In this subroutine one wave function given as the argument \c psout
!!is stepped forward in time by the interval \c dt. 
!>
!> @details
!!The method used is the expansion of the exponential time-development operator
!!cut off at the power of \c mxp, which in practice is usually around 6. 
!!Suppressing the argument of \f$ \hat h \f$ for brevity, we can write
!!\f[ \hat U(t,t+\Delta t)\,\phi\approx \sum_{n=0}^m \phi^{(n)} \f]
!!with
!!\f[ \phi^{(0)}=\phi,\qquad \phi^{(k+1)}=\frac{-\I\,\Delta t}
!! {\hbar c k}\,\hat h\,\phi^{(k)},\quad k=0,\ldots,m-1. \f]
!!
!!Thus the application of the polynomial to a single-particle wave
!!function can be evaluated simply in a loop applying \f$ \hat h\phi_k \f$
!!repeatedly and accumulating the results in wave function \c psout.
!!
!!The argument \c iq is only necessary because \c hpsi needs
!!information about the isospin of the wave function.
!>
!> @param[in] iq
!> INTEGER, takes the isospin.
!> @param[in] mxp
!> REAL(db), takes the power to which the time-development operator is expanded.
!> @param[in,out] psout
!> REAL(db), array, takes the wave function and returns evolved wave function.
!--------------------------------------------------------------------------- 
  SUBROUTINE tstep(iq,mxp,psout)
    INTEGER,INTENT(IN) :: iq,mxp
    COMPLEX(db),INTENT(INOUT) :: psout(:,:,:,:)
    INTEGER :: m,is,ix,iy,iz
    COMPLEX(db) :: ps1(nx,ny,nz,2),ps2(nx,ny,nz,2)  
    REAL(db) :: fmd
    ps1=psout
    !***********************************************************************
    !        compute exp(-I*dt*h) by power series expansion              *
    !***********************************************************************
    DO m=1,mxp
       fmd=-dt/(hbc*m)  
       CALL hpsi(iq,esf,ps1,ps2)
       DO is=1,2  
          DO iz=1,nz  
             DO iy=1,ny  
                DO ix=1,nx  
                   ps1(ix,iy,iz,is)=& 
                        CMPLX(-fmd*AIMAG(ps2(ix,iy,iz,is)),   &
                        fmd*REAL(ps2(ix,iy,iz,is)),db)
                   psout(ix,iy,iz,is)=psout(ix,iy,iz,is)+ps1(ix,iy,iz,is)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    END DO
  END SUBROUTINE tstep
  !---------------------------------------------------------------------------  
! DESCRIPTION: tinfo
!> @brief
!!This subroutine is used to output various pieces of information
!!relevant especially to the dynamic mode of the code. 
!>
!> @details
!!It is called at
!!the end of the full time step and consists of the following steps:
!!  - <b> Step 1: initialization: </b> the flag \c printnow is
!!    calculated to keep track of whether this is the proper time step for
!!    a full printout (determined by \c mprint).
!!  - <b> Step 2: twobody analysis </b> the twobody analysis is
!!    performed, but only if the calculation started as a twobody
!!    scenario.
!!  - <b> Step 3: moments: </b> the moments of the distribution are
!!    calculated using subroutine \c moments. This includes total mass,
!!    momenta, and angular momenta. They are printed out if indicated by
!!    \c printnow, both on the large output and in the specialized
!!    files \c dipolesfile, \c momentafile, and
!!    \c spinfile. If there is an external field, the routine
!!    \c print_extfield is called to print the current expectation
!!    value of the external field.
!!
!!    <b>Note that the moments need to be calculated every time step,
!!    because some of the logic may depend on them, especially the
!!    twobody-analysis, which needs the correct c.m., for example and is
!!    calculated at every time step and on every node.</b>
!! 
!!  - <b> Step 4: single-particle quantities:</b> the single-particle
!!    energies are calculated straightforwardly as expectation values 
!!    of the Hamiltonian.
!!    The routine \c sp_properties is then called to obtain the other
!!    single-particle properties like angular momenta. They are
!!    communicated between the processors. The angular momenta are written
!!    to \c spinfile.
!!  - <b> Step 5: total energies:</b> The integrated energy \c ehfint
!!    and its contributions are calculated in \c integ_energy. The
!!    subroutine \c sum_energy is called to calculate the three-body
!!    energy and the single-particle based total energy \c ehf. The
!!    collective kinetic energy \c ecoll is computed directly
!!    here, because it is needed only in the dynamic calculations and only
!!    for output. It is defined as
!!    \f[ E_{\rm coll}=\frac{\hbar^2}{2m}\int \D^3r \frac{{\vec \jmath\,}^2}{\rho} \f]
!!    The energies are protocolled in \c energiesfile and on standard
!!    output.
!!  - <b> Step 6: density output:</b> at intervals of \c mprint or in
!!    the first time step the density printer plot is generated using 
!!    \c plot_densities and the binary densities are written onto <tt> *.tdd </tt>
!!    files using \c write_densities.
!!  - <b> Step 7: other output: </b> in the proper \c mprint interval
!!    the two-body analysis results, the single-particle state
!!    information, and the moments are printed on standard output, using
!!    also the routines \c twobody_print and \c moment_print. 
!!
!!    It is important to note that when \c tinfo is called before the 
!!    time step iteration starts, the
!!    two-body analysis cannot work because the fragment centers of mass
!!    from the previous time step  are either not known yet (restart) or
!!    identical to the present ones (initialization). The subroutine \c twobody_case 
!!    is still called to find the fragment properties,
!!    but the derived kinetic energy etc. will be incorrect. Therefore
!!    the call to \c twobody_print is suppressed in this case. The
!!    logical variable \c initialcall is used to recognize this case.
!!
!!  - <b> Step 8: check for final separation: </b> for the twobody case
!!    it is checked whether the separation found between the two fragments
!!    is larger than the input quantity \c rsep with positive time
!!    derivative \c rdot of the separation distance, in which case the
!!    program is terminated with an appropriate message.
!--------------------------------------------------------------------------- 
  SUBROUTINE tinfo
    REAL(db), DIMENSION(2) :: ecoll     ! storage for collective-flow energy
    INTEGER :: il
    CHARACTER(*),PARAMETER :: &
         header='  #    v**2    Norm     Ekin    Energy &
         &    Lx      Ly      Lz     Sx     Sy     Sz  '   
    !***********************************************************************
    !     calculates dynamic observables for printout                      *
    !***********************************************************************
    INTEGER :: iq,nst
    COMPLEX(db),ALLOCATABLE :: ps1(:,:,:,:)
    LOGICAL,SAVE :: initialcall=.TRUE.
    ALLOCATE(ps1(nx,ny,nz,2))
    ! Step 1
    printnow=mprint>0.AND.MOD(iter,mprint)==0
    ! Step 2: twobody analysis
    IF(nof/=2) THEN  
       istwobody=.FALSE.
    ELSE
       CALL twobody_case(dt)
    ENDIF
    ! Step 3: moments calculated
    CALL moments
    IF(printnow.AND.wflag) THEN
       OPEN(unit=scratch,file=dipolesfile,POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,6E14.4)') iter,cmtot,cm(:,2)-cm(:,1)
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=momentafile, POSITION='APPEND')  
       WRITE(scratch,'(1x,f10.2,3g14.6)') time,pcm(:,1)+pcm(:,2)
       CLOSE(unit=scratch)
       CALL moment_shortprint
       IF(texternal) CALL print_extfield()
    ENDIF
    ! Step 4: single-particle properties
    IF(printnow) THEN  
       sp_energy=0.0D0
       sp_norm=0.0D0
       DO nst=1,nstloc
          iq=isospin(globalindex(nst))
          CALL hpsi(iq,esf,psi(:,:,:,:,nst),ps1)
          sp_energy(globalindex(nst))=overlap(psi(:,:,:,:,nst),ps1)
          sp_norm(globalindex(nst))=&
               overlap(psi(:,:,:,:,nst),psi(:,:,:,:,nst))
       ENDDO
       CALL sp_properties
       IF(tmpi) THEN
          CALL collect_sp_properties
       ENDIF
    END IF
    ! Step 5: total energies & angular momenta
    IF(printnow.AND.wflag) THEN
       CALL integ_energy
       CALL sum_energy
       DO iq=1,2
          ecoll(iq)=f%h2m(iq)*SUM( &
               (current(:,:,:,1,iq)**2+current(:,:,:,2,iq)**2&
               +current(:,:,:,3,iq)**2)/rho(:,:,:,iq) )
       ENDDO
       OPEN(unit=scratch,file=spinfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,9F10.4)') iter,orbital,spin,total_angmom 
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=energiesfile,POSITION='APPEND')  
       WRITE(scratch,'(F10.2,2F8.3,2F15.7,3(1PG13.5))') &
            time,pnr,ehf,ehfint,tke,ecoll
       CLOSE(unit=scratch)
       WRITE(*,'(/A,I7,A,F12.4,A/2(A,F12.4),A)') &
            ' ***** Time step #',iter,' at time ',time,' fm/c *****', &
            ' Total energy: ',ehf,' MeV  Total kinetic energy: ', &
            tke,' MeV.'
       WRITE(*,'(/A)') ' Energies integrated from density functional:'
       WRITE(*,'(4(A,1PE14.6),A/26X,3(A,1PE14.6),A)') &
            ' Total:',ehfint,' MeV. t0 part:',ehf0,' MeV. t1 part:', &
            ehf1,' MeV. t2 part:',ehf2,' MeV',' t3 part:',ehf3, &
            ' MeV. t4 part:',ehfls,' MeV. Coulomb:',ehfc,' MeV.'
    ENDIF
    ! Step 6: density output
    IF(mplot/=0) THEN
       IF(wflag.AND.MOD(iter,mplot)==0) THEN
          CALL plot_density
          CALL write_densities
       ENDIF
    ENDIF
    ! Step 7: print other output
    IF(printnow.AND.wflag) THEN
       IF(nof==2.AND.istwobody.AND..NOT.initialcall) &
            CALL twobody_print
       WRITE(*,'(/A)') ' Neutron Single Particle States:',header
       DO il=1,nstmax
          IF(il==npmin(2)) THEN
             WRITE(*,'(/A)') ' Proton Single Particle States:',header  
          END IF
          WRITE(*,'(1X,I3,F8.5,F9.6,F8.3,F10.3,3F8.3,3F7.3)') &
               il,wocc(il),sp_norm(il),sp_kinetic(il), &
               sp_energy(il),sp_orbital(:,il),sp_spin(:,il)
       ENDDO
       CALL moment_print
    ENDIF
    DEALLOCATE(ps1)
    ! Step 8: check whether final distance is reached in twobody case
    IF(istwobody.AND.roft>rsep.AND.rdot>=0.D0) THEN  
       CALL twobody_print
       CALL write_wavefunctions
       IF(wflag) WRITE(*,*) ' Final separation distance reached'
       STOP ' Final distance reached'  
    ENDIF
    initialcall=.FALSE.
  END SUBROUTINE tinfo
!---------------------------------------------------------------------------  
! DESCRIPTION: resetcm
!> @brief
!!This subroutine resets the center-of-mass velocity to zero. 
!>
!> @details
!!The velocity, or rather the corresponding wave vector, is calculated from
!!the current density using
!!\f[
!!\rho(\vec r)  \vec k=\frac{1}{2\I}\sum_{\alpha\in q}w_\alpha^2
!!    \sum_s\left(\psi_\alpha^*(\vec r,s)\nabla\psi_\alpha(\vec
!!      r,s)-\psi_\alpha(\vec r,s) \nabla\psi_\alpha^*(\vec r,s)\right).
!!\f]
!!The wave functions are then multiplied by a common plane-wave phase
!!factor \f$ \exp(-\I\vec k\cdot \vec r) \f$ to give a counter boost.
!--------------------------------------------------------------------------- 
  SUBROUTINE resetcm
    INTEGER :: ix,iy,iz,is,nst
    REAL(db) :: akf(3),totmass
    ! compute c.m. momenta and prepare corrective phase factor
    totmass=wxyz*SUM(rho)
    DO is=1,3
       akf(is)=-wxyz*SUM(current(:,:,:,is,:))/totmass
    END DO
    ! now apply correction to each wavefunction
    DO nst=1,nstloc
       FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)             
          psi(ix,iy,iz,is,nst)=&
               psi(ix,iy,iz,is,nst) * EXP( &
               CMPLX(0.D0,akf(1)*x(ix)+akf(2)*y(iy)+akf(3)*z(iz),db))
       END FORALL
    ENDDO
  END SUBROUTINE resetcm
END MODULE DYNAMIC
