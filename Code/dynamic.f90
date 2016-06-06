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
  USE abso_bc
  IMPLICIT NONE
  SAVE
  INTEGER :: nt        ! number of time steps
  REAL(db) :: dt
  INTEGER :: mxpact=6 ! number of terms in expon. expansion
  INTEGER :: mrescm=0 ! frequency of c.m. motion correction
  INTEGER :: nabsorb=0  !  number of absorbing points aling each direction
  REAL(db) :: rsep     ! distance where the calculation is stopped 
  LOGICAL :: texternal=.FALSE.  ! must be true if an external field is present
  LOGICAL :: text_timedep ! true for time-dependent external field
  REAL(db),PARAMETER :: esf=0.0D0  
  REAL(db) :: pnrold(2)=(/0D0,0D0/),ehfold=0D0,eintold=0D0,ekinold=0D0,&
              ecollold(2)=0D0
CONTAINS
  !*************************************************************************
  SUBROUTINE getin_dynamic
    NAMELIST /dynamic/ nt,dt,mxpact,mrescm,rsep,texternal,nabsorb
    READ(5,dynamic)  
    IF(wflag) THEN
       WRITE(*,*) '***** Parameters for the dynamic calculation *****'
       WRITE(*,"(A,I8,A,F10.5,A)") " Number of time steps:",nt, &
            ' Time step size: ',dt,' fm/c'
       WRITE(*,'(A,F7.2,A)') ' The calculation stops at ',rsep, &
            ' fm fragment separation'
       WRITE(*,'(A,I3)') ' Power limit in operator expansion:',mxpact
       WRITE(*,'(A,L7)') ' External field is invoked:',texternal
       WRITE(*,'(A,I3)') ' Number of absorbing points:',nabsorb
    ENDIF
    IF(texternal) CALL getin_external
  END SUBROUTINE getin_dynamic
  !*************************************************************************
  SUBROUTINE dynamichf
    INTEGER                   :: nst,istart,number_threads
    COMPLEX(db),ALLOCATABLE   :: ps4(:,:,:,:)
    INTEGER,  EXTERNAL        :: omp_get_num_threads 
    ALLOCATE(ps4(nx,ny,nz,2))
    ! Step 1: Preparation phase
    number_threads=1
    !$OMP PARALLEL
    !$ number_threads=OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    WRITE(*,*)'number of threads= ',number_threads
    IF(.NOT.trestart) THEN
       iter=0
       time=0D0  
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
            &      Ecoll(n)     Ecoll(p)  E_ext')
       CALL start_protocol(diffenergiesfile, &
            '#    Time N_n(0)-N_n(t) N_p(0)-N_p(t) diff-E(sum)  diff-E(integ)    diff-Ekin &
            &  diff-Ecoll(n)  diff-Ecoll(p)')
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
    CALL skyrme(.FALSE.,'N')
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
       CALL skyrme(.FALSE.,'N')
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
       IF(nabsorb > 0) CALL absbc(nabsorb,iter,nt,time)
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
       CALL skyrme(.FALSE.,'N')
       IF(text_timedep) CALL extfld(time+dt)
       IF(MOD(iter,mrest)==0) THEN  
          CALL write_wavefunctions
          IF(wflag) WRITE(*,*) ' Wrote restart file at end of  iter=',iter
       ENDIF
    END DO Timestepping
    DEALLOCATE(ps4)
  END SUBROUTINE dynamichf
  !***********************************************************************
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
  !***********************************************************************
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
       WRITE(scratch,'(1x,i5,6(1pg14.4))') iter,cmtot,cm(:,2)-cm(:,1)
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=momentafile, POSITION='APPEND')  
       WRITE(scratch,'(1x,f10.2,6(1pg14.6))') time,pcm(:,1)+pcm(:,2),&
            (pnr(2)*pcm(:,2)-pnr(1)*pcm(:,1))/(pnr(1)+pnr(2)) 
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
       IF(time==0D0) THEN
          pnrold=pnr
          ehfold=ehf
          eintold=ehfint
          ekinold=tke
          ecollold=ecoll
       END IF
       OPEN(unit=scratch,file=spinfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,9F10.4)') iter,orbital,spin,total_angmom 
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=energiesfile,POSITION='APPEND')  
       IF(ipulse==0) THEN  
         WRITE(scratch,'(F10.2,2F8.3,2F15.7,3(1PG13.5))') &
              time,pnr,ehf,ehfint,tke,ecoll
       ELSE
         WRITE(scratch,'(F10.2,2F8.3,2F15.7,4(1PG13.5))') &
              time,pnr,ehf,ehfint,tke,ecoll,e_extern
       END IF
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=diffenergiesfile,POSITION='APPEND')  
       WRITE(scratch,'(F10.2,9(1pg13.5))') &
            time,pnr-pnrold,ehf-ehfold,ehfint-eintold,&
            tke-ekinold,ecoll-ecollold
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
  !***********************************************************************
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
