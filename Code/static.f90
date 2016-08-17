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
  USE Pairs, ONLY: pair,epair,avdelt,avg,eferm
  USE Parallel, ONLY:ttabc, tabc_av
  IMPLICIT NONE
  LOGICAL :: tdiag=.FALSE.
  LOGICAL :: tlarge=.FALSE.
  LOGICAL :: tvaryx_0=.FALSE.
  INTEGER :: maxiter,outerpot=0
  REAL(db) :: radinx,radiny,radinz, &
       serr,delesum,x0dmp=0.2D0,e0dmp=100.D0,x0dmpmin=0.2d0,ttime(10,20)=0.0d0
  CHARACTER(1) :: outertype='N'
CONTAINS
  !*************************************************************************
  SUBROUTINE getin_static
    NAMELIST/static/ tdiag,tlarge,maxiter, &
         radinx,radiny,radinz,serr,x0dmp,e0dmp,nneut,nprot,npsi,tvaryx_0,&
         outerpot,outertype
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
    END IF
  END SUBROUTINE getin_static
  !*************************************************************************
  SUBROUTINE init_static
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
            '# Iter    N(n)    N(p)       E(sum)         E(integ)       Ekin         &
            &E_Coul         ehfCrho0       ehfCrho1       ehfCdrho0      ehfCdrh     & 
            &ehfCtau0       ehfCtau1       ehfCdJ0        ehfCdJ1')
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
  END SUBROUTINE init_static
  !*************************************************************************
  SUBROUTINE statichf
    LOGICAL, PARAMETER   :: taddnew=.TRUE. ! mix old and new densities
    INTEGER              :: iq,nst,firstiter,number_threads
    REAL(db)             :: sumflu,denerg
    REAL(db) , PARAMETER :: addnew=0.2D0,addco=1.0D0-addnew  
    INTEGER,  EXTERNAL   :: omp_get_num_threads 
    !****************************************************   
    ! Step 1: initialization
    !****************************************************  
    number_threads=1
    !$OMP PARALLEL
    !$ number_threads=OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    IF(wflag)WRITE(*,*)'number of threads= ',number_threads
    IF(wflag)WRITE(*,*)
    ALLOCATE(hmatr(nstmax,nstmax))
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
       WRITE(*,'(A25)',advance="no") 'Initial schmid... '
       CALL schmid
       WRITE(*,*)'DONE'
    END IF
    !****************************************************  
    ! Step 2: calculate densities and mean field
    !****************************************************  
    rho=0.0D0
    tau=0.0D0
    current=0.0D0
    sdens=0.0D0
    sodens=0.0D0
    WRITE(*,'(A25)',advance="no")'Initial add_density... '
    DO nst=1,nstmax
       CALL add_density(isospin(nst),wocc(nst),psi(:,:,:,:,nst), &
            rho,tau,current,sdens,sodens)  
    ENDDO
    WRITE(*,*) 'DONE'
    WRITE(*,'(A25)',advance="no")'Initial skyrme... '
    CALL skyrme(iter<=outerpot,outertype)
    WRITE(*,*) 'DONE'
    !****************************************************  
    ! Step 3: initial gradient step
    !****************************************************  
    delesum=0.0D0  
    sumflu=0.0D0  
    WRITE(*,'(A25)',advance="no")'Initial grstep... '
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
    !$OMP SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
    DO nst=1,nstmax
       CALL grstep(nst,isospin(nst),sp_energy(nst),denerg, &
            psi(:,:,:,:,nst))
       sumflu=sumflu+wocc(nst)*sp_efluct1(nst)  
       delesum=delesum+wocc(nst)*denerg  
    ENDDO
    !$OMP END PARALLEL DO
    WRITE(*,*) 'DONE'
    ! pairing and orthogonalization
    IF(ipair/=0) CALL pair
    WRITE(*,'(A25)',advance="no") 'Initial schmid... '
    CALL schmid
    WRITE(*,*) 'DONE'
    ! produce and print detailed information
    CALL sp_properties
    CALL sinfo(.TRUE.)
    !set x0dmp to 3* its value to get faster convergence
    IF(tvaryx_0) x0dmp=3.0d0*x0dmp
    !****************************************************  
    ! step 4: start static iteration loop
    !****************************************************  
    Iteration: DO iter=firstiter,maxiter
       IF(wflag)WRITE(*,'(a,i6)') ' Static Iteration No.',iter
       !****************************************************  
       ! Step 5: gradient step
       !****************************************************  
       delesum=0.0D0  
       sumflu=0.0D0
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,denerg) &
       !$OMP SCHEDULE(STATIC) REDUCTION(+: sumflu , delesum)
       DO nst=1,nstmax
          CALL grstep(nst,isospin(nst),sp_energy(nst),denerg, &
               psi(:,:,:,:,nst))
          sumflu=sumflu+wocc(nst)*sp_efluct1(nst)  
          delesum=delesum+wocc(nst)*denerg  
       ENDDO
       !$OMP END PARALLEL DO
       !****************************************************
       ! Step 6: diagonalize and orthonormalize
       !****************************************************
       DO iq=1,2
          CALL diagstep(iq,npsi(iq)-npmin(iq)+1)
       ENDDO
       !****************************************************
       ! Step 7: do pairing
       !****************************************************
       IF(ipair/=0) CALL pair
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
       CALL skyrme(iter<=outerpot,outertype)
       ! calculate and print information
       CALL sp_properties
       CALL sinfo(mprint>0.AND.MOD(iter,mprint)==0)
       !****************************************************
       ! Step 9: check for convergence, saving wave functions, update stepsize
       !****************************************************
       IF(sumflu/nstmax<serr.AND.iter>1.AND..NOT.ttabc) THEN
          CALL write_wavefunctions
          EXIT Iteration  
       END IF
       IF(MOD(iter,mrest)==0) THEN  
          CALL write_wavefunctions
       ENDIF
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
    DEALLOCATE(hmatr)
  END SUBROUTINE statichf
  !*************************************************************************
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
    IF(mprint>0.AND.MOD(iter,mprint)==0) THEN  
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
  !*************************************************************************
  SUBROUTINE diagstep(iq,nlin)
    !***********************************************************************
    !                                                                      *
    !     diagstep= diagonalize Hamiltonian matrix of active shells      *
    !                                                                      *
    !***********************************************************************
    USE Trivial, ONLY: overlap,rpsnorm
    INTEGER,INTENT(IN) :: iq,nlin
    INTEGER :: nst,nst2,noffset,i,j,ix,iy,iz,is
    INTEGER :: infoconv
    REAL(db) :: eigen_h(nlin)
    COMPLEX(db) :: unitary_h(nlin,nlin),unitary_rho(nlin,nlin),unitary(nlin,nlin)
    COMPLEX(db) :: pstemp(nlin),pstemp1(nlin)
    COMPLEX(db) :: cmplxone=CMPLX(1.0,0.0),cmplxzero=CMPLX(0.0,0.0)
    COMPLEX(db) :: rhomatr_lin(nlin,nlin)
    COMPLEX(db) :: cwork(2*nlin*nlin)
    REAL(db)    :: rwork(2*nlin*nlin+5*nlin+1)
    INTEGER     :: iwork(5*nlin+3)
    EXTERNAL    :: zgemv,zheevd,zgemm,zheev
    ! Step 1: copy matrix into symmetric storage mode, then diagonalize
    noffset=npmin(iq)-1
    unitary_h=0.0d0
    rhomatr_lin=0.0d0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nst,nst2) &
    !$OMP SCHEDULE(STATIC)
    DO i=nlin/2+1,nlin
       DO j=1,2*(nlin/2)+1
          IF(j>i) THEN
             nst2=noffset+2*(nlin/2)-i+1
             nst=noffset+2*(nlin/2)+2-j
          ELSE
             nst2=noffset+i
             nst=noffset+j
          END IF
          unitary_h(nst2-noffset,nst-noffset)=&
               0.5D0*(CONJG(hmatr(nst,nst2))+hmatr(nst2,nst))
          rhomatr_lin(nst-noffset,nst2-noffset)=overlap(psi(:,:,:,:,nst),psi(:,:,:,:,nst2))
          IF(nst==nst2) THEN
             sp_norm(nst)=rhomatr_lin(nst-noffset,nst2-noffset)
          ELSE
             rhomatr_lin(nst2-noffset,nst-noffset)=CONJG(rhomatr_lin(nst-noffset,nst2-noffset))
          END IF
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
! Diagonalization
    CALL zheevd('V','L',nlin,unitary_h,nlin,eigen_h,cwork,nlin*nlin*2,rwork,&
                2*nlin*nlin+5*nlin+1,iwork,5*nlin+3,infoconv)
! Loewdin Orthonormalization
    CALL zheevd('V','L',nlin,rhomatr_lin,nlin,eigen_h,cwork,nlin*nlin*2,rwork,&
                2*nlin*nlin+5*nlin+1,iwork,5*nlin+3,infoconv)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO i=1,nlin
      rhomatr_lin(:,i)=rhomatr_lin(:,i)*sqrt(eigen_h(i))
    END DO
    CALL zgemm('N','C',nlin,nlin,nlin,cmplxone,rhomatr_lin,nlin,rhomatr_lin,nlin,&
               cmplxzero,unitary_rho,nlin)  
! Recombination
    CALL zgemm('N','N',nlin,nlin,nlin,cmplxone,&
               unitary_rho,nlin,unitary_h,nlin,cmplxzero,unitary,nlin)
    noffset=npmin(iq)-1
    !$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2) PRIVATE(pstemp,pstemp1)&
    !$OMP SCHEDULE(STATIC) 
    DO ix=1,nx
       DO iy=1,ny
          DO iz=1,nz
             DO is=1,2
                pstemp=psi(ix,iy,iz,is,npmin(iq):npsi(iq))
                CALL zgemv('T',nlin,nlin,cmplxone,unitary,nlin,&
                     pstemp,1,cmplxzero,pstemp1,1)
                psi(ix,iy,iz,is,npmin(iq):npsi(iq))=pstemp1
             END DO
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE diagstep
  !*************************************************************************
  SUBROUTINE sinfo(printing)
    INTEGER :: il
    LOGICAL :: printing
    REAL(db):: tabc_energy, tabc_ekin, tabc_ecoul, tabc_eskyrme
    CHARACTER(*),PARAMETER :: &
         header='  #  Par   v**2   var_h1   var_h2    Norm     Ekin    Energy &
         &    Lx      Ly      Lz     Sx     Sy     Sz  '   
    ! calculate static observables for printout                       *
    CALL moments
    CALL integ_energy
    CALL sum_energy
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
       WRITE(scratch,'(1x,i5,2F8.3,12F15.7)') &
            iter,pnr,ehf,ehfint,tke,ehfc,ehfCrho0,ehfCrho1,ehfCdrho0,ehfCdrho1,ehfCtau0,&
            ehfCtau1,ehfCdJ0,ehfCdJ1
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=converfile,POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,f9.2,3(1pg11.3),2(0pf8.3),f6.1,f10.7)') &
            iter,ehf,delesum/pnrtot,efluct1,efluct2,rmstot,beta,gamma,x0dmp
       CLOSE(scratch)
       OPEN(unit=scratch,file=dipolesfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,6E14.4)') iter,cmtot,cm(:,2)-cm(:,1)
       CLOSE(unit=scratch)
       OPEN(unit=scratch,file=spinfile, POSITION='APPEND')  
       WRITE(scratch,'(1x,i5,9F10.4)') iter,orbital,spin,total_angmom 
       CLOSE(unit=scratch)
       WRITE(*,'(/,A,I7,A/2(A,F12.4),A/(3(A,E12.5),A))') &
            ' ***** Iteration ',iter,' *************************************************&
            &***********************************',' Total energy: ',ehf,&
            ' MeV  Total kinetic energy: ', tke,' MeV',' de/e:      ',delesum,&
            '      h**2  fluct.:    ',efluct1,' MeV, h*hfluc.:    ',efluct2,' MeV', &
            ' MeV. Rearrangement E: ',e3corr,' MeV. Coul.Rearr.: ', &
            ecorc,' MeV'
       ! detail printout
       WRITE(*,'(/A)') ' Energies integrated from density functional:********************&
                  &********************************************'
       WRITE(*,'(4(A,1PE14.6),A/26X,3(A,1PE14.6),A)') &
            ' Total:',ehfint,' MeV. t0 part:',ehf0,' MeV. t1 part:',ehf1, &
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
         WRITE(*,'(a)') '   e_ferm      e_pair     aver_gap    aver_force '
         DO il=1,2  
            WRITE(*,'(a,i2,a,4(1pg12.4))') 'iq=',il,': ',eferm(il) , &
               epair(il) ,avdelt(il),avg(il)
         ENDDO
         WRITE(*,*)'**********************************************************************&
                  &**************************************'
       END IF
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
          WRITE(*,'(1X,I3,F4.0,F8.5,2F9.5,F9.6,F8.3,F10.3,3F8.3,3F7.3)') &
               il,sp_parity(il),wocc(il),sp_efluct1(il),sp_efluct2(il), &
               sp_norm(il),sp_kinetic(il),sp_energy(il), &
               sp_orbital(:,il),sp_spin(:,il)
       ENDDO
       CALL moment_print
    END IF
  END SUBROUTINE sinfo
  !*************************************************************************
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
