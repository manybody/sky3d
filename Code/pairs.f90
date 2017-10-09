!------------------------------------------------------------------------------
! MODULE: Pairs
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!The principal part of this module is the subroutine \c pair, which
!!computes the pairing solution based on the BCS model. It is the only
!!public part of this module. The other subroutines are helper routines
!!that directly use the single-particle properties defined in module
!!\c Levels. The module variable \c iq controls whether the
!!solution is sought for neutrons (<tt> iq=1 </tt> or protons <tt> iq=2 </tt>
!!and accordingly the single-particle levels from <tt> npmin(iq) </tt> to
!!<tt> npsi(iq) </tt> are affected.
!>
!>@details
!!The principal procedure followed is to first calculate the pairing gap
!!for each single-particle state. This determines the occupation numbers
!!\c wocc, which of course are used throughout the program.  Then the
!!Fermi energy for the given isospin is determined such that the correct
!!particle number results.
!!
!!<b> Note that there is a factor of one half in many formulas compared
!!to what is usually found in textbooks. This is because here the sum
!!over states \f$ \sum_k\ldots \f$ runs over all states, while in textbooks
!!the sum is over pairs, giving half of that result.</b>
!------------------------------------------------------------------------------
MODULE Pairs
  USE Params, ONLY: db,iter,printnow,wflag,hbc
  USE Forces, ONLY: ipair,p,pair_reg,delta_fit,pair_cutoff,cutoff_factor,ecut_stab
  USE Grids, ONLY: nx,ny,nz,wxyz
  USE Densities, ONLY:rho
  USE Levels
  USE Parallel, ONLY:globalindex,collect_density,collect_sp_property,tmpi,wflag
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pair,epair,avdelt,avg,eferm,avdeltv2,eferm_cutoff,partnum_cutoff
  INTEGER :: iq                          !<Index labeling the isospin.
  REAL(db),SAVE :: eferm(2)=(/0D0,0D0/)  !<Fermi energy in MeV for the two isospins.
  REAL(db),SAVE :: epair(2)              !<Pairing energy in MeV for the two isospins. It is given by
  !!\f[ E_{\rm pair}=\frac1{2}\sum_k \Delta_k u_k v_k. \f]
  !!This is a public variable and the sum of the two values is subtracted from the
  !total energies in module \c Energies.  
  REAL(db),SAVE :: avdelt(2)             !<Average gap in MeV for the two isospins. It is
  !!given by \f[ \frac{\sum_k\Delta_ku_kv_k}{\sum_k u_kv_k} \f] where the
  !!sum is over states with the given isospin.
  REAL(db),SAVE :: avdeltv2(2)             !<Average gap in MeV for the two isospins with weight
  !!\f$ v^2 \f$. It is
  !!given by \f[ \frac{\sum_k\Delta_kv_k^2}{\sum_k u_kv_k} \f] where the
  !!sum is over states with the given isospin.
  REAL(db),SAVE :: avg(2)                !<The average pairing force for each isospin, given by
  !!\f[ \frac{E_{\rm pair}}{\sum_k u_k v_k/2}.\f]
  REAL(db),SAVE,ALLOCATABLE :: deltaf(:) !<single-particle gap in MeV for each
  !!single-particle state.
  REAL(db),SAVE,ALLOCATABLE,PUBLIC :: pairwg(:) !<Weights for pairing phase space
  !!in case of soft cutoff function often denoted \f$ f_n \f$.
  REAL(db),SAVE :: delesmooth            !<The width of the Fermi distribution defining pairing space.
  REAL(db),SAVE :: eferm_cutoff(2)=(/0D0,0D0/)   !<The Fermi energy for pairing cutoff.
  REAL(db),SAVE :: partnum_cutoff(2)     !<Particle number in pairing phase space.
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: pair
!> @brief
!!This is the only routine visible from outside the module. It solves
!!the pairing problem and prints out summary information. The principal
!!results used in the rest of the code are the BCS occupation numbers
!!\f$ v_k^2\rightarrow {\tt wocc} \f$ and the pairing energies \c epair.
!>
!> @details
!!The subroutine is structured straightforwardly: it first calculates
!!the pairing gaps \c deltaf by calling \c pairgap. Then for the
!!two isospin values \c pairdn is called with the correct particle
!!number as argument. This does the real work of solving the equations.
!!Finally summary information is printed.
!--------------------------------------------------------------------------- 
  SUBROUTINE pair  
    REAL(db) :: particle_number
    ! prepare phase space weight for subroutines
    IF(.NOT.ALLOCATED(deltaf)) ALLOCATE(deltaf(nstmax))
    IF(.NOT.ALLOCATED(pairwg)) THEN
      ALLOCATE(pairwg(nstmax))
      pairwg=1D0                        ! default if no soft cutoff is invoked
    END IF
    ! calculate gaps
    CALL pairgap
    ! solve pairing problem for each isospin iq
    DO iq=1,2  
       IF(iq==2) THEN
          particle_number=charge_number
       ELSE  
          particle_number=mass_number-charge_number
       ENDIF
       CALL pairdn(particle_number)
       IF(pair_cutoff(iq)>0.0d0.AND.sp_energy(npsi(iq))-eferm(iq)<pair_cutoff(iq)) &
       WRITE(*,*) '**Warning: Not enough particles for requested cutoff. Actual cutoff is: ',&
                  sp_energy(npsi(iq))-eferm(iq), 'iq= ',iq
    ENDDO
    ! print pairing information
    IF(printnow.AND.wflag) THEN  
       WRITE(*,'(/7x,a)') '   e_ferm      e_pair     <uv delta>   <v2 delta>    aver_force '
       DO iq=1,2  
          WRITE(*,'(a,i2,a,5(1pg12.4))') 'iq=',iq,': ',eferm(iq) , &
               epair(iq) ,avdelt(iq), avdeltv2(iq),avg(iq)
       ENDDO
!       IF(cutoff_factor>0D0) THEN
       WRITE(*,'(/7x,a)') '  e_ferm_cut    pairing-band   pairing space '
       DO iq=1,2  
          WRITE(*,'(a,i2,a,5(1pg12.4))') 'iq=',iq,': ',eferm_cutoff(iq) , &
               eferm_cutoff(iq)-eferm(iq),partnum_cutoff(iq)
       ENDDO
!       END IF
    ENDIF
    CALL flush(6)
  END SUBROUTINE pair
!---------------------------------------------------------------------------  
! DESCRIPTION: pairgap
!> @brief
!!This subroutine calculates the pairing gaps \f$ \Delta_k \f$ stored in the
!!array \c deltaf for all single-particle states.
!>
!> @details
!!First a simplified version is returned if <tt> ipair=1 </tt> or for the
!!initial <tt> itrsin </tt> (at present set to 10) iterations of a static
!!calculation. All gaps are set equal to \f$ 11.2\,{\rm MeV}/\sqrt{A} \f$ in
!!this case.
!!
!!In the general case, there is a loop over the two isospin values 
!!\c iq. The pairing density is obtained by evaluation of
!!\f[ {\tt work}(\vec r)=\sum_k u_kv_k\left|\phi_k(\vec r)\right|^2 \f]
!!where the simple conversion
!!\f[ u_kv_k=v_k\sqrt{1-v_k^2}=\sqrt{v_k^2(1-v_k^2)}=\sqrt{{\tt wocc}-{\tt wocc}^2} \f]
!!is used.
!!
!!The isospin-dependent pairing strength \c v0act is obtained from
!!the force definition. The pairing field \f$ V_P(\vec r) \f$ is then given by
!!two different expressions: for \c VDI pairing (<tt> ipair=5 </tt>),
!!the pairing density is simply multiplied by \c v0act, while for
!!\c DDDI pairing (<tt> ipair=6</tt>) it is
!!\f[ V_P(\vec r)={\tt v0act}\cdot{\tt work}\cdot (1-\rho(\vec r))/{\tt rho0pr} \f]
!!involving the <em> total </em> density \f$ \rho \f$ and the parameter 
!!\f[ rho0pr \f] from the pairing force definition.
!!
!!In the final step the gaps are computed as the expectation values of
!!the pairing field,
!!\f[ \Delta_k=\int\D^3r V_P(\vec r)\left|\phi_k(\vec r)\right|^2. \f]
!--------------------------------------------------------------------------- 
  SUBROUTINE pairgap
    !     computes the pair density
    INTEGER,PARAMETER :: itrsin=10
    ! iqq is different from module variable iq
    INTEGER :: iqq,nst,is,nstind
    REAL(db),PARAMETER :: smallp=0.000001D0
    REAL(db) :: v0act,v2,vol,sumduv,edif
    REAL(db),ALLOCATABLE :: work(:,:,:,:)  
    ! constant gap in early stages of iteration

    ALLOCATE(work(nx,ny,nz,2))
    IF(iter<=itrsin) THEN  
       deltaf=11.2/SQRT(mass_number)
       RETURN  
    ENDIF
    ! now the detailed gaps:
    work=0.D0
    DO iqq=1,2  
       ! accumulate new pair-density
       !DO nst=npmin(iqq),npsi(iqq)  
       DO is=1,2 
         DO nst=1,nstloc 
           nstind=globalindex(nst)
           IF(isospin(nstind)==iqq)&
            work(:,:,:,iqq)=work(:,:,:,iqq)+ pairwg(nstind)* &
              SQRT(MAX(wocc(nstind)-wocc(nstind)**2,smallp))*0.5* &
             (REAL(psi(:,:,:,is,nst))**2+AIMAG(psi(:,:,:,is,nst))**2)
         ENDDO
       ENDDO
    END DO
    IF(tmpi) CALL collect_density(work)
    IF(wflag) WRITE(*,*)1,SUM(work(:,:,:,1))
    IF(wflag) WRITE(*,*)2,SUM(work(:,:,:,2))
       ! determine pairing strength
    deltaf=0.0d0
    DO iqq=1,2
       IF(iqq==2) THEN  
          v0act=p%v0prot  
       ELSE  
          v0act=p%v0neut  
       ENDIF
       ! now multiply with strength to obtain local pair-potential
       IF(ipair==6) THEN
         IF(pair_reg) THEN
           work(:,:,:,iqq)=g_eff(v0act,iqq)*work(:,:,:,iqq)*(1D0-(rho(:,:,:,1)+rho(:,:,:,2))/p%rho0pr)
         ELSE
           work(:,:,:,iqq)=v0act*work(:,:,:,iqq)*(1D0-(rho(:,:,:,1)+rho(:,:,:,2))/p%rho0pr)
         END IF
       ELSE
          work(:,:,:,iqq)=v0act*work(:,:,:,iqq)  
       END IF
       ! finally compute the actual gaps as s.p. expectation values with
       ! the pair potential
       DO nst=1,nstloc
          nstind=globalindex(nst)
          IF(isospin(nstind)==iqq) THEN
             DO is=1,2  
                deltaf(nstind)=deltaf(nstind)+wxyz*SUM(work(:,:,:,iqq)* &
                pairwg(nstind)* &
                (REAL(psi(:,:,:,is,nst))**2+AIMAG(psi(:,:,:,is,nst))**2))
             ENDDO 
          END IF
       ENDDO
    ENDDO

    DEALLOCATE(work)

!!! PGR
    IF(ecut_stab.NE.0D0) THEN
       ! compute pairing energy from scratch
       sumduv=0.0D0 
       DO iqq=1,2
          DO nst=npmin(iqq),npsi(iqq)
             edif=sp_energy(nst)-eferm(iqq)  
             v2=0.5D0-0.5D0*edif/SQRT(edif*edif+(deltaf(nst)*pairwg(nst))**2)  
             vol=0.5D0*SQRT(MAX(v2-v2*v2,1D-20))  
             sumduv=vol*deltaf(nst)+sumduv  
          ENDDO
          epair(iqq)=sumduv  
       END DO
       ! now modify gap
       DO nst=1,nstloc
          nstind=globalindex(nst)
          IF(isospin(nstind)==1) THEN
            deltaf(nstind) = deltaf(nstind)*(1D0+(ecut_stab/epair(1))**2)
          ELSE
            deltaf(nstind) = deltaf(nstind)*(1D0+(ecut_stab/epair(2))**2)
          END IF
          deltaf(nstind) = min(deltaf(nstind),5D0)       
       END DO


    END IF
!!! PGR

    IF(tmpi) CALL collect_sp_property(deltaf)
    
  END SUBROUTINE pairgap
!---------------------------------------------------------------------------  
! DESCRIPTION: pairdn
!> @brief
!!The subroutine \c pairdn determines the pairing solution by using
!!\c rbrent to find the correct Fermi energy for the given particle
!!number. After that a few averaged or integral quantities are
!!calculated.
!>
!> @details
!!As the starting value for the Fermi energy the one for gap zero is
!!used, i. e., the average of the first unfilled and last filled
!!single-particle energies. Then \c rbrent is called to calculate the
!!correct solution, after which there is a loop for a straightforward
!!evaluation of the module variables \c epair, \c avdelt,
!!and \c avg. 
!> @param[in] particle_number
!> REAL(db), takes the particle number. 
!--------------------------------------------------------------------------- 
  SUBROUTINE pairdn(particle_number)
    REAL(db),INTENT(IN) :: particle_number
    REAL(db),PARAMETER :: xsmall=1.d-20
    INTEGER :: it,na
    REAL(db) :: sumuv,sumduv,sumepa,edif,equasi,v2,vol,sumv2,sumdv2


    ! start with non-pairing value of Fermi energy        !!! only first time
    IF(eferm(iq)==0D0) THEN
      it=npmin(iq)+NINT(particle_number)-1  
      eferm(iq)=0.5D0*(sp_energy(it)+sp_energy(it+1))  
      wocc(npmin(iq):it)=1.D0
      wocc(it+1:npsi(iq))=0.D0
    END IF

    ! Determine soft cutoff for pairing space
    IF(cutoff_factor>0D0) THEN
      ! In first call, determine guess for Fermi energy of cutoff 
      partnum_cutoff(iq)=particle_number+&
                         cutoff_factor*particle_number**0.6666666667D0
      IF(eferm_cutoff(iq)==0D0) THEN
!        it=npmin(iq)+NINT(partnum_cutoff(iq))-1
        it=NINT(partnum_cutoff(iq))-1
        IF(it+NINT(particle_number/10)>npsi(iq)) THEN
          WRITE(*,*) 'particle_number etc:',particle_number,npmin(iq), &
            it,it+NINT(particle_number/10),npsi(iq)
          STOP "not enough states to cover pairing space"
        END IF
        it=npmin(iq)+it
        eferm_cutoff(iq)=0.5D0*(sp_energy(it)+sp_energy(it+1))  
        pairwg(npmin(iq):it)=1.D0
        pairwg(it+1:npsi(iq))=0.D0
      END IF
      delesmooth=(eferm_cutoff(iq)-eferm(iq))/1D1
      ! determine 'pairgw' within function 'fermi_partnum'
      eferm_cutoff(iq)=rbrent(partnum_cutoff(iq),fermi_partnum)
      IF(eferm_cutoff(iq)<eferm(iq)) THEN
        WRITE(*,*) 'eferm_cutoff(1:2),eferm(1:2)=',eferm_cutoff,eferm
        STOP "cutoff Fermi energy below pairing Fermi energy"
      END IF
    END IF

    ! determine true fermi energy for given 'delta' and 'e'
    eferm(iq)=rbrent(particle_number,bcs_partnum)
    ! finally compute average gap and pairing energy
    sumuv=0.0D0  
    sumduv=0.0D0 
    sumepa=0.0D0
    sumv2=0.0D0
    sumdv2=0.0D0
    DO na=npmin(iq),npsi(iq)
       edif=sp_energy(na)-eferm(iq)  
       equasi=SQRT(edif*edif+deltaf(na)**2)  
       v2=0.5D0-0.5D0*edif/equasi  
       vol=0.5D0*SQRT(MAX(v2-v2*v2,xsmall))  
       sumuv=vol+sumuv  
       sumduv=vol*deltaf(na)+sumduv  
       sumv2=sumv2+v2
       sumdv2=sumdv2+deltaf(na)*v2
       sumepa=0.5D0*deltaf(na)**2/equasi+sumepa  
    ENDDO
    sumuv=MAX(sumuv,xsmall)  
    avdelt(iq)=sumduv/sumuv
    avdeltv2(iq)=sumdv2/sumv2  
    epair(iq)=sumduv              !!! PGR: to be checked ???
    avg(iq)=epair(iq)/sumuv**2  
!    IF(iq==1) WRITE(*,*) avdelt(1),delta_fit(iq)+1d-2,delta_fit(iq)-1d-2
    IF(delta_fit(iq)>1.0d-5.AND.avdeltv2(iq)>delta_fit(iq)+1d-3.AND.iq==1) THEN
      p%v0neut=p%v0neut-MAX(ABS(avdeltv2(iq)-delta_fit(iq)),0.1d-1)
      WRITE(*,*) ' V0_neut= ',p%v0neut
    END IF
    IF(delta_fit(iq)>1.0d-5.AND.avdeltv2(iq)<delta_fit(iq)-1d-3.AND.iq==1) THEN
      p%v0neut=p%v0neut+MAX(ABS(avdeltv2(iq)-delta_fit(iq)),0.1d-1)
      WRITE(*,*) ' V0_neut= ',p%v0neut
    END IF
    IF(delta_fit(iq)>1.0d-5.AND.avdeltv2(iq)>delta_fit(iq)+1d-3.AND.iq==2) THEN
      p%v0prot=p%v0prot-MAX(ABS(avdeltv2(iq)-delta_fit(iq)),0.1d-1)
      WRITE(*,*) ' V0_prot= ',p%v0prot
    END IF
    IF(delta_fit(iq)>1.0d-5.AND.avdeltv2(iq)<delta_fit(iq)-1d-3.AND.iq==2) THEN
      p%v0prot=p%v0prot+MAX(ABS(avdeltv2(iq)-delta_fit(iq)),0.1d-1)
      WRITE(*,*) ' V0_prot= ',p%v0prot
    END IF
  END SUBROUTINE pairdn
!---------------------------------------------------------------------------  
! DESCRIPTION: rbrent
!> @brief
!!This subroutine is an adapted version of the <em> Van
!!  Wijngaarden-Dekker-Brent</em> method for finding the root of a function
!!(see \cite Pre92aB ). Given the desired particle number as an
!!argument, it searches for the value of the Fermi energy that makes
!!this particle number agree with that returned by \c bcs_partnum.
!>
!> @details
!!It is clear that this subroutine is in a very antiquated style of
!!Fortran; it will be replaced at some time in the future.
!>
!> @param[in] particle_number
!> REAL(db), takes the particle number. 
!--------------------------------------------------------------------------- 
  FUNCTION rbrent(particle_number,funk) RESULT(res)
    !      data for wijngaarden-dekker-brent method for root
    !
    REAL(db),INTENT(in) :: particle_number
    REAL(db) :: res
    REAL(db),PARAMETER :: ra0=-100.D0,rb0=100.D0,reps=1.d-14,rtol=1.d-16
    INTEGER,PARAMETER :: iitmax=100
    REAL(db) :: ra,rb,rfa,rfb,rfc,rc,rd,re,rtol1,rxm,rp,rq,rr,rs, &
         rmin1,rmin2
    INTEGER :: ii
    REAL(db), EXTERNAL :: funk
!   INTERFACE
!     REAL(8) FUNCTION funk(x)
!     REAL(8),INTENT(in) :: x
!     END FUNCTION funk
!   END INTERFACE
    ra=ra0  
    rb=rb0
    rfa=particle_number-funk(ra)
    rfb=particle_number-funk(rb)
    rfc=rfb  
    rc=rb  
    rd=rb-ra  
    re=rd  
    IF((rfa>0.D0) .AND.(rfb>0.D0)) THEN  
       IF(rfa>rfb) THEN  
          rb=1.0  
       ELSE  
          rb=0.D0  
       ENDIF
    ELSEIF((rfa<0.D0) .AND.(rfb<0.D0)) THEN  
       IF(rfa>rfb) THEN  
          rb=0.D0  
       ELSE  
          rb=1.0  
       ENDIF
    ELSE  
       iteration: DO ii=1,iitmax  
          IF(((rfb>0.D0) .AND.(rfc>0.D0)).OR.((rfb<0.D0) &
               .AND.(rfc<0.D0))) THEN
             rc=ra  
             rfc=rfa  
             rd=rb-ra  
             re=rd  
          ENDIF
          IF(ABS(rfc)<ABS(rfb)) THEN  
             ra=rb  
             rb=rc  
             rc=ra  
             rfa=rfb  
             rfb=rfc  
             rfc=rfa  
          ENDIF
          !     convergence check
          rtol1=2.0*reps*ABS(rb)+0.5*rtol  
          rxm=0.5 *(rc-rb)  
          IF((ABS(rxm) <=rtol1).OR.((ABS(rfb)==0.D0))) THEN
             res=rb
             rfa = funk(res)
             RETURN  
          ENDIF
          IF((ABS(re)>=rtol1).OR.(ABS(rfa)>ABS(rfb))) &
               THEN
             rs=rfb/rfa  
             IF(ra==rc) THEN  
                rp=2.0*rxm*rs  
                rq=1.0-rs  
             ELSE  
                rq=rfa/rfc  
                rr=rfb/rfc  
                rp=rs *(2.0*rxm*rq *(rq-rr)-(rb-ra) *(rr- &
                     1.0))
                rq=(rq-1.0) *(rr-1.0) *(rs-1.0)  
             ENDIF
             IF(rp>0.D0) THEN  
                rq=-rq  
             ENDIF
             rp=ABS(rp)  
             rmin1=3.0*rxm*rq-ABS(rtol1*rq)  
             rmin2=ABS(rq*re)  
             IF(2.0*rp<MIN(rmin1,rmin2)) THEN  
                re=rd  
                rd=rp/rq  
             ELSE  
                rd=rxm  
                re=rd  
             ENDIF
          ELSE  
             rd=rxm  
             re=rd  
          ENDIF
          ra=rb  
          rfa=rfb  
          IF(ABS(rd) >rtol1) THEN  
             rb=rb+rd  
          ELSE  
             rb=rb+SIGN(rtol1,rxm)  
          ENDIF
          rfb=particle_number-funk(rb)
       ENDDO iteration
       STOP 'No solution found in pairing iterations'  
    ENDIF
  END FUNCTION rbrent
!---------------------------------------------------------------------------  
! DESCRIPTION: bcs_partnum
!> @brief
!!For a given Fermi energy \f$ \epsilon_F \f$ passed as argument \c efermi,
!!this subroutine evaluates the particle number that would result with
!!such a Fermi energy and returns it as function value \c bcs_partnum. 
!>
!> @details
!!The isospin is controlled by module variable \c iq. First the occupation 
!!probabilities are calculated using the standard BCS expression
!!\f[ v_k^2=\frac1{2}\left(1-\frac{\epsilon_k-\epsilon_F}
!!    {\sqrt{(\epsilon_k-\epsilon_F)^2+\Delta_k^2}}\right). \f]
!!They are stored in \c wocc. A small correction is added, so
!!that they are not exactly identical to 1 or0. The particle number
!!is finally obtained as \f$ N=\sum_k v_k^2 \f$.
!>
!> @param[in] efermi
!> REALD(db), takes the fermi energy.
!> @param[out] bcs_partnum
!> REAL(db), returns the particle number
!--------------------------------------------------------------------------- 
  REAL(db) FUNCTION bcs_partnum(efermi)
    REAL(db),PARAMETER :: smal=1.0d-10  
    REAL(db),INTENT(in) :: efermi
!    REAL(db),INTENT(OUT) :: bcs_partnum
    INTEGER :: k
    REAL(db) :: edif,bcs_accum
    bcs_accum=0D0
    DO k=npmin(iq),npsi(iq)  
       edif=sp_energy(k)-efermi  
       IF(pair_cutoff(iq)>0.0d0.AND.edif>pair_cutoff(iq)) THEN
         wocc(k)=smal
       ELSE
         wocc(k)=0.5D0 *(1.0D0-edif/SQRT(edif**2+(deltaf(k)*pairwg(k))**2))  
         wocc(k)=MIN(MAX(wocc(k),smal),1.D0-smal)  
         bcs_accum=bcs_accum+wocc(k)
       END IF
    ENDDO
    bcs_partnum=bcs_accum
  END FUNCTION bcs_partnum
!---------------------------------------------------------------------------  
! DESCRIPTION: fermi_partnum
!> @brief
!!For a given Fermi energy \f$ \epsilon_F \f$ passed as argument \c efermi,
!!this subroutine evaluates the particle number that would result for
!!a Fermi distribution with 
!!such a Fermi energy and returns it as function value \c fermi_partnum. 
!>
!> @details
!!The isospin is controlled by module variable \c iq. First the occupation 
!!probabilities are calculated using the standard BCS expression
!!\f[ v_k^2=\frac1{2}\left(1-\frac{\epsilon_k-\epsilon_F}
!!    {\sqrt{(\epsilon_k-\epsilon_F)^2+\Delta_k^2}}\right). \f]
!!They are stored in \c wocc. A small correction is added, so
!!that they are not exactly identical to 1 or0. The particle number
!!is finally obtained as \f$ N=\sum_k v_k^2 \f$.
!>
!> @param[in] efermi
!> REALD(db), takes the fermi energy.
!> @param[out] bcs_partnum
!> REAL(db), returns the particle number
!--------------------------------------------------------------------------- 
  REAL(db) FUNCTION fermi_partnum(efermi)
    REAL(db),PARAMETER :: smal=1.0d-20  
    REAL(db),INTENT(in) :: efermi
!    REAL(db),INTENT(OUT) :: bcs_partnum
    INTEGER :: k
    REAL(db) :: edif,fermi_accum
    fermi_accum=0D0
    DO k=npmin(iq),npsi(iq)  
       edif=sp_energy(k)-efermi  
       pairwg(k)=1D0/(1.0D0+EXP(edif/delesmooth))
       pairwg(k)=MIN(MAX(pairwg(k),smal),1.D0-smal)  
       fermi_accum=fermi_accum+pairwg(k)
    ENDDO
    fermi_partnum=fermi_accum
  END FUNCTION fermi_partnum
!---------------------------------------------------------------------------  
! DESCRIPTION: g_eff
!> @brief
!!g_eff is the regulated pairing strength.
!>
!> @param[in] g
!> REALD(db), takes the unregulated pairing strength.
!> @param[in] iq
!> REAL(db), takes the isospin.
!--------------------------------------------------------------------------- 
  FUNCTION g_eff(g,iq)
    USE Levels    , ONLY :  sp_energy
    REAL(db),INTENT(IN) :: g
    INTEGER ,INTENT(IN) :: iq
    REAL(db)            :: g_eff,delta
    REAL(db)            :: e_l,e_u
    e_l=sp_energy(npmin(iq))
    IF(pair_cutoff(iq)>0.0d0)THEN
      e_u=MIN(sp_energy(npsi(iq)),eferm(iq)+pair_cutoff(iq))
    ELSE
      e_u=sp_energy(npsi(iq))
    END IF
    IF (avdeltv2(iq)<0.001) avdeltv2(iq)=1.0d0
    delta=MIN(2.0d0,MAX(0.5,avdeltv2(iq)))
    g_eff=g/log((e_u-e_l)/delta)
    WRITE(*,*)iq,g_eff,g,e_l,e_u,delta,log((e_u-e_l)/delta)
  END FUNCTION g_eff
END MODULE Pairs
