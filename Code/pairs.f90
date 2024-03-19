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
  USE Params, ONLY: db,iter,printnow,wflag
  USE Forces, ONLY: ipair,p
  USE Grids, ONLY: nx,ny,nz,wxyz
  USE Densities, ONLY:rho
  USE Levels
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pair,epair,avdelt
  INTEGER :: iq                          !<Index labeling the isospin.
  REAL(db),SAVE :: eferm(2)              !<Fermi energy in MeV for the two isospins.
  REAL(db),SAVE :: epair(2)              !<Pairing energy in MeV for the two isospins. It is given by
  !!\f[ E_{\rm pair}=\frac1{2}\sum_k \Delta_k u_k v_k. \f]
  !!This is a public variable and the sum of the two values is subtracted from the
  !total energies in module \c Energies.  
  REAL(db),SAVE :: avdelt(2)             !<Average gap in MeV for the two isospins. It is
  !!given by \f[ \frac{\sum_k\Delta_ku_kv_k}{\sum_k u_kv_k} \f] where the
  !!sum is over states with the given isospin.
  REAL(db),SAVE :: avg(2)                !<The average pairing force for each isospin, given by
  !!\f[ \frac{E_{\rm pair}}{\sum_k u_k v_k/2}.\f]
  REAL(db),SAVE,ALLOCATABLE :: deltaf(:) !<single-particle gap in MeV for each
  !!single-particle state.
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
    ENDDO
    ! print pairing information
    IF(printnow.AND.wflag) THEN  
       WRITE(*,'(/7x,a)') '   e_ferm      e_pair     aver_gap    aver_force '
       DO iq=1,2  
          WRITE(*,'(a,i2,a,4(1pg12.4))') 'iq=',iq,': ',eferm(iq) , &
               epair(iq) ,avdelt(iq),avg(iq)
       ENDDO
    ENDIF
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
!!\f[ V_P(\vec r)={\tt v0act}\cdot{\tt work}\cdot (1-x\rho(\vec r)/{\tt rho0pr}) \f]
!!involving the <em> total </em> density \f$ \rho \f$ and the parameter 
!!\f[ rho0pr \f] from the pairing force definition and the parameter \f$x\f$ is used to create 
!! mix pairing when its value is set to 0.5, it is read from the pairing definition in the 
!! forces.data file.
!!
!!In the final step the gaps are computed as the expectation values of
!!the pairing field,
!!\f[ \Delta_k=\int\D^3r V_P(\vec r)\left|\phi_k(\vec r)\right|^2. \f]
!--------------------------------------------------------------------------- 
  SUBROUTINE pairgap
    !     computes the pair density
    INTEGER,PARAMETER :: itrsin=10
    ! iqq is different from module variable iq
    INTEGER :: iqq,nst,is
    REAL(db),PARAMETER :: smallp=0.000001  
    REAL(db) :: v0act
    REAL(db) :: work(nx,ny,nz)  
    ! constant gap in early stages of iteration
    IF(iter<=itrsin) THEN  
       deltaf=11.2/SQRT(mass_number)
       RETURN  
    ENDIF
    ! now the detailed gaps:
    DO iqq=1,2  
       work=0.D0
       ! accumulate new pair-density
       DO nst=npmin(iqq),npsi(iqq)  
          DO is=1,2  
             work=work+SQRT(MAX(wocc(nst)-wocc(nst)**2,smallp))*0.5* &
                  (REAL(psi(:,:,:,is,nst))**2+AIMAG(psi(:,:,:,is,nst))**2)
          ENDDO
       ENDDO
       ! determine pairing strength
       IF(iqq==2) THEN  
          v0act=p%v0prot  
       ELSE  
          v0act=p%v0neut  
       ENDIF
       ! now multiply with strength to obtain local pair-potential
       IF(ipair==6) THEN
          work=v0act*work*(1D0-p%mixture*(rho(:,:,:,1)+rho(:,:,:,2))/p%rho0pr)
       ELSE
          work=v0act*work  
       END IF
       ! finally compute the actual gaps as s.p. expectation values with
       ! the pair potential
       DO nst=npmin(iqq),npsi(iqq)  
          deltaf(nst)=0.0D0  
          DO is=1,2  
             deltaf(nst)=deltaf(nst)+wxyz*SUM(work* &
                  (REAL(psi(:,:,:,is,nst))**2+AIMAG(psi(:,:,:,is,nst))**2))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE pairgap
!---------------------------------------------------------------------------  
! DESCRIPTION: Routinename
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
    REAL(db) :: sumuv,sumduv,sumepa,edif,equasi,v2,vol
    ! start with non-pairing value of Fermi energy
    it=npmin(iq)+NINT(particle_number)-1  
    eferm(iq)=0.5D0*(sp_energy(it)+sp_energy(it+1))  
    wocc(npmin(iq):it)=1.D0
    wocc(it+1:npsi(iq))=0.D0
    ! determine true fermi energy for given 'delta' and 'e'
    eferm(iq)=rbrent(particle_number)
    ! finally compute average gap and pairing energy
    sumuv=0.0D0  
    sumduv=0.0D0 
    sumepa=0.0D0
    DO na=npmin(iq),npsi(iq)
       edif=sp_energy(na)-eferm(iq)  
       equasi=SQRT(edif*edif+deltaf(na)**2)  
       v2=0.5D0-0.5D0*edif/equasi  
       vol=0.5D0*SQRT(MAX(v2-v2*v2,xsmall))  
       sumuv=vol+sumuv  
       sumduv=vol*deltaf(na)+sumduv  
       sumepa=0.5D0*deltaf(na)**2/equasi+sumepa  
    ENDDO
    sumuv=MAX(sumuv,xsmall)  
    avdelt(iq)=sumduv/sumuv  
    epair(iq)=sumduv  
    avg(iq)=epair(iq)/sumuv**2  
  END SUBROUTINE pairdn
!---------------------------------------------------------------------------  
! DESCRIPTION: rbrent
!> @brief
!!This subroutine is an adapted version of the <em> Van
!!  Wijngaarden-Dekker-Brent</em> method for finding the root of a function
!!(see \cite Pre92aB ). Given the desired particle number as an
!!argument, it searches for the value of the Fermi energy that makes
!!this particle number agree with that returned by \c bcs_occupation.
!>
!> @details
!!It is clear that this subroutine is in a very antiquated style of
!!Fortran; it will be replaced at some time in the future.
!>
!> @param[in] particle_number
!> REAL(db), takes the particle number. 
!--------------------------------------------------------------------------- 
  FUNCTION rbrent(particle_number) RESULT(res)
    REAL(db),INTENT(in) :: particle_number
    REAL(db) :: res
    REAL(db),PARAMETER :: ra0=-100.D0,rb0=100.D0,reps=1.d-14,rtol=1.d-16
    INTEGER,PARAMETER :: iitmax=100
    REAL(db) :: ra,rb,rfa,rfb,rfc,rc,rd,re,rtol1,rxm,rp,rq,rr,rs, &
         rmin1,rmin2
    INTEGER :: ii
    ra=ra0  
    rb=rb0
    CALL bcs_occupation(ra,rfa)
    rfa=particle_number-rfa
    CALL bcs_occupation(rb,rfb)
    rfb=particle_number-rfb
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
             CALL bcs_occupation(res,rfa)
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
          CALL bcs_occupation(rb,rfb)
          rfb=particle_number-rfb
       ENDDO iteration
       STOP 'No solution found in pairing iterations'  
    ENDIF
  END FUNCTION rbrent
!---------------------------------------------------------------------------  
! DESCRIPTION: Routinename
!> @brief
!!For a given Fermi energy \f$ \epsilon_F \f$ passed as argument \c efermi,
!!this subroutine evaluates the particle number that would result with
!!such a Fermi energy and returns it as its second argument, \c bcs_partnum. 
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
!> REAL(db), returns the oarticle number
!--------------------------------------------------------------------------- 
  SUBROUTINE bcs_occupation(efermi,bcs_partnum)
    REAL(db),PARAMETER :: smal=1.0d-10  
    REAL(db),INTENT(in) :: efermi
    REAL(db),INTENT(OUT) :: bcs_partnum
    INTEGER :: k
    REAL(db) :: edif
    bcs_partnum=0.0  
    DO k=npmin(iq),npsi(iq)  
       edif=sp_energy(k)-efermi  
       wocc(k)=0.5D0 *(1.0D0-edif/SQRT(edif**2+deltaf(k)**2))  
       wocc(k)=MIN(MAX(wocc(k),smal),1.D0-smal)  
       bcs_partnum=bcs_partnum+wocc(k)
    ENDDO
  END SUBROUTINE bcs_occupation
END MODULE Pairs
