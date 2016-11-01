MODULE Pairs
  USE Params, ONLY: db,iter,printnow,wflag,hbc
  USE Forces, ONLY: ipair,p,pair_reg,delta_fit,pair_cutoff
  USE Grids, ONLY: nx,ny,nz,wxyz
  USE Densities, ONLY:rho
  USE Levels
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pair,epair,avdelt,avg,eferm,avdeltv2
  INTEGER :: iq
  REAL(db),SAVE :: eferm(2),epair(2),avdelt(2),avdeltv2(2),avg(2)  
  REAL(db),SAVE,ALLOCATABLE :: deltaf(:)
CONTAINS
  !***********************************************************************
  ! Determine pairing solution
  !***********************************************************************
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
    ENDIF
  END SUBROUTINE pair
  !***********************************************************************
  ! Calculate pairing gaps for s.p. levels
  !***********************************************************************
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
         IF(pair_reg) THEN
           work=g_eff(v0act,iqq)*work*(1D0-(rho(:,:,:,1)+rho(:,:,:,2))/p%rho0pr)
         ELSE
           work=v0act*work*(1D0-(rho(:,:,:,1)+rho(:,:,:,2))/p%rho0pr)
         END IF
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
  !***********************************************************************
  ! Calculate the Fermi energy for isospin iq such that the particle
  ! number condition is fulfilled. 
  !***********************************************************************
  SUBROUTINE pairdn(particle_number)
    REAL(db),INTENT(IN) :: particle_number
    REAL(db),PARAMETER :: xsmall=1.d-20
    INTEGER :: it,na
    REAL(db) :: sumuv,sumduv,sumepa,edif,equasi,v2,vol,sumv2,sumdv2
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
    epair(iq)=sumduv  
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
  !***********************************************************************
  ! Search the solution for efermi where the particle number is correct
  !***********************************************************************
  FUNCTION rbrent(particle_number) RESULT(res)
    !      data for wijngaarden-dekker-brent method for root
    !
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
  !***********************************************************************
  ! Calc. occupation numbers and particle number for given Fermi energy
  !***********************************************************************
  SUBROUTINE bcs_occupation(efermi,bcs_partnum)
    REAL(db),PARAMETER :: smal=1.0d-10  
    REAL(db),INTENT(in) :: efermi
    REAL(db),INTENT(OUT) :: bcs_partnum
    INTEGER :: k
    REAL(db) :: edif
    bcs_partnum=0.0  
    DO k=npmin(iq),npsi(iq)  
       edif=sp_energy(k)-efermi  
       IF(pair_cutoff(iq)>0.0d0.AND.edif>pair_cutoff(iq)) THEN
         wocc(k)=smal
       ELSE
         wocc(k)=0.5D0 *(1.0D0-edif/SQRT(edif**2+deltaf(k)**2))  
         wocc(k)=MIN(MAX(wocc(k),smal),1.D0-smal)  
         bcs_partnum=bcs_partnum+wocc(k)
       END IF
    ENDDO
  END SUBROUTINE bcs_occupation
  !***********************************************************************
  ! Calc. pairing regulator
  !***********************************************************************
  FUNCTION g_eff(g,iq)
    USE Levels    , ONLY :  sp_energy
    REAL(db),INTENT(IN) :: g
    INTEGER ,INTENT(IN) :: iq
    REAL(db)            :: g_eff,delta
    REAL(db)            :: e_l,e_u,efermapprox,particle_number
    INTEGER             :: it
    e_l=sp_energy(npmin(iq))
    IF(pair_cutoff(iq)>0.0d0)THEN
      IF(iq==2) THEN
        particle_number=charge_number
      ELSE  
        particle_number=mass_number-charge_number
      ENDIF
      it=npmin(iq)+NINT(particle_number)-1  
      efermapprox=0.5D0*(sp_energy(it)+sp_energy(it+1))  
      e_u=MIN(sp_energy(npsi(iq)),efermapprox+pair_cutoff(iq))
    ELSE
      e_u=sp_energy(npsi(iq))
    END IF
    IF (avdeltv2(iq)<0.001) avdeltv2(iq)=1.0d0
    delta=12.0d0/sqrt(mass_number)
    g_eff=g/log((e_u-e_l)/delta)
    WRITE(*,*)iq,g_eff,g,e_l,e_u,delta,log((e_u-e_l)/delta)
  END FUNCTION g_eff
END MODULE Pairs
