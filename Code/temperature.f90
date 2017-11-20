!------------------------------------------------------------------------------
! MODULE: Temperature
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!! Calculates the Fermi occupation distribution for a given temperature
!------------------------------------------------------------------------------
MODULE Temperature

USE Levels, ONLY: wocc,charge_number,mass_number,sp_energy,npmin,npsi
USE Params, ONLY: db

IMPLICIT NONE
REAL(db) :: kbT,chem_pot(2)
INTEGER :: iq

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: temp_dist
!> @brief
!!This routine sets the occupation probabilities for a given temperature
!--------------------------------------------------------------------------- 
  SUBROUTINE temp_dist
    IMPLICIT NONE
    REAL(db) :: particle_number
    DO iq=1,2
      IF(iq==1) particle_number=mass_number-charge_number
      IF(iq==2) particle_number=charge_number
      chem_pot(iq)=rbrent(particle_number,fermi_dist)
      WRITE(*,*)'chem. potential',iq,chem_pot(iq)
    END DO
  END SUBROUTINE temp_dist
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
! DESCRIPTION: fermi_dist
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
  REAL(db) FUNCTION fermi_dist(efermi)
    REAL(db),PARAMETER :: smal=1.0d-20  
    REAL(db),INTENT(in) :: efermi
    INTEGER :: k
    REAL(db) :: edif,fermi_accum
    fermi_accum=0D0
    DO k=npmin(iq),npsi(iq)  
       edif=sp_energy(k)-efermi  
       wocc(k)=1D0/(1.0D0+EXP(edif/kbT))
       wocc(k)=MIN(MAX(wocc(k),smal),1.D0-smal)  
       fermi_accum=fermi_accum+wocc(k)
    ENDDO
    fermi_dist=fermi_accum
  END FUNCTION fermi_dist
  
END MODULE Temperature
