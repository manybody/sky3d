!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module computes the total energy and the various contributions to
!!it in two ways. The first method evaluates the density functional, 
!!by direct integration to compute the <em> integrated
!!energy</em> \c ehfint. The second method uses the sum of
!!single-particle energies plus the rearrangement energies, \f$ E_{3,\rm
!!corr} \f$ for the density dependent part and \f$ E_{C,\rm corr} \f$ for
!!Coulomb exchange.
!>
!>@details 
!!The two ways of calculating the energy are assigned to the subroutines
!!\c integ_energy, which also calculates the rearrangement energies,
!!and \c sum_energy. Note that since \c integ_energy is always
!!called briefly before \c sum_energy, the rearrangement energies
!!are correctly available.
!!
!!In addition subroutine \c sum_energies also calculates the summed
!!spin, orbital, and total angular momenta.
!!
!!TODO: SHOULD WE ADD THE EQUATIONS FOR THE ENERGIES HERE?
!------------------------------------------------------------------------------
MODULE Energies
  USE Params, ONLY: db,tcoul
  USE Forces
  USE Densities
  USE Levels
  USE Pairs, ONLY: epair
  IMPLICIT NONE
  SAVE
  ! total energies calculated in subroutine "energy" and "hfenergy"
  REAL(db) :: ehft            !< kinetic energy
  REAL(db) :: ehf0            !< t0 contribution
  REAL(db) :: ehf1            !< b1 contribution (current part)
  REAL(db) :: ehf2            !< b2 contribution (Laplacian part)
  REAL(db) :: ehf3            !< t3 contribution which models density dependence
  REAL(db) :: ehfls           !< spin-orbit contribution (time even)
  REAL(db) :: ehflsodd        !< spin-orbit contribution (odd-odd)
  REAL(db) :: ehfc            !< Coulomb contribution
  REAL(db) :: ecorc           !< Slater & Koopman exchange
  REAL(db) :: ehfint          !< integrated total energy
  REAL(db) :: efluct1         !< energyfluctuation \f$ \hat{h}^2 \f$
  REAL(db) :: efluct1prev     !< energyfluctuation \f$ \hat{h}^2 \f$ of previous iteration
  REAL(db) :: efluct2         !< fluctuation \f$ \hat{h}\cdot {\tt efluct}\f$
  REAL(db) :: efluct2prev     !< fluctuation \f$ \hat{h}\cdot {\tt efluct}\f$ of previous iteration
  REAL(db) :: tke             !< kinetic energy summed
  REAL(db) :: ehf             !< Hartree-Fock energy from s.p. levels
  REAL(db) :: ehfprev         !< Hartree-Fock energy from s.p. levels of previous iteration
  REAL(db) :: e3corr          !< rearrangement energy
  REAL(db) :: orbital(3)      !< the three components of the total orbital
                              !!angular momentum in units of \f$ \hbar \f$.
  REAL(db) :: spin(3)         !< the three components of the total spin in
                              !!units of \f$ \hbar \f$.
  REAL(db) :: total_angmom(3) !< the three components of the total
                              !!angular momentum in units of \f$ \hbar \f$.
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: integ_energy
!> @brief
!!The purpose of this subroutine is to calculate the integrated energy.
!>
!> @details
!!This is implemented pretty straightforwardly. The only
!!programming technique worth noting is that intermediate variables such
!!as \c rhot for the total density are used to avoid repeating the
!!lengthy index lists. Compilers will eliminate these by optimization.
!!
!!In principle the integration loops in the subroutine could be
!!combined, but some space is saved by using the array \c worka for
!!different purposes in different loops.
!!
!!The calculation proceeds in the following steps:
!!  - <b> Step 1:</b> the Laplacian of the densities is calculated in
!!    \c worka, then the integrals for 
!!    \c ehf0,\c ehf2, and \c ehf3 are performed.
!!    After the loop the result for \c ehf3 is also used to calculate
!!    \c e3corr.
!!  - <b> Step 2:</b> the integral for \c ehf1 is evaluated
!!    using \c worka for the \f$ \vec\jmath_q{}^2 \f$ term.
!!  - <b> Step 3:</b> the spin-orbit contribution of 
!!    \c ehfls is calculated using \c worka as storage for
!!    \f$ \nabla\cdot\vec J_q \f$.
!!  - <b> Step 4:</b> the Coulomb energy \c ehfc is evaluated 
!!    with the Slater correction taken into account if the
!!    force's \c ex is nonzero. At the same time the Coulomb correction
!!    for the summed energy is calculated
!!    and stored in \c ecorc. It will be used in
!!    the subroutine \c sum_energy.
!!  - <b> Step 5: </b> the kinetic energy is integrated for 
!!    \c ehft. Note that only at this point the correct prefactor
!!    \f$ \hbar^2/2m \f$ is added; the use of \c tau in other expressions
!!    assumes its absence.
!!  - <b> Step 6: </b> Finally all terms are added to produce the total
!!  energy, \c efundet, from which the pairing energies are subtracted.
!--------------------------------------------------------------------------- 
  SUBROUTINE integ_energy
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    USE Grids, ONLY: wxyz,der1x,der2x,der1y,der2y,der1z,der2z
    USE Coulomb, ONLY: wcoul
    INTEGER :: ix,iy,iz,iq
    REAL(db) :: rhot,rhon,rhop,d2rho,d2rhon,d2rhop,sc
    REAL(db) :: worka(nx,ny,nz,2)
    REAL(db) :: workb(nx,ny,nz,3,2)
    ! Step 1: compute laplacian of densities, then ehf0, ehf2, and ehf3
    DO iq=1,2
       CALL rmulx(der2x,rho(:,:,:,iq),worka(:,:,:,iq),0)
       CALL rmuly(der2y,rho(:,:,:,iq),worka(:,:,:,iq),1)
       CALL rmulz(der2z,rho(:,:,:,iq),worka(:,:,:,iq),1)
    ENDDO
    ehf0=0.0D0
    ehf3=0.0D0
    ehf2=0.0D0
    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             rhot=rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
             rhon=rho(ix,iy,iz,1)
             rhop=rho(ix,iy,iz,2)
             ehf0=ehf0+wxyz*(b0*rhot**2-b0p*(rhop**2+rhon**2))/2.D0
             ehf3=ehf3+wxyz*rhot**f%power*(b3*rhot**2 &
                  -b3p*(rhop**2+rhon**2))/3.D0
             d2rho=worka(ix,iy,iz,1)+worka(ix,iy,iz,2)
             d2rhon=worka(ix,iy,iz,1)
             d2rhop=worka(ix,iy,iz,2)
             ehf2=ehf2+wxyz*(-b2*rhot*d2rho+b2p*(rhop* &
                  d2rhop+rhon*d2rhon))/2.D0
          ENDDO
       ENDDO
    ENDDO
    e3corr=-f%power*ehf3/2.D0
    ! Step 2: b1-contribution. worka=square of the current vector
    worka=current(:,:,:,1,:)**2+current(:,:,:,2,:)**2+current(:,:,:,3,:)**2
    ehf1=wxyz*SUM(b1*((rho(:,:,:,1)+rho(:,:,:,2))*(tau(:,:,:,1)+tau(:,:,:,2)) &
         -((current(:,:,:,1,1)+current(:,:,:,1,2))**2 &
         + (current(:,:,:,2,1)+current(:,:,:,2,2))**2 &
         + (current(:,:,:,3,1)+current(:,:,:,3,2))**2))&
         -b1p*(rho(:,:,:,1)*tau(:,:,:,1)+rho(:,:,:,2)*tau(:,:,:,2)&
         -worka(:,:,:,1)-worka(:,:,:,2)))
    ! Step 3: the spin-orbit contribution
    !         (a) Time even part: worka=div J
    DO iq=1,2
       CALL rmulx(der1x,sodens(:,:,:,1,iq),worka(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),worka(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),worka(:,:,:,iq),1)
    ENDDO
    ehfls=wxyz*SUM(-b4*(rho(:,:,:,1)+rho(:,:,:,2)) &
         *(worka(:,:,:,1)+worka(:,:,:,2)) &
         -b4p *(rho(:,:,:,2)*worka(:,:,:,2)+rho(:,:,:,1)*worka(:,:,:,1)))
    ! Step 3: the spin-orbit contribution
    !         (b) odd-odd part: workb = curl J
    DO iq = 1,2
       CALL rmuly(der1y,current(:,:,:,3,iq),workb(:,:,:,1,iq),0)
       CALL rmulz(der1z,current(:,:,:,2,iq),workb(:,:,:,1,iq),-1)
       CALL rmulz(der1z,current(:,:,:,1,iq),workb(:,:,:,2,iq),0)
       CALL rmulx(der1x,current(:,:,:,3,iq),workb(:,:,:,2,iq),-1)
       CALL rmulx(der1x,current(:,:,:,2,iq),workb(:,:,:,3,iq),0)
       CALL rmuly(der1y,current(:,:,:,1,iq),workb(:,:,:,3,iq),-1)
    ENDDO
    ehflsodd=wxyz*SUM(-b4*(sdens(:,:,:,:,1)+sdens(:,:,:,:,2))* &
         (workb(:,:,:,:,1)+workb(:,:,:,:,2))-b4p*(sdens(:,:,:,:,1)* &
         workb(:,:,:,:,1)+sdens(:,:,:,:,2)*workb(:,:,:,:,2)))
    ehfls = ehfls + ehflsodd
    !
    ! Step 4: Coulomb energy with Slater term, also correction
    ! term to Koopman formula
    ehfc=0.0D0
    ecorc=0.0D0
    IF(tcoul) THEN
       IF(f%ex/=0) THEN
          sc=-3.0D0/4.0D0*slate
       ELSE
          sc=0.D0
       END IF
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx
                rhop=rho(ix,iy,iz,2)
                ehfc=ehfc+wxyz *(0.5D0*rhop*wcoul(ix,iy,iz) &
                     +sc*rhop**(4.0D0/3.0D0))
                ecorc=ecorc+wxyz*sc/3.0D0*rhop**(4.0D0/3.0D0)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! Step 5: kinetic energy contribution
    ehft=wxyz*SUM(f%h2m(1)*tau(:,:,:,1)+f%h2m(2)*tau(:,:,:,2))
    ! Step 6: form total energy
    ehfint=ehft+ehf0+ehf1+ehf2+ehf3+ehfls+ehfc-epair(1)-epair(2)
  END SUBROUTINE integ_energy
!---------------------------------------------------------------------------  
! DESCRIPTION: sum_energy
!> @brief
!!This subroutine mainly computes the Koopman sum, but also
!!sums up a number of other single-particle properties. 
!>
!> @details
!!For systematics, the latter should be done in a different
!!place, but at present is left here.
!!
!!The summation of the total energy uses \c spenerg to compute 
!!\f[ \sum_k (\epsilon_k-\tfrac1{2}v_k)=\tfrac1{2}\sum_k(2t_k+v_k)=
!!\tfrac1{2}\sum_k(t_k+\epsilon_k). \f]
!!The last sum is calculated, the rearrangement corrections 
!!are added and the pairing energies subtracted.
!!
!!The subroutine then sums up the single-particle energy fluctuation
!!\c sp_efluct1 and \c sp_efluct2, dividing them by the nucleon
!!number.  Finally the orbital and spin angular
!!momentum components are summed to form the total ones. 
!--------------------------------------------------------------------------- 
  SUBROUTINE sum_energy
    USE Moment, ONLY: pnrtot
    INTEGER :: i
    ehf=SUM(wocc*(sp_kinetic+sp_energy))/2.D0+e3corr+ecorc &
         -epair(1)-epair(2)
    tke=SUM(wocc*sp_kinetic)
    efluct1=SUM(wocc*sp_efluct1)/pnrtot
    efluct2=SUM(wocc*sp_efluct2)/pnrtot
    DO i=1,3
       orbital(i)=SUM(wocc*sp_orbital(i,:))
       spin(i)=SUM(wocc*sp_spin(i,:))
       total_angmom(i)=orbital(i)+spin(i)
    END DO
  END SUBROUTINE sum_energy
  !***************************************************************************
END MODULE Energies

