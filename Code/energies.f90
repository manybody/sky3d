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
  USE Grids, ONLY: wxyz
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
  REAL(db) :: ecmcorr=0D0     !< energy of c.m. correction (a posteriori)
  REAL(db) :: orbital(3)      !< the three components of the total orbital
                              !!angular momentum in units of \f$ \hbar \f$.
  REAL(db) :: spin(3)         !< the three components of the total spin in
                              !!units of \f$ \hbar \f$.
  REAL(db) :: total_angmom(3) !< the three components of the total
                              !!angular momentum in units of \f$ \hbar \f$.
  REAL(db) :: e_extern=0D0    !< energy from external field (accumulator)
  REAL(db) :: ehfCrho0        !<\f$ C^\rho_0 \f$ contribution
  REAL(db) :: ehfCrho1        !<\f$ C^\rho_1 \f$ contribution
  REAL(db) :: ehfCdrho0       !<\f$ C^{\nabla\rho}_0 \f$ contribution
  REAL(db) :: ehfCdrho1       !<\f$ C^{\nabla\rho}_1 \f$ contribution
  REAL(db) :: ehfCtau0        !<\f$ C^{\tau}_0 \f$ contribution
  REAL(db) :: ehfCtau1        !<\f$ C^{\tau}_1 \f$ contribution
  REAL(db) :: ehfCdJ0         !<\f$ C^{\nabla\vec{J}}_0 \f$ contribution
  REAL(db) :: ehfCdJ1         !<\f$ C^{\nabla\vec{J}}_1 \f$ contribution
  REAL(db) :: ehfCj0          !<\f$ C^{\vec{j}}_0 \f$ contribution
  REAL(db) :: ehfCj1          !<\f$ C^{\vec{j}}_1 \f$ contribution
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
    REAL(db) :: rho0,rho1,d2rho0,d2rho1,tau0,tau1
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
    ehfCrho0=0.0d0
    ehfCrho1=0.0d0
    ehfCdrho0=0.0d0
    ehfCdrho1=0.0d0
    ehfCtau0=0.0d0
    ehfCtau1=0.0d0
    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             rhot=rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
             rho0=rhot
             rho1=-rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
             rhon=rho(ix,iy,iz,1)
             rhop=rho(ix,iy,iz,2)
             d2rhon=worka(ix,iy,iz,1)
             d2rhop=worka(ix,iy,iz,2)
             d2rho=worka(ix,iy,iz,1)+worka(ix,iy,iz,2)
             d2rho0=d2rho
             d2rho1=-worka(ix,iy,iz,1)+worka(ix,iy,iz,2)
             tau0=tau(ix,iy,iz,1)+tau(ix,iy,iz,2)
             tau1=-tau(ix,iy,iz,1)+tau(ix,iy,iz,2)
             ehf0=ehf0+wxyz*(b0*rhot**2-b0p*(rhop**2+rhon**2))/2.D0
             ehfCrho0=ehfCrho0+wxyz*(Crho0+Crho0D*rho0**f%power)*rho0**2
             ehfCrho1=ehfCrho1+wxyz*(Crho1+Crho1D*rho0**f%power)*rho1**2
             ehf3=ehf3+wxyz*rhot**f%power*(b3*rhot**2 &
                  -b3p*(rhop**2+rhon**2))/3.D0
             ehf2=ehf2+wxyz*(-b2*rhot*d2rho+b2p*(rhop* &
                  d2rhop+rhon*d2rhon))/2.D0
             ehfCdrho0=ehfCdrho0+wxyz*Cdrho0*rho0*d2rho0
             ehfCdrho1=ehfCdrho1+wxyz*Cdrho1*rho1*d2rho1
             ehfCtau0=ehfCtau0+wxyz*Ctau0*tau0*rho0
             ehfCtau1=ehfCtau1+wxyz*Ctau1*tau1*rho1
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
    ehf1=wxyz*SUM(b1*((rho(:,:,:,1)+rho(:,:,:,2))*(tau(:,:,:,1)+tau(:,:,:,2))) &
         -b1p*(rho(:,:,:,1)*tau(:,:,:,1)+rho(:,:,:,2)*tau(:,:,:,2)))
    ehfCj0=-wxyz*SUM(Ctau0*((current(:,:,:,1,1)+current(:,:,:,1,2))**2 &
         + (current(:,:,:,2,1)+current(:,:,:,2,2))**2 &
         + (current(:,:,:,3,1)+current(:,:,:,3,2))**2))
    ehfCj1=-wxyz*SUM(Ctau1*((-current(:,:,:,1,1)+current(:,:,:,1,2))**2 &
         + (-current(:,:,:,2,1)+current(:,:,:,2,2))**2 &
         + (-current(:,:,:,3,1)+current(:,:,:,3,2))**2))
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
    ehfCdJ0=wxyz*sum(CdJ0*(rho(:,:,:,1)+rho(:,:,:,2))*(worka(:,:,:,1)+worka(:,:,:,2)))
    ehfCdJ1=wxyz*sum(CdJ1*(-rho(:,:,:,1)+rho(:,:,:,2))*(-worka(:,:,:,1)+worka(:,:,:,2)))
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
    ! Step 6: form total energy (without c.m. correction)
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
!---------------------------------------------------------------------------  
! DESCRIPTION: cm_correction
!> @brief
!!This subroutine computes the center-of-mass energy from the 
!!variance of the momentum operator. 
!>
!> @details
!!The c.m. correction from variance of total momentum is an
!!approximation to center-of-mass projection. It reads 
!! \f$ E_\mathrm{cm}=\langle\hat{P}_\mathrm{cm}^2\rangle/(2mA) $\f.
!!The computation is a bit expensive. The correction is thus 
!!evaluated only in the few analyzing steps (modulus 'mprint'). 
!!The concept is developed and tested for static states. An 
!!extension to dynamical situations is presently "off label" use.
!!
  SUBROUTINE cm_correction()
  USE Parallel, ONLY: globalindex
    COMPLEX(db),ALLOCATABLE :: ps1(:,:,:,:)  
    INTEGER :: nst,nst2,iso1,ix,iy,iz
    REAL(db) :: accvar(2),v2a,uva,v2b
    COMPLEX(db) :: offdig
    LOGICAL,PARAMETER :: testprint=.FALSE.

    ALLOCATE(ps1(nx,ny,nz,2))
    accvar=0D0
    DO nst=1,nstloc
       iso1=isospin(globalindex(nst))
       v2a=wocc(globalindex(nst))
       uva=SQRT(v2a-v2a*v2a)

       !********************************************************************
       ! x-derivatives
       !********************************************************************
       IF(TFFT) THEN
          CALL cdervx(psi(:,:,:,:,nst),ps1)  
       ELSE
          CALL cmulx(der1x,psi(:,:,:,:,nst),ps1,0)  
       ENDIF
       accvar(iso1) = accvar(iso1) + wxyz*SUM(REAL(CONJG(ps1)*ps1,db))*v2a
       DO nst2=1,nstloc
          IF(iso1==isospin(globalindex(nst2))) THEN
             v2b=wocc(globalindex(nst2))
             offdig=wxyz*SUM(CONJG(psi(:,:,:,:,nst2))*ps1)
             accvar(iso1) = accvar(iso1) - &
               REAL(CONJG(offdig)*offdig,db)*(v2a*v2b+uva*SQRT(v2b-v2b*v2b))
          END IF
       ENDDO

       !********************************************************************
       ! y-derivatives
       !********************************************************************
       IF(TFFT) THEN
          CALL cdervy(psi(:,:,:,:,nst),ps1)  
       ELSE
          CALL cmuly(der1y,psi(:,:,:,:,nst),ps1,0)  
       ENDIF
       accvar(iso1) = accvar(iso1) + wxyz*SUM(REAL(CONJG(ps1)*ps1,db))*v2a
       DO nst2=1,nstloc
          IF(iso1==isospin(globalindex(nst2))) THEN
             v2b=wocc(globalindex(nst2))
             offdig=wxyz*SUM(CONJG(psi(:,:,:,:,nst2))*ps1)
             accvar(iso1) = accvar(iso1) - &
               REAL(CONJG(offdig)*offdig,db)*(v2a*v2b+uva*SQRT(v2b-v2b*v2b))
          END IF
       ENDDO

       !********************************************************************
       ! z-derivatives
       !********************************************************************
       IF(TFFT) THEN
          CALL cdervz(psi(:,:,:,:,nst),ps1)  
       ELSE
          CALL cmulz(der1z,psi(:,:,:,:,nst),ps1,0)  
       ENDIF
       accvar(iso1) = accvar(iso1) + wxyz*SUM(REAL(CONJG(ps1)*ps1,db))*v2a
       DO nst2=1,nstloc
          IF(iso1==isospin(globalindex(nst2))) THEN
             v2b=wocc(globalindex(nst2))
             offdig=wxyz*SUM(CONJG(psi(:,:,:,:,nst2))*ps1)
             accvar(iso1) = accvar(iso1) - &
               REAL(CONJG(offdig)*offdig,db)*(v2a*v2b+uva*SQRT(v2b-v2b*v2b))
          END IF
       ENDDO

    ENDDO
    DEALLOCATE(ps1)

    ecmcorr=(accvar(1)+accvar(2))*(f%h2m(1)+f%h2m(2))/2D0/mass_number
    IF(testprint) WRITE(*,*) ' final check c.m. correction: <P**2>,Ecm=',&
            accvar,ecmcorr,mass_number

    RETURN
  END SUBROUTINE cm_correction

END MODULE Energies

