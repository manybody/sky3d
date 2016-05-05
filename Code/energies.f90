MODULE Energies
  USE Params, ONLY: db,tcoul
  USE Forces
  USE Densities
  USE Levels
  USE Pairs, ONLY: epair
  IMPLICIT NONE
  SAVE
  ! total energies calculated in subroutine "energy" and "hfenergy"
  REAL(db) :: ehft      ! kinetic energy
  REAL(db) :: ehf0      ! t0 contribution
  REAL(db) :: ehf1      ! b1 contribution (current part)
  REAL(db) :: ehf2      ! b2 contribution (Laplacian part)
  REAL(db) :: ehf3      ! t3 contribution
  REAL(db) :: ehfls     ! spin-orbit contribution (time even)
  REAL(db) :: ehflsodd  ! spin-orbit contribution (odd-odd)
  REAL(db) :: ehfc      ! Coulomb contribution
  REAL(db) :: ecorc     ! Slater & Koopman exchange
  REAL(db) :: ehfint    ! integrated total energy
  REAL(db) :: efluct1   ! energyfluctuation h**2
  REAL(db) :: efluct2   ! fluctuation h*efluct
  REAL(db) :: tke       ! kinetic energy summed
  REAL(db) :: ehf       ! Hartree-Fock energy from s.p. levels
  REAL(db) :: e3corr    ! rearrangement energy
  REAL(db) :: orbital(3),spin(3),total_angmom(3)
CONTAINS
  !***************************************************************************
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
  !***************************************************************************
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

