MODULE Energies
  USE Params, ONLY: db,tcoul,pi,e2
  USE Forces
  USE Densities
  USE Levels
  USE Pairs, ONLY: epair
  IMPLICIT NONE
  SAVE
  ! total energies calculated in subroutine "energy" and "hfenergy"
  REAL(db) :: ehfs0	! S^2 part of t0 contribution
  REAL(db) :: ehf12	! t1 & t2 contribution (Laplacian part)
  REAL(db) :: ehf22	! t1 & t2 contribution (current part)
  REAL(db) :: ehfcurr	! t1 & t2 j^2 contribution
  REAL(db) :: ehf32	! t1 & t2 contribution (invariant part - S.T-J^2)
  REAL(db) :: ehf52	! t1 & t2 contribution (invariant part - S.F-J^2)
  REAL(db) :: ehfst	! t1 & t2 contribution (time odd part - S.T)
  REAL(db) :: ehfsf	! S.F contribution (time odd)
  REAL(db) :: ehfds	! divS contribution (time odd)
  REAL(db) :: ehf42	! t1 & t2 contribution (time odd part - Lap S)
  REAL(db) :: ehf3	! t3 contribution
  REAL(db) :: ehfs3	! S^2 part of t3 contribution
  REAL(db) :: ehfls	! spin-orbit contribution
  REAL(db) :: ehflsodd	! time-odd spin-orbit contribution
  REAL(db) :: ehfj	! J_ij^2 contribution - (Lesinski sum_ij J_ijJ_ij)
  REAL(db) :: ehfjuu	! J_ii^2 tensor contribution - (Lesinski sum_i J_iiJ_ii)
  REAL(db) :: ehfjvu	! J_ji^2 tensor contribution - (Lesinski sum_i J_ijJ_ji)
  REAL(db) :: ehfjvu1	! J_ji^2 tensor contribution - (Lesinski sum_i J_ijJ_ji)
  REAL(db) :: ehfjvu2	! J_ji^2 tensor contribution - (Lesinski sum_i J_ijJ_ji)
  REAL(db) :: ehfj0	! J_0^2 contribution - (Perlinska (sum_i J_ii)^2)
  REAL(db) :: ehfj0f	! J_0^2 tensor contribution - (Perlinska (sum_i J_ii)^2)
  REAL(db) :: ehfj1f	! J_1^2 tensor contribution - (Perlinska (sum_ij E_kij J_ij)^2)
  REAL(db) :: ehfj1	! J_1^2 contribution - (Perlinska (sum_ij E_kij J_ij)^2)
  REAL(db) :: ehfj2	! J_2^2 contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj2f	! J_2^2 tensor contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj20	! J_2^2 contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj20f	! J_2^2 tensor contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj21	! J_2^2 contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj21f	! J_2^2 tensor contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj22a	! J_2^2 contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj22af	! J_2^2 tensor contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj22b	! J_2^2 contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfj22bf	! J_2^2 tensor contribution - (Perlinska sum_ij J_ij(2)J_ij(2))
  REAL(db) :: ehfjp	! J^2 contribution - total Perlinska contribution
  REAL(db) :: ehfjl	! J^2 contribution - total Lesinski contribution
  REAL(db) :: ehfint	! integrated total energy
  REAL(db) :: tkerr	! error in kinetic energy
  REAL(db) :: ecorr	! rearrangement energy
  REAL(db) :: ecorrs	! rearrangement energy from time-odd
  REAL(db) :: efluct1   ! energyfluctuation h**2
  REAL(db) :: efluct1prev  
  REAL(db) :: efluct2   ! fluctuation h*efluct
  REAL(db) :: efluct2prev
  REAL(db) :: ehft      ! kinetic energy
  REAL(db) :: ehf0      ! t0 contribution
  REAL(db) :: ehfc      ! Coulomb contribution
  REAL(db) :: ecorc     ! Slater & Koopman exchange
  REAL(db) :: tke       ! kinetic energy summed
  REAL(db) :: ehf       ! Hartree-Fock energy from s.p. levels
  REAL(db) :: ehfprev
  REAL(db) :: orbital(3),spin(3),total_angmom(3)
  REAL(db) :: e_extern=0D0  ! energy from external field (accumulator)
CONTAINS
!***************************************************************************
  SUBROUTINE integ_energy
    !
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    USE Grids, ONLY: wxyz,der1x,der2x,der1y,der2y,der1z,der2z
    USE Coulomb, ONLY: wcoul
    IMPLICIT NONE
    !
    INTEGER :: ix,iy,iz,iq,icomp,p
    REAL(db) :: rhot,rhon,rhop,sdenst,sdensp,sdensn,d2rho,d2rhon,d2rhop, &
         t0a,t0c,t3a,t3b,t1a,t1b,tj1,tj2,turnx, &
         scurrent2,scurrentp,scurrentn,j1t,j1n,j1p,t0ja, &
         t0jb,t0jaf,t0jbf,jot,jon,jop,t1ja,t1jb,j21t,j21n,j21p, &
         j122n,j122p,t0d,t0f,sdent,t1jaf,t1jbf,t2jaf,t2jbf, &
         sdenn,sdenp,t3d,t3f,tsta,tstb,tlsa,tlsb,stden,stdenn, &
         stdenp,t2ja,t2jb,tjuua,tjuub,tlsc,tlsd,t20jaf,t20jbf, &
         t20ja,t20jb,tstat,tstbt,sc,b3,b3p,b2,b2p
    INTEGER :: il
    REAL(db) :: weight
    !***********************************************************************
    !calculates the hf energy by integrating the energy functional         *
    !***********************************************************************
    REAL(db) :: worka(nx,ny,nz,2)
    REAL(db) :: workb(nx,ny,nz,3,2)
    REAL(db) :: workc(nx,ny,nz,3)
    b3=f%t3*(1.D0+0.5D0*f%x3)/4.D0
    b3p=f%t3*(0.5D0+f%x3)/4.D0 
    b2=f%t3*f%x3/8.D0
    b2p=f%t3/8.D0
    !***********************************************************************
    !calculate the correction term to the hf energy due to the 3-body force*
    !***********************************************************************
    ecorr=0.0d0
    ecorrs=0.0d0
    ecorr=-sum(f%power*wxyz*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2))**f%power*&
              (b3*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2))**2 &
               -b3p*(rho(ix,iy,iz,2)**2+ rho(ix,iy,iz,1)**2))/6.D0)

    IF(s2on) THEN
      ecorrs=-sum(f%power*wxyz*(rho(:,:,:,1)+rho(:,:,:,2))**f%power*(b2*(sdens(:,:,:,:,1)+sdens(:,:,:,:,2))**2 &
              -b2p*(sdens(:,:,:,:,1)**2+sdens(:,:,:,:,2)**2))/6.D0)
    ENDIF
    !***********************************************************************
    !        t0 and t3 terms                                               *
    !***********************************************************************
    t0a = f%t0/2.0d0 *(1.0d0+0.5d0*f%x0)
    t0c = f%t0/2.0d0 *(0.5d0+f%x0)
    t3a = f%t3/12.0d0 *(1.0d0+0.5d0*f%x3)
    t3b = f%t3/12.0d0 *(0.5d0+f%x3)
    ehf0 = 0.0d0
    ehf3 = 0.0d0
    DO iz = 1,nz
       DO iy = 1,ny
          DO ix = 1,nx
             rhot = rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
             rhon = rho(ix,iy,iz,1)
             rhop = rho(ix,iy,iz,2)
             ehf0 = ehf0+wxyz *(t0a*rhot**2-t0c *(rhop**2+rhon**2))
             ehf3 = ehf3+wxyz*rhot**f%power *(t3a*rhot**2-t3b* &
                  (rhop**2+rhon**2))
          ENDDO
       ENDDO
    ENDDO
    !***********************************************************************
    !        time-odd t0 and t3 terms S^2                                  *
    !***********************************************************************
    ehfs0 = 0.0d0
    ehfs3 = 0.0d0
    t0d = (f%t0*f%x0)/4.0d0
    t0f = -f%t0/4.0d0
    t3d = (f%t3*f%x3)/24.0d0
    t3f = -f%t3/24.0d0
    !
    IF (s2on) THEN
      DO iz = 1,nz
         DO iy = 1,ny
            DO ix = 1,nx
               rhot = rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
               sdent =  (sdens(ix,iy,iz,1,1) + sdens(ix,iy,iz,1,2))**2&
                      + (sdens(ix,iy,iz,2,1) + sdens(ix,iy,iz,2,2))**2&
                      + (sdens(ix,iy,iz,3,1) + sdens(ix,iy,iz,3,2))**2
               sdenn = sdens(ix,iy,iz,1,1)**2+sdens(ix,iy,iz,2,1)**2+&
                       sdens(ix,iy,iz,3,1)**2
               sdenp = sdens(ix,iy,iz,1,2)**2+sdens(ix,iy,iz,2,2)**2+&
                       sdens(ix,iy,iz,3,2)**2
               ehfs0 = ehfs0 + wxyz*(t0d*sdent+t0f*(sdenp+sdenn))
               ehfs3 = ehfs3 + wxyz*rhot**f%power*(t3d*sdent+t3f*(sdenp+sdenn))
            ENDDO
         ENDDO
      ENDDO
    ENDIF
    ehf0 = ehf0 + ehfs0
    ehf3 = ehf3 + ehfs3
    !***********************************************************************
    !        terms depending on the laplacian of density(t1,t2)            *
    !***********************************************************************
    ehf12 = 0.0d0
    DO iq=1,2
       CALL rmulx(der2x,rho(:,:,:,iq),worka(:,:,:,iq),0)
       CALL rmuly(der2y,rho(:,:,:,iq),worka(:,:,:,iq),1)
       CALL rmulz(der2z,rho(:,:,:,iq),worka(:,:,:,iq),1)
    ENDDO
    t1a =(f%t2+0.5d0*f%t2*f%x2-3.0d0*f%t1-1.5d0*f%t1*f%x1)/16.0d0
    t1b =(1.5d0*f%t1+3.0d0*f%t1*f%x1+0.5d0*f%t2+f%t2*f%x2)/16.0d0
    !
    DO iz = 1,nz
       DO iy = 1,ny
          DO ix = 1,nx
             rhon = rho(ix,iy,iz,1)
             rhop = rho(ix,iy,iz,2)
             d2rhon = worka(ix,iy,iz,1)
             d2rhop = worka(ix,iy,iz,2)
             ehf12 = ehf12+wxyz *(t1a*(rhon+rhop)*(d2rhon+d2rhop)+t1b *(rhop* &
                  d2rhop+rhon*d2rhon))
          ENDDO
       ENDDO
    ENDDO
    !***********************************************************************
    !        square of the current(dot product of the current vector)      *
    !        calculate the other t1,t2 term(rho*tau-j**2)                  *
    !***********************************************************************
    worka=current(:,:,:,1,:)**2+current(:,:,:,2,:)**2+current(:,:,:,3,:)**2
    workc=(current(:,:,:,:,1)+current(:,:,:,:,2))**2
    t1a =(f%t1+0.5d0*f%t1*f%x1+f%t2+0.5d0*f%t2*f%x2)/4.0d0
    t1b =(0.5d0*f%t2+f%t2*f%x2-0.5d0*f%t1-f%t1*f%x1)/4.0d0
    ehf22=wxyz*SUM(t1a*((rho(:,:,:,1)+rho(:,:,:,2))*(tau(:,:,:,1)+tau(:,:,:,2)) &
         -workc(:,:,:,1)-workc(:,:,:,2)-workc(:,:,:,3)) &
         +t1b*(rho(:,:,:,1)*tau(:,:,:,1)+rho(:,:,:,2)*tau(:,:,:,2) &
         -worka(:,:,:,1)-worka(:,:,:,2)))
    ! Work out j^2 component separately
    ehfcurr=wxyz*SUM(-t1a*(workc(:,:,:,1)+workc(:,:,:,2)+workc(:,:,:,3)) &
                -t1b*(worka(:,:,:,1)+worka(:,:,:,2)))
    !***********************************************************************
    !        S.T  term                                                     *
    !***********************************************************************
    ehf32=0.0d0
    ehfst=0.0d0
    tsta = (f%t1*f%x1+f%t2*f%x2)/8.0d0
    tstb = (f%t2-f%t1)/8.0d0
    tstat = -(te+toten)/4.0d0
    tstbt = (te-toten)/4.0d0
    IF(ston) THEN
      ehfst=wxyz*SUM((tsta+tstat)*(sdens(:,:,:,:,1)+sdens(:,:,:,:,2))* &
              (tdens(:,:,:,:,1)+tdens(:,:,:,:,2)) &
              +(tstb+tstbt)*(sdens(:,:,:,:,1)*tdens(:,:,:,:,1)+ &
              sdens(:,:,:,:,2)*tdens(:,:,:,:,2)))
    ENDIF
    !***********************************************************************
    !        J_ij^2 terms (Lesinski notation)                              *
    !***********************************************************************
    ehfj=0.0d0
    ehfjl=0.0d0
    ehfjp=0.0d0
    ehfj21=0.0d0
    ehfj21f=0.0d0
    tj1 =-(f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/8.0d0
    tj2 =(f%t1-f%t2+2.0d0*(toten-te))/8.0d0
    t2ja = -(f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/8.0d0
    t2jb = (f%t1-f%t2-2.0d0*(te-toten))/8.0d0
    t2jaf = -3.0d0*(te+toten)/8.0d0
    t2jbf = 3.0d0*(te-toten)/8.0d0
    scurrent2=SUM((scurrentx(:,:,:,:,1)+scurrentx(:,:,:,:,2))**2 &
           +(scurrenty(:,:,:,:,1)+scurrenty(:,:,:,:,2))**2 &
           +(scurrentz(:,:,:,:,1)+scurrentz(:,:,:,:,2))**2)
    scurrentn=SUM(scurrentx(:,:,:,:,1)**2+scurrenty(:,:,:,:,1)**2 &
           +scurrentz(:,:,:,:,1)**2)
    scurrentp=SUM(scurrentx(:,:,:,:,2)**2+scurrenty(:,:,:,:,2)**2 &
           +scurrentz(:,:,:,:,2)**2)
    IF(j2on) THEN
      ehfj = wxyz*(tj1*scurrent2 + tj2*(scurrentn+scurrentp))
      ehfj21 = wxyz*(0.5d0*t2ja*scurrent2 + 0.5d0*t2jb*(scurrentn+scurrentp))
    ENDIF
    IF(jfon) THEN
      ehfj21f = wxyz*(0.5d0*t2jaf*scurrent2 + 0.5d0*t2jbf*(scurrentn+scurrentp))
    ENDIF
    !***********************************************************************
    !        J_0^2 terms                                                   *
    !***********************************************************************
    ehfj0=0.0d0
    ehfj0f=0.0d0
    ehfj2=0.0d0
    ehfj2f=0.0d0
    ehfjuu=0.0d0
    ehfj20=0.0d0
    ehfj20f=0.0d0
    t0ja = -(f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/24.0d0
    t0jb = (f%t1-f%t2-2.0d0*(te-toten))/24.0d0
    t0jaf = -(te+toten)/2.0d0
    t0jbf = (te-toten)/2.0d0
    tjuua = -3.0d0*(te+toten)/8.0d0
    tjuub = 3.0d0*(te-toten)/8.0d0
    t20ja = (f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/24.0d0
    t20jb = -(f%t1-f%t2-2.0d0*(te-toten))/24.0d0
    t20jaf = 3.0d0*(te+toten)/24.0d0
    t20jbf = -3.0d0*(te-toten)/24.0d0
    jot=SUM(((scurrentx(:,:,:,1,1)+scurrentx(:,:,:,1,2)) &
         +(scurrenty(:,:,:,2,1)+scurrenty(:,:,:,2,2)) &
         +(scurrentz(:,:,:,3,1)+scurrentz(:,:,:,3,2)))**2)
    jon=SUM((scurrentx(:,:,:,1,1)+scurrenty(:,:,:,2,1) &
         +scurrentz(:,:,:,3,1))**2)
    jop=SUM((scurrentx(:,:,:,1,2)+scurrenty(:,:,:,2,2) &
         +scurrentz(:,:,:,3,2))**2)
    IF(j2on) THEN
      ehfj0 = wxyz*(t0ja*jot + t0jb*(jon+jop))
      ehfj20 = wxyz*(t20ja*jot + t20jb*(jon+jop))
    ENDIF
    IF(jfon) THEN
      ehfjuu = wxyz*(tjuua*jot + tjuub*(jon+jop))
      ehfj0f = wxyz*(t0jaf*jot + t0jbf*(jon+jop))
      ehfj20f = wxyz*(t20jaf*jot + t20jbf*(jon+jop))
    ENDIF
    !***********************************************************************
    !        J_1^2 terms                                                   *
    !***********************************************************************
    ehfj1=0.0d0
    t1ja = -(f%t1*f%x1+f%t2*f%x2-2.0d0*(te+toten))/16.0d0
    t1jb = (f%t1-f%t2+2.0d0*(toten-te))/16.0d0
    t1jaf = 3.0d0*(te+toten)/16.0d0
    t1jbf = 3.0d0*(toten-te)/16.0d0
    j1t=SUM((sodens(:,:,:,1,1)+sodens(:,:,:,1,2))**2 &
         +(sodens(:,:,:,2,1)+sodens(:,:,:,2,2))**2 &
         +(sodens(:,:,:,3,1)+sodens(:,:,:,3,2))**2)
    j1n=SUM(sodens(:,:,:,1,1)**2+sodens(:,:,:,2,1)**2 &
         +sodens(:,:,:,3,1)**2)
    j1p=SUM(sodens(:,:,:,1,2)**2+sodens(:,:,:,2,2)**2 &
         +sodens(:,:,:,3,2)**2)
    IF(j2on) THEN
      ehfj1 = wxyz*(t1ja*j1t + t1jb*(j1n+j1p))
    ENDIF
    IF(jfon) THEN
      ehfj1f = wxyz*(t1jaf*j1t + t1jbf*(j1n+j1p))
    ENDIF
    !***********************************************************************
    !        J_2^2 terms                                                   *
    !***********************************************************************
    ehfjvu1=0.0d0
    ehfjvu2=0.0d0
    ehfjvu=0.0d0
    ehfj22a=0.0d0
    ehfj22b=0.0d0
    j21t=SUM((scurrentx(:,:,:,1,1)+scurrentx(:,:,:,1,2))**2 &
         +(scurrenty(:,:,:,2,1)+scurrenty(:,:,:,2,2))**2 &
         +(scurrentz(:,:,:,3,1)+scurrentz(:,:,:,3,2))**2)
    j21n=SUM(scurrentx(:,:,:,1,1)**2+scurrenty(:,:,:,2,1)**2 &
         +scurrentz(:,:,:,3,1))**2
    j21p=SUM(scurrentx(:,:,:,1,2)**2+scurrenty(:,:,:,2,2)**2 &
         +scurrentz(:,:,:,3,2))**2
    IF(j2on) THEN
      ehfj22a = wxyz*(0.5d0*t2ja*j21t + 0.5d0*t2jb*(j21n+j21p))
    ENDIF
    IF(jfon) THEN
      ehfjvu1 = wxyz*(tjuua*j21t + tjuub*(j21n+j21p))
      ehfj22af = wxyz*(0.5d0*t2jaf*j21t + 0.5d0*t2jbf*(j21n+j21p))
    ENDIF
    j122n=SUM(scurrentx(:,:,:,2,1)*scurrenty(:,:,:,1,1)+ &
          scurrentx(:,:,:,3,1)*scurrentz(:,:,:,1,1)+ &
          scurrenty(:,:,:,3,1)*scurrentz(:,:,:,2,1))
    j122p=SUM(scurrentx(:,:,:,2,2)*scurrenty(:,:,:,1,2)+ &
          scurrentx(:,:,:,3,2)*scurrentz(:,:,:,1,2)+ &
          scurrenty(:,:,:,3,2)*scurrentz(:,:,:,2,2))
    IF(j2on) THEN
      ehfj22b = wxyz*(t2ja*(j122n+j122p) + t2jb*(j122n+j122p))
    ENDIF
    IF(jfon) THEN
      ehfjvu2 = wxyz*(2.0d0*tjuua*(j122n+j122p) + 2.0d0*tjuub*(j122n+j122p))
      ehfj22bf = wxyz*(t2jaf*(j122n+j122p) + t2jbf*(j122n+j122p))
    ENDIF
    ehfjvu=ehfjvu1+ehfjvu2
    ehfj2=ehfj21+ehfj22a+ehfj22b+ehfj20
    ehfj2f=ehfj21f+ehfj22af+ehfj22bf+ehfj20f
    ehfjp=ehfj0+ehfj0f+ehfj1+ehfj1f+ehfj2+ehfj2f
    ehfjl=ehfj+ehfjuu+ehfjvu
    ehf32=ehfst+ehfj
    !***********************************************************************
    !        S.F  term                                                     *
    !***********************************************************************
    ehfsf=0.0d0
    IF(sfon) THEN
      ehfsf=wxyz*SUM(3.0d0*(te+toten)/4.0d0* &
                (sdens(:,:,:,:,1)+sdens(:,:,:,:,2))* &
                (fdens(:,:,:,:,1)+fdens(:,:,:,:,2)) &
             +3.0d0*(toten-te)/4.0d0*(sdens(:,:,:,:,1)*fdens(:,:,:,:,1)+ &
              sdens(:,:,:,:,2)*fdens(:,:,:,:,2)))
    ENDIF
    ehf52=ehfsf+ehfjuu+ehfjvu
    !***********************************************************************
    !        spin-orbit contribution to hf energy                          *
    !***********************************************************************
    DO iq = 1,2
       CALL rmulx(der1x,sodens(:,:,:,1,iq),worka(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),worka(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),worka(:,:,:,iq),1)
    ENDDO
    b4 = f%t4/2.0d0

    ehfls=wxyz*SUM(-b4*(rho(:,:,:,1)+ &
         rho(:,:,:,2))*(worka(:,:,:,1)+worka(:,:,:,2))-&
         f%b4p *(rho(:,:,:,2)*worka(:,:,:,2)+rho(:,:,:,1)*worka(:,:,:,1)))

    !***********************************************************************
    !        s.(nabla x j)                                                 *
    !***********************************************************************
    DO iq = 1,2
        CALL rmuly(der1y,current(:,:,:,3,iq),workb(:,:,:,1,iq),0)
        CALL rmulz(der1z,current(:,:,:,2,iq),workb(:,:,:,1,iq),-1)
        CALL rmulz(der1z,current(:,:,:,1,iq),workb(:,:,:,2,iq),0)
        CALL rmulx(der1x,current(:,:,:,3,iq),workb(:,:,:,2,iq),-1)
        CALL rmulx(der1x,current(:,:,:,2,iq),workb(:,:,:,3,iq),0)
        CALL rmuly(der1y,current(:,:,:,1,iq),workb(:,:,:,3,iq),-1)
    ENDDO
    ehflsodd=wxyz*SUM(-b4*(sdens(:,:,:,:,1)+sdens(:,:,:,:,2))* &
             (workb(:,:,:,:,1)+workb(:,:,:,:,2))-f%b4p*(sdens(:,:,:,:,1)* &
             workb(:,:,:,:,1)+sdens(:,:,:,:,2)*workb(:,:,:,:,2)))

    ehfls=ehfls+ehflsodd

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
            ehfc= wyxz*sum((0.5D0*rho(:,:,:,2)*wcoul(:,:,:) &
                 +sc*rho(:,:,:,2)**(4.0D0/3.0D0)))
            ecorc= wxyz*sum(sc/3.0D0*(rho(:,:,:,2)**(4.0D0/3.0D0)))
    ENDIF
    !***********************************************************************
    !        kinetic energy contribution                                   *
    !***********************************************************************
    ehft=wxyz*SUM(f%h2m(1)*tau(:,:,:,1)+f%h2m(2)*tau(:,:,:,2))
    ehfint = ehft+ehf0+ehf12+ehf22+ehf32+ehf42+ehf52+ehfds+ehf3+ehfls+ehfc-&
             epair(1)-epair(2)
  END SUBROUTINE integ_energy
!***************************************************************************
  SUBROUTINE sum_energy
    USE Moment, ONLY: pnrtot
    INTEGER :: i
    ehf = 0.0d0
    tke = 0.0d0
    efluct1=0.0d0
    efluct2=0.0d0
    ehf=SUM(wocc*(sp_kinetic+sp_energy))/2.D0 +ecorr+ecorc +ecorrs - &
        epair(1) - epair(2)
    tke=SUM(wocc*sp_kinetic)
    efluct1=SUM(wocc*sp_efluct1)/pnrtot
    efluct2=SUM(wocc*sp_efluct2)/pnrtot
    DO i=1,3
      orbital(i)=SUM(wocc*sp_orbital(i,:))
      spin(i)=SUM(wocc*sp_spin(i,:))
      total_angmom(i)=orbital(i)+spin(i)
    END DO
    tkerr = ehft-tke
  END SUBROUTINE sum_energy
!***************************************************************************
END MODULE Energies
