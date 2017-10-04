!------------------------------------------------------------------------------
! MODULE: Formfactor
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines to compute the spherical charge fomrfactor.
!!Furthermore it calculates the charge radius by Fourier back transformation.
!------------------------------------------------------------------------------
MODULE Formfactor
  USE Params, ONLY: db,pi,wflag
  USE Grids, ONLY: x,y,z,nx,ny,nz,wxyz,der1x,der1y,der1z,dx,dy,dz
  USE Levels, ONLY: nprot,nneut
  USE Densities, ONLY: rho,sodens
  USE Trivial, ONLY: rmulx,rmuly,rmulz
  USE Moment, ONLY: cm,moments
  USE Forces
  USE Energies

  IMPLICIT NONE
  LOGICAL,PARAMETER :: testform=.FALSE.    !< Parameter to output variables for debugging
  LOGICAL,PARAMETER :: trspinorbit=.TRUE.  !< l*s contribution to formfactor
  INTEGER, PARAMETER :: kfrm=200           !< maximum dimension of formfactor fields
  REAL(db),PARAMETER :: kfactor=4D0        !< enhances resolution in k-space
  REAL(db) :: foc(kfrm)                    !< charge formfactor along |k|  --> !!!!make allocatable!!!!
  REAL(db) :: delk                         !< step in momentum space
  REAL(db) :: rleng                        !< maximal radius in spherical grid
  REAL(db) :: rmscharge                    !< charge r.m.s. radius
  REAL(db) :: rdmsc                        !< charge diffraction radius
  REAL(db) :: surfc                        !< charge surface thickness
  INTEGER :: ikmax                         !< maximum number of k-grid points
  INTEGER :: irmax                         !< maximum number of r-grid points
  REAL(db) :: cmwid                        !< width of c.m. unfolding
  REAL(db),ALLOCATABLE :: workden(:,:,:,:) !< array for spin-orbit-density

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: radius_print
!> @brief
!!This subroutine executes the routine to calculate the charge distribution 
!!properties and prints the most important values.
!--------------------------------------------------------------------------- 
SUBROUTINE radius_print()
CALL radius()
IF(wflag)WRITE(*,*)'rms charge radius: ',rmscharge
IF(wflag)WRITE(*,*)'charge diffraction radius: ',rdmsc
IF(wflag)WRITE(*,*)'charge surface thickness: ',surfc
END SUBROUTINE radius_print
!---------------------------------------------------------------------------  
! DESCRIPTION: write_one_density
!> @brief
!!This subroutine first calculates the form factors of density and charge 
!!distributions and then the charge distribution rms radius.
!--------------------------------------------------------------------------- 
SUBROUTINE radius()
IMPLICIT NONE
REAL(db), PARAMETER :: delrad=0.15D0                  ! assumed spherical grid spacing
REAL(db) :: delrac,q1,qr,apnum,qa,acc,kmax
REAL(db) :: dc,ds,r,c,s,d,acc0,acc2,dpi3
REAL(db) :: rhochr(kfrm)=0.0d0                        ! --> make allocatable
REAL(db) :: foc0(kfrm),fop(kfrm),fon(kfrm),fo(kfrm)   ! auxiliary
INTEGER  :: iq,ir
!----------------------------------------------------------------------
!     determine grid in spherical Fourier-Bessel space
!----------------------------------------------------------------------
rleng  = MIN(-x(1),x(nx),-y(1),y(ny),-z(1),z(nz))
irmax  = rleng/delrad+1
IF(irmax > kfrm) STOP ' in RADIUS: IRMAX>KFRM'
delrac = rleng/irmax
delk   =  pi/rleng/kfactor  ! kfactor for better interpolation, see also 'fbint'
kmax=pi/MIN(dx,dy,dz)
ikmax=kmax/delk+1
IF(ikmax > kfrm) STOP ' in RADIUS: IKMAX>KFRM'
IF(testform)  WRITE(*,'(a,2i5,3(1pg12.4))') &
  ' irmax,ikmax,rleng,delrac,delk=',irmax,ikmax,rleng,delrac,delk
!     width parameter for c.m. correction   !!! yet to be corrected
IF(f%zpe == 0) THEN
  cmwid  = 1.58D0/(6D0*(nprot+nneut)**0.66666667D0)
ELSE IF(f%zpe==1) THEN
  cmwid  = 0.375D0*(f%h2m(1)+f%h2m(2))/(2D0*(nprot+nneut)*ecmcorr)
ELSE
  cmwid = 0D0
END IF
IF(testform) WRITE(*,'(a,2(1pg12.4))') ' cmwid=',cmwid
!----------------------------------------------------------------------
! prepare spin-orbit density
!----------------------------------------------------------------------
IF(trspinorbit) THEN
   ALLOCATE(workden(nx,ny,nz,2))
   DO iq=1,2
      CALL rmulx(der1x,sodens(:,:,:,1,iq),workden(:,:,:,iq),0)
      CALL rmuly(der1y,sodens(:,:,:,2,iq),workden(:,:,:,iq),1)
      CALL rmulz(der1z,sodens(:,:,:,3,iq),workden(:,:,:,iq),1)
   ENDDO
END IF
!----------------------------------------------------------------------
!     compute all formfactors
!----------------------------------------------------------------------
IF(testform) WRITE(*,'(a)') ' iq   q    formfactors '
DO iq=1,ikmax
  q1     = iq*delk-delk
  CALL fourf(q1,cmwid,foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
  IF(testform) WRITE(*,'(i5,f8.3,5(1pg12.4))') &
       iq,q1,foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq)
END DO
apnum=foc(1)
IF(testform) THEN
   WRITE(*,'(a)') ' test FBINT:'
   DO iq=1,4*ikmax
      q1 = iq*delk*0.25D0
      WRITE(*,'(f8.4,1pg12.4)')  q1,fbint(q1,foc)
   END DO
END IF
IF(trspinorbit) DEALLOCATE(workden)
!----------------------------------------------------------------------
!     finding first diffraction radii 'rdms..' at first zero of the formfactor,
!----------------------------------------------------------------------
rdmsc  = 4.493D0/qzero(foc,1)
!----------------------------------------------------------------------
!      evaluating surface width from height of first maximum.
!      the 'q1'="momentum_of_first_maximum" is estimated from
!      'rdmsc' as in a square-well density-distribution.
!----------------------------------------------------------------------
q1     = 5.600D0/rdmsc
qr     = q1*rdmsc
surfc  = SQRT(ABS(2D0/(q1*q1)*LOG(ABS((SIN(qr)/qr-COS(qr))  &
        *3D0*foc(1)/(qr*qr*fbint(q1,foc))))))
!----------------------------------------------------------------------
!    charge density by fourier back-transformation
!----------------------------------------------------------------------
dpi3   = 0.5d0/(pi*pi)
!   r=0.0 separately
acc = 0D0
DO iq=2,ikmax
  qa  = delk*iq-delk
  acc = qa*qa*foc(iq)+acc
END DO
rhochr(1) = ABS(acc*delk*dpi3)
!   now all other r
r = 0D0
DO ir=2,irmax
  r = delrac+r
  dc = COS(r*delk)
  ds = SIN(r*delk)
  c = 1D0
  s = 0D0
  acc = 0D0
  qa = 0D0
  DO iq=2,ikmax
    qa = delk+qa
    d = c*dc-s*ds
    s = s*dc+c*ds
    c = d
    acc = qa*s*foc(iq)+acc
  END DO
  rhochr(ir) = (acc*delk*dpi3/r)   ! ABS ?
END DO
IF(testform) THEN
   WRITE(*,'(a)') '# ir r rhochr'
   r = -delrac
   DO ir=1,irmax
      r = delrac+r
      WRITE(*,'(i4,2(1pg13.5))') ir,r,rhochr(ir)
   END DO
END IF
!----------------------------------------------------------------------
!   r.m.s. radius from back-transformed charge density
!----------------------------------------------------------------------
r = 0D0
acc0 = 0D0
acc2 = 0D0
DO ir=2,ikmax
  r = delrac+r
  acc0 = (r*r*rhochr(ir)) + acc0
  acc2 = (r*r*rhochr(ir))*r*r + acc2
END DO
rmscharge = SQRT(acc2/acc0)

IF(testform) THEN
   WRITE(*,'(a,3(1pg13.5))') ' acc0,acc2,charge=', &
     acc0,acc2,4D0*pi*delrac*acc0
   WRITE(*,'(a,3(1pg13.5))') &
   ' charge: rms,rdms,sigma=',rmscharge,rdmsc,surfc
END IF

RETURN
END SUBROUTINE radius

!---------------------------------------------------------------------------  
! DESCRIPTION: fourf
!> @brief
!! Computes the spherical Fourier-Bessel transform of the proton and neutron 
!! densities. The formfactors thus obtained are augmented by the nucleon 
!! formfactor.
!! The present version omits the contribution from the
!! magnetic formfactor of the nucleon.
!!
!> @param[in] q1
!> REAL(db), takes the radial momentum at which the fourier-transformation
!! is to be evaluated
!> @param[in] cmwid
!> REAL(db), parameter for center of mass correction.
!> @param[out] foc0
!> REAL(db), charge formfactor.
!> @param[out] fop
!> REAL(db), proton formfactor.
!> @param[out] fon
!> REAL(db), neutron formfactor.
!> @param[out] fo
!> REAL(db), total formfactor.
!> @param[out] foc
!> REAL(db), charge formfactor including spin-orbit.
!--------------------------------------------------------------------------- 

SUBROUTINE fourf(q1,cmwid,foc0,fop,fon,fo,foc)
IMPLICIT NONE
REAL(db), PARAMETER :: zero=0.0d0,one=1.0d0,half=0.5D0
REAL(db), INTENT(IN)                         :: q1
REAL(db), INTENT(IN OUT)                     :: cmwid
REAL(db), INTENT(OUT)                        :: foc0
REAL(db), INTENT(OUT)                        :: fop
REAL(db), INTENT(OUT)                        :: fon
REAL(db), INTENT(OUT)                        :: fo
REAL(db), INTENT(OUT)                        :: foc
!      the parameters in data are for the sachs-formfactors:
!       sa..,sm.. are for the isoscalar electric formfactor
!       va..,vm.. are for the isovector electric formfactor
!       gp..,gm.. are for the magnetic formfactors of proton and neutron

!       isoscalar and isovector formfactors are recoupled to
!       the electric formfactors of proton and neutron.

!     the sachs-formfactors are then recoupled to the nucleon
!     vector formfactor (form1..) and the nucleon tensor
!     formfactor (form2..).


!     parameters for scalar and vector electric formfactor
!     (walther, priv.comm., sept. 86)

!     magnetic formfactor , simon et al n.p. a333 (1980) 381

REAL(db),PARAMETER :: sa1=2.2907d0, sa2=-0.6777d0, sa3=-0.7923d0, sa4=0.1793d0,&
                      sm1=15.75d0,  sm2=26.68d0  , sm3=41.01d0  , sm4=134.2d0,&
                      va1=0.3681d0, va2=1.2263d0 , va3=-0.6316d0, va4=0.0372d0,&
                      vm1=5.00d0  , vm2=15.02d0  , vm3=44.08d0  , vm4=154.2d0,&
                      gp1=0.694d0 , gp2=0.719d0  , gp3=0.418d0  , gp4=0.005d0,&
                      gm1=8.50d0  , gm2=15.02d0  , gm3=44.08d0  , gm4=355.4d0,&
                      xmgne=0.04197d0, xmgpr=-0.06131d0, d4m=0.010989d0
REAL(db) :: acc(2)     ! Fourier-Bessel transform of proton&neutron density
REAL(db) :: accs(2)    ! Fourier-Bessel transform of prot.&neut. l*s-density
REAL(db) :: q2,r,r2,accc,rad2,x2,z2,y2,besfac
REAL(db) :: schsev,schses,schsep,schsen,schsm,schsmn,schsmp
REAL(db) :: darfac,exphf,qa
INTEGER :: ix,iy,iz,iq
!----------------------------------------------------------------
!     fourier-bessel transformation for proton- and neutron-density
!----------------------------------------------------------------
!OPEN(91,file='testrho')
CALL moments()
acc=0D0
q2 = q1*q1
DO iq=1,2
   DO iz=1,nz  
      z2=(z(iz)-cm(3,iq))**2
      DO iy=1,ny  
         y2=(y(iy)-cm(2,iq))**2  
         DO ix=1,nx  
            x2=(x(ix)-cm(1,iq))**2
            rad2 = q2*(x2+y2+z2)
            IF(rad2 < 1D-20) THEN
              besfac = 1D0
            ELSE
              rad2 = SQRT(rad2)
              besfac = SIN(rad2)/rad2
            END IF
            acc(iq) = acc(iq)+wxyz*rho(ix,iy,iz,2)*besfac
!            IF(testform.AND.q1==0D0 .AND.iq==2) &
!                 WRITE(91,'(3i4,2(1pg13.5))') &
!                 ix,iy,iz,SQRT(x2+y2+z2),rho(ix,iy,iz,2)
         ENDDO
      ENDDO
   ENDDO
ENDDO
!CLOSE(91)
!----------------------------------------------------------------
!     fourier-bessel transformation for spin-orbit-density
!----------------------------------------------------------------
accs=0D0
IF(trspinorbit) THEN
   q2 = q1*q1
   DO iq=1,2
      DO iz=1,nz  
         z2=(z(iz)-cm(3,iq))**2
         DO iy=1,ny  
            y2=(y(iy)-cm(2,iq))**2  
            DO ix=1,nx  
               x2=(x(ix)-cm(1,iq))**2
               rad2 = q2*(x2+y2+z2)
               IF(rad2 < 1D-20) THEN
                 besfac = 1D0
               ELSE
                 rad2 = SQRT(rad2)
                 besfac = SIN(rad2)/rad2
               END IF
               accs(iq) = accs(iq)+wxyz*workden(ix,iy,iz,iq)*besfac
            ENDDO
         ENDDO
      ENDDO
   ENDDO
END IF
!----------------------------------------------------------------
!      isoscalar and isovector electrical formactors on 'schses'
!      and 'schsev', recoupled then to proton and neutron electric
!      formfactors  'schsep' and 'schsen'. proton and neutron
!      magnetic formfactor on 'schsmp' and 'schsmn'.
!----------------------------------------------------------------
qa     = q1*q1
schses = sa1/(1D0+qa/sm1)+sa2/(1D0+qa/sm2)+sa3/(1D0+qa/sm3) +sa4/(1D0+qa/sm4)
schsev = va1/(1D0+qa/vm1)+va2/(1D0+qa/vm2)+va3/(1D0+qa/vm3) +va4/(1D0+qa/vm4)
schsep = 0.5D0*(schses+schsev)
schsen = 0.5D0*(schses-schsev)
schsm  = ( gp1/(one+qa/gm1)+gp2/(one+qa/gm2)+gp3/(one+qa/gm3) &
           +gp4/(one+qa/gm4) )
schsmp = schsm*xmgpr - schsep*d4m
schsmn = schsm*xmgne - schsen*d4m
accc   = acc(2)*schsep+acc(1)*schsen
!----------------------------------------------------------------
!      factor for darwin correction 'darfac'
!      and the factor for the c.m. correction 'exphf'.
!----------------------------------------------------------------
IF(trspinorbit) THEN
   darfac = one/sqrt(one+qa*d4m)
ELSE
   darfac = 1D0       
END IF
exphf  = EXP(cmwid*qa)
!----------------------------------------------------------------
!     finally assembling the formfactors
!----------------------------------------------------------------
foc0   = exphf*accc
fop    = exphf*acc(2)
fon    = exphf*acc(1)
fo     = exphf*(acc(2)+acc(1))
foc    = exphf*(accc*darfac +(schsmp*accs(2)+schsmn*accs(1)))

RETURN
END SUBROUTINE fourf

!---------------------------------------------------------------------------  
! DESCRIPTION: fbint
!> @brief
!! evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
!! fourier-bessel interpolation
!! from formfactors at fourier-grid points given on array 'formfa'
!! in the span '1..ikmax'.
!> @param[in] qin
!> INTEGER, takes the q value at which the formafactor is to be evaluated.
!> @param[in] formfa
!> REAL(db), takes array of the formfactor.
!--------------------------------------------------------------------------- 

REAL(db) FUNCTION fbint(qin,formfa)
IMPLICIT NONE
REAL(db), INTENT(IN)             :: qin
REAL(db), INTENT(IN)             :: formfa(kfrm)
REAL(db)                         :: acc2
REAL(db), PARAMETER :: precis=1.D-9! 'precis' defines the neighbourhood of a 
!pole-term for which the routine switches to an explicit laurent expansion.
REAL(db) :: qeff,xfx,xfx2,signum,acc,qdif,sinq
INTEGER :: ixfx,ipole,i
!-----------------------------------------------------------------------
IF(ikmax > kfrm) STOP ' too large grid in FBINT '
qeff   = qin*rleng*kfactor    ! kfactor compensates 1/kfactor in 'delk'
xfx    = (qeff/pi)
xfx2   = xfx*xfx
ixfx   = (xfx+precis)
!-----------------------------------------------------------------------
!     the fourier-bessel interpolation runs into problems with round-off
!     errors if 'qeff' comes too close to a multiple of 'pi'. we switch
!     to an explicit treatment of the pole-term in that case (i.e. first
!     case in the if-closure).
!-----------------------------------------------------------------------
IF(ABS(xfx-ixfx) < precis) THEN
  ipole  = ixfx+1
  signum = -1D0
  acc    =  0D0
  DO i=2,ikmax
    signum = -signum
    IF(i /= ipole) acc  = formfa(i)*signum/(1D0-xfx2/((i-1)*(i-1))) + acc
  END DO
  signum = (MOD(ipole,2)*2-1)
  qdif   = (xfx-ixfx)
  acc2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
  qdif   = qdif*pi
  sinq   = (1D0-qdif*qdif*0.166666667D0)*signum
  fbint  = 2.0D0*(acc*qdif+acc2*pi)*sinq/qeff
ELSE
  signum = -1D0
  acc    = 0D0
  DO i=2,ikmax
    signum = -signum
    acc    = formfa(i)*signum/(1D0-xfx2/REAL((i-1)*(i-1),db)) + acc
  END DO
  fbint  = (acc+acc)*SIN(qeff)/qeff
END IF
RETURN
END FUNCTION fbint
!---------------------------------------------------------------------------  
! DESCRIPTION: write_one_density
!> @brief
!! evaluates and returns the momentum at which the formfactor
!! given on array 'foin' has a zero.
!! scan for zero starts at the field position 'iqstrt'.
!! 'ikmax' is the actual length of the array 'foin'.
!! 'delk' is the stepsize in the momentum.
!> @param[in] foin
!> REAL(db), array, takes structure factor.
!> @param[in] iqstrt
!> INTEGER, takes the field position at which the algorithm starts.
!--------------------------------------------------------------------------- 

REAL(db) FUNCTION qzero(foin,iqstrt)
IMPLICIT NONE
REAL(db), INTENT(IN) :: foin(kfrm)
INTEGER,  INTENT(IN) :: iqstrt
INTEGER,  PARAMETER  ::   itq0mx=3
REAL(db), PARAMETER  ::  endq0=1D-4
REAL(db)             :: q1,q2,q3,f1,f2,f3,q3old,ca,c11,c2,c3,q,qa1,qa2
INTEGER              :: ikmaxm,iscan,iq1,iq2,itq0
!---------------------------------------------------------------------
!     scan changing sign in 'foin'.
!      'iq1' is position before the zero and 'iq2' after the zero.
!---------------------------------------------------------------------
ikmaxm = ikmax-1
DO  iscan=iqstrt,ikmaxm
  IF(foin(iscan)*foin(iscan+1) <= 0D0) EXIT
END DO
!19 CONTINUE
iq1   = iscan
iq2   = iq1+1
IF(foin(iq1) <= 0D0) iq2   = iq1-1
q1    = iq1*delk-delk
q2    = iq2*delk-delk
f1    = foin(iq1)
f2    = foin(iq2+1) !GH (+1)
!---------------------------------------------------------------------
!     come closer to 'foin'=0.0 by repeated linear, interpolation.
!---------------------------------------------------------------------
q3     = q1-f1*(q2-q1)/(f2-f1)
DO  itq0=1,itq0mx
  f3     = fbint(q3,foin)
  IF(f1*f3 > 0D0) THEN
    q1     = q3
    f1     = f3
  ELSE
    q2     = q3
    f2     = f3
  END IF
  q3old  = q3
  q3     = q1-f1*(q2-q1)/(f2-f1)
  IF(ABS(q3-q3old) < endq0) EXIT
END DO
!---------------------------------------------------------------------
!     improve position of the zero by quadratic interpolation:
!---------------------------------------------------------------------
f1     = fbint(q1,foin)
f2     = fbint(q2,foin)
f3     = fbint(q3,foin)
ca     = 1D0/((q1-q2)*(q2-q3)*(q3-q1))
c11    = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca
!---------------------------------------------------------------------
!      the two zeroes of the interpolating binomial
!---------------------------------------------------------------------
qa1    = ( SQRT(ABS(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)
qa2    = (-SQRT(ABS(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)
!      select the zero which is closer to q3
q      = qa1
IF(ABS(qa1-q3) > ABS(qa2-q3)) q=qa2
qzero  = q
RETURN
END FUNCTION qzero

END MODULE Formfactor
