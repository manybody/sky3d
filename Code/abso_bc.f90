!------------------------------------------------------------------------------
! MODULE: abso_bc
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!Package containing absorbing boundary conditions and
!!subsequent analysis of observables from electron emission.
!!One may also consider to perform the cumulation of total
!!absorbed density within the absorbing loop.
!------------------------------------------------------------------------------
MODULE abso_bc
  USE Params
  USE Grids, ONLY: dx,dy,dz,x,y,z,wxyz
  USE Levels, ONLY: nstmax,psi,npsi,npmin,wocc,isospin,nprot,nneut
  USE Meanfield
  USE Parallel
  IMPLICIT NONE
  SAVE
  !***********************************************************************
  !                                                                      *
  !        General variables for the absorbing bounds                    *
  !                                                                      *
  !***********************************************************************
  INTEGER,PARAMETER :: nocc=1    !< degeneracy of states (outdated)
  INTEGER  :: ispherabso=0       !< switch to spherical absorbing  bounds
  INTEGER  :: iangabso=0
  INTEGER  :: ipes=0
  INTEGER  :: ifabsorbital=0     !< switch to orbit-wise accumulation
  INTEGER  :: nangtheta=0        !< number of reference points for outgoing spectra
  INTEGER  :: nangphi=0
  INTEGER  :: jescmaskorb=0
  INTEGER  :: jescmask=10
  REAL(db) :: powabs=0.0375D0    !< power of absorbing mask

  REAL(db),PRIVATE :: xango
  REAL(db),PRIVATE :: yango
  REAL(db),PRIVATE :: zango
  INTEGER,PRIVATE :: ix
  INTEGER,PRIVATE :: iy
  INTEGER,PRIVATE :: iz
  INTEGER,PRIVATE :: nst
  INTEGER,PRIVATE :: iq
  INTEGER,PRIVATE :: nta
  INTEGER,PRIVATE :: ita
  INTEGER,PRIVATE :: nabso=0     !<number of absorbing grid points in each direction
  !!0 signals no absorbing bounds
  REAL(db),PRIVATE :: timea
  REAL(db),PRIVATE,ALLOCATABLE :: rhoabso(:,:,:,:)
  REAL(db),PRIVATE,ALLOCATABLE :: absomask(:,:,:)
  REAL(db),PRIVATE,ALLOCATABLE :: spherloss(:,:,:)
  REAL(db),PRIVATE,ALLOCATABLE :: rhoescmaskorb(:,:,:,:)
  LOGICAL,PRIVATE,ALLOCATABLE :: tgridabso(:,:,:)
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: absbc
!> @brief
!!Apply absorbing bounds, optionally accumulate absorbed density.
!--------------------------------------------------------------------------- 
  SUBROUTINE absbc(nabsorb,it,nt,time)
    INTEGER,INTENT(in) :: it,nt,nabsorb
    REAL(db),INTENT(in) :: time
    REAL(db) :: weight
    INTEGER :: iqa

    LOGICAL,SAVE :: firstcall = .TRUE.
    REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: rhoabso_all

    !------------------------------------------------------------

    nabso = nabsorb
    IF(nabsorb <= 0) RETURN

    ! set calling parameters to internal common variables
    ita = it
    nta = nt
    timea = time

    ! initializations

    IF (firstcall) THEN
       IF(nabsorb>=NX/2.OR.nabsorb>=NY/2 .OR. nabsorb>=NZ/2) &
            STOP 'ABSBC: problem with nabsorb'
       CALL init_absbc()
       firstcall = .FALSE.
    END IF

    ! applying absorbing masks and accumulating absorption




    ! apply mask function (and accumulate absorption per state)
    ! and compute new density

    DO nst = 1,nstloc
       iqa = isospin(globalindex(nst))
       weight = nocc*wocc(nst)  
       DO iz=1,NZ; DO iy=1,NY; DO ix=1,NX
          IF(tgridabso(ix,iy,iz)) THEN
             rhoabso(ix,iy,iz,iqa) =  rhoabso(ix,iy,iz,iqa) + &
                  weight*spherloss(ix,iy,iz)*(  &
                  abs2(psi(ix,iy,iz,1,nst))+abs2(psi(ix,iy,iz,2,nst)) )
             psi(ix,iy,iz,1,nst)=absomask(ix,iy,iz)*psi(ix,iy,iz,1,nst)
             psi(ix,iy,iz,2,nst)=absomask(ix,iy,iz)*psi(ix,iy,iz,2,nst)
          END IF
       END DO;END DO;END DO
    END DO

    IF(ifabsorbital==1) THEN
       DO nst = 1,nstloc
          DO iz=1,NZ; DO iy=1,NY; DO ix=1,NX
             IF(tgridabso(ix,iy,iz)) THEN
                rhoescmaskorb(ix,iy,iz,nst) = rhoescmaskorb(ix,iy,iz,nst) &
                     + spherloss(ix,iy,iz)*(  &
                     abs2(psi(ix,iy,iz,1,nst))+abs2(psi(ix,iy,iz,2,nst)) )
             END IF
          END DO;END DO;END DO
       END DO
    END IF



    ! print escaped nucleons

    IF(tmpi .AND. (MOD(ita,mprint)==1 .OR. MOD(ita,mrest)==0 .OR. &
         MOD(ita,jescmask) == 0 .OR. ita == nta) ) THEN
       ALLOCATE(rhoabso_all(NX,NY,NZ,2))
       CALL mpi_barrier (mpi_comm_world, mpi_ierror)
       CALL mpi_allreduce(rhoabso,rhoabso_all,2*NX*NY*NZ,        &
            mpi_double_precision,mpi_sum,          &
            mpi_comm_world,mpi_ierror)
       IF(MOD(ita,mprint)==1) CALL nescape(rhoabso_all)
       CALL escmask(rhoabso_all)
    ELSE
       IF(MOD(ita,mprint)==1) CALL nescape(rhoabso)
       CALL escmask(rhoabso)
    END IF

    ! prints wavefunction on measuring points on file

    !  IF (ipes /= 0) CALL evalMP()      ! not yet implemented

    IF(MOD(ita,mrest)==0) THEN  
       IF(tmpi) THEN
          CALL wrtabso(ita,rhoabso_all)
       ELSE
          CALL wrtabso(ita,rhoabso)
       END IF
       WRITE(6,*) ' Augment restart file at it=',ita
    ENDIF
    IF(ALLOCATED(rhoabso_all)) DEALLOCATE(rhoabso_all)

    RETURN
  END SUBROUTINE absbc
!---------------------------------------------------------------------------  
! DESCRIPTION: init_absbc
!> @brief
!!Initializes absorbing boundary conditions.
!---------------------------------------------------------------------------

  SUBROUTINE init_absbc()

    NAMELIST /absobc/ ispherabso,iangabso,ipes,ifabsorbital, &
         nangtheta,nangphi,jescmaskorb,jescmask,powabs

    READ(5,absobc)

    ALLOCATE(tgridabso(NX,NY,NZ))
    ALLOCATE(absomask(NX,NY,NZ))
    ALLOCATE(spherloss(NX,NY,NZ))
    IF(ispherabso /= 0) THEN
       CALL init_spherabso()
    ELSE
       CALL init_abso()
    END IF
    IF(ifabsorbital==1) THEN
       ALLOCATE(rhoescmaskorb(NX,NY,NZ,nstloc))
       IF(.NOT.trestart) rhoescmaskorb = 0D0
    END IF
    ALLOCATE(rhoabso(NX,NY,NZ,2))
    IF(.NOT.trestart) rhoabso = 0D0
    !  IF (ipes /= 0) CALL init_MP()      ! not yet implemented
    IF(trestart) THEN
       CALL readabso()
       WRITE(6,*) ' Reading abso-info from restart file'
    ENDIF

  END SUBROUTINE init_absbc




  REAL(db) FUNCTION dist_min()

    ! determines minimum radius for last non-absorbing point

    REAL(db) :: rmin

    xango=0D0
    yango=0D0
    zango=0D0

    rmin = abs(x(NX-nabso)-xango)
    rmin = min(rmin,abs(x(nabso+1)-xango))
    rmin = min(rmin,abs(y(NY-nabso)-yango))
    rmin = min(rmin,abs(y(nabso+1)-yango))
    rmin = min(rmin,abs(z(NZ-nabso)-zango))
    rmin = min(rmin,abs(z(nabso+1)-zango))

    dist_min = rmin


    RETURN  
  END FUNCTION dist_min
!---------------------------------------------------------------------------  
! DESCRIPTION: init_sphereabso
!> @brief
!!Initializes sphericalabsorbing boundary conditions.
!---------------------------------------------------------------------------
  SUBROUTINE init_spherabso()

    !     Initializes mask function for spherical boundary conditions

    REAL(db) :: dmin1,dmin2,bcrad,dmin12,dmin22,dist,dist2,cosact

    !------------------------------------------------------------

    bcrad = nabso*dx
    dmin1 = dist_min()
    IF (DMIN1 < 0.) STOP 'Error in abso: dmin1<0'
    dmin2 = dmin1+bcrad
    dmin22=dmin2**2
    dmin12=DMIN1**2

    DO iz=1,NZ; DO iy=1,NY; DO ix=1,NX
       dist2 = x(ix)**2+y(iy)**2+z(iz)**2
       IF (dist2 <= dmin12) THEN
          absomask(ix,iy,iz)=1D0
          tgridabso(ix,iy,iz)=.false.
       ELSE IF (dist2 > dmin22) THEN
          absomask(ix,iy,iz)=0D0
          tgridabso(ix,iy,iz)=.true.
       ELSE
          dist = MAX(1D-20,SQRT(dist2))
          cosact = COS((dist-dmin1)*0.5D0*pi/bcrad)
          IF(cosact > 0D0) THEN
             absomask(ix,iy,iz)= cosact**powabs
          ELSE
             absomask(ix,iy,iz)= 0D0
          END IF
          tgridabso(ix,iy,iz)=.true.
       END IF
       spherloss(ix,iy,iz) = 1D0-absomask(ix,iy,iz)**2
    END DO;END DO;END DO

    RETURN
  END SUBROUTINE init_spherabso
!---------------------------------------------------------------------------  
! DESCRIPTION: init_abso
!> @brief
!!Initializes normal absorbing boundary conditions.
!---------------------------------------------------------------------------
  SUBROUTINE init_abso()
    USE Parallel

    !     Initializes mask for rectangular absorbing boundaries conditions

    REAL(db),ALLOCATABLE,DIMENSION(:) :: xmask,ymask,zmask

    LOGICAL :: wflagabs=.true.

    ALLOCATE(xmask(NX))
    ALLOCATE(ymask(NY))
    ALLOCATE(zmask(NZ))

    !     prepare mask functions in each direction separately

    zmask = 1D0
    DO iz = 1,nabso
       zmask(iz) = COS(pi*0.5D0*(nabso+1.0-iz)/nabso)**powabs
       zmask(NZ+1-iz) = zmask(iz)
    END DO

    ymask = 1D0
    DO iy = 1,nabso
       ymask(iy) = COS(pi*0.5D0*(nabso+1.0-iy)/nabso)**powabs
       ymask(NY+1-iy) = ymask(iy)
    END DO

    xmask = 1D0
    DO ix = 1,nabso
       xmask(ix) = COS(pi*0.5D0*(nabso+1.0D0-ix)/nabso)**powabs
       xmask(NX+1-ix) = xmask(ix)
    END DO

    IF(mpi_myproc == 0 .AND. wflagabs) THEN
       WRITE(6,'(a)') ' ZMASK:'
       WRITE(6,'(1x,5(1pg12.4))') zmask
       WRITE(6,'(a)') ' YMASK:'
       WRITE(6,'(1x,5(1pg12.4))') ymask
       WRITE(6,'(a)') ' XMASK:'
       WRITE(6,'(1x,5(1pg12.4))') xmask
    END IF


    !     compose to one mask function on all grid points

    FORALL(ix=1:NX,iy=1:NY,iz=1:NZ)
       absomask(ix,iy,iz)=xmask(ix)*ymask(iy)*zmask(iz)
       tgridabso(ix,iy,iz)= absomask(ix,iy,iz) < 0.999999999999D0
       spherloss(ix,iy,iz) = 1D0-absomask(ix,iy,iz)**2
    END FORALL

    RETURN
  END SUBROUTINE init_abso



  !-----escmask---------------------------------------------------------

  SUBROUTINE escmask(rhoabsoprint)
    USE Parallel

    ! print collected information on escaping nucleons

    REAL(db),INTENT(IN),DIMENSION(:,:,:,:) :: rhoabsoprint

    CHARACTER(LEN=5) :: scratchline
    CHARACTER(LEN=5) :: readstring
    CHARACTER(LEN=4) :: str

    !--------------------------------------------------------------------




    IF(ifabsorbital==1) THEN
       IF (MOD(ita,jescmaskorb) == 0 .OR. ita==nta) THEN
          DO nst=1,nstmax
             IF(nstmax > 9999)  STOP 'ERROR: Too many states for escmaskOrb'
             IF(mpi_myproc == node(nst)) THEN
                WRITE(scratchline,'(i5)') 10000+nst
                READ(scratchline,'(a)') readstring
                str = readstring(2:5)
                OPEN(588,STATUS='unknown', FILE='escmaskOrb.'//str)
                WRITE(588,'(a,f12.4,a,i5,/a,3i5/a)') &
                     '# distribution of escaped nucleon at time=',timea,&
                     ' for state=',nst, &
                     '#nx,ny,nz:',NX,NY,NZ,                                    &
                     '#      x         y          z         loss of density'

                DO iz=1,NZ; DO iy=1,NY; DO ix=1,NX
                   WRITE(588,'(3f12.4,1e17.7)') &
                        x(ix),y(iy),z(iz),rhoescmaskorb(ix,iy,iz,node(nst))
                END DO; END DO; END DO

                CLOSE(588)
             END IF
          END DO
       END IF
    END IF


    IF(mpi_myproc == 0) THEN
       IF (MOD(ita,jescmask) == 0 .OR. ita == nta) THEN
          DO iq=1,2
             IF(iq == 1)THEN
                OPEN(589,STATUS='unknown',FILE='escmask.1')
                WRITE(589,'(a,f12.4/a,3i5/a)') &
                     '# total distribution of escaped protons at time=',timea, &
                     '#nx,ny,nz:',NX,NY,NZ,                                    &
                     '#      x         y          z         loss of density'
             ELSE
                OPEN(589,STATUS='unknown',FILE='escmask.2')
                WRITE(589,'(a,f12.4/a,3i5/a)') &
                     '# total distribution of escaped neutrons at time=',timea, &
                     '#nx,ny,nz:',NX,NY,NZ,                                     &
                     '#      x         y          z         loss of density'
             END IF

             DO iz=1,NZ; DO iy=1,NY; DO ix=1,NX
                WRITE(589,'(3f12.4,1e17.7)') &
                     x(ix),y(iy),z(iz),rhoabsoprint(ix,iy,iz,iq)
             END DO; END DO; END DO

             CLOSE(589)
          END DO
       END IF
    END IF


    RETURN
  END SUBROUTINE escmask




  SUBROUTINE nescape(rhoabsoprint)

    !  Compute and print total number of escaped nucleons

    REAL(db),INTENT(IN),DIMENSION(:,:,:,:) :: rhoabsoprint

    LOGICAL,SAVE :: firstcall = .TRUE.
    LOGICAL :: texist
    REAL(db) :: aprot,aneut

    !------------------------------------------------------------

    IF(firstcall) THEN
       INQUIRE(file='nescape.res',exist=texist)
       IF(.NOT.texist) THEN
          OPEN(unit=scratch,file='nescape.res')
          WRITE(scratch,'(a)') &
               '# number of escaped nucleons ', &
               '#   time    lost protons    lost neutrons '
          CLOSE(unit=scratch)
       END IF
       firstcall = .FALSE.
    END IF

    aprot = 1D0*nprot   !SUM(fcharg(1:2))
    aneut = 1D0*nneut   !SUM(fmass(1:2))-aprot
    IF(mpi_myproc == 0) THEN
       OPEN(unit=scratch,file='nescape.res',POSITION='APPEND')
       WRITE(scratch,'(1x,f10.4,4(1pg20.10))') timea, &
            aprot-wxyz*sum(rho(:,:,:,1)),aneut-wxyz*sum(rho(:,:,:,2)), &
            wxyz*sum(rhoabsoprint(:,:,:,1)),wxyz*sum(rhoabsoprint(:,:,:,2))
       CLOSE(unit=scratch)
    END IF

    RETURN
  END SUBROUTINE nescape

  !------------------------------------------------------------


  SUBROUTINE wrtabso(iteration,rhoabsoprint)
    !
    USE Parallel, ONLY: mpi_myproc,mpi_nprocs,nstloc
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: iteration
    INTEGER :: nst,lenstr
    CHARACTER(120) :: rsfp
    REAL(db),INTENT(IN),DIMENSION(:,:,:,:) :: rhoabsoprint
    !
    !***********************************************************************
    !                                                                      *
    !        complement wffile by info on accumulated absorption      *
    !                                                                      *
    !***********************************************************************
    !
    ! first write general information only with processor 0
    !
    IF(wffile=='none'.OR.wffile=='NONE') THEN
       WRITE(6,*) " wffile='NONE'  --> no file written "
       RETURN
    ENDIF
    lenstr = LEN_TRIM(wffile)
    IF(mpi_myproc==0) THEN
       OPEN(UNIT=scratch2,FILE=wffile(1:lenstr)//'.a',STATUS='REPLACE',   &
            FORM='UNFORMATTED')
       WRITE(scratch2) rhoabsoprint
    ENDIF
    ! now if we are running with one processor, also write the wave functions
    ! onto that file; otherwise open new individualized file
    !
    IF(ifabsorbital==1) THEN
       IF(mpi_nprocs>1) THEN
          IF(mpi_myproc/=0) CLOSE(UNIT=scratch2,STATUS='KEEP')
          WRITE(rsfp,'(I3.3,''.'',A)') mpi_myproc,wffile(1:lenstr)//'.a'
          OPEN(UNIT=scratch2,FILE=rsfp,STATUS='REPLACE',  &
               FORM='UNFORMATTED')
       ENDIF
       !
       DO nst = 1,nstloc
          WRITE(scratch2) rhoescmaskOrb(:,:,:,nst)
       ENDDO
    END IF
    CLOSE(UNIT=scratch2,STATUS='KEEP')
    WRITE(6,*) 'Restart file augmented at iteration ',iteration
  END SUBROUTINE wrtabso

  SUBROUTINE readabso()
    !
    USE Parallel, ONLY: mpi_myproc,mpi_nprocs,nstloc
    IMPLICIT NONE
    INTEGER :: nst,lenstr
    CHARACTER(120) :: rsfp
    !
    !***********************************************************************
    !                                                                      *
    !        complement wffile by info on accumulated absorption      *
    !                                                                      *
    !***********************************************************************
    !
    ! first write general information only with processor 0
    !
    IF(wffile=='none'.OR.wffile=='NONE') THEN
       WRITE(6,*) " wffile='NONE'  --> no file read "
       RETURN
    ENDIF
    lenstr = LEN_TRIM(wffile)
    IF(mpi_myproc==0) THEN
       OPEN(UNIT=scratch2,FILE=wffile(1:lenstr)//'.a',STATUS='OLD',   &
            FORM='UNFORMATTED')
       READ(scratch2) rhoabso
    ENDIF
    ! now if we are running with one processor, also write the wave functions
    ! onto that file; otherwise open new individualized file
    !
    IF(ifabsorbital==1) THEN
       IF(mpi_nprocs>1) THEN
          IF(mpi_myproc/=0) CLOSE(UNIT=scratch2,STATUS='KEEP')
          WRITE(rsfp,'(I3.3,''.'',A)') mpi_myproc,wffile(1:lenstr)//'.a'
          OPEN(UNIT=scratch2,FILE=rsfp,STATUS='OLD',  &
               FORM='UNFORMATTED')
       ENDIF
       !
       DO nst = 1,nstloc
          READ(scratch2) rhoescmaskOrb(:,:,:,nst)
       ENDDO
    END IF
    CLOSE(UNIT=scratch2,STATUS='KEEP')
    WRITE(6,*) 'Absorption data read from restart file'
  END SUBROUTINE readabso

  PURE REAL(db) FUNCTION abs2(c)
    COMPLEX(db),INTENT(IN) :: c
    abs2 = REAL(c)**2 + AIMAG(c)**2
    RETURN
  END FUNCTION abs2


END MODULE abso_bc
