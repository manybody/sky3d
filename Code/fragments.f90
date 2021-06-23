MODULE Fragments
  USE Params
  USE Grids 
  USE Forces, ONLY: f,nucleon_mass
  USE Levels
  IMPLICIT NONE 
  SAVE
  PRIVATE
  LOGICAL :: fix_boost
  REAL(db) :: ecm      ! center-of-mass energy
  REAL(db) :: b        ! impact parameter
  CHARACTER(64) :: filename(mnof) ! names of fragment wf files
  REAL(db) :: fcent(3,mnof)   ! c.m. vectors for the fragments
  REAL(db) :: fboost(3,mnof)  ! boost vectors for the fragments
  REAL(db) :: fcmtot(3,mnof) ! center of mass in data file
  REAL(db) :: fmass(mnof),fcharge(mnof) ! A and Z of fragments
  INTEGER :: fnneut(mnof),fnprot(mnof),fnstmax(mnof),fnpmin(2,mnof), &
       fnpsi(2,mnof),fnumber(2,mnof),fnewnpmin(2,mnof),fnewnpsi(2,mnof)
  INTEGER :: fnx(mnof),fny(mnof),fnz(mnof)
  PUBLIC :: getin_fragments, read_fragments
  PRIVATE :: locate,phases
  REAL(db) :: start_fragments, finish_fragments
CONTAINS
  !*******************************************************************
  SUBROUTINE getin_fragments
    USE Twobody, ONLY: centerx,centerz
    INTEGER :: i
    REAL(db) :: fdx,fdy,fdz,fwxyz
    CHARACTER(8) :: forcename
    NAMELIST /fragments/ filename,fcent,fboost,ecm,b,fix_boost
    IF(trestart) THEN
       filename(1)=wffile
       fcent(:,1)=0.D0
       fboost(:,1)=0.D0
       fix_boost=.TRUE.
       IF(wflag) WRITE(*,*) 'Reading restart as one fragment from file ',wffile
    ELSE
       fix_boost=.FALSE.
       READ(5,fragments)
    END IF
    IF(nof/=2.AND..NOT.fix_boost) THEN  
       IF(wflag) WRITE(*,"(//,a,//)") "Non-fixed boost only for nof=2"  
       STOP  
    ENDIF
    IF(wflag) THEN
       IF(fix_boost) THEN
          WRITE(*,*) '***** Fixed boost values used *****'
       ELSE
          WRITE(*,*) '***** Boost values computed from twobody kinematics *****'
       END IF
    END IF
    fnewnpmin(:,1)=1
    DO i=1,nof
       OPEN(UNIT=scratch,FILE=filename(i),STATUS='old',FORM=&
            'unformatted')
       READ(scratch) iter,time,forcename,fnstmax(i),fnneut(i),fnprot(i), &
            fnumber(:,i),fnpsi(:,i),fcharge(i),fmass(i),fcmtot(:,i)
       IF(trestart) fcmtot(:,i)=0.D0
       IF(forcename/=f%name) THEN
          IF(wflag) WRITE(*,*) 'Forces do not agree: ',forcename, &
               ' : ',f%name,' Fragment:',filename(i)
          STOP
       END IF
       READ(scratch) fnx(i),fny(i),fnz(i),fdx,fdy,fdz,fwxyz
       CLOSE(UNIT=scratch)
       IF(fnx(i)>nx.OR.fny(i)>ny.OR.fnz(i)>nz) &
            CALL errormsg('Fragment dimensioning larger than present', &
            filename(i))
       IF(MAX(ABS(fdx-dx),ABS(fdy-dy),ABS(fdz-dz)) &
            >0.01D0) CALL errormsg('Grid spacing different for', &
            filename(i))
       IF(wflag) THEN
          WRITE(*,"(A,I2,2A)") " ***** Data for fragment #",i,' from file ',&
               filename(i)
          WRITE(*,"(A,1P,3E12.4)") " Location(x,y,z):  ",fcent(:,i)
          WRITE(*,"(2(A,F9.4))") " Mass number: ",fmass(i), &
               ' Charge number: ',fcharge(i)   
          IF(fix_boost) &
               WRITE(*,"(A,1P,3E12.4)") " Boost(x,y,z):     ",fboost(:,i)
       END IF
       ! for static restart read also unoccupied wave functions
       IF(tstatic.AND.nof==1) THEN
          fnumber(1,i)=fnpsi(1,i)
          fnumber(2,i)=fnpsi(2,i)-fnpsi(1,i)
       ENDIF
       ! get isospin starting points in fragment file
       fnpmin(1,i)=1
       fnpmin(2,i)=fnpsi(1,i)+1
       IF(i<nof) &
            fnewnpmin(:,i+1)=fnewnpmin(:,i)+fnumber(:,i)
       fnewnpsi(:,i)=fnewnpmin(:,i)+fnumber(:,i)-1
    END DO
    ! determine total wave function numbers
    nneut=SUM(fnneut)
    nprot=SUM(fnprot)
    npsi(1)=SUM(fnumber(1,:))
    npsi(2)=SUM(fnumber(2,:))
    nstmax=npsi(1)+npsi(2)
    npmin(1)=1
    npmin(2)=npsi(1)+1
    npsi(2)=nstmax
    fnewnpmin(2,:)=fnewnpmin(2,:)+npsi(1)
    fnewnpsi(2,:)=fnewnpsi(2,:)+npmin(2)-1
    charge_number=nprot
    mass_number=nneut+nprot
    IF(.NOT.fix_boost.AND.nof==2.AND.tdynamic) CALL twobody_init
    IF(nof==2) THEN
       ! Give fragment positions to 2-body analysis as initial guess
       centerx=fcent(1,1:2)
       centerz=fcent(3,1:2)
    END IF
  END SUBROUTINE getin_fragments
  !*******************************************************************
  SUBROUTINE read_fragments
    INTEGER :: iff
    IF(wflag) &
         WRITE(*,*) '***** Input of fragment wave functions *****'
    DO iff=1,nof  
       CALL read_one_fragment(iff)
       CALL boost_fragment(iff)  
    ENDDO
    IF(wflag) &
         WRITE(*,*) '***** All fragments loaded and boosted *****'
  END SUBROUTINE read_fragments
  !*******************************************************************
  SUBROUTINE read_one_fragment(iff)
    USE Parallel, ONLY: node,mpi_myproc,localindex
    USE Fourier
#ifdef CUDA
    USE cudafor
    USE cufft_m
#endif
    INTEGER,INTENT(IN) :: iff
    LOGICAL :: multifile
    INTEGER :: ipn
    COMPLEX(db) :: ps1(nx,ny,nz,2),akx(nx),aky(ny),akz(nz)
#ifdef CUDA
    INTEGER :: istat
    COMPLEX(db), ALLOCATABLE, DEVICE :: ps1_d(:,:,:,:)
#endif
    INTEGER :: iq,is,nst,oldnst,newnst,ix,iy,iz,iold,inew
    REAL(db) :: cmi(3)
    REAL(db) :: fx(fnx(iff)),fy(fny(iff)),fz(fnz(iff))
    REAL(db),DIMENSION(fnstmax(iff)) :: fwocc,fsp_energy,fsp_parity, &
         fsp_norm,fsp_kinetic,fsp_efluct1
    INTEGER,DIMENSION(fnstmax(iff)) :: fnode,flocalindex
#ifdef CUDA
    ALLOCATE(ps1_d(nx,ny,nz,2))
#endif
    OPEN(UNIT=scratch,FILE=filename(iff),STATUS='old',FORM=&
         'unformatted')
    READ(scratch) 
    READ(scratch) 
    READ(scratch) fx,fy,fz
    ! read the s.p.quantities, reinserting them into the correct positions
    READ(scratch) fwocc,fsp_energy,fsp_parity,fsp_norm,fsp_kinetic,fsp_efluct1
    ps1=0.0D0
    IF(wflag) THEN
       WRITE(*,'(A,I3)') ' ***** S.p. levels for fragment #',iff
       WRITE(*,'(A,3I5,A,2I5)') '   Iso   #new  #old  Occup.     E_sp        &
            &Parity       Norm       E_kin      E_fluct' 
    END IF
    DO iq=1,2
       DO inew=fnewnpmin(iq,iff),fnewnpsi(iq,iff)
          iold=fnpmin(iq,iff)+inew-fnewnpmin(iq,iff)
          wocc(inew)=fwocc(iold)
          sp_energy(inew)=fsp_energy(iold)
          sp_parity(inew)=fsp_parity(iold)
          sp_norm(inew)=fsp_norm(iold)
          sp_kinetic(inew)=fsp_kinetic(iold)
          sp_efluct1(inew)=fsp_efluct1(iold)
          isospin(inew)=iq
          IF(wflag) THEN
             WRITE(*,'(3I6,6G12.5)') &
                  iq,inew,iold,wocc(inew),sp_energy(inew),sp_parity(inew), &
                  sp_norm(inew),sp_kinetic(inew),sp_efluct1(inew)
          ENDIF
       END DO
    ENDDO
    cmi(1)=(fcent(1,iff)-fcmtot(1,iff)-x(1)+fx(1))/(nx*dx)
    cmi(2)=(fcent(2,iff)-fcmtot(2,iff)-y(1)+fy(1))/(ny*dy)
    cmi(3)=(fcent(3,iff)-fcmtot(3,iff)-z(1)+fz(1))/(nz*dz)
    IF(nof==2.AND.wflag) WRITE(*,'(A,I3,3F10.5)') &
         ' Translation for fragment ',iff,nx*cmi(1),ny*cmi(2),nz*cmi(3)
    CALL phases(nx,akx,cmi(1))
    CALL phases(ny,aky,cmi(2))
    CALL phases(nz,akz,cmi(3))
    ! Code for multiple files: find out which file each w.f. is on and at which index
    READ(scratch) fnode,flocalindex
    multifile=ANY(fnode/=0)
    ! start reading wave functions
    DO iq=1,2
       ! skip unoccupied neutron states in one-file case
       IF(iq==2.AND..NOT.multifile) THEN
          DO nst=fnumber(1,iff)+1,fnpsi(1,iff)
             READ(scratch)
          END DO
       END IF
       ! Read occupied states
       DO nst=1,fnumber(iq,iff)
          newnst=fnewnpmin(iq,iff)+nst-1
          oldnst=fnpmin(iq,iff)+nst-1 ! includes empty states
          IF(node(newnst)==mpi_myproc) THEN           
             ipn=localindex(newnst)
             ! open correct file and position for multifile case
             IF(multifile) CALL locate(fnode(oldnst),flocalindex(oldnst))
             ps1=0.D0
             READ(scratch) ps1(1:fnx(iff),1:fny(iff),1:fnz(iff),:)
             DO is=1,2
#ifdef CUDA
                ! Copy to device, perform FFT, copy back to host
                istat = cudaMemcpy(ps1_d(1,1,1,is), ps1(1,1,1,is), nx*ny*nz)
                CALL cufftExecZ2Z(pforward,ps1_d(:,:,:,is),ps1_d(:,:,:,is),CUFFT_FORWARD)
                istat = cudaMemcpy(ps1(1,1,1,is), ps1_d(1,1,1,is), nx*ny*nz)
#else
                CALL dfftw_execute_dft(pforward,ps1(:,:,:,is),ps1(:,:,:,is))
#endif
                FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
                   ps1(ix,iy,iz,is)=ps1(ix,iy,iz,is)*akx(ix)*aky(iy)*akz(iz) &
                        /DBLE(nx*ny*nz)
                END FORALL
#ifdef CUDA
                istat = cudaMemcpy(ps1_d(1,1,1,is), ps1(1,1,1,is), nx*ny*nz)
                CALL cufftExecZ2Z(pbackward,ps1_d(:,:,:,is),ps1_d(:,:,:,is),CUFFT_INVERSE)
                istat = cudaMemcpy(ps1(1,1,1,is), ps1_d(1,1,1,is), nx*ny*nz)
#else
                CALL dfftw_execute_dft(pbackward,ps1(:,:,:,is),ps1(:,:,:,is))
#endif
                WHERE(ABS(ps1)==0.D0)
                   psi(:,:,:,:,ipn)=1.0D-20
                ELSEWHERE
                   psi(:,:,:,:,ipn)=ps1
                ENDWHERE
             ENDDO
          ELSE
             IF(.NOT.multifile) READ(scratch)
          ENDIF
       ENDDO
    ENDDO
    CLOSE(unit=scratch)
#ifdef CUDA
    DEALLOCATE(ps1_d)
#endif
  END SUBROUTINE read_one_fragment
  !*******************************************************************
  ! position scratch file in correct location
  SUBROUTINE locate(fileno,pos)
    INTEGER,INTENT(IN) :: fileno,pos
    INTEGER,SAVE :: presentfile=-1
    INTEGER :: i
    CHARACTER(120) :: rsfp
    IF(fileno/=presentfile) THEN
       CLOSE(scratch)
       WRITE(rsfp,'(I3.3,''.'',A)') fileno,wffile
       OPEN(scratch,FORM='unformatted',FILE=rsfp,STATUS='OLD')
       DO i=1,pos-1
          READ(scratch)
       END DO
    END IF
  END SUBROUTINE locate
  !*******************************************************************
  ! phases calculates the phases for the translation
  PURE SUBROUTINE phases(n,a,c)
    INTEGER :: n,i,k
    COMPLEX(db),INTENT(OUT) :: a(n)
    REAL(db) :: c
    INTENT(IN) :: n,c
    DO i=1,n
       k=i-1
       IF(i>(n+1)/2) k=k-n
       a(i)=EXP(CMPLX(0.D0,-2.D0*pi*k*c,db))
    END DO
  END SUBROUTINE phases
  !*******************************************************************
  SUBROUTINE twobody_init
    REAL(db) :: vrel,vrel_d,b_d,v1,v2,ec,dix,diz,totmass,sint,cost, &
         xli,xmu,roft
    INTEGER :: i
    IF(ABS(fcent(2,1)-fcent(2,2))>1.D-5) &
         STOP 'Two fragments must be in the x-z plane'
    fboost=0.D0
    totmass=SUM(fmass(1:2))
    xmu=fmass(1)*fmass(2)/totmass*nucleon_mass
    vrel=SQRT(2.D0*ecm/xmu)
    xli=xmu*vrel*b/hbc
    dix=fcent(1,1)-fcent(1,2)
    diz=fcent(3,1)-fcent(3,2)
    roft=SQRT(dix**2+diz**2)
    dix=dix/roft; diz=diz/roft
    ec=e2*fcharge(1)*fcharge(2)/roft
    IF(ec>ecm) STOP 'Not enough energy to reach this distance'
    vrel_d=SQRT(2.D0*(ecm-ec)/xmu)
    v1=fmass(2)/totmass*vrel_d
    v2=fmass(1)/totmass*vrel_d
    b_d=xli*hbc/(xmu*vrel_d)
    sint=b_d/roft
    cost=SQRT(1.D0-sint**2)
    fboost(1,1)=-v1*(dix*cost-diz*sint)
    fboost(3,1)=-v1*(dix*sint+diz*cost)
    fboost(1,2)=v2*(dix*cost-diz*sint)
    fboost(3,2)=v2*(dix*sint+diz*cost)
    IF(wflag) THEN
       WRITE(*,*) '***** Two-body initialization *****'
       WRITE(*,"(2(A,F12.4),A)") " c. m. Energy:",ecm, &
            ' MeV. Impact parameter b: ',b,' fm'
       WRITE(*,"(A,F12.4)") " xli/hbar    :",xli  
       WRITE(*,*) 'Computed boost velocities in units of c'
       WRITE(*,'(A,I2,3G15.6)') (' Fragment #',i,fboost(:,i),i=1,2)
    END IF
    !convert into kinetic energy for that direction
    FORALL(i=1:2)
       fboost(:,i)=SIGN(0.5D0*fmass(i)*fboost(:,i)**2*nucleon_mass, &
            fboost(:,i))
    END FORALL
    IF(wflag) THEN
       WRITE(*,*) 'Total translational kinetic energies in units of MeV'
       WRITE(*,'(A,I2,3G15.6)') (' Fragment #',i,fboost(:,i),i=1,2)
    END IF
  END SUBROUTINE twobody_init
  !*******************************************************************
  SUBROUTINE boost_fragment(iff)  
    USE Parallel, ONLY: node,localindex,mpi_myproc
    INTEGER,INTENT(IN) :: iff
    REAL(db) :: sb(3),tb(3),akf(3)
    INTEGER :: iq,nst,is,ix,iy,iz
    sb=0.0D0  
    tb=fboost(:,iff)  
    WHERE(tb/=0.0D0) 
       sb=tb/ABS(tb)  
    END WHERE
    DO iq=1,2  
       akf=sb*SQRT(1.0D0/f%h2m(iq)*ABS(tb)/fmass(iff))  
       DO nst=fnewnpmin(iq,iff),fnewnpsi(iq,iff)
          IF(node(nst)==mpi_myproc) THEN
             FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)             
                psi(ix,iy,iz,is,localindex(nst))=&
                     psi(ix,iy,iz,is,localindex(nst)) * EXP( &
                     CMPLX(0.D0,akf(1)*x(ix)+akf(2)*y(iy)+akf(3)*z(iz),db))
             END FORALL
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE boost_fragment
  !*******************************************************************
  SUBROUTINE errormsg(msg1,msg2)
    CHARACTER(*),INTENT(IN) :: msg1,msg2
    IF(wflag) WRITE(*,*) msg1,msg2
    STOP
  END SUBROUTINE errormsg
END MODULE Fragments
