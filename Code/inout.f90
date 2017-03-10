MODULE Inout
  ! This module contains routines related to producing output files
  ! containing wave functions, densities, etc.
  USE Params
  USE Parallel, ONLY: node,localindex,mpi_myproc,mpi_nprocs,tabc_dens,tabc_vec_dens
  USE Grids
  USE Forces, ONLY:f
  USE Moment, ONLY: cm,cmtot
  USE Densities, ONLY: rho,tau,current,sdens,sodens,localization
  USE Meanfield, ONLY: upot
  USE Coulomb, ONLY: wcoul
  USE Levels
  IMPLICIT NONE
  PRIVATE :: write_one_density,write_vec_density
CONTAINS
  !**************************************************************************
  SUBROUTINE write_wavefunctions
    USE Parallel, ONLY: mpi_myproc,mpi_nprocs,nstloc,node,localindex
    INTEGER :: nst,iq,number(2)
    CHARACTER(120) :: rsfp
    IF(wflag.AND.(wffile=='none'.OR.wffile=='NONE')) THEN
       WRITE(*,*) " wffile='NONE'  --> no file written "
       RETURN
    ENDIF
    ! Determine number of states with non-zero occupation
    DO iq=1,2
       number(iq)=COUNT(wocc(npmin(iq):npsi(iq))>0.D0)
    END DO
    IF(mpi_myproc==0) THEN
       OPEN(UNIT=scratch2,FILE=wffile,STATUS='REPLACE',FORM='UNFORMATTED')
       WRITE(scratch2) iter,time,f%name,nstmax,nneut,nprot,number,npsi, &
            charge_number,mass_number,cm
       WRITE(scratch2) nx,ny,nz,dx,dy,dz,wxyz
       WRITE(scratch2) x,y,z
       WRITE(scratch2) wocc,sp_energy,sp_parity,sp_norm,sp_kinetic, &
            sp_efluct1
       WRITE(scratch2) node,localindex
       IF(mpi_nprocs>1) CLOSE(UNIT=scratch2,STATUS='KEEP')
    ENDIF
    IF(mpi_nprocs>1) THEN
       WRITE(rsfp,'(I3.3,''.'',A)') mpi_myproc,wffile
       OPEN(UNIT=scratch2,FILE=rsfp,STATUS='REPLACE',FORM='UNFORMATTED')
    ENDIF
    DO nst=1,nstloc
       WRITE(scratch2) psi(:,:,:,:,nst)
    ENDDO
    CLOSE(UNIT=scratch2,STATUS='KEEP')
  END SUBROUTINE write_wavefunctions
  !**************************************************************************
  SUBROUTINE write_densities
    CHARACTER(10) :: filename
    CHARACTER(1) :: c
    INTEGER :: i
    IF(.NOT.wflag) RETURN
    WRITE(filename,'(I6.6,A4)') iter,'.tdd'
    IF(.NOT.tdynamic) time=0.D0
    IF(tabc_myid==0) OPEN(UNIT=scratch,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE')
    WRITE(scratch) iter,time,nx,ny,nz
    WRITE(scratch) dx,dy,dz,wxyz,x,y,z
    DO i=1,nselect
       c=writeselect(i:i)
       SELECT CASE(c)
       CASE('r','R')
          CALL write_one_density('Rho',tabc_dens(rho))
       CASE('t','T')
          CALL write_one_density('Tau',tabc_dens(tau))
       CASE('u','U')
          CALL write_one_density('Upot',tabc_dens(upot))
       CASE('l','L')
          IF(write_isospin.EQV..FALSE.) STOP 'Please set <write_isospin=T> for localization plots'
          CALL write_one_density('Loc',tabc_dens(localization))
       CASE('w','W')
          WRITE(scratch) 'Wcoul     ',.FALSE.,.FALSE.
          WRITE(scratch) wcoul
       CASE('c','C')
          CALL write_vec_density('Current',tabc_vec_dens(current))
       CASE('s','S')
          CALL write_vec_density('Spindens',tabc_vec_dens(sdens))
       CASE('o','O')
          CALL write_vec_density('s-o-Dens',tabc_vec_dens(sodens))
       END SELECT
    END DO
    IF(tabc_myid==0) CLOSE(UNIT=scratch)
  END SUBROUTINE write_densities
  !**************************************************************************
  SUBROUTINE write_one_density(name,values)
    CHARACTER(*),INTENT(IN) :: name
    REAL(db),INTENT(IN) :: values(nx,ny,nz,2)
    CHARACTER(10) :: stored_name
    REAL(db) a(nx,ny,nz)
    stored_name=name
    WRITE(scratch) stored_name,.FALSE.,write_isospin
    IF(write_isospin) THEN
       IF(tabc_myid==0) WRITE(scratch) values
    ELSE
       a=values(:,:,:,1)+values(:,:,:,2)
       IF(tabc_myid==0) WRITE(scratch) a
    END IF
  END SUBROUTINE write_one_density
  !**************************************************************************
  SUBROUTINE write_vec_density(name,values)
    CHARACTER(*),INTENT(IN) :: name
    REAL(db),INTENT(IN) :: values(nx,ny,nz,3,2)
    CHARACTER(10) :: stored_name
    REAL(db) a(nx,ny,nz,3)
    stored_name=name
    WRITE(scratch) stored_name,.TRUE.,write_isospin
    IF(write_isospin) THEN
       IF(tabc_myid==0) WRITE(scratch) values
    ELSE
       a=values(:,:,:,:,1)+values(:,:,:,:,2)
       IF(tabc_myid==0) WRITE(scratch) a
    END IF
  END SUBROUTINE write_vec_density
  !**************************************************************************
  SUBROUTINE plot_density
    REAL(db),PARAMETER :: density_scale=0.14D0
    INTEGER,PARAMETER :: ixsc=10,izsc=6
    REAL(db) :: rhoplt(nx,nz),xperi
    CHARACTER(1) :: ifun(121),ibor(121)
    CHARACTER(1),PARAMETER :: icar(25)=(/' ','1',' ','3',' ', &
         '5',' ','7',' ','9',' ','b',' ','d',' ','f',' ','h', &
         ' ','j',' ','l',' ','n','*' /)
    REAL(db) ::  xco(nx+1),zco(nz+1),dimx,dimz,dxp,dzp, &
         xcu,zcu
    INTEGER :: nhor,nver,ntkx,ntkz,jcar,ntz,i,j
    IF(.NOT.wflag) RETURN
    rhoplt=rho(:,ny/2,:,1)+rho(:,ny/2,:,2)
    xperi=(x(nx)-x(1))/12.0D0
    ibor='-'
    ibor(1:121:10)='+'
    dimx=x(nx)-x(1); dimz=z(nz)-z(1)
    dxp=xperi/ixsc; dzp=xperi/izsc
    nhor=INT(dimx/dxp+1); nver=INT(dimz/dzp+1)
    nver=MIN(nver,121)
    ntkx=(nhor+ixsc-1)/ixsc; ntkz=(nver+izsc-1)/izsc
    WRITE(*,'(/,A,F12.4,A)') ' Contour 5=',density_scale,' n/fm**3'
    FORALL(i=1:ntkx) xco(i)=x(1)+(i-1)*ixsc*dxp
    FORALL(i=1:ntkz) zco(i)=z(nz)-(i-1)*izsc*dzp
    WRITE(*,'(A)') '  z '
    WRITE(*,'(1X,F6.2,1X,121A)') zco(1),ibor(1:nhor)
    DO i=2,nver-1
       zcu=z(nz)-(i-1)*dzp
       DO j=2,nhor-1
          xcu=x(1)+(j-1)*dxp
          jcar=1+INT(0.5D0+ &
               bplina(nx,nz,x,z,rhoplt,xcu,zcu) &
               *5.D0/density_scale)
          jcar=MIN(jcar,25)
          ifun(j)=icar(jcar)
       END DO
       IF(MOD(i-1,izsc)==0) THEN
          ntz=(i+izsc-1)/izsc
          ifun(nhor)='+'
          WRITE(*,'(1X,F6.2,1X,A,120A)') zco(ntz),'+',ifun(2:nhor)
       ELSE
          ifun(nhor)='i'
          WRITE(*,'(8X,A,120A)') 'i',ifun(2:nhor)
       END IF
    END DO
    IF(MOD(nver-1,izsc)==0) THEN
       WRITE(*,'(1X,F6.2,1X,121A)') zco(ntkz),ibor(1:nhor)
    ELSE
       WRITE(*,'(8X,A,120A)') '-',ibor(2:nhor)
    END IF
    WRITE(*,'(A,12(F6.2,4X),F6.2)') '  x= ',(xco(i),i=1,ntkx)
  END SUBROUTINE plot_density
  !**************************************************************************
  PURE FUNCTION bplina(n,m,xar,zar,fun,xcu,zcu) RESULT(ff)
    ! to do a bilinear inter. of fun(nx,nz), on xar(nx) and zar(nz) 
    INTEGER,INTENT(IN) :: n,m
    REAL(db),INTENT(IN) :: xar(n),zar(m),fun(n,m),xcu,zcu
    REAL(db) :: ff,dxf,dzf
    INTEGER :: i,icu,jcu
    !
    ff=0.D0
    DO i=1,n-1
       IF(xar(i)<=xcu.AND.xcu<=xar(i+1)) GO TO 10
    END DO
    RETURN
10  icu=i
    DO i=1,m-1
       IF(zar(i)<=zcu.AND.zcu<=zar(i+1)) GO TO 20
    END DO
    RETURN
20  jcu=i
    dxf=(xcu-xar(icu))/(xar(icu+1)-xar(icu))
    dzf=(zcu-zar(jcu))/(zar(jcu+1)-zar(jcu))
    ff=dxf*(fun(icu+1,jcu+1)*dzf+fun(icu+1,jcu)*(1.0D0-dzf)) &
         + (1.D0-dxf)*(fun(icu,jcu+1)*dzf+fun(icu,jcu)*(1.0D0-dzf))
  END FUNCTION bplina
  !**************************************************************************
  SUBROUTINE sp_properties
    USE Trivial, ONLY: cmulx,cmuly,cmulz
    INTEGER :: nst,ix,iy,iz,is,ixx,iyy,izz
    COMPLEX(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: pst,psx,psy,psz,psw,ps2
    REAL(db) ::rp,ip,xx(nx),yy(ny),zz(nz),cc(3),ss(3),kin,xpar
    ALLOCATE(pst(nx,ny,nz,2),psx(nx,ny,nz,2),psy(nx,ny,nz,2),psz(nx,ny,nz,2),&
             psw(nx,ny,nz,2),ps2(nx,ny,nz,2))
    sp_orbital=0.D0
    sp_spin=0.D0
    sp_kinetic=0.0D0  
    sp_parity=0.0D0

    xx=x-cmtot(1)
    yy=y-cmtot(2)
    zz=z-cmtot(3)
    !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)&
    !$OMP PRIVATE(nst,ix,iy,iz,is,ixx,iyy,izz,pst,psx,psy,psz,psw,ps2,rp,ip,cc,ss,kin,xpar)
    DO nst=1,nstmax
       IF(node(nst)/=mpi_myproc) CYCLE
       pst=psi(:,:,:,:,localindex(nst))
       IF(TFFT) THEN
        CALL cdervx(psi(:,:,:,:,localindex(nst)),psx,psw)  
        CALL cdervy(psi(:,:,:,:,localindex(nst)),psy,ps2)
        psw=psw+ps2  
        CALL cdervz(psi(:,:,:,:,localindex(nst)),psz,ps2)  
        psw=psw+ps2  
       ELSE  
          CALL cmulx(der1x,pst,psx,0)  
          CALL cmuly(der1y,pst,psy,0)  
          CALL cmulz(der1z,pst,psz,0)  
          CALL cmulx(der2x,pst,psw,0)  
          CALL cmuly(der2y,pst,psw,1)  
          CALL cmulz(der2z,pst,psw,1)  
       ENDIF
       cc=0.D0
       ss=0.D0
       kin=0.D0
       xpar=0.D0
       DO iz=1,nz  
          izz=nz-iz+1  
          DO iy=1,ny  
             iyy=ny-iy+1  
             DO ix=1,nx  
                ixx=nx-ix+1  
                DO is=1,2  
                   rp=REAL(pst(ix,iy,iz,is))
                   ip=AIMAG(pst(ix,iy,iz,is))
                   cc(1)=cc(1)+ &
                        rp*(yy(iy)*AIMAG(psz(ix,iy,iz,is))-zz(iz)*AIMAG(psy(ix,iy,iz,is))) &
                        +ip*(zz(iz)*REAL(psy(ix,iy,iz,is))-yy(iy)*REAL(psz(ix,iy,iz,is)))
                   cc(2)=cc(2)+ &
                        rp*(zz(iz)*AIMAG(psx(ix,iy,iz,is))-xx(ix)*AIMAG(psz(ix,iy,iz,is))) &
                        +ip*(xx(ix)*REAL(psz(ix,iy,iz,is))-zz(iz)*REAL(psx(ix,iy,iz,is)))
                   cc(3)=cc(3)+ &
                        rp*(xx(ix)*AIMAG(psy(ix,iy,iz,is))-yy(iy)*AIMAG(psx(ix,iy,iz,is))) &
                        +ip*(yy(iy)*REAL(psx(ix,iy,iz,is))-xx(ix)*REAL(psy(ix,iy,iz,is)))
                   kin=kin-rp*REAL(psw(ix,iy,iz,is))-ip*AIMAG(psw(ix,iy,iz,is))
                   xpar=xpar+REAL(pst(ix,iy,iz,is))*REAL(pst(ixx,iyy,izz,is)) &
                        +AIMAG(pst(ix,iy,iz,is))*AIMAG(pst(ixx,iyy,izz,is))
                END DO
                ss(1)=ss(1) + REAL(CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,2)) &
                     + REAL(CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,1))
                ss(2)=ss(2) + REAL(CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,2)*CMPLX(0.D0,-1.D0,db)) &
                     + REAL(CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,1)*CMPLX(0.D0,1.D0,db))
                ss(3)=ss(3) + REAL(CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,1)) &
                     - REAL(CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,2))
             ENDDO
          ENDDO
       END DO
       sp_spin(:,nst)=0.5D0*wxyz*ss(:)
       sp_orbital(:,nst)=wxyz*cc
       sp_kinetic(nst)=wxyz*f%h2m(isospin(nst))*kin
       sp_parity(nst)=wxyz*xpar
    END DO
    !$OMP END PARALLEL DO
    DEALLOCATE(pst,psx,psy,psz,psw)
  END SUBROUTINE sp_properties
  !************************************************************
  SUBROUTINE start_protocol(filename,header)
    ! if the protocol file exists, do nothing, since later writes will be
    ! appended. Otherwise write the title and header lines into the new file
    CHARACTER(*),INTENT(IN) :: filename,header
    LOGICAL :: exists
    INQUIRE(FILE=filename,EXIST=exists)
    IF(exists.AND.trestart) RETURN
    OPEN(UNIT=scratch,FILE=filename,FORM='FORMATTED',STATUS='REPLACE')
    WRITE(scratch,'(A)') header
    CLOSE(UNIT=scratch)
  END SUBROUTINE start_protocol
END MODULE Inout
