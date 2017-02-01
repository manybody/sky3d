!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the procedures for binary input and output 
!!of the larger fields. 
!>
!>@details 
!!There are two variants, both of which are written
!!at regular intervals: the wave function file \c wffile and the
!!files containing the densities and currents. Since the former is
!!extremely space-consuming, each output normally overwrites the
!!previous one. These files are intended to be used for a restart, as
!!initialization input (static solution for one fragment) for another
!!run, or for a final analysis of the wave functions.
!!
!!The densities, on the other hand, are written on a series of file <tt>
!!nnnnnn.tdd</tt>, where \c nnnnnn is the number of the time step or
!!iteration. This is useful for later graphical or other types of
!!analysis.
!!
!!\b Note: where the variable name \c wffile is used inside file
!!names in the following, it should not be taken literally but is
!!replaced by the character string it contains.
!!
!!In addition the routine for printer plots, \c plot_density, is
!!included in this module, as well as the subroutines \c sp_properties
!!and \c start_protocol, which do not completely
!!match the purpose of this module but are placed here for convenience.
!------------------------------------------------------------------------------
MODULE Inout
  USE Params
  USE Parallel, ONLY: node,localindex,mpi_myproc
  USE Grids
  USE Forces, ONLY:f
  USE Moment, ONLY: cm,cmtot
  USE Densities, ONLY: rho,tau,current,sdens,sodens
  USE Meanfield, ONLY: upot
  USE Coulomb, ONLY: wcoul
  USE Levels
  IMPLICIT NONE
  PRIVATE :: write_one_density,write_vec_density
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: write_wavefunctions
!> @brief
!!This subroutine writes the wave functions to the disk and has two modes of
!!operation depending on whether the code runs on distributed-memory
!!systems in \c MPI mode or on a shared-memory or single-processor
!!machine. 
!>
!> @details
!!In both cases it first determines the number of filled
!!single-particle states <tt> number(iq)</tt>, which need not be the same as
!!either the number of particles or the number of states, since pairing
!!may lead to partial occupation and in addition there can be empty
!!states.
!!
!!  - <b> Sequential operation:</b>
!!    this case is recognized recognized by <tt> mpi_nprocs==1</tt>. Open 
!!    \c wffile, then write four records  containing general information.
!!    - <em>Record 1:</em> <tt> iter, nstmax, nneut, nprot, number, npsi, 
!!      charge_number, mass_number, cm.</tt>
!!    - <em>Record 2:</em> <tt> nx, ny, nz, dx, dy, dz, wxyz.</tt>
!!    - <em>Record 3:</em> <tt> x, y, z.</tt>
!!    - <em>Record 4:</em> <tt> wocc, sp\_energy, sp_parity,
!!      sp_norm, sp_kinetic, sp_efluct1.</tt>
!!    .
!!    These are followed by one record containing information for the \c
!!    MPI case, which is included here only for compatibility:
!!    \c node, \c localindex.  This is then followed by a series of
!!    \c nstloc records (in the sequential case, \c nstloc equals
!!    \c nstmax), containing the array of <tt> nx*ny*nz*2</tt> wave
!!    function values for each single-particle state (including spin).
!!
!!  - <b> MPI operation:</b> in this case processor #0 writes the same
!!    general data as in the sequential case onto file \c wffile, which
!!    is then closed. The purpose of record 5 in this case is to record
!!    for each wave function (in global index space) which node it is in
!!    and what the index on that node is. Since each node produces a
!!    separate output file with only its wave functions, this allows reading
!!    any wave function correctly from the set of files.
!!
!!    Each processor thus only writes the wave function data for its
!!    locally stored set of \c nstloc wave functions onto files with
!!    the names composed (in variable \c rsf}) of the number of the
!!    processor and \c wffile in the form \c nnn.wffile. For
!!    example, if \c wffile has the value 'Ca40', these files will be
!!    \c 000.Ca40, \c 001.Ca40, \c 002.Ca40, etc. up to the number of processors.
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: write_densities
!> @brief
!!This subroutine produces a file \c iter.tdd with density and
!!current data for the present time step or iteration with number \c iter. 
!!In the file name \c iter is given in 6 decimal digits. 
!>
!> @details
!!The record structure is as follows:
!!
!!  - <em>Record 1:</em> this contains the variables <tt> iter, time, nx, 
!!    ny, and nz to define the dimensions of the fields.</tt>
!!  - <em>Record 2:</em> contains the variables <tt> dx, dy, dz, wxyz, x, 
!!    y, and z </tt> to allow proper labelling of axes in plots, etc.
!!  - <em>Further records:</em> for each field to be written, a record is
!!    produced with the following information:
!!    -# Name of the field with up to 10 characters
!!    -# Logical value \c scalar to indicate whether it is a scalar
!!       (\c .FALSE.) or a vector field (\c .TRUE.).
!!    -# Logical value \c write_isospin to indicate whether the
!!       field is summed over protons and neutrons ((\c .FALSE. or not
!!       (\c .TRUE.). In the latter case the field has a last index
!!       running from 1 to 2 for neutrons and protons, respectively.  <b>
!!       This selection applies to all fields equally (except the Coulomb
!!       potential)</b>.
!!    .
!!    After this identification record, the corresponding field itself is
!!    written. 
!!    The dimension varies in the following way: 
!!<table>
!!<caption id="multi_row">Dimensions of arrays</caption>
!!<tr><th>scalar<th>write_isospin<th>dimension
!!<tr><td>.FALSE.<td>.FALSE.<td>(nx,ny,nz)
!!<tr><td>.FALSE.<td>.TRUE.<td>(nx,ny,nz,2)
!!<tr><td>.TRUE.<td>.FALSE.<td>(nx,ny,nz,3)
!!<tr><td>.TRUE.<td>.TRUE.<td>(nx,ny,nz,3,2)
!!</table>
!!The <b> selection of fields to be output </b> is handled through variable
!!\c writeselect consisting of \c nselect characters. Each field
!!is selected by a one-character code, where both lower and upper case
!!are acceptable. At present the choices are:
!!  -<b>R</b>: density \c rho (scalar). Name \c Rho.
!!  -<b>T</b>: kinetic energy density \c tau (scalar). Name \c Tau
!!  -<b>U</b>: local mean field \c upot. Name \c Upot.
!!  -<b>W</b>: Coulomb potential \c wcoul (scalar). This has to be
!!   handled specially, since it has no isospin index. Name \c Wcoul.
!!  -<b>C</b>: current density \c current (vector). Name \c Current.
!!  -<b>S</b>: spin density \c sdens (vector). Name \c Spindens.
!!  -<b>O</b>: spin-orbit density \c sodens. Name \c s-o-Dens.
!!
!!This system is set up to be easily modified for writing additional fields.
!--------------------------------------------------------------------------- 
  SUBROUTINE write_densities
    CHARACTER(10) :: filename
    CHARACTER(1) :: c
    INTEGER :: i
    IF(.NOT.wflag) RETURN
    WRITE(filename,'(I6.6,A4)') iter,'.tdd'
    IF(.NOT.tdynamic) time=0.D0
    OPEN(UNIT=scratch,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE')
    WRITE(scratch) iter,time,nx,ny,nz
    WRITE(scratch) dx,dy,dz,wxyz,x,y,z
    DO i=1,nselect
       c=writeselect(i:i)
       SELECT CASE(c)
       CASE('r','R')
          CALL write_one_density('Rho',rho)
       CASE('t','T')
          CALL write_one_density('Tau',tau)
       CASE('u','U')
          CALL write_one_density('Upot',upot)
       CASE('w','W')
          WRITE(scratch) 'Wcoul     ',.FALSE.,.FALSE.
          WRITE(scratch) wcoul
       CASE('c','C')
          CALL write_vec_density('Current',current)
       CASE('s','S')
          CALL write_vec_density('Spindens',sdens)
       CASE('o','O')
          CALL write_vec_density('s-o-Dens',sodens)
       END SELECT
    END DO
    CLOSE(UNIT=scratch)
  END SUBROUTINE write_densities
!---------------------------------------------------------------------------  
! DESCRIPTION: write_one_density
!> @brief
!!This subroutines does the actual output for \c write_densities in
!!the case of a scalar field. Its functioning should be clear from the
!!description above.
!> @param[in] name
!> INTEGER, takes the name of the density.
!> @param[in] values
!> REAL(db), takes the density.
!--------------------------------------------------------------------------- 
  SUBROUTINE write_one_density(name,values)
    CHARACTER(*),INTENT(IN) :: name
    REAL(db),INTENT(IN) :: values(nx,ny,nz,2)
    CHARACTER(10) :: stored_name
    REAL(db) a(nx,ny,nz)
    stored_name=name
    WRITE(scratch) stored_name,.FALSE.,write_isospin
    IF(write_isospin) THEN
       WRITE(scratch) values
    ELSE
       a=values(:,:,:,1)+values(:,:,:,2)
       WRITE(scratch) a
    END IF
  END SUBROUTINE write_one_density
!---------------------------------------------------------------------------  
! DESCRIPTION: write_vec_density
!> @brief
!!This also does the actual output for subroutines \c write_densities 
!!for the case of a vector field. Its functioning should be clear from the
!!description above.
!> @param[in] name
!> INTEGER, takes the name of the density.
!> @param[in] values
!> REAL(db), takes the vector density.
!--------------------------------------------------------------------------- 
  SUBROUTINE write_vec_density(name,values)
    CHARACTER(*),INTENT(IN) :: name
    REAL(db),INTENT(IN) :: values(nx,ny,nz,3,2)
    CHARACTER(10) :: stored_name
    REAL(db) a(nx,ny,nz,3)
    stored_name=name
    WRITE(scratch) stored_name,.TRUE.,write_isospin
    IF(write_isospin) THEN
       WRITE(scratch) values
    ELSE
       a=values(:,:,:,:,1)+values(:,:,:,:,2)
       WRITE(scratch) a
    END IF
  END SUBROUTINE write_vec_density
!---------------------------------------------------------------------------  
! DESCRIPTION: plot_density
!> @brief
!!Produces a simple printer plot of the density distribution in the
!!reaction plane. This is not supposed to replace better plotting codes,
!!but simply allows a quick glance at what is happening in the code,
!!even while it is running.
!>
!> @details
!!It is based on a very old routine found at ORNL and was translated
!!into modern Fortran. It uses helper function \c bplina for
!!interpolation.  
!--------------------------------------------------------------------------- 
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
    nhor=dimx/dxp+1; nver=dimz/dzp+1
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
!---------------------------------------------------------------------------  
! DESCRIPTION: bplina
!> @brief
!!Does a bilinear interpolation of fun(nx,nz), on xar(n) and zar(m)
!>
!> @param[in] n
!> INTEGER, takes the number of grid points in x-direction.
!> @param[in] m
!> INTEGER, takes the number of grid points in z-direction.
!> @param[in] xar
!> REAL(db), array, takes coordinates in x direction.
!> @param[in] zar
!> REAL(db), array, takes coordinates in z direction.
!> @param[in] fun
!> REAL(db), array, takes the function.
!> @param[in] xcu
!> REAL(db), takes the x-value at which the interpolation is performed
!> @param[in] zcu
!> REAL(db), takes the z-value at which the interpolation is performed  
!--------------------------------------------------------------------------- 
  PURE FUNCTION bplina(n,m,xar,zar,fun,xcu,zcu) RESULT(ff)
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
!---------------------------------------------------------------------------  
! DESCRIPTION: sp_properties
!> @brief
!!In this routine the kinetic energy, orbital and spin angular momenta
!!expectation values, \c sp_kinetic, \c sp_orbital} and
!!\c sp_spin of the single-particle states are calculated. The
!!latter are both three-dimensional vectors.
!>
!> @details
!!Note that the single-particle energy \c sp_energy itself is not
!!calculated here but in the main static and dynamic routines, since it
!!is obtained by applying the single-particle Hamiltonian, which is
!!done more conveniently there.
!!
!!The procedure is quite simple: in a loop over wave functions the
!!active one is copied into \c pst for convenience. Then its three
!!directional derivatives \c psx, \c psy, and \c psz and
!!Laplacian \c psw are calculated. In the big loop over the grid they
!!are combined to the desired matrix elements; the only technical point
!!to remark is that since the result must be real, efficiency can be
!!achieved by formulating the complex products in an explicit way. Then
!!\c kin contains the kinetic energy (without the \f$ \hbar^2/2m \f$), 
!!\c cc the orbital and \c ss then spin matrix elements.
!!
!!Finally only the volume element, the factor of one half for the
!!spin ad the prefactor of the kinetic energy are added. 
!--------------------------------------------------------------------------- 
  SUBROUTINE sp_properties
    USE Trivial, ONLY: cmulx,cmuly,cmulz
    INTEGER :: nst,ix,iy,iz,is,ixx,iyy,izz
    COMPLEX(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: pst,psx,psy,psz,psw
    REAL(db) ::rp,ip,xx(nx),yy(ny),zz(nz),cc(3),ss(3),kin,xpar
    ALLOCATE(pst(nx,ny,nz,2),psx(nx,ny,nz,2),psy(nx,ny,nz,2),psz(nx,ny,nz,2),psw(nx,ny,nz,2))
    sp_orbital=0.D0
    sp_spin=0.D0
    sp_kinetic=0.0D0  
    sp_parity=0.0D0

    xx=x-cmtot(1)
    yy=y-cmtot(2)
    zz=z-cmtot(3)
    DO nst=1,nstmax
       IF(node(nst)/=mpi_myproc) CYCLE
       pst=psi(:,:,:,:,localindex(nst))
       IF(TFFT) THEN
          CALL cdervx(pst,psx)  
          CALL cdervy(pst,psy)  

          CALL cdervz(pst,psz)  
          CALL laplace(pst,psw)  
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
                ss(1)=ss(1) + CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,2) &
                     + CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,1)
                ss(2)=ss(2) + CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,2)*CMPLX(0.D0,-1.D0,db) &
                     + CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,1)*CMPLX(0.D0,1.D0,db)
                ss(3)=ss(3) + CONJG(pst(ix,iy,iz,1))*pst(ix,iy,iz,1) &
                     - CONJG(pst(ix,iy,iz,2))*pst(ix,iy,iz,2)
             ENDDO
          ENDDO
       END DO
       sp_spin(:,nst)=0.5D0*wxyz*ss(:)
       sp_orbital(:,nst)=wxyz*cc
       sp_kinetic(nst)=wxyz*f%h2m(isospin(nst))*kin
       sp_parity(nst)=wxyz*xpar
    END DO
    DEALLOCATE(pst,psx,psy,psz,psw)
  END SUBROUTINE sp_properties
!---------------------------------------------------------------------------  
! DESCRIPTION: Routinename
!> @brief
!!This is given a file name and a character
!!string for a header line to start the file contents. It is used for
!!the <tt>*.res</tt> files.  If the file already exists, nothing is done,
!!since this probably a restart job and output should just be added at
!!the end of the file. 
!> @param[in] filename
!> CHARACTER, array, takes the filename.
!> @param[in] header
!> CHARACTER, array, takes the intended header for the file.
!--------------------------------------------------------------------------- 
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
