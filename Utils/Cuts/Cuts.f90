PROGRAM Tdhf2Silo
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  ! define data to be read from files
  INTEGER :: iter,nx,ny,nz
  LOGICAL :: vector,isospin
  REAL(db) :: time,dx,dy,dz,wxyz
  REAL(db),ALLOCATABLE :: x(:),y(:),z(:),rho(:,:,:,:),vec(:,:,:,:,:)
  CHARACTER(13) :: filename
  CHARACTER :: stored_name*10
  ! Store the length of each string
1 READ(*,'(A10)',END=4) filename
  WRITE(*,*) 'Converting file ',filename
  OPEN(UNIT=10,FILE=filename,FORM='UNFORMATTED')
  filename(10:13)='.txt'
  !
  ! read general information about grid, time, etc., allocate arrays
  !
  READ(10) iter,time,nx,ny,nz
  ALLOCATE(x(nx),y(ny),z(nz),rho(nx,ny,nz,2),vec(nx,ny,nz,3,2))
  READ(10) dx,dy,dz,wxyz,x,y,z
  !
  ! now start reading various densities
2 READ(10,END=3) stored_name,vector,isospin
  IF(TRIM(stored_name)=='Rho') THEN
     IF(isospin) THEN
        READ(10) rho
        rho(:,:,:,1)=rho(:,:,:,1)+rho(:,:,:,2)
     ELSE
        READ(10) rho(:,:,:,1)
     END IF
     CALL writeit(0.5D0*(rho(:,:,nz/2,1)+rho(:,:,nz/2+1,1)),x,y,nx,ny,'rxy')
     CALL writeit(0.5D0*(rho(:,ny/2,:,1)+rho(:,ny/2+1,:,1)),x,z,nx,nz,'rxz')
     CALL writeit(0.5D0*(rho(nx/2,:,:,1)+rho(nx/2+1,:,:,1)),y,z,ny,nz,'ryz')
  END IF
  GOTO 2
3 DEALLOCATE(x,y,z,rho,vec)
  GOTO 1
4 STOP
  !*******************************************************************
CONTAINS
  SUBROUTINE writeit(a,b,c,n,m,name)
    REAL(db) :: a(:,:),b(:),c(:)
    INTEGER :: n,m
    CHARACTER(3) :: name
    INTENT(IN) :: a,b,c,n,m,name
    INTEGER :: i,j
    filename(7:9)=name
    WRITE(*,*) 'Writing ',filename
    OPEN(UNIT=11,FILE=filename,FORM='FORMATTED',STATUS='REPLACE')
    DO j=1,m
       DO i=1,n 
          WRITE(11,'(2F8.3,E14.6)') b(i),c(j),a(i,j)
       END DO
       WRITE(11,'(1X)')
    ENDDO
    CLOSE(UNIT=11)
  END SUBROUTINE writeit
END PROGRAM Tdhf2Silo

