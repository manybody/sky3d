PROGRAM Inertia
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  INTEGER :: iter,nx,ny,nz,i,j,k
  LOGICAL :: vector,isospin
  REAL(db) :: time,dx,dy,dz,wxyz,mass,cm(3),tensor(3,3)
  REAL(db),ALLOCATABLE :: x(:),y(:),z(:),rho(:,:,:,:)
  CHARACTER(11) :: filename
  CHARACTER :: stored_name*10
1 READ(*,'(A10)',END=4) filename
  OPEN(UNIT=10,FILE=filename,FORM='UNFORMATTED')
  READ(10) iter,time,nx,ny,nz
  ALLOCATE(x(nx),y(ny),z(nz),rho(nx,ny,nz,2))
  READ(10) dx,dy,dz,wxyz,x,y,z
  WRITE(*,'(A,I9,A,F10.2)') 'Iteration number ',iter,' Physical time ',time
! search for a field named "rho" and not a vector
2 READ(10,END=1) stored_name,vector,isospin
  IF(stored_name/='Rho'.OR.vector) GOTO 2
! If p and n are stored separately, add up to toal density
  IF(isospin) THEN
    READ(10) rho
    rho(:,:,:,1)=rho(:,:,:,1)+rho(:,:,:,2)
  ELSE
    READ(10) rho(:,:,:,1)
  END IF
! calculate mass and c.m. vector
  mass=SUM(rho(:,:,:,1))*wxyz
  cm=0.D0
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        cm(1)=cm(1)+rho(i,j,k,1)*x(i)
        cm(2)=cm(2)+rho(i,j,k,1)*y(j)
        cm(3)=cm(3)+rho(i,j,k,1)*z(k)
      END DO
    END DO
  END DO
  cm=cm/mass
! correct coordinates
  x=x-cm(1)
  y=y-cm(2)
  z=z-cm(3)
! now do the real calculation
  tensor=0.D0
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tensor(1,1)=tensor(1,1)+rho(i,j,k,1)*(y(j)**2+z(k)**2)
        tensor(1,2)=tensor(1,2)+rho(i,j,k,1)*x(i)*y(j)
        tensor(1,3)=tensor(1,3)+rho(i,j,k,1)*x(i)*z(k)
        tensor(2,2)=tensor(2,2)+rho(i,j,k,1)*(x(i)**2+z(k)**2)
        tensor(2,3)=tensor(2,3)+rho(i,j,k,1)*y(j)*z(k)
        tensor(3,3)=tensor(3,3)+rho(i,j,k,1)*(x(i)**2+y(j)**2)
      END DO
    END DO
  END DO
  tensor(2,1)=tensor(1,2)
  tensor(3,1)=tensor(1,3)
  tensor(3,2)=tensor(2,3)
  tensor=tensor*wxyz
  WRITE(*,'(("(",F12.3,",",F12.3,",",F12.3,")"))') tensor
  DEALLOCATE(x,y,z,rho)
  CLOSE(10)
  GOTO 1
4 STOP
5 WRITE(*,*) 'Could not open file ',filename
END PROGRAM Inertia
