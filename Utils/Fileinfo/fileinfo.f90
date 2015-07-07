PROGRAM Fileinfo
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
  INTEGER :: iter,nx,ny,nz,nstmax,nneut,nprot,number(2),npsi(2)
  LOGICAL :: vector,isospin
  REAL(db) :: time,dx,dy,dz,charge_number,mass_number,cm(3),wxyz
  CHARACTER :: filename*11,forcename*8
  CHARACTER :: stored_name*10,vecind*6,isoind*17
1 WRITE(*,*) 'Please enter a file name-'
  READ(*,'(A80)',END=4) filename
  WRITE(*,*) 'Opening file ',filename
  OPEN(UNIT=10,FILE=filename,FORM='UNFORMATTED')
! Try to read as wave function file first
  READ(10,ERR=3,END=3) iter,time,forcename,nstmax,nneut,nprot,number,npsi, &
         charge_number,mass_number,cm  
  READ(10) nx,ny,nz,dx,dy,dz,wxyz
  WRITE(*,'(3A,I9,A,F10.2)') 'File ',filename,'is a wave function file at iteration', &
       iter,' or time ',time
  WRITE(*,'(3(A,I6))') 'Total s.p. states: ',nstmax,' Neutrons:',nneut,' Protons:',nprot
  WRITE(*,'(2(A,2I6))') 'Occupied (n/p)',number,' npsi:',npsi
  WRITE(*,*) 'Skyrme force used: ',forcename
  WRITE(*,'(A,3F10.4)') 'Center-of-mass vector:',cm
  CLOSE(10)
  GOTO 1
3 REWIND(10)
  READ(10) iter,time,nx,ny,nz
  READ(10) dx,dy,dz
  WRITE(*,'(A,I9,A,F10.2)') 'Iteration number ',iter,' Physical time ',time
  WRITE(*,'(A,3I5,A,3F8.3)') ' Dimensions ',nx,ny,nz,' Spacings ',dx,dy,dz
  WRITE(*,*) 'Fields stored in this file:'
2 READ(10,END=6) stored_name,vector,isospin
  IF(vector) THEN
    vecind='vector'
  ELSE
    vecind='scalar'
  END IF
  IF(isospin) THEN
    isoind='isospin separated'
  ELSE
    isoind='isospin summed   '
  END IF
  WRITE(*,'(A,"(",A,",",A,")")') stored_name,vecind,isoind
  READ(10) ! ignore actual data
  GOTO 2
6 CLOSE(10)
  GOTO 1
4 STOP
5 WRITE(*,*) 'Could not open file ',filename
END PROGRAM Fileinfo

