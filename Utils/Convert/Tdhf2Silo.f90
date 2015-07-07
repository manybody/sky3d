PROGRAM Tdhf2Silo
  IMPLICIT NONE
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100), &
       maxdefs=20
  ! define data to be read from files
  INTEGER :: iter,nx,ny,nz
  LOGICAL :: vector,isospin
  REAL(db) :: time,dx,dy,dz,wxyz
  REAL(db),ALLOCATABLE :: x(:),y(:),z(:),rho(:,:,:,:),vec(:,:,:,:,:)
  INCLUDE 'silo.inc'
  CHARACTER(11) :: filename
  CHARACTER :: stored_name*10
  INTEGER :: is,ret,iret,dims(3),dbid,op,i,ndefs
  INTEGER err, ierr, types(maxdefs), lnames(maxdefs), ldefs(maxdefs)
  INTEGER oldlen
  ! Initialize some 20 character length strings
  CHARACTER(40) :: names(maxdefs),defs(maxdefs)
  ! Store the length of each string
1 READ(*,'(A10)',END=4) filename
  WRITE(*,*) 'Converting file ',filename
  OPEN(UNIT=10,FILE=filename,FORM='UNFORMATTED')
  filename(7:11)='.silo'
  ret=dbcreate(filename,11,DB_CLOBBER,DB_LOCAL,'TDHF DATA',9&
       &,DB_PDB,dbid)
  CALL errchk(ret,' could not open')
  !
  ! read general information about grid, time, etc., allocate arrays
  !
  READ(10) iter,time,nx,ny,nz
  ALLOCATE(x(nx),y(ny),z(nz),rho(nx,ny,nz,2),vec(nx,ny,nz,3,2))
  READ(10) dx,dy,dz,wxyz,x,y,z
  !
  ! put this into silo file
  !
  ! set dimensions for following zone-centered arrays
  dims=(/nx,ny,nz/)
  ! make options to add cycle number and time
  ret=dbmkoptlist(5,op)
  CALL errchk(ret,' could not make optlist')
  ret=dbadddopt(op,DBOPT_DTIME,time)
  CALL errchk(ret,' could not add time parameter')
  ret=dbaddiopt(op,DBOPT_CYCLE,iter)
  CALL errchk(ret,' could not add ncyc parameter')
  ret=dbaddcopt(op,DBOPT_XUNITS,'fm',2)
  CALL errchk(ret,' could not add xunits parameter')
  ret=dbaddcopt(op,DBOPT_YUNITS,'fm',2)
  CALL errchk(ret,' could not add yunits parameter')
  ret=dbaddcopt(op,DBOPT_ZUNITS,'fm',2)
  CALL errchk(ret,' could not add zunits parameter')
  ! write regular mesh
  ret=dbputqm(dbid,'Mesh',4,'x',1,'y',1,'z',1, x,y,z,dims,3,DB_DOUBLE&
       &,DB_COLLINEAR,op,iret)
  !
  ! now start reading various densities
  ndefs=0
2 READ(10,END=3) stored_name,vector,isospin
  IF(vector) THEN
    IF(isospin) THEN
      READ(10) vec
      CALL write_density(TRIM(stored_name) // 'nx',vec(:,:,:,1,1))
      CALL write_density(TRIM(stored_name) // 'ny',vec(:,:,:,2,1))
      CALL write_density(TRIM(stored_name) // 'nz',vec(:,:,:,3,1))
      CALL write_density(TRIM(stored_name) // 'px',vec(:,:,:,1,2))
      CALL write_density(TRIM(stored_name) // 'py',vec(:,:,:,2,2))
      CALL write_density(TRIM(stored_name) // 'pz',vec(:,:,:,3,2))
    ELSE
      READ(10) vec(:,:,:,:,1)
      CALL write_density(TRIM(stored_name) // 'x',vec(:,:,:,1,1))
      CALL write_density(TRIM(stored_name) // 'y',vec(:,:,:,2,1))
      CALL write_density(TRIM(stored_name) // 'z',vec(:,:,:,3,1))
      ndefs=ndefs+1
      IF(ndefs>maxdefs) STOP 'Too many vector quantities'
      names(ndefs)=stored_name
      lnames(ndefs)=LEN(stored_name)
      WRITE(defs(ndefs),'("{",A12,",",A12,",",A12,"}")') &
           TRIM(stored_name) // 'x',TRIM(stored_name) // 'y' &
           ,TRIM(stored_name) // 'z'
      ldefs(ndefs)=40
      types(ndefs)=DB_VARTYPE_VECTOR
    END IF
  ELSE
    IF(isospin) THEN
      READ(10) rho
      CALL write_density(TRIM(stored_name) // 'n',rho(:,:,:,1))
      CALL write_density(TRIM(stored_name) // 'p',rho(:,:,:,2))
    ELSE
      READ(10) rho(:,:,:,1)
      CALL write_density(TRIM(stored_name),rho(:,:,:,1))
    END IF
  END IF
  GOTO 2
3 IF(ndefs>0) THEN
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(40)
    err = dbputdefvars(dbid, "defvars", 7, ndefs, &
         names, lnames, types, defs, ldefs, DB_F77NULL, iret)
    CALL errchk(ret,'could not write defvars')
    err = dbset2dstrlen(oldlen)
  ENDIF
  ret=dbfreeoptlist(op)
  CALL errchk(ret,' could not free optlist')
  ret=dbclose(dbid)
  CALL errchk(ret,' could not close')
  DEALLOCATE(x,y,z,rho,vec)
  GOTO 1
4 STOP
  !*******************************************************************
CONTAINS
  SUBROUTINE write_density(rname,data)
    CHARACTER*(*),INTENT(IN) :: rname
    REAL(8),INTENT(IN) :: data(*)
    WRITE(*,*) 'Writing ',rname,'  ',LEN(rname)
    ret=dbputqv1(dbid,rname,LEN(rname),'Mesh',4,data,dims,3, &
         DB_F77NULL,0,DB_DOUBLE,DB_NODECENT,op,iret)
    CALL errchk(ret,'could not write '//rname//' in '//filename)
  END SUBROUTINE write_density
  SUBROUTINE errchk(ret,msg)
    INTEGER,INTENT(IN) :: ret
    CHARACTER(*),INTENT(IN) :: msg
    IF(ret==0) RETURN
    WRITE(6,*) 'Error in SILO: ',msg
    STOP
  END SUBROUTINE errchk
END PROGRAM Tdhf2Silo

