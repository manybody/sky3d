PROGRAM tdhf3d  
  USE Params
  USE Fourier
  USE Forces, ONLY: read_force
  USE Densities, ONLY: alloc_densities
  USE Meanfield, ONLY: alloc_fields
  USE Levels
  USE Grids, ONLY: init_grid
  USE Fragments
  USE Parallel
  USE Dynamic, ONLY: getin_dynamic,dynamichf
  USE Static, ONLY: getin_static,init_static,statichf,harmosc
  USE Coulomb, ONLY: coulinit
  USE User
  IMPLICIT NONE
  INTEGER :: imode,nofsave
  !***********************************************************************
  NAMELIST /files/ wffile,converfile,monopolesfile,dipolesfile, &
       momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile
  NAMELIST /main/ tcoul,mprint,mplot,trestart, &
       writeselect,write_isospin,mrest,imode,tfft,nof,r0
  !********************************************************************
  ! Step 1: filename definitions
  !********************************************************************
  CALL init_all_mpi
  OPEN(unit=05,file='for005',status='old',form='formatted')
  READ(5,files)
  !********************************************************************
  ! Step 2: read force definition and determine force
  !********************************************************************
  CALL read_force
  !********************************************************************
  ! Step 3: read principal parameters
  !********************************************************************
  READ(5,main)
  tstatic=imode==1
  tdynamic=imode==2
  IF(tmpi.AND.tstatic) STOP 'MPI not implemented for static mode'
  IF(.NOT.(tstatic.OR.tdynamic)) THEN
     IF(wflag) WRITE(*,*) 'Illegal value for imode:',imode
     STOP
  END IF
  IF(wflag) THEN
     WRITE(*,*) '***** Main parameter input *****'
     IF(tstatic) WRITE(*,*) 'This is a static calculation'
     IF(tdynamic) WRITE(*,*) 'This is a dynamic calculation'
     WRITE(*,'(3(A16,I5))') 'Print interval:',mprint, &
          'Plot interval:',mplot,'Save interval:',mrest
     WRITE(*,'(A,F7.4)') 'Radius parameter r0=',r0
     IF(trestart) THEN
        WRITE(*,*) 'The calculation is being restarted from file ',wffile
     ELSE
        SELECT CASE(nof)
        CASE(0)
           IF(tdynamic) STOP &
                'Harmonic oscillator initialization not allowed for dynamic case'
           WRITE(*,*) 'Harmonic oscillator initialization'
        CASE(1:)
           WRITE(*,'(A,I4,A)') ' Initialization from ',nof,' fragments'
        CASE DEFAULT
           WRITE(*,*) 'User initialization'
        END SELECT
     END IF
     IF(tcoul) THEN
        WRITE(*,*) "Coulomb field is included"
     ELSE
        WRITE(*,*) " Coulomb field *not* included"
     END IF
     WRITE(*,*) "Field output selection: ",writeselect
  END IF
  IF(trestart) THEN
     nofsave=nof ! save whether initial was 2-body
     nof=1 ! restart simulated as fragment input
  END IF
  IF(.NOT.(tstatic.OR.tdynamic)) STOP 'Illegal value of imode'
  !********************************************************************
  ! Step 4: initialize grid, density and mean-field arrays, and FFTW
  !********************************************************************
  CALL init_grid
  CALL alloc_densities
  CALL alloc_fields
  CALL init_fft
  !********************************************************************
  ! Step 5: get input for static or dynamic case
  !********************************************************************
  IF(tstatic) THEN
     CALL getin_static
  ELSE
     CALL getin_dynamic
  END IF
  !********************************************************************
  ! Step 6: get input for static or dynamic case
  !********************************************************************
  ! Determine wave function numbers etc.
  IF(nof>mnof) THEN
     STOP 'nof > mnof'
  ELSE IF(nof>0) THEN
     CALL getin_fragments
  ELSE
     IF(nof==0.AND.tdynamic) STOP 'Dynamic case with nof==0 not allowed'
     npmin(1)=1
     npmin(2)=npsi(1)+1
     npsi(2)=npsi(2)+npsi(1)
     nstmax =npsi(2)
  END IF
  !********************************************************************
  ! Step 7: allocate wave functions
  !********************************************************************
  CALL alloc_nodes
  CALL associate_nodes
  CALL alloc_levels
  !********************************************************************
  ! Step 8: initialize wave functions
  !********************************************************************
  IF(nof>0) THEN
     CALL read_fragments
     IF(.NOT.(tmpi.OR.trestart)) THEN
        CALL schmid
        WRITE(*,*) 'Reorthogonalization complete'
     END IF
  ELSEIF(nof==0) THEN    
     CALL harmosc
  ELSE
     CALL init_user
  END IF
  CLOSE(5)
  !********************************************************************
  ! Step 9: Coulomb initialization
  !********************************************************************
  IF(tcoul) CALL coulinit
  !********************************************************************
  ! Step 10: static or dynamic  calculation performed
  !********************************************************************
  IF(tstatic) THEN
     IF(tmpi .AND. wflag) STOP ' static should not run parallel'
     CALL init_static
     CALL statichf
  ELSE
     !*************************************************************************
     ! Dynamic branch
     !*************************************************************************
     IF(tmpi) CALL mpi_barrier (mpi_comm_world, mpi_ierror)
     IF(trestart) nof=nofsave ! restore 2-body status so analysis is done
     CALL dynamichf
  ENDIF
  CALL finish_mpi
END PROGRAM tdhf3d
