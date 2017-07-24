!------------------------------------------------------------------------------
! Programm: tdhf3d
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This is the main program that organizes the reading of the input and
!!funnels the calculation into the correct subroutines.
!>
!>@details
!!It consists of a number of simple steps. They are:
!! -# initialize \c MPI in case of an \c MPI parallel job.
!!    Start reading from standard input  beginning with 
!!    namelist files for any changes in the file names.
!! -# read the definition of the force to be used and set
!!    it up.
!! -# read namelist main, which contains the
!!    overall controlling parameters. The choice of \c imode is used to
!!    set up \c tstatic and \c tdynamic as logical variables. In
!!    addition, if this is a restart, the number of fragments is set to
!!    one. The properties of this one fragment are then filled in in
!!    subroutine \c getin_fragments.
!! -# read the grid definition.  With the dimensions known, all grid-based
!!    arrays (but not the wave functions) can be allocated and the
!!    initialization of the \c FFTW system can be done.
!! -# the appropriate namelist is read for the static or
!!    dynamic case, defining some parameters in those modules.
!! -# determine wave function numbers. The way this is
!!    done depends on the choice of \c nof. For positive \c nof, the
!!    fragment files are consulted to find out the properties of each and
!!    add up the numbers. For \c nof=0 the numbers are taken from
!!    namelist \c static, which was read before. If \c nof<0, they
!!    are also given in namelist \c static, but the wave functions will
!!    be replaced by user-calculated ones.
!! -# now that the numbers are known, the wave function
!!    distribution over nodes is computed or set to trivial for a
!!    sequential calculation, and the the arrays related with wave
!!    functions are allocated.
!! -# the initial values of the wave functions are
!!    calculated. For \c nof>0 the wave functions are read from the
!!    fragment files and inserted into the proper positions. For \c nof=0 the routine
!!    \c harmosc in module \c Static is called to calculate
!!    harmonic-oscillator wave function, and for \c nof<0 the routine
!!    \c init_user does an arbitrary user initialization.
!! -# the mean-field arrays are zeroed and the Coulomb
!!    solver is initialized.
!! -# the calculation branches into either static or dynamic mode.
!!
!!    In the static calculation, \c init_static just prints some
!!    information and sets up damping  and
!!    \c statichf does the real calculation.
!!
!!    In the dynamic case, the subroutine \c dynamichf does all the
!!    work.
!! -# finally the \c MPI system is terminated.
!------------------------------------------------------------------------------
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
  USE Static, ONLY: getin_static,init_static,statichf,harmosc,planewaves
  USE Coulomb, ONLY: coulinit
  USE User
  IMPLICIT NONE
  INTEGER :: imode,nofsave=0
  !***********************************************************************
  NAMELIST /files/ wffile,converfile,monopolesfile,dipolesfile, &
       momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile,&
       diffenergiesfile
  NAMELIST /main/ tcoul,mprint,mplot,trestart, &
       writeselect,write_isospin,mrest,imode,tfft,nof,r0
  !********************************************************************
  ! Step 1: filename definitions
  !********************************************************************
  CALL init_all_mpi
  OPEN(unit=05,file='for005',status='old',form='formatted')
  READ(5,files)
  CALL mpi_init_filename
  !********************************************************************
  ! Step 2: read force definition and determine force
  !********************************************************************
  CALL read_force
  !********************************************************************
  ! Step 3: read principal parameters
  !********************************************************************
  READ(5,main)
  CALL mpi_init_filename
  tstatic=imode==1
  tdynamic=imode==2
  IF(.NOT.(tstatic.OR.tdynamic)) THEN
     IF(wflag) WRITE(*,*) 'Illegal value for imode:',imode
     STOP
  END IF
  IF(wflag) THEN
     WRITE(*,*)
     WRITE(*,*) '***** Main parameter input *****'
     IF(tstatic) WRITE(*,*) 'This is a static calculation'
     IF(tdynamic) WRITE(*,*) 'This is a dynamic calculation'
     WRITE(*,'(3(A16,I5))') 'Print interval:',mprint, &
          'Plot interval:',mplot,'Save interval:',mrest
     WRITE(*,'(X,A,F7.4)') 'Radius parameter r0=',r0
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
  IF(tmpi.AND.tstatic) CALL init_mpi_2d
  CALL alloc_levels
  !********************************************************************
  ! Step 8: initialize wave functions
  !********************************************************************
  IF(nof>0) THEN
     CALL read_fragments
  ELSEIF(nof==0) THEN    
     CALL harmosc
  ELSEIF(nof==-1) THEN
     CALL planewaves
  ELSE
     CALL init_user
  END IF
!  CLOSE(5)
  !********************************************************************
  ! Step 9: Coulomb initialization
  !********************************************************************
  IF(tcoul) CALL coulinit
  !********************************************************************
  ! Step 10: static or dynamic  calculation performed
  !********************************************************************
  IF(tstatic) THEN
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
!  CALL init_user
  CALL finish_mpi
END PROGRAM tdhf3d
