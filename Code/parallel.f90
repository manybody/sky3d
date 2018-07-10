!------------------------------------------------------------------------------
! MODULE: Parallel
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module organizes the execution on distributed-memory machines
!!using \c MPI. Its interface is made such that on sequential
!!machines \c parallel.f90 can simply be replaced by a special
!!version \c sequential.f90 which contains a definition of this
!!module generating trivial data.
!> 
!> @details
!!\c MPI parallelization is
!!based on distributing the wave functions onto the different nodes.
!!Thus all operations acting directly on the wave functions can be done
!!in parallel, not only the time evolution by the application of the
!!single-particle Hamiltonian, but also the summing up of the densities
!!and currents over those wave function stored on the node. This means
!!that only the final summation of the densities and the calculations
!!done with them have to be communicated across the nodes.
!!
!!It is important that the single-particle properties also defined in
!!\c Levels, e.g. \c sp_energy are not split up up for the nodes
!!but the full set is present on each node. The values are communicated
!!by summing from all nodes with zeroes in those index positions not
!!present on a particular one. This method of handling them avoids having to
!!communicate many small arrays.
!------------------------------------------------------------------------------

MODULE Parallel
  USE Params, ONLY : wflag
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  !>@name Block parameters. 
  !>@{
  INTEGER, PARAMETER   :: NB=2
  INTEGER, PARAMETER   :: MB=2
  INTEGER, PARAMETER   :: NB_psi =2
  !>@}
  LOGICAL,PARAMETER :: tmpi=.TRUE.       !<a logical variable set to true if \c MPI
  !!parallelization is activated. It is used to turn the calling of all
  !!the \c MPI routines in the code on or off.
  LOGICAL, PARAMETER   :: ttabc=.FALSE.
  !>@name Variables characterizing the different distributions of wave functions. 
  !>@{
  INTEGER, ALLOCATABLE :: node(:)        !<For the single-particle state
  !!with index \c i the wave function is stored on computing node \c node(i).
  INTEGER, ALLOCATABLE :: localindex(:)  !<For the single-particle state
  !!with index \c i the wave function is stored on computing node
  !!\c node(i) and its index on that node is \c localindex(i).
  INTEGER, ALLOCATABLE :: globalindex(:) !<tells the index of the single-particle
  !!state in the whole array of \c nstmax states in 1d distribution. (it could be
  !!dimensioned \c nstloc but is dimensioned as \c nstmax 
  !!to make its allocation simpler). So for wave function index \c i 
  !!<em> on the local node</em>, <tt> i=1..nstloc</tt>, the
  !!single-particle energy must be obtained using <tt> sp_energy(globalindex(i))</tt>. 
  INTEGER, ALLOCATABLE :: globalindex_x(:)!<tells the index of the single-particle
  !!state in the whole array of \c nstmax states in 2d distribution for rows.
  INTEGER, ALLOCATABLE :: globalindex_y(:)!<tells the index of the single-particle
  !!state in the whole array of \c nstmax states in 2d distribution for columns.
  INTEGER, ALLOCATABLE :: globalindex_diag_x(:)!<tells the index of the single-particle
  !!state in the whole array of \c nstmax states in splitted 2d distribution for rows.
  INTEGER, ALLOCATABLE :: globalindex_diag_y(:)!<tells the index of the single-particle
  !!state in the whole array of \c nstmax states in splitted 2d distribution for columns.
  INTEGER              :: nstloc_x!<Size of local array for matrices in 2d distribution for rows.
  INTEGER              :: nstloc_y!<Size of local array for matrices in 2d distribution for columns.
  INTEGER              :: psiloc_x!<Size of local array for wave functions in 2d distribution (rows).
  INTEGER              :: psiloc_y!<Size of local array for wave functions in 2d distribution (rows).
  INTEGER              :: nstloc_diag_x!<Size of local array for splitted matrices in 2d distribution for rows.
  INTEGER              :: nstloc_diag_y!<Size of local array for splitted matrices in 2d distribution for columns.
  !>@}
  !>@name Variables that specify computation grids.
  !>@{ 
  INTEGER              :: mpi_nprocs            !<number of MPI processes.
  INTEGER              :: mpi_ierror            !<varable for error output of MPI routines.
  INTEGER              :: mpi_myproc            !<the number of the local MPI process
  INTEGER              :: mpi_nprocs_iso(2)     !<number of MPI processes for each isospin.
  INTEGER              :: my_iso                !<the local isospin that is calculated.
  INTEGER              :: mpi_myproc_iso        !<local MPI process number in the isospin group
  INTEGER              :: mpi_size_x            !<number of MPI ranks in row direction.
  INTEGER              :: mpi_size_y            !<number of MPI ranks in column direction.
  INTEGER              :: mpi_rank_x            !<local MPI process number in row direction.
  INTEGER              :: mpi_rank_y            !<local MPI process number in column direction.
  INTEGER              :: nprd                  !<number of local rows for diagonalization/orthonormalization.
  INTEGER              :: npcd                  !<number of local columns for diagonalization/orthonormalization.
  INTEGER              :: myprowd               !<local row rank for diagonalization/orthonormalization.
  INTEGER              :: mypcold               !<local column rank for diagonalization/orthonormalization.
  INTEGER              :: mpi_dims(2)           !<dimensions for 2d computaion grid.
  INTEGER              :: mpi_mycoords(2)       !<local location in 2d computaion grid
  INTEGER              :: mpi_nprocs_diag(4)    !<number of processes in each diagonalization/orthonormalization group.
  INTEGER              :: my_diag               !<determines the local group for diagonalization/orthonormalization.
  INTEGER              :: NPROCS                !<dummy variable for number of processes.
  INTEGER              :: IAM                   !<dummy variable for local process number.
  INTEGER              :: NPROW                 !<number of local rows for 2d distribution.
  INTEGER              :: NPCOL                 !<number of local columns for 2d distribution.
  INTEGER              :: MYPROW                !<local row rank for 2d distribution.
  INTEGER              :: MYPCOL                !<local column rank for 2d distribution.
  !>@}
  !>@name Communicators 
  !>@{  
  INTEGER              :: comm2d     !<communicator of 2d distribution.
  INTEGER              :: comm_iso   !<communicator for local isospin group.
  INTEGER              :: contxt_d   !<BLACS context for diagonalization group.
  INTEGER              :: contxt_o   !<BLACS context for orthonormalization group.
  INTEGER              :: CONTXT     !<BLACS context for 2d distribution.    
  INTEGER              :: CONTXT1D   !<BLACS context for 1d distribution.
  INTEGER              :: CONTXT_DO  !<BLACS context for local either diagonalization or orthonormalization group.
  !>@}
   
  INTEGER, EXTERNAL    :: NUMROC,INDXL2G,INDXG2L,INDXG2P
  REAL(db)             :: timer(20)
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_nodes
!> @brief This subroutine merely allocates the internal arrays of module \c Parallel.
!---------------------------------------------------------------------------   
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             globalindex_x(nstmax),globalindex_y(nstmax),&
             globalindex_diag_x(nstmax),globalindex_diag_y(nstmax))
  END SUBROUTINE alloc_nodes
!---------------------------------------------------------------------------  
! DESCRIPTION: init_all_mpi
!> @brief
!!This subroutine initializes \c MPI and
!!finds out the number of processors \c mpi_nprocs as well as the
!!index of the current one \c mpi_myproc. The flag \c wflag is
!!set to true only for the processor numbered 0.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_all_mpi
    CALL mpi_init(mpi_ierror)
    CALL mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)
    CALL mpi_comm_rank(mpi_comm_world,mpi_myproc,mpi_ierror)
    wflag=mpi_myproc==0
    IF(wflag) WRITE(*,'(a,i5)') ' number of nodes=',mpi_nprocs
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE init_all_mpi
!---------------------------------------------------------------------------  
! DESCRIPTION: mpi_init_filename
!> @brief
!!If calculated with TABC it sets the filenames for output files.
!---------------------------------------------------------------------------  
  SUBROUTINE mpi_init_filename
    CONTINUE
  END SUBROUTINE mpi_init_filename   
!---------------------------------------------------------------------------  
! DESCRIPTION: tabc_av
!> @brief
!!Averages a value if TABC are used
!---------------------------------------------------------------------------  
  FUNCTION tabc_av(val)
    REAL(db),INTENT(IN) :: val
    REAL(db)            :: tabc_av
    tabc_av=val
  END FUNCTION tabc_av
!---------------------------------------------------------------------------  
! DESCRIPTION: tabc_dens
!> @brief
!!Averages a density if TABC are used
!--------------------------------------------------------------------------- 
  FUNCTION tabc_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,2)
    REAL(db)            :: tabc_dens(nx,ny,nz,2)
    tabc_dens=density
  END FUNCTION tabc_dens
!---------------------------------------------------------------------------  
! DESCRIPTION: tabc_vec_dens
!> @brief
!!Averages a vector density if TABC are used
!---------------------------------------------------------------------------
  FUNCTION tabc_vec_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,3,2)
    REAL(db)            :: tabc_vec_dens(nx,ny,nz,3,2)
    tabc_vec_dens=density
  END FUNCTION tabc_vec_dens
!---------------------------------------------------------------------------  
! DESCRIPTION: tabc_filename
!> @brief
!!Changes a filename for the use of TABC.
!--------------------------------------------------------------------------- 
  FUNCTION tabc_filename(filename)
    CHARACTER(64),INTENT(IN) :: filename
    CHARACTER(64)            :: tabc_filename
    tabc_filename=filename
  END FUNCTION tabc_filename
!---------------------------------------------------------------------------  
! DESCRIPTION: init_mpi_2d
!> @brief
!!This subroutine sets up the 2d caomputation grids with MPI routines. It 
!!determines the best distributions for each group and determines \c globalindex
!!and \localindex for the 2d case.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_mpi_2d
    !***********************************************************************
    ! calculates 2d wave function distribution over nodes
    !***********************************************************************
    LOGICAL :: isperiodic(2),reorder
    INTEGER :: i
    isperiodic=.FALSE.
    reorder=.FALSE.
    
    !Calculate best dimensions for the given number of processes
    CALL mpi_dims_create(mpi_nprocs_iso(my_iso),2,mpi_dims,mpi_ierror)
    !Create 2-dimensional grid of processes to calculate the matrices in diagstep
    CALL mpi_cart_create(comm_iso,2,mpi_dims,isperiodic,reorder,comm2d,mpi_ierror)

    !get my coordinates on the grid
    CALL mpi_cart_get(comm2d,2,mpi_dims,isperiodic,mpi_mycoords,mpi_ierror)
    IF(mpi_myproc_iso==0) WRITE(*,'(A36,I4,A4,I4,A9,I2)')&
                                'Initialized 2d-grid with dimensions ',mpi_dims(1),&
                                ' and',mpi_dims(2),' for iso=',my_iso
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
!
    CALL init_blacs
    nstloc_x = NUMROC(npsi(my_iso)-npmin(my_iso)+1,NB,MYPROW,0,NPROW)
    nstloc_y = NUMROC(npsi(my_iso)-npmin(my_iso)+1,MB,MYPCOL,0,NPCOL)
    psiloc_x = NUMROC(nx*ny*nz*2,NB,MYPROW,0,NPROW)
    psiloc_y = nstloc_y
    nstloc_diag_x = NUMROC(npsi(my_iso)-npmin(my_iso)+1,NB,myprowd,0,nprd)
    nstloc_diag_y = NUMROC(npsi(my_iso)-npmin(my_iso)+1,NB,mypcold,0,npcd)

    DO i=1,nstloc_x
      globalindex_x(i)=INDXL2G(i, NB, MYPROW, 0, NPROW )
    END DO
    DO i=1,nstloc_y
      globalindex_y(i)=INDXL2G(i, MB, MYPCOL, 0, NPCOL )
    END DO
    IF(my_iso==2) THEN
      globalindex_x=globalindex_x+npsi(1)
      globalindex_y=globalindex_y+npsi(1)
    END IF

    DO i=1,nstloc_diag_x
      globalindex_diag_x(i)=INDXL2G(i, NB, myprowd, 0, nprd )
    END DO
    DO i=1,nstloc_diag_y
      globalindex_diag_y(i)=INDXL2G(i, MB, mypcold, 0, npcd )
    END DO
    IF(my_iso==2) THEN
      globalindex_diag_x=globalindex_diag_x+npsi(1)
      globalindex_diag_y=globalindex_diag_y+npsi(1)
    END IF
  END SUBROUTINE init_mpi_2d
!---------------------------------------------------------------------------  
! DESCRIPTION: init_blacs
!> @brief
!!This routine initializes all contexts and handlers for ScaLAPACK with BLACS.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_blacs
    INTEGER             :: i,j,k,is,dims(2),CTXTTMP,NPROW1d,NPCOL1d,MYPROW1d,MYPCOL1d,&
                           diad_dims(2)
    INTEGER,ALLOCATABLE :: IMAP(:,:)
      CALL BLACS_PINFO(IAM,NPROCS)
      IF (NPROCS.LT.1) THEN
        CALL BLACS_SETUP(IAM,NPROCS)
      END IF
      DO is=1,2
        dims=0
        CALL mpi_dims_create(mpi_nprocs_iso(is),2,dims,mpi_ierror)
        NPROW=dims(1)
        NPCOL=dims(2)
        CALL BLACS_GET(0,0,CTXTTMP)
        ALLOCATE(IMAP(NPROW,NPCOL))
        K=0
        IF(is==2) K=mpi_nprocs_iso(1)
!        WRITE(*,*),'is = ',is,'K = ',K
        DO I = 1, NPROW
          DO J = 1, NPCOL
            IMAP(I, J) = K
            K = K + 1
          END DO
        END DO
!        WRITE(*,*),'imap = ',IMAP      
        CALL BLACS_GRIDMAP( CTXTTMP, IMAP, NPROW, NPROW, NPCOL )
        DEALLOCATE(IMAP)
        IF(is==my_iso) CONTXT=CTXTTMP

      END DO

      CALL BLACS_GRIDINFO(CONTXT,NPROW,NPCOL,MYPROW,MYPCOL)
      CALL BLACS_GET(CONTXT,10,CONTXT1D)
      CALL BLACS_GRIDINIT(CONTXT1D,'Row',1,mpi_nprocs_iso(my_iso))
      CALL BLACS_GRIDINFO(CONTXT1D,NPROW1d,NPCOL1d,MYPROW1d,MYPCOL1d)
      WRITE(*,*),'proc = ',mpi_myproc,nprow1d,npcol1d,myprow1d,mypcol1d
      
      DO is=1,4
        diad_dims=0
        CALL mpi_dims_create(mpi_nprocs_diag(is),2,diad_dims,mpi_ierror)
        nprd=diad_dims(1)
        npcd=diad_dims(2)
        CALL BLACS_GET(0,0,CTXTTMP)
        ALLOCATE(IMAP(nprd,npcd))
        K=0

        IF(is==2) K = mpi_nprocs_iso(1)/2
        IF(is==3) K = mpi_nprocs_iso(1)
        IF(is==4) K = mpi_nprocs_iso(1)+mpi_nprocs_iso(2)/2
!        WRITE(*,*),'is = ',is,'K = ',K
        DO I = 1, nprd
          DO J = 1, npcd
            IMAP(I, J) = K
            K = K + 1
          END DO
        END DO
!        WRITE(*,*),'imap = ',IMAP      
        CALL BLACS_GRIDMAP( CTXTTMP, IMAP, nprd, nprd, npcd )
        DEALLOCATE(IMAP)
        IF(is==my_diag) &
          CONTXT_DO=CTXTTMP
        IF((my_diag==1.OR.my_diag==2).AND.is==1)CONTXT_D=CTXTTMP
        IF((my_diag==1.OR.my_diag==2).AND.is==2)CONTXT_O=CTXTTMP
        IF((my_diag==3.OR.my_diag==4).AND.is==3)CONTXT_D=CTXTTMP
        IF((my_diag==3.OR.my_diag==4).AND.is==4)CONTXT_O=CTXTTMP

      END DO
      IF(my_diag==1.OR.my_diag==3)THEN
        CALL BLACS_GRIDINFO(CONTXT_D,nprd,npcd,myprowd,mypcold)
      ELSE IF(my_diag==2.OR.my_diag==4)THEN
        CALL BLACS_GRIDINFO(CONTXT_O,nprd,npcd,myprowd,mypcold)
      ENDIF
!      WRITE(*,*),'proc = ',mpi_myproc,nprd,npcd,myprowd,mypcold


      CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE init_blacs
!---------------------------------------------------------------------------  
! DESCRIPTION: associate_nodes
!> @brief
!!This subroutine determines the number of processes for neutrons and protons
!!and associates MPI ranks with wave functions. It sets the number of processes 
!!in each group and determines \c globalindex and \c localindex.
!--------------------------------------------------------------------------- 
  SUBROUTINE associate_nodes
    !***********************************************************************
    ! calculates 1d wave function distribution over nodes
    !***********************************************************************
    INTEGER :: ncount,nst,ip,iloc,iq,mpi_ierror
    REAL(db) :: N,Z,sq_N,sq_Z
    ncount=0
    globalindex=0
    node=0
    Z=npsi(2)-npsi(1)
    N = npsi(1)
    sq_N = N*N
    sq_Z = Z*Z
!    mpi_nprocs_iso(1)=NINT(REAL(npsi(1))/REAL(npsi(2))/2.0*mpi_nprocs)*2
    IF(mpi_nprocs <128) THEN
     mpi_nprocs_iso(1)=NINT(sq_N/(sq_N+sq_Z)/2.0*mpi_nprocs)*2
     mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
    ELSEIF(mpi_nprocs >= 128 .AND. mpi_nprocs < 256)THEN
     mpi_nprocs_iso(1)=NINT(sq_N/(sq_N+sq_Z)/2.0*mpi_nprocs)*2
     mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
    ELSEIF(mpi_nprocs >= 256 .AND. mpi_nprocs < 512) THEN
     mpi_nprocs_iso(1)=NINT(sq_N/(sq_N+sq_Z)/16.0*mpi_nprocs)*16
     mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
    ELSEIF(mpi_nprocs >=512) THEN
     mpi_nprocs_iso(1)=NINT(sq_N/(sq_N+sq_Z)/32.0*mpi_nprocs)*32
     mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
    ENDIF
    IF(mpi_nprocs_iso(1)<4) THEN
     mpi_nprocs_iso(1)=4
     mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
    END IF
    IF(mpi_nprocs_iso(2)<4) THEN
     mpi_nprocs_iso(2)=4
     mpi_nprocs_iso(1)=mpi_nprocs-mpi_nprocs_iso(2)
    END IF
    my_iso=1
    IF(mpi_myproc>=mpi_nprocs_iso(1)) my_iso=2
    CALL mpi_comm_split(mpi_comm_world,my_iso,mpi_myproc,comm_iso,mpi_ierror)
    CALL mpi_comm_rank(comm_iso,mpi_myproc_iso,mpi_ierror)
    
    IF(mpi_myproc < mpi_nprocs_iso(1)/2) my_diag=1
    IF(mpi_myproc >= mpi_nprocs_iso(1)/2.AND.mpi_myproc < mpi_nprocs_iso(1))my_diag = 2
      
    IF(mpi_myproc >= mpi_nprocs_iso(1).AND.mpi_myproc < mpi_nprocs_iso(1)+mpi_nprocs_iso(2)/2)my_diag = 3

    IF(mpi_myproc >= (mpi_nprocs_iso(1)+mpi_nprocs_iso(2)/2)) my_diag = 4
 
    mpi_nprocs_diag(1) = mpi_nprocs_iso(1)/2
    mpi_nprocs_diag(2) = mpi_nprocs_iso(1)/2
    mpi_nprocs_diag(3) = mpi_nprocs_iso(2)/2
    mpi_nprocs_diag(4) = mpi_nprocs_iso(2)/2
     
      
!      WRITE(*,*),'diag procs = ',mpi_nprocs_diag
!    WRITE(*,*),'proc = ',mpi_myproc,'iso = ',my_iso,'my diag = ',my_diag

    DO iq=1,2
      DO nst=npmin(iq),npsi(iq)
        node(nst)=MOD((nst-npmin(iq))/nb_psi,mpi_nprocs_iso(iq))
        IF(iq==2) node(nst)=node(nst)+mpi_nprocs_iso(1)
      END DO
    END DO
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    nstloc=0
    DO nst=1,nstmax
       IF(node(nst)==mpi_myproc) THEN
          nstloc=nstloc+1
          globalindex(nstloc)=nst
       ENDIF
    ENDDO
    DO ip=0,mpi_nprocs-1
      iloc=0
      DO nst=1,nstmax
        IF(node(nst)==ip) THEN
          iloc=iloc+1
          localindex(nst)=iloc
        ENDIF
      ENDDO
    ENDDO
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE associate_nodes
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_densities
!> @brief
!!This subroutine uses the \c MPI routine \c mpi_allreduce to sum
!!up the partial densities from the different nodes.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_densities
    USE Densities, ONLY : rho,tau,current,sodens,sdens
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    CALL collect_density(rho)
    CALL collect_density(tau)
    CALL collect_vecdens(current)
    CALL collect_vecdens(sodens)
    CALL collect_vecdens(sdens)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    RETURN
  END SUBROUTINE collect_densities
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_density
!> @brief
!!This subroutine uses the \c MPI routine \c mpi_allreduce to sum
!!up one partial density from the different nodes.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_density(dens)
    REAL(db), INTENT(INOUT)  :: dens(:,:,:,:)
    INTEGER                  :: cnt
    cnt=SIZE(dens)/2
    IF(mpi_myproc_iso==0) THEN
      CALL mpi_reduce(MPI_IN_PLACE,dens(:,:,:,my_iso),cnt,&
                      mpi_double_precision,mpi_sum,0,comm_iso,mpi_ierror)
    ELSE
      CALL mpi_reduce(dens(:,:,:,my_iso),dens(:,:,:,my_iso),cnt,&
                      mpi_double_precision,mpi_sum,0,comm_iso,mpi_ierror)
    END IF
    CALL mpi_bcast(dens(:,:,:,1),cnt,mpi_double_precision,&
                   0,mpi_comm_world,mpi_ierror)
    CALL mpi_bcast(dens(:,:,:,2),cnt,mpi_double_precision,&
                   mpi_nprocs_iso(1),mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_density
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_vecdens
!> @brief
!!This subroutine uses the \c MPI routine \c mpi_allreduce to sum
!!up one partial vector density from the different nodes.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_vecdens(dens)
    REAL(db), INTENT(INOUT)  :: dens(:,:,:,:,:)
    INTEGER                  :: cnt
    cnt=SIZE(dens)/2
    IF(mpi_myproc_iso==0) THEN
      CALL mpi_reduce(MPI_IN_PLACE,dens(:,:,:,:,my_iso),cnt,&
                      mpi_double_precision,mpi_sum,0,comm_iso,mpi_ierror)
    ELSE
      CALL mpi_reduce(dens(:,:,:,:,my_iso),dens(:,:,:,:,my_iso),cnt,&
                      mpi_double_precision,mpi_sum,0,comm_iso,mpi_ierror)
    END IF
    CALL mpi_bcast(dens(:,:,:,:,1),cnt,mpi_double_precision,&
                   0,mpi_comm_world,mpi_ierror)
    CALL mpi_bcast(dens(:,:,:,:,2),cnt,mpi_double_precision,&
                   mpi_nprocs_iso(1),mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_vecdens
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_sp_property
!> @brief
!!This subroutine collects one single-particle property calculated
!!from the wave functions and thus available only for the local wave
!!functions on each node. 
!>
!> @details
!!It uses a simple trick: the arrays like \c sp_energy 
!!are defined for the full set of indices but set to zero
!!before the calculation of these properties. On each node then the
!!local values are calculated but inserted at the proper index for the
!!full set of wave functions. In this subroutine the results from all
!!the nodes are added up using \c mpi_reduce, so that effectively
!!for each index one node contributes the correct value and the others
!!zeroes. 
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_sp_property(sp_prop)
    REAL(db),INTENT(INOUT) :: sp_prop(:)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_prop,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_sp_property
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_sp_properties
!> @brief
!!This subroutine collects the single-particle properties calculated
!!from the wave functions and thus available only for the local wave
!!functions on each node. 
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_sp_properties
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_kinetic,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_orbital,3*nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_spin,3*nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_energy,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_efluct1,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_efluct2,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_norm,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sp_parity,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_sp_properties
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_energies
!> @brief
!!This subroutine sums up three energies which are partly calculated on each node.
!> @param[in,out] delesum
!> REAL(db), takes and returns energy differences.
!> @param[in,out] sumflu
!> REAL(db), takes and returns fluctuations.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_energies(delesum,sumflu)
    !***********************************************************************
    ! collects s.p. energies, fluctuation measure and energy differences
    !***********************************************************************
    REAL(db),INTENT(INOUT) :: delesum,sumflu
    INTEGER  :: nst
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    DO nst=1,nstmax
      IF(node(nst)/=mpi_myproc) THEN
        sp_energy(nst)=0.0d0
        sp_efluct1(nst)=0.0d0
        sp_efluct2(nst)=0.0d0
      END IF
    ENDDO
    CALL mpi_allreduce(MPI_IN_PLACE,sp_energy,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sumflu,1,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,delesum,1,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_energies
!---------------------------------------------------------------------------  
! DESCRIPTION: finish_mpi
!> @brief
!!This is just a wrapper for the \c MPI finalization call.
!--------------------------------------------------------------------------- 
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
!---------------------------------------------------------------------------  
! DESCRIPTION: mpi_start_timer
!> @brief
!!This subroutine starts an MPI timer with index \c index for the global group.
!> @param[in] index
!> INTEGER, takes the index from the array where the start time will be stored.
!--------------------------------------------------------------------------- 
  SUBROUTINE mpi_start_timer(index)
    INTEGER,INTENT(IN) :: index
    INTEGER :: ierr
    CALL mpi_barrier (mpi_comm_world,ierr)
    timer(index)=mpi_wtime()
  END SUBROUTINE mpi_start_timer
!---------------------------------------------------------------------------  
! DESCRIPTION: mpi_stop_timer
!> @brief
!!This subroutine stops an MPI timer with index \c index for the global group.
!> @param[in] index
!> INTEGER, takes the index from the array where the start time was stored.
!> @param[in] textline
!> INTEGER, takes a text which is outputted in front of the time.
!--------------------------------------------------------------------------- 
  SUBROUTINE mpi_stop_timer(index,textline)
    INTEGER,INTENT(IN) :: index
    CHARACTER(*),INTENT(IN) :: textline
    INTEGER :: ierr
    CALL mpi_barrier (mpi_comm_world,ierr)
    IF(wflag)WRITE(*,'(A20,F10.4)')textline,mpi_wtime()-timer(index)
  END SUBROUTINE mpi_stop_timer
!---------------------------------------------------------------------------  
! DESCRIPTION: mpi_start_timer_iq
!> @brief
!!This subroutine starts an MPI timer with index \c index for the isospin group.
!> @param[in] index
!> INTEGER, takes the index from the array where the start time will be stored.
!--------------------------------------------------------------------------- 
  SUBROUTINE mpi_start_timer_iq(index)
    INTEGER,INTENT(IN) :: index
    INTEGER :: ierr
!    CALL mpi_barrier (comm_iso,ierr)
    CALL mpi_barrier(mpi_comm_world,ierr)
    timer(index)=mpi_wtime()
  END SUBROUTINE mpi_start_timer_iq
!---------------------------------------------------------------------------  
! DESCRIPTION: mpi_stop_timer_iq
!> @brief
!!This subroutine stops an MPI timer with index \c index for the isospin group.
!> @param[in] index
!> INTEGER, takes the index from the array where the start time was stored.
!> @param[in] textline
!> INTEGER, takes a text which is outputted in front of the time.
!--------------------------------------------------------------------------- 
  SUBROUTINE mpi_stop_timer_iq(index,textline)
    INTEGER,INTENT(IN) :: index
    CHARACTER(*),INTENT(IN) :: textline
    INTEGER :: ierr
    CALL mpi_barrier (comm_iso,ierr)
    IF(mpi_myproc_iso==0)WRITE(*,'(A20,F10.4,A4,I4)')textline,mpi_wtime()-timer(index),' iq = ',my_iso
  END SUBROUTINE mpi_stop_timer_iq
END MODULE Parallel
