MODULE Parallel
  USE Params, ONLY: wflag,db,tabc_nprocs,tabc_myid,converfile,wffile,dipolesfile,&
  momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile,diffenergiesfile
  USE Levels, ONLY: nstmax,npsi,nstloc,npmin
  USE Grids, ONLY: nx,ny,nz
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.FALSE.       !<a logical variable set to true if \c MPI
  !!parallelization is activated. It is used to turn the calling of all
  !!the \c MPI routines in the code on or off.
  LOGICAL, PARAMETER   :: ttabc=.TRUE.
  !>@name Variables characterizing the different distributions of wave functions. 
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
  INTEGER              :: mpi_nprocs            !<number of MPI processes.
  INTEGER              :: mpi_ierror            !<varable for error output of MPI routines.
  INTEGER              :: mpi_myproc            !<the number of the local MPI process
  INTEGER              :: my_diag               !<determines the local group for diagonalization/orthonormalization.
  INTEGER              :: processor_name,proc_namelen
  INTEGER              :: my_iso                !<the local isospin that is calculated.
CONTAINS     !  all dummy subroutines to run on a sequential machine
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_nodes
!> @brief This subroutine merely allocates the internal arrays of module \c Parallel.
!---------------------------------------------------------------------------  
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             globalindex_x(nstmax),globalindex_y(nstmax),&
             globalindex_diag_x(nstmax),globalindex_diag_y(nstmax))
  END SUBROUTINE alloc_nodes
  !************************************************************************
  SUBROUTINE init_all_mpi
    WRITE(*,*) '***** Running TABC version *****'
    CALL mpi_init(mpi_ierror)
    CALL mpi_comm_size(mpi_comm_world,tabc_nprocs,mpi_ierror)
    CALL mpi_comm_rank(mpi_comm_world,tabc_myid,mpi_ierror)
    mpi_nprocs=1
    mpi_myproc=0
    wflag=.TRUE.
  END SUBROUTINE init_all_mpi
  !************************************************************************
  SUBROUTINE mpi_init_filename
    WRITE(wffile,'(A3,I0.3)')'wf-',tabc_myid
    WRITE(converfile,'(A5,I0.3,A4)')'conver',tabc_myid,'.res'
    WRITE(dipolesfile,'(A6,I0.3,A4)')'dipole',tabc_myid,'.res'
    WRITE(momentafile,'(A3,I0.3,A4)')'mom',tabc_myid,'.res'
    WRITE(energiesfile,'(A5,I0.3,A4)')'energ',tabc_myid,'.res'
    WRITE(diffenergiesfile,'(A9,I0.3,A4)')'diffenerg',tabc_myid,'.res'
    WRITE(quadrupolesfile,'(A4,I0.3,A4)')'quad',tabc_myid,'.res'
    WRITE(spinfile,'(A4,I0.3,A4)')'spin',tabc_myid,'.res'
    WRITE(extfieldfile,'(A3,I0.3,A4)')'ext',tabc_myid,'.res'  
    WRITE(*,*)'Number of processes: ', tabc_nprocs,'Process number: ',tabc_myid
  END SUBROUTINE mpi_init_filename
  !************************************************************************
  FUNCTION tabc_filename(filename)
    CHARACTER(64),INTENT(IN) :: filename
    CHARACTER(64)            :: tabc_filename,myid
    WRITE(myid,'(A1,I0.3)')'-',tabc_myid
    tabc_filename=TRIM(filename)//TRIM(myid)
  END FUNCTION tabc_filename
  !************************************************************************
  FUNCTION tabc_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,2)
    REAL(db)            :: tabc_dens(nx,ny,nz,2)
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    WRITE(*,*) 'tabc_dens'
    CALL mpi_allreduce(density,tabc_dens,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tabc_dens=tabc_dens/REAL(tabc_nprocs)
  END FUNCTION tabc_dens
  !************************************************************************
  FUNCTION tabc_av(val)
    REAL(db),INTENT(IN) :: val
    REAL(db)            :: tabc_av
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    WRITE(*,*) 'tabc_av'
    CALL mpi_allreduce(val,tabc_av,1,mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tabc_av=tabc_av/REAL(tabc_nprocs)
  END FUNCTION tabc_av
  !************************************************************************
  FUNCTION tabc_vec_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,3,2)
    REAL(db)            :: tabc_vec_dens(nx,ny,nz,3,2)
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    WRITE(*,*) 'tabc_vec_dens'
    CALL mpi_allreduce(density,tabc_vec_dens,2*3*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tabc_vec_dens=tabc_vec_dens/REAL(tabc_nprocs)
  END FUNCTION tabc_vec_dens
  !************************************************************************
  SUBROUTINE mpi_get_processor_name(processor_name,proc_namelen,ierror)
    INTEGER :: ierror, processor_name, proc_namelen
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_get_processor_name
  !************************************************************************
!---------------------------------------------------------------------------  
! DESCRIPTION: associate_nodes
!> @brief
!!This subroutine determines the number of processes for neutrons and protons
!!and associates MPI ranks with wave functions. It sets the number of processes 
!!in each group and determines \c globalindex and \c localindex.
!--------------------------------------------------------------------------- 
  SUBROUTINE associate_nodes
    INTEGER :: i,is,noffset
    node=0
    nstloc=nstmax
    FORALL(i=1:nstmax) 
      globalindex(i)=i
      localindex(i)=i
    END FORALL
  END SUBROUTINE associate_nodes
  !************************************************************************
  SUBROUTINE collect_densities
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_densities
  !************************************************************************
  SUBROUTINE collect_density(dens)
    REAL(db):: dens(:,:,:,:)
    STOP ' COLLECT_DENSITY: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_density
  !************************************************************************
  SUBROUTINE collect_sp_properties
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_properties
  !************************************************************************
  SUBROUTINE collect_sp_property(dens)
    REAL(db):: dens(:)
    STOP ' COLLECT_SP_PROPERTY: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_property
  !************************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
END MODULE Parallel
