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
  USE Params, ONLY: wflag,db
  USE Levels, ONLY: nstmax,npsi,nstloc,npmin
  USE Grids, ONLY: nx,ny,nz
  IMPLICIT NONE
  SAVE
  LOGICAL,PARAMETER :: tmpi=.FALSE.       !<a logical variable set to true if \c MPI
  !!parallelization is activated. It is used to turn the calling of all
  !!the \c MPI routines in the code on or off.
  LOGICAL, PARAMETER   :: ttabc=.FALSE.
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
  INTEGER              :: mpi_comm_world,mpi_sum,mpi_double_precision
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
!---------------------------------------------------------------------------  
! DESCRIPTION: init_all_mpi
!> @brief
!!This subroutine initializes \c MPI and
!!finds out the number of processors \c mpi_nprocs as well as the
!!index of the current one \c mpi_myproc. The flag \c wflag is
!!set to true only for the processor numbered 0.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_all_mpi
    WRITE(*,*) '***** Running sequential version *****'
    mpi_nprocs=1
    mpi_myproc=0
    wflag=.TRUE.
    my_diag=0
    my_iso=0
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
  !************************************************************************
  SUBROUTINE mpi_init(ierror)
    INTEGER :: ierror
    STOP ' MPI_INIT: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_init
  !************************************************************************
  SUBROUTINE mpi_comm_size(comm_world,nprocs,ierror)
    INTEGER :: ierror, nprocs, comm_world
    STOP ' MPI_COMM_SIZE: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_comm_size
  !************************************************************************
  SUBROUTINE mpi_comm_rank(comm_world,myproc,ierror)
    INTEGER :: ierror, myproc, comm_world
    STOP ' MPI_COMM_RANK: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_comm_rank
  !************************************************************************
  SUBROUTINE mpi_get_processor_name(processor_name,proc_namelen,ierror)
    INTEGER :: ierror, processor_name, proc_namelen
    STOP ' MPI_GET_PROCESSOR_NAME: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_get_processor_name
  !************************************************************************
  SUBROUTINE mpi_barrier (comm_world, ierror)
    INTEGER :: ierror, comm_world
    STOP ' MPI_BARRIER: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_barrier
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
  SUBROUTINE mpi_allreduce(rho,tmp_rho,length,        &
       i_double_precision,sum,  &
       comm_world,ierror)
    INTEGER :: ierror, comm_world, i_double_precision, length, sum
    REAL(db), DIMENSION(*), INTENT(IN) :: rho,tmp_rho
    STOP ' MPI_ALLREDUCE: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_allreduce
  !************************************************************************
  SUBROUTINE collect_densities
    STOP ' COLLECT_DENSITIES: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_densities
  !************************************************************************
  SUBROUTINE collect_density(dens)
    REAL(db):: dens(:,:,:,:)
    STOP ' COLLECT_DENSITY: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_density
  !************************************************************************
  SUBROUTINE collect_sp_property(dens)
    REAL(db):: dens(:)
    STOP ' COLLECT_SP_PROPERTY: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_property
  !************************************************************************
  SUBROUTINE collect_sp_properties
    STOP ' COLLECT_SP_PROPERTIES: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_properties
  !************************************************************************
  SUBROUTINE collect_energies(delesum,sumflu)
    REAL(db), INTENT(IN) :: delesum,sumflu
    STOP ' COLLECT_ENERGIES: parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_energies
  !************************************************************************
  SUBROUTINE finish_mpi
  END SUBROUTINE finish_mpi
END MODULE Parallel
