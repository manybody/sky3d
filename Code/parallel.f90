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
!!
!!Since an efficient way of dealing with Gram-Schmidt orthogonalization
!!in this case was yet not found, at present the code can be run in \c MPI 
!!parallel mode only for the dynamic case. 
!------------------------------------------------------------------------------
MODULE Parallel
  USE Params, ONLY : wflag
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.TRUE.       !<a logical variable set to true if \c MPI
  !!parallelization is activated. It is used to turn the calling of all
  !!the \c MPI routines in the code on or off.
  INTEGER, ALLOCATABLE :: node(:)        !<For the single-particle state
  !!with index \c i the wave function is stored on computing node \c node(i).
  INTEGER, ALLOCATABLE :: localindex(:)  !<For the single-particle state
  !!with index \c i the wave function is stored on computing node
  !!\c node(i) and its index on that node is \c localindex(i).
  INTEGER, ALLOCATABLE :: globalindex(:) !<tells the index of the single-particle
  !!state in the whole array of \c nstmax states (it could be
  !!dimensioned \c nstloc but is dimensioned as \c nstmax 
  !!to make its allocation simpler). So for wave function index \c i 
  !!<em> on the local node</em>, <tt> i=1..nstloc</tt>, the
  !!single-particle energy must be obtained using <tt> sp_energy(globalindex(i))</tt>.  
  INTEGER :: mpi_nprocs                  !<number of MPI processes.
  INTEGER :: mpi_ierror                  !<varable for error output of MPI routines.
  INTEGER :: mpi_myproc                  !<the number of the local MPI process
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_nodes
!> @brief This subroutine merely allocates the internal arrays of module \c Parallel.

!--------------------------------------------------------------------------- 
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax))
  END SUBROUTINE alloc_nodes
!---------------------------------------------------------------------------  
! DESCRIPTION: init_all_mpi
!> @brief
!!\c init_all_mpi This subroutine initializes \c MPI and
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
! DESCRIPTION: associate_nodes
!> @brief
!!The first loop in this subroutine distributes the wave functions over
!!the nodes. This is done by looping over the wave functions and
!!assigning one to each processor in turn. When the number of processors
!!has been reached, it restarts from processor 0. This way of allocation
!!is to some extent arbitrary and can be changed.
!!
!!The second loop then calculates which wave functions are present on
!!the local node and records their index \c gobalindex in the
!!complete sequence. The third loop sets up the reverse pointers
!!\c localindex, which has to be done in a loop over all processors to
!!set up the information for the proper global indices.
!--------------------------------------------------------------------------- 
  SUBROUTINE associate_nodes
    INTEGER :: ncount,nst,ip,iloc
    ncount=0
    DO nst=1,nstmax  
       ncount=MOD(ncount,mpi_nprocs)
       node(nst)=ncount
       ncount=ncount+1
    ENDDO
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
    IF(wflag) THEN
       WRITE(*,'(A/(1X,20I4))')   &
            ' sorting of wfs on nodes:',node(1:nstmax)
    ENDIF
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE associate_nodes
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_densities
!> @brief
!!This subroutine uses the \c MPI routine \c mpi_allreduce to sum
!!up the partial densities from the different nodes, using temporary
!!arrays \c tmp_rho and \c tmp_current (depending on whether it
!!is a scalar or vector field) in the process.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_densities
    USE Densities, ONLY : rho,tau,current,sodens,sdens
    REAL(db) :: tmp_rho(nx,ny,nz,2),tmp_current(nx,ny,nz,3,2)
    REAL(db) :: rsum
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    CALL mpi_allreduce(rho,tmp_rho,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    rho=tmp_rho
    CALL mpi_allreduce(tau,tmp_rho,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tau=tmp_rho
    CALL mpi_allreduce(current,tmp_current,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    current=tmp_current
    CALL mpi_allreduce(sodens,tmp_current,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    sodens=tmp_current
    !     collect spin densities
    CALL mpi_allreduce(sdens,tmp_current,3*2*nx*ny*nz,   &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    sdens=tmp_current
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    RETURN
  END SUBROUTINE collect_densities
!---------------------------------------------------------------------------  
! DESCRIPTION: collect_sp_properties
!> @brief
!!This subroutine collects the single-particle properties calculated
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
!!zeroes. This process sounds inefficient but considering the small size
!!of the arrays that does not matter.
!--------------------------------------------------------------------------- 
  SUBROUTINE collect_sp_properties
    REAL(db) :: tmpgat(nstmax),tmpgat3(3,nstmax)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(sp_kinetic,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_kinetic=tmpgat
    CALL mpi_allreduce(sp_orbital,tmpgat3,3*nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_orbital=tmpgat3
    CALL mpi_allreduce(sp_spin,tmpgat3,3*nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_spin=tmpgat3
    CALL mpi_allreduce(sp_energy,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_energy=tmpgat
    CALL mpi_allreduce(sp_efluct1,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_efluct1=tmpgat
    CALL mpi_allreduce(sp_efluct2,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_efluct2=tmpgat
    CALL mpi_allreduce(sp_norm,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_norm=tmpgat
    CALL mpi_allreduce(sp_parity,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_parity=tmpgat
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_sp_properties
!---------------------------------------------------------------------------  
! DESCRIPTION: finish_mpi
!> @brief
!!This is just a wrapper for the \c MPI finalization call.
!--------------------------------------------------------------------------- 
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
END MODULE Parallel
