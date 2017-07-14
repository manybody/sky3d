MODULE Parallel
  USE Params, ONLY: wflag,db
  USE Levels, ONLY: nstmax,npsi,nstloc,npmin
  USE Grids, ONLY: nx,ny,nz
  IMPLICIT NONE
  SAVE
  LOGICAL,PARAMETER    :: tmpi=.FALSE.,ttabc=.FALSE.
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:),&
                          globalindex_x(:),globalindex_y(:),&
                          globalindex_diag_x(:),globalindex_diag_y(:)
  INTEGER              :: mpi_nprocs,mpi_ierror,mpi_myproc,my_diag,psiloc_x,psiloc_y, &
                          processor_name,proc_namelen,nstloc_x,nstloc_y,&
                          nstloc_diag_x,nstloc_diag_y,my_iso
  INTEGER              :: mpi_comm_world,mpi_sum,mpi_double_precision
CONTAINS     !  all dummy subroutines to run on a sequential machine
  !************************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             globalindex_x(nstmax),globalindex_y(nstmax),&
             globalindex_diag_x(nstmax),globalindex_diag_y(nstmax))
  END SUBROUTINE alloc_nodes
  !************************************************************************
  SUBROUTINE init_all_mpi
    WRITE(*,*) '***** Running sequential version *****'
    mpi_nprocs=1
    mpi_myproc=0
    wflag=.TRUE.
    my_diag=0
    my_iso=0
  END SUBROUTINE init_all_mpi
  !************************************************************************  
  SUBROUTINE mpi_init_filename
    CONTINUE
  END SUBROUTINE mpi_init_filename   
  !************************************************************************  
  FUNCTION tabc_av(val)
    REAL(db),INTENT(IN) :: val
    REAL(db)            :: tabc_av
    tabc_av=val
  END FUNCTION tabc_av
  !************************************************************************
  FUNCTION tabc_filename(filename)
    CHARACTER(64),INTENT(IN) :: filename
    CHARACTER(64)            :: tabc_filename
    tabc_filename=filename
  END FUNCTION tabc_filename
  !************************************************************************
  FUNCTION tabc_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,2)
    REAL(db)            :: tabc_dens(nx,ny,nz,2)
    tabc_dens=density
  END FUNCTION tabc_dens
  !************************************************************************
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
  !************************************************************************
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
