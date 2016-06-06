MODULE Parallel
  USE Params, ONLY: wflag,db,tabc_nprocs,tabc_myid,converfile,wffile,dipolesfile,&
  momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile
  USE Levels, ONLY: nstmax,npsi,nstloc
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.FALSE.,ttabc=.TRUE.
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:)
  INTEGER :: mpi_nprocs,mpi_ierror,mpi_myproc, &
       processor_name,proc_namelen
CONTAINS     !  all dummy subroutines to run tabc on a parallel machine
  !************************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax))
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
    WRITE(wffile,'(A2,I0.3,A5)')'wf',tabc_myid,'.data'
    WRITE(converfile,'(A5,I0.3,A4)')'conver',tabc_myid,'.res'
    WRITE(dipolesfile,'(A6,I0.3,A4)')'dipole',tabc_myid,'.res'
    WRITE(momentafile,'(A3,I0.3,A4)')'mom',tabc_myid,'.res'
    WRITE(energiesfile,'(A5,I0.3,A4)')'energ',tabc_myid,'.res'
    WRITE(quadrupolesfile,'(A4,I0.3,A4)')'quad',tabc_myid,'.res'
    WRITE(spinfile,'(A4,I0.3,A4)')'spin',tabc_myid,'.res'
    WRITE(extfieldfile,'(A3,I0.3,A4)')'ext',tabc_myid,'.res'  
    WRITE(*,*)'Number of processes: ', tabc_nprocs,'Process number: ',tabc_myid
  END SUBROUTINE mpi_init_filename
  !************************************************************************
  SUBROUTINE mpi_filename(filename)
    CHARACTER(64),INTENT(INOUT) :: filename
    WRITE(filename,'(A2,I0.3)')'wf',mpi_myproc
  END SUBROUTINE mpi_filename
  !************************************************************************
  SUBROUTINE mpi_get_processor_name(processor_name,proc_namelen,ierror)
    INTEGER :: ierror, processor_name, proc_namelen
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_get_processor_name
  !************************************************************************
  SUBROUTINE mpi_barrier (comm_world, ierror)
    INTEGER :: ierror, comm_world
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_barrier
  !************************************************************************
  SUBROUTINE associate_nodes
    INTEGER :: i
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
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_allreduce
  !************************************************************************
  SUBROUTINE collect_densities
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_densities
  !************************************************************************
  SUBROUTINE collect_sp_properties
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_properties
  !************************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
END MODULE Parallel
