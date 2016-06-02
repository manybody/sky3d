MODULE Parallel
  USE Params, ONLY : wflag,nprocs,myid,converfile,wffile,dipolesfile,&
       momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.TRUE. 
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:)
  INTEGER :: mpi_nprocs,mpi_ierror,mpi_myproc
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax))
  END SUBROUTINE alloc_nodes
  !***********************************************************************
  SUBROUTINE init_all_mpi
    CALL mpi_init(mpi_ierror)
    CALL mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)
    CALL mpi_comm_rank(mpi_comm_world,mpi_myproc,mpi_ierror)
    nprocs=mpi_nprocs
    myid=mpi_myproc
    wflag=.TRUE.
    IF(wflag) WRITE(*,'(a,i5)') ' number of nodes=',mpi_nprocs
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE init_all_mpi
  !***********************************************************************
  SUBROUTINE mpi_init_filename
    WRITE(wffile,'(A2,I0.3,A5)')'wf',mpi_myproc,'.data'
    WRITE(converfile,'(A5,I0.3,A4)')'conver',mpi_myproc,'.res'
    WRITE(dipolesfile,'(A6,I0.3,A4)')'dipole',mpi_myproc,'.res'
    WRITE(momentafile,'(A3,I0.3,A4)')'mom',mpi_myproc,'.res'
    WRITE(energiesfile,'(A5,I0.3,A4)')'energ',mpi_myproc,'.res'
    WRITE(quadrupolesfile,'(A4,I0.3,A4)')'quad',mpi_myproc,'.res'
    WRITE(spinfile,'(A4,I0.3,A4)')'spin',mpi_myproc,'.res'
    WRITE(extfieldfile,'(A3,I0.3,A4)')'ext',mpi_myproc,'.res'  
    WRITE(*,*)'Number of processes: ', mpi_nprocs,'Process number: ', mpi_myproc
  END SUBROUTINE mpi_init_filename
  !***********************************************************************
  SUBROUTINE associate_nodes
    INTEGER :: i
    node=mpi_myproc
    nstloc=nstmax
    FORALL(i=1:nstmax) 
       globalindex(i)=i
       localindex(i)=i
    END FORALL
  END SUBROUTINE associate_nodes
  !***********************************************************************
  SUBROUTINE collect_densities
    USE Densities, ONLY : rho,tau,current,sodens,sdens
    REAL(db) :: tmp_rho(nx,ny,nz,2),tmp_current(nx,ny,nz,3,2)
    REAL(db) :: rsum
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    CALL mpi_allreduce(rho,tmp_rho,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    rho=tmp_rho/mpi_nprocs
    CALL mpi_allreduce(tau,tmp_rho,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tau=tmp_rho/mpi_nprocs
    CALL mpi_allreduce(current,tmp_current,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    current=tmp_current/mpi_nprocs
    CALL mpi_allreduce(sodens,tmp_current,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    sodens=tmp_current/mpi_nprocs
    CALL mpi_allreduce(sdens,tmp_current,3*2*nx*ny*nz,   &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    sdens=tmp_current/mpi_nprocs
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    RETURN
  END SUBROUTINE collect_densities
  !***********************************************************************
  SUBROUTINE collect_sp_properties
    REAL(db) :: tmpgat(nstmax),tmpgat3(3,nstmax)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_sp_properties
  !***********************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
END MODULE Parallel
