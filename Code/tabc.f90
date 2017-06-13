MODULE Parallel
  USE Params, ONLY: wflag,db,tabc_nprocs,tabc_myid,converfile,wffile,dipolesfile,&
  momentafile,energiesfile,quadrupolesfile,spinfile,extfieldfile,diffenergiesfile
  USE Grids, ONLY: nx,ny,nz
  USE Levels, ONLY: nstmax,npsi,nstloc,npmin
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.FALSE.,ttabc=.TRUE.
  INTEGER, ALLOCATABLE ::node(:),localindex(:),globalindex(:),node_x(:),node_y(:),&
                         localindex_x(:),localindex_y(:),globalindex_x(:,:),globalindex_y(:,:)
  INTEGER :: mpi_nprocs,mpi_ierror,mpi_myproc,processor_name,proc_namelen,nstloc_x(2),nstloc_y(2)
CONTAINS     !  all dummy subroutines to run tabc on a parallel machine
  !************************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             localindex_x(nstmax),localindex_y(nstmax),&
             globalindex_x(nstmax,2),globalindex_y(nstmax,2),&
             node_x(nstmax),node_y(nstmax))
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
    CALL mpi_allreduce(density,tabc_dens,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tabc_dens=tabc_dens/REAL(tabc_nprocs)
  END FUNCTION tabc_dens
  !************************************************************************
  FUNCTION tabc_av(val)
    REAL(db),INTENT(IN) :: val
    REAL(db)            :: tabc_av
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    WRITE(*,*)'++++++++++++++',val,tabc_av,tabc_myid
    CALL mpi_allreduce(val,tabc_av,1,mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    tabc_av=tabc_av/REAL(tabc_nprocs)
  END FUNCTION tabc_av
  !************************************************************************
  FUNCTION tabc_vec_dens(density)
    REAL(db),INTENT(IN) :: density(nx,ny,nz,3,2)
    REAL(db)            :: tabc_vec_dens(nx,ny,nz,3,2)
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
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
  SUBROUTINE associate_nodes
    INTEGER :: i,is,noffset
    node=0
    node_x=0
    node_y=0
    nstloc=nstmax
    FORALL(i=1:nstmax)
       globalindex(i)=i
       localindex(i)=i
    END FORALL
    DO is=1,2
      noffset=npmin(is)-1
      nstloc_x(is)=npsi(is)-npmin(is)+1
      nstloc_y(is)=npsi(is)-npmin(is)+1
      DO i=npmin(is),npsi(is)
        localindex_x(i) = i-noffset
        localindex_y(i) = i-noffset
        globalindex_x(localindex_x(i),is)=i
        globalindex_y(localindex_y(i),is)=i
      END DO
    END DO
  END SUBROUTINE associate_nodes
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
