MODULE Parallel
  USE Params, ONLY : wflag
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  LOGICAL,PARAMETER :: tmpi=.TRUE.,ttabc=.FALSE.
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
    wflag=mpi_myproc==0
    IF(wflag) WRITE(*,'(a,i5)') ' number of nodes=',mpi_nprocs
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  END SUBROUTINE init_all_mpi
  !************************************************************************     
  SUBROUTINE mpi_init_filename
    CONTINUE
  END SUBROUTINE mpi_init_filename
  !***********************************************************************
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
  !***********************************************************************
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
  !***********************************************************************
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
  !***********************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
END MODULE Parallel