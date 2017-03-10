MODULE Parallel
  USE Params, ONLY : wflag
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  INTEGER, PARAMETER   :: NB=32,MB=32,NB_psi = 32
  LOGICAL, PARAMETER   :: tmpi=.TRUE.,ttabc=.FALSE.
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:),&
                          node_x(:),node_y(:),localindex_x(:),localindex_y(:),&
                          globalindex_x(:,:),globalindex_y(:,:)
  INTEGER, ALLOCATABLE :: recvcounts(:,:),displs(:,:)
  INTEGER, ALLOCATABLE :: node_2dto1d(:,:),node_1dto2d_x(:),node_1dto2d_y(:)
  INTEGER              :: mpi_nprocs,mpi_ierror,mpi_myproc
  INTEGER              :: comm2d,mpi_dims(2),mpi_mycoords(2),nstloc_iso(2),nstloc_x(2),&
                          nstloc_y(2),first(2),npmin_loc(2),npsi_loc(2)
  INTEGER              :: comm2d_x,comm2d_y,mpi_size_x,mpi_size_y,mpi_rank_x,mpi_rank_y
  INTEGER              :: NPROCS,NPROW,NPCOL,MYPROW,MYPCOL,CONTXT,IAM
  INTEGER, EXTERNAL    :: NUMROC,INDXL2G,INDXG2L,INDXG2P
  REAL(db)             :: timer(20)
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             localindex_x(nstmax),localindex_y(nstmax))
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
  !************************************************************************
  FUNCTION tabc_av(val)
    REAL(db),INTENT(IN) :: val
    REAL(db)            :: tabc_av
    tabc_av=val
  END FUNCTION tabc_av
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
  FUNCTION tabc_filename(filename)
    CHARACTER(64),INTENT(IN) :: filename
    CHARACTER(64)            :: tabc_filename
    tabc_filename=filename
  END FUNCTION tabc_filename
  !***********************************************************************
  SUBROUTINE init_mpi_2d
    !***********************************************************************
    ! calculates 2d wave function distribution over nodes
    !***********************************************************************
    LOGICAL :: isperiodic(2),reorder
    INTEGER :: i,is,noffset
    isperiodic=.FALSE.
    reorder=.FALSE.
    
    !Calculate best dimensions for the given number of processes
    CALL mpi_dims_create(mpi_nprocs,2,mpi_dims,mpi_ierror)

    !Create 2-dimensional grid of processes to calculate the matrices in diagstep
    CALL mpi_cart_create(mpi_comm_world,2,mpi_dims,isperiodic,reorder,comm2d,mpi_ierror)

    !get my coordinates on the grid
    CALL mpi_cart_get(comm2d,2,mpi_dims,isperiodic,mpi_mycoords,mpi_ierror)
    IF(wflag)WRITE(*,*)'Initialized 2d-grid with dimensions ',mpi_dims(1),' and',mpi_dims(2)

    !Create communicators for x- and y-directions
    CALL mpi_comm_split(comm2d,mpi_mycoords(2),mpi_mycoords(1),comm2d_x,mpi_ierror)
    CALL mpi_comm_split(comm2d,mpi_mycoords(1),mpi_mycoords(2),comm2d_y,mpi_ierror)

    !determine their sizes and ranks (is already known, just to make sure)
    CALL mpi_comm_size(comm2d_x,mpi_size_x,mpi_ierror)
    CALL mpi_comm_rank(comm2d_x,mpi_rank_x,mpi_ierror)
    CALL mpi_comm_size(comm2d_y,mpi_size_y,mpi_ierror)
    CALL mpi_comm_rank(comm2d_y,mpi_rank_y,mpi_ierror)

    CALL init_blacs

    DO is=1,2
      nstloc_x(is) = NUMROC(npsi(is)-npmin(is)+1,NB,MYPROW,0,NPROW)
      nstloc_y(is) = NUMROC(npsi(is)-npmin(is)+1,MB,MYPCOL,0,NPCOL)
    END DO

    ALLOCATE(node_x(nstmax),node_y(nstmax),globalindex_x(nstmax,2),globalindex_y(nstmax,2))

    globalindex_x=0
    globalindex_y=0

    DO is=1,2
      noffset=npmin(is)-1
      DO i=npmin(is),npsi(is)
        localindex_x(i) = INDXG2L(i-noffset, NB, MYPROW, 0, NPROW)
        localindex_y(i) = INDXG2L(i-noffset, MB, MYPCOL, 0, NPCOL)
        node_x(i)       = INDXG2P(i-noffset, NB, MYPROW, 0, NPROW)
        node_y(i)       = INDXG2P(i-noffset, MB, MYPCOL, 0, NPCOL)
        IF(node_x(i)==mpi_rank_x) globalindex_x(localindex_x(i),is)=i
        IF(node_y(i)==mpi_rank_y) globalindex_y(localindex_y(i),is)=i
      END DO
    END DO

  END SUBROUTINE init_mpi_2d
!***************************************************************************
  SUBROUTINE init_blacs
        CALL BLACS_PINFO(IAM,NPROCS)
        IF (NPROCS.LT.1) THEN
          CALL BLACS_SETUP(IAM,NPROCS)
        END IF
        NPROW=mpi_dims(1)
        NPCOL=mpi_dims(2)
        CALL BLACS_GET(0,0,CONTXT)
        CALL BLACS_GRIDINIT(CONTXT,'Row',NPROW,NPCOL)
        CALL BLACS_GRIDINFO(CONTXT,NPROW,NPCOL,MYPROW,MYPCOL)
        IF(mpi_rank_x/=MYPROW.OR.mpi_rank_y/=MYPCOL) STOP 'BLACS and MPI init is different'
  END SUBROUTINE init_blacs
  !***********************************************************************
  SUBROUTINE associate_nodes
    !***********************************************************************
    ! calculates 1d wave function distribution over nodes
    !***********************************************************************
    INTEGER :: ncount,nst,ip,iloc,iq,mpi_ierror
    ncount=0
    globalindex=0
    node=0
    DO iq=1,2
      DO nst=npmin(iq),npsi(iq)
        node(nst)=MOD((nst-npmin(iq))/nb_psi,mpi_nprocs)
        IF(node(nst)==mpi_myproc) nstloc_iso(iq)=nstloc_iso(iq)+1
      END DO
    END DO
    nstloc=0
    DO nst=1,nstmax
       IF(node(nst)==mpi_myproc) THEN
          nstloc=nstloc+1
          globalindex(nstloc)=nst
       ENDIF
    ENDDO
    npmin_loc(1)=1
    npmin_loc(2)=nstloc_iso(1)+1
    npsi_loc(1)=nstloc_iso(1)
    npsi_loc(2)=nstloc
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
    CALL mpi_relate_comm
  END SUBROUTINE associate_nodes
!***********************************************************************
  SUBROUTINE mpi_relate_comm
    ALLOCATE(node_1dto2d_x(0:mpi_nprocs-1),node_1dto2d_y(0:mpi_nprocs-1),&
             node_2dto1d(0:mpi_size_x-1,0:mpi_size_y-1))
    node_1dto2d_x=0
    node_1dto2d_y=0
    node_2dto1d=0
    node_1dto2d_x(mpi_myproc)=mpi_rank_x
    node_1dto2d_y(mpi_myproc)=mpi_rank_y
    node_2dto1d(mpi_rank_x,mpi_rank_y)=mpi_myproc 
    CALL mpi_allreduce(MPI_IN_PLACE,node_1dto2d_x,mpi_nprocs,mpi_integer,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,node_1dto2d_y,mpi_nprocs,mpi_integer,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,node_2dto1d,mpi_nprocs,mpi_integer,mpi_sum,mpi_comm_world,mpi_ierror)
  END SUBROUTINE mpi_relate_comm
  !***********************************************************************
  SUBROUTINE collect_densities
    USE Densities, ONLY : rho,tau,current,sodens,sdens
    REAL(db) :: tmp_rho(nx,ny,nz,2),tmp_current(nx,ny,nz,3,2)
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
  SUBROUTINE collect_energies(delesum,sumflu)
    !***********************************************************************
    ! collects s.p. energies, fluctuation measure and energy differences
    !***********************************************************************
    REAL(db),INTENT(INOUT) :: delesum,sumflu
    REAL(db) :: tmpgat(nstmax),tmpgat3(3,nstmax)
    INTEGER  :: nst
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    DO nst=1,nstmax
      IF(node(nst)/=mpi_myproc) THEN
        sp_energy(nst)=0.0d0
        sp_efluct1(nst)=0.0d0
        sp_efluct2(nst)=0.0d0
      END IF
    ENDDO
    CALL mpi_allreduce(sp_energy,tmpgat,nstmax,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sp_energy=tmpgat
    CALL mpi_allreduce(sumflu,tmpgat(1),1,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    sumflu=tmpgat(1)
    CALL mpi_allreduce(delesum,tmpgat(1),1,mpi_double_precision,mpi_sum,  &
         mpi_comm_world,mpi_ierror)
    delesum=tmpgat(1)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
  END SUBROUTINE collect_energies
!***********************************************************************i
  SUBROUTINE mpi_wf_1d2x(psi,psi_x,iq)
    USE Trivial
    COMPLEX(db),INTENT(IN)    :: psi(:,:,:,:,:)
    COMPLEX(db),INTENT(OUT)   :: psi_x(:,:)
    INTEGER,    INTENT(IN)    :: iq
    INTEGER                   :: nst,ierr,is
    IF(.NOT.ALLOCATED(recvcounts)) THEN
      ALLOCATE(recvcounts(0:mpi_size_y-1,2),displs(0:mpi_size_y-1,2))
      recvcounts=0
      displs=0
      first(1)=1
      DO is=1,2
        DO nst=1,nstloc_x(is)
          recvcounts(node_1dto2d_y(node(globalindex_x(nst,is))),is)=&
          recvcounts(node_1dto2d_y(node(globalindex_x(nst,is))),is)+1
        END DO
      END DO
      first(2)=recvcounts(mpi_rank_y,1)+1
      recvcounts=recvcounts*size(psi(:,:,:,:,1))
      DO is=1,2
        DO nst=1,mpi_size_y-1
          displs(nst,is)=SUM(recvcounts(0:nst-1,is))
        END DO
      END DO
    END IF
    psi_x=0.0d0
    CALL mpi_allgatherv(psi(:,:,:,:,first(iq)),recvcounts(mpi_rank_y,iq),mpi_double_complex,&
                        psi_x,recvcounts(:,iq),displs(:,iq),mpi_double_complex,comm2d_y,ierr)
  END SUBROUTINE mpi_wf_1d2x
!***********************************************************************
  SUBROUTINE mpi_wf_x2y(psi_x,psi_y,iq)
    COMPLEX(db),INTENT(IN)    :: psi_x(:,:)
    COMPLEX(db),INTENT(OUT)   :: psi_y(:,:)
    INTEGER,    INTENT(IN)    :: iq
    INTEGER                   :: nst,ierr,is,lastnode,ip,first,last
    INTEGER,ALLOCATABLE,SAVE  :: rootnode(:,:),firstwf(:,:),nwf(:,:)
    IF(.NOT.ALLOCATED(rootnode)) THEN
      nst=MAX(nstloc_y(1)+1,nstloc_y(2)+1)
      ALLOCATE(rootnode(nst,2),firstwf(nst,2),nwf(nst,2))
      rootnode=0
      firstwf=0
      nwf=0
      DO is=1,2
        lastnode=-1
        ip=0
        DO nst=1,nstloc_y(is)
          IF(lastnode==node_x(globalindex_y(nst,is))) THEN
            nwf(ip,is)=nwf(ip,is)+1
          ELSE
            ip=ip+1
            firstwf(ip,is)=nst
            nwf(ip,is)=1
            lastnode=node_x(globalindex_y(nst,is))
            rootnode(ip,is)=lastnode
          END IF
        END DO
      END DO
    END IF
    psi_y=0.0d0
    ip=1
    DO WHILE (nwf(ip,iq)>0)
      first=firstwf(ip,iq)
      last=firstwf(ip,iq)+nwf(ip,iq)-1
      IF(rootnode(ip,iq)==mpi_rank_x) &
        psi_y(:,first:last)=&
          psi_x(:,localindex_x(globalindex_y(first,iq)):localindex_x(globalindex_y(last,iq)))
      CALL mpi_bcast(psi_y(:,first:last),size(psi_y(:,first:last)),mpi_double_complex,&
                     rootnode(ip,iq),comm2d_x,ierr)
      ip=ip+1
    END DO
  END SUBROUTINE mpi_wf_x2y
!***********************************************************************
  SUBROUTINE collect_wf_1d_x(psi,psi_x,iq)
  !*********************************************************************************
  ! adds all |psi> together and copies them to 1d distribution over nodes. 
  ! Copies each wave function at a time
  !*********************************************************************************
    USE Trivial
    INTEGER,INTENT(IN) :: iq
    COMPLEX(db),INTENT(INOUT) :: psi_x(:,:)
    COMPLEX(db),INTENT(OUT)   :: psi(:,:,:,:,:)
    INTEGER :: ierr
     CALL mpi_reduce_scatter(psi_x,psi(:,:,:,:,first(iq)),recvcounts(:,iq),mpi_double_complex,&
          mpi_sum,comm2d_y,ierr)
  END SUBROUTINE collect_wf_1d_x
  !***********************************************************************
  SUBROUTINE finish_mpi
    INTEGER :: ierr    
    CALL mpi_finalize(ierr)
  END SUBROUTINE finish_mpi
!***********************************************************************
  SUBROUTINE mpi_start_timer(index)
    INTEGER,INTENT(IN) :: index
    INTEGER :: ierr
    CALL mpi_barrier (mpi_comm_world,ierr)
    timer(index)=mpi_wtime()
  END SUBROUTINE mpi_start_timer
!***********************************************************************
  SUBROUTINE mpi_stop_timer(index,textline)
    INTEGER,INTENT(IN) :: index
    CHARACTER(*),INTENT(IN) :: textline
    INTEGER :: ierr
    CALL mpi_barrier (mpi_comm_world,ierr)
    IF(wflag)WRITE(*,'(A20,F10.4)')textline,mpi_wtime()-timer(index)
  END SUBROUTINE mpi_stop_timer
END MODULE Parallel
