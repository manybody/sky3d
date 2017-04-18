MODULE Parallel
  USE Params, ONLY : wflag
  USE Grids, ONLY: nx,ny,nz
  USE Levels
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  SAVE
  INTEGER, PARAMETER   :: NB=2,MB=2,NB_psi =2
  LOGICAL, PARAMETER   :: tmpi=.TRUE.,ttabc=.FALSE.
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:)
  INTEGER, ALLOCATABLE :: recvcounts(:,:),displs(:,:),globalindex_x(:),globalindex_y(:),&
                          globalindex_diag_x(:),globalindex_diag_y(:)
  INTEGER              :: mpi_nprocs,mpi_ierror,mpi_myproc,mpi_nprocs_iso(2),my_iso,mpi_myproc_iso
  INTEGER              :: comm2d,comm_iso,mpi_dims(2),mpi_mycoords(2),nstloc_x,&
                          nstloc_y,psiloc_x,psiloc_y,mpi_nprocs_diag(4),my_diag,&
                          nstloc_diag_x,nstloc_diag_y
  INTEGER              :: comm2d_x,comm2d_y,mpi_size_x,mpi_size_y,mpi_rank_x,mpi_rank_y,&
                          nprd,npcd,myprowd,mypcold,contxt_d,contxt_o
  INTEGER              :: NPROCS,NPROW,NPCOL,MYPROW,MYPCOL,CONTXT,IAM,CONTXT1D,CONTXT_DO
  
!  INTEGER              :: CONTXT_DP,CONTXT_OP,CONTXT_DN,CONTXT_ON,NPROW_DP,NPCOL-DP,&
                          
   
  INTEGER, EXTERNAL    :: NUMROC,INDXL2G,INDXG2L,INDXG2P
  REAL(db)             :: timer(20)
CONTAINS
  !***********************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax),&
             globalindex_x(nstmax),globalindex_y(nstmax),&
             globalindex_diag_x(nstmax),globalindex_diag_y(nstmax))
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
!      WRITE(*,*),'proc = ',mpi_myproc,my_iso,npsi(my_iso),npmin(my_iso),nstloc_x,nstloc_y
      nstloc_diag_x = NUMROC(npsi(my_iso)-npmin(my_iso)+1,NB,myprowd,0,nprd)
      nstloc_diag_y = NUMROC(npsi(my_iso)-npmin(my_iso)+1,NB,mypcold,0,npcd)
!      WRITE(*,*),'proc = ',mpi_myproc,'diag dim = ',nstloc_diag_x,nstloc_diag_y
   

!
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
!***************************************************************************
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

!      WRITE(*,*),'proc = ',mpi_myproc,'contxt = ',CONTXT_DO


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
      WRITE(*,*),'proc = ',mpi_myproc,nprd,npcd,myprowd,mypcold


      CALL mpi_barrier (mpi_comm_world, mpi_ierror)
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
    mpi_nprocs_iso(1)=NINT(REAL(npsi(1))/REAL(npsi(2))/2.0*mpi_nprocs)*2
    mpi_nprocs_iso(2)=mpi_nprocs-mpi_nprocs_iso(1)
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
    WRITE(*,*),'proc = ',mpi_myproc,'iso = ',my_iso,'my diag = ',my_diag

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
!***********************************************************************
  SUBROUTINE collect_densities
    USE Densities, ONLY : rho,tau,current,sodens,sdens
    CALL mpi_barrier (mpi_comm_world, mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,rho,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,tau,2*nx*ny*nz,        &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,current,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sodens,3*2*nx*ny*nz,  &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_allreduce(MPI_IN_PLACE,sdens,3*2*nx*ny*nz,   &
         mpi_double_precision,mpi_sum,mpi_comm_world,mpi_ierror)
    CALL mpi_barrier (mpi_comm_world,mpi_ierror)
    RETURN
  END SUBROUTINE collect_densities
  !***********************************************************************
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
!***********************************************************************
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
!***********************************************************************
  SUBROUTINE mpi_start_timer_iq(index)
    INTEGER,INTENT(IN) :: index
    INTEGER :: ierr
    CALL mpi_barrier (comm_iso,ierr)
    timer(index)=mpi_wtime()
  END SUBROUTINE mpi_start_timer_iq
!***********************************************************************
  SUBROUTINE mpi_stop_timer_iq(index,textline)
    INTEGER,INTENT(IN) :: index
    CHARACTER(*),INTENT(IN) :: textline
    INTEGER :: ierr
    CALL mpi_barrier (comm_iso,ierr)
    IF(mpi_myproc_iso==0)WRITE(*,'(A20,F10.4,A4,I4)')textline,mpi_wtime()-timer(index),' iq = ',my_iso
  END SUBROUTINE mpi_stop_timer_iq
END MODULE Parallel
