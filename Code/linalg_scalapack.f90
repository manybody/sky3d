MODULE LINALG
  USE Params,   ONLY: db,cmplxzero,cmplxone
  USE Levels
  USE Parallel, ONLY: nb,mb,contxt,contxt1d,nstloc_x,nstloc_y,globalindex_x,&
                      globalindex_y,nb_psi,psiloc_x,psiloc_y,my_iso,my_diag,&
                      nstloc_diag_x,nstloc_diag_y,contxt_do,mpi_myproc,contxt_d,contxt_o,&
                      globalindex_diag_x,globalindex_diag_y,my_diag
  USE Grids,    ONLY: wxyz
!
  IMPLICIT NONE
  INTEGER                 :: nlin,gridsize
  INTEGER                 :: desca(10),descz(10),descc(10),desc_to(10),&
                             desc_psi1d(10),desc_psi2d(10),desc_diag(10),desc_ortho(10),&
                             work_t_size,iwork_t_size,rwork_t_size,desc_d(10),desc_o(10)
!
  REAL(db)                :: VL,VU
  INTEGER                 :: IL,IU,M,NZ_t
  REAL(db)   ,ALLOCATABLE :: rwork_t(:),evals(:)  
  COMPLEX(db),ALLOCATABLE :: work_t(:),matr_lin(:,:),unitary(:,:),matr_lin_d(:,:),unitary_d(:,:)
  INTEGER    ,ALLOCATABLE :: iwork_t(:)
!
  CONTAINS
  !************************************************************
  SUBROUTINE init_linalg
    INTEGER                 :: infoconv
    nlin=npsi(my_iso)-npmin(my_iso)+1
    gridsize=nx*ny*nz*2
    CALL DESCINIT(DESCA(1:10),nlin,nlin,NB,MB,0,0,CONTXT,nstloc_x,infoconv)
    CALL DESCINIT(DESCZ(1:10),nlin,nlin,NB,MB,0,0,CONTXT,nstloc_x,infoconv)
    CALL DESCINIT(DESCC(1:10),nlin,nlin,NB,MB,0,0,CONTXT,nstloc_x,infoconv)
    CALL DESCINIT(DESC_psi1D(1:10),gridsize,nlin,gridsize,nb_psi,0,0,CONTXT1D,gridsize,infoconv)
    CALL DESCINIT(DESC_TO(1:10),nlin,nlin,nlin,nb_psi,0,0,CONTXT1D,nlin,infoconv)
    CALL DESCINIT(DESC_psi2D(1:10),gridsize,nlin,NB,MB,0,0,CONTXT,psiloc_x,infoconv)
    CALL DESCINIT(DESC_D(1:10),nlin,nlin,NB,MB,0,0,CONTXT_D,nstloc_diag_x,infoconv)
    CALL DESCINIT(DESC_O(1:10),nlin,nlin,NB,MB,0,0,CONTXT_O,nstloc_diag_x,infoconv)
    work_t_size  = -1
    iwork_t_size = -1
    rwork_t_size = -1
    ALLOCATE(work_t(1),iwork_t(1),rwork_t(1),matr_lin(nstloc_x,nstloc_y),&
             unitary(nstloc_x,nstloc_y),evals(nlin))
    ALLOCATE(matr_lin_d(nstloc_diag_x,nstloc_diag_y),unitary_d(nstloc_diag_x,nstloc_diag_y))


    IF(my_diag==1.OR.my_diag==3)THEN
      CALL PZHEEVR('V','A','L',nlin,matr_lin_d,1,1,DESC_D(1:10),VL,VU,IL,IU,&
                M,NZ_t,evals,unitary_d,1,1,DESC_D(1:10),work_t,work_t_size,rwork_t,&
                 rwork_t_size,iwork_t,iwork_t_size,infoconv)
    ELSE IF(my_diag==2.OR.my_diag==4)THEN
      CALL PZHEEVR('V','A','L',nlin,matr_lin_d,1,1,DESC_O(1:10),VL,VU,IL,IU,&
                M,NZ_t,evals,unitary_d,1,1,DESC_O(1:10),work_t,work_t_size,rwork_t,&
                 rwork_t_size,iwork_t,iwork_t_size,infoconv)
    ENDIF
    work_t_size = INT(ABS(work_t(1)))
    iwork_t_size = INT(ABS(iwork_t(1)))
    rwork_t_size = INT(ABS(rwork_t(1)))
    DEALLOCATE(work_t,iwork_t,rwork_t,matr_lin,unitary,evals,matr_lin_d,unitary_d)

  END SUBROUTINE init_linalg
  !************************************************************
  SUBROUTINE wf_1dto2d(psi_1d,psi_2d)
    COMPLEX(db), INTENT(IN)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: psi_2d(:,:)
    CALL PZGEMR2D(gridsize,nlin,psi_1d,1,1,desc_psi1d(1:10),psi_2d,&
                  1,1,desc_psi2d(1:10),contxt)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE wf_2dto1d(psi_2d,psi_1d)
    COMPLEX(db), INTENT(OUT)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(IN)   :: psi_2d(:,:)
    CALL PZGEMR2D(gridsize,nlin,psi_2d,1,1,desc_psi2d(1:10),psi_1d,&
                  1,1,desc_psi1d(1:10),contxt)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE matrix_split(rhomatr_lin,hmatr_lin,rhomatr_lin_d,hmatr_lin_d)
    COMPLEX(db), INTENT(OUT)  :: rhomatr_lin_d(:,:),hmatr_lin_d(:,:)
    COMPLEX(db), INTENT(IN)   :: rhomatr_lin(:,:),hmatr_lin(:,:)
    CALL PZGEMR2D(nlin,nlin,rhomatr_lin,1,1,desca,rhomatr_lin_d,1,1,desc_o,contxt)
    CALL PZGEMR2D(nlin,nlin,hmatr_lin,1,1,desca,hmatr_lin_d,1,1,desc_d,contxt)    
  END SUBROUTINE matrix_split
  !************************************************************
  SUBROUTINE matrix_gather(unitary_rho,unitary_h,unitary_rho_d,unitary_h_d)
    COMPLEX(db), INTENT(IN)  :: unitary_rho_d(:,:),unitary_h_d(:,:)
    COMPLEX(db), INTENT(OUT) :: unitary_rho(:,:),unitary_h(:,:)
    CALL PZGEMR2D(nlin,nlin,unitary_rho_d,1,1,desc_o,unitary_rho,1,1,desca,contxt)
    CALL PZGEMR2D(nlin,nlin,unitary_h_d,1,1,desc_d,unitary_h,1,1,desca,contxt)
  END SUBROUTINE matrix_gather
  !************************************************************
  SUBROUTINE calc_matrix(psi_1,psi_2,matrix)
    COMPLEX(db),INTENT(IN)  :: psi_1(:,:),psi_2(:,:)
    COMPLEX(db), INTENT(OUT):: matrix(:,:)
    CALL PZGEMM('C','N',nlin,nlin,gridsize,cmplxone,psi_1,1,1,&
           desc_psi2d(1:10),psi_2,1,1,desc_psi2d(1:10),cmplxzero,matrix,1,1,desca(1:10))
    matrix=matrix*wxyz
  END SUBROUTINE calc_matrix
  !************************************************************
  SUBROUTINE eigenvecs(matr_in,evecs,evals_out)
    COMPLEX(db), INTENT(IN)           :: matr_in(:,:)
    COMPLEX(db), INTENT(OUT)          :: evecs(:,:)
    REAL(db),    INTENT(OUT),OPTIONAL :: evals_out(:)
    COMPLEX(db), ALLOCATABLE          :: matr_in_d(:,:),evecs_d(:,:)
    INTEGER                           :: infoconv
    ALLOCATE(work_t(work_t_size),rwork_t(rwork_t_size),&
             iwork_t(iwork_t_size),evals(nlin))
    IF(PRESENT(evals_out)) THEN
     CALL PZHEEVR('V','A','L',nlin,matr_in,1,1,DESC_O(1:10),VL,VU,IL,IU,M,NZ_t,evals,evecs,1,1,DESC_O(1:10),&
                  work_t,work_t_size,rwork_t,rwork_t_size,iwork_t,iwork_t_size,infoconv)
    ELSE
      CALL PZHEEVR('V','A','L',nlin,matr_in,1,1,DESC_D(1:10),VL,VU,IL,IU,M,NZ_t,evals,evecs,1,1,DESC_D(1:10),&
                  work_t,work_t_size,rwork_t,rwork_t_size,iwork_t,iwork_t_size,infoconv)
    ENDIF  
    IF (PRESENT(evals_out))evals_out=evals
    DEALLOCATE(evals,work_t,iwork_t,rwork_t)
  END SUBROUTINE eigenvecs
  !************************************************************
  SUBROUTINE loewdin(imatr,smatr)
    COMPLEX(db),INTENT(IN)  :: imatr(:,:)
    COMPLEX(db),INTENT(OUT) :: smatr(:,:)
    COMPLEX(db),ALLOCATABLE :: tmatr(:,:)
    REAL(db)   ,ALLOCATABLE :: eigen_h(:)
    INTEGER                 :: it,jt,iy
    ALLOCATE(tmatr(nstloc_diag_x,nstloc_diag_y),eigen_h(nlin))
    CALL eigenvecs(imatr,tmatr,eigen_h)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO it = 1,nstloc_diag_x
      DO jt = 1,nstloc_diag_y
        iy = globalindex_diag_y(jt)-npmin(my_iso)+1
        tmatr(it,jt) = tmatr(it,jt)*sqrt(eigen_h(iy))
      ENDDO
    ENDDO 
    CALL PZGEMM('N','C',nlin,nlin,nlin,cmplxone,tmatr,1,1,DESC_O(1:10),&
                tmatr,1,1,DESC_O(1:10),cmplxzero,smatr,1,1,DESC_O(1:10))
    DEALLOCATE(tmatr,eigen_h)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE comb_orthodiag(unitary_1,unitary_2,unitary)
    COMPLEX(db),INTENT(IN)  :: unitary_1(:,:),unitary_2(:,:)
    COMPLEX(db),INTENT(OUT) :: unitary(:,:)  
    CALL PZGEMM('T','T',nlin,nlin,nlin,cmplxone,unitary_1,1,1,DESCA(1:10),unitary_2,&
                    1,1,DESCZ(1:10),cmplxzero,unitary,1,1,DESCC(1:10))
  END SUBROUTINE
  !************************************************************
  SUBROUTINE recombine(psi_in,matrix,psi_out)
    COMPLEX(db),INTENT(IN)  :: psi_in(:,:),matrix(:,:)
    COMPLEX(db),INTENT(OUT) :: psi_out(:,:)
    CALL PZGEMM('N','T',gridsize,nlin,nlin,cmplxone,psi_in,1,1,desc_psi2d(1:10),&
           matrix,1,1,desca(1:10),cmplxzero,psi_out,1,1,desc_psi2d(1:10))
  END SUBROUTINE
END MODULE LINALG
