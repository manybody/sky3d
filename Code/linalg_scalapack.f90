MODULE LINALG
  USE Params,   ONLY: db,cmplxzero,cmplxone
  USE Levels
  USE Parallel, ONLY: nb,mb,contxt,contxt1d,nstloc_x,nstloc_y,globalindex_x,&
                      globalindex_y,nb_psi,psiloc_x,psiloc_y,my_iso
  USE Grids,    ONLY: wxyz
!
  IMPLICIT NONE
  INTEGER                 :: nlin,gridsize
  INTEGER                 :: desca(10),descz(10),descc(10),desc_to(10),&
                             desc_psi1d(10),desc_psi2d(10),&
                             work_t_size,iwork_t_size,rwork_t_size
!
  REAL(db)   ,ALLOCATABLE :: rwork_t(:),evals(:)  
  COMPLEX(db),ALLOCATABLE :: work_t(:),matr_lin(:,:),unitary(:,:)
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
    work_t_size  = -1
    iwork_t_size = -1
    rwork_t_size = -1
    ALLOCATE(work_t(1),iwork_t(1),rwork_t(1),matr_lin(nstloc_x,nstloc_y),&
             unitary(nstloc_x,nstloc_y),evals(nlin))
    CALL PZHEEVD('V','L',nlin,matr_lin,1,1,DESCA(1:10),evals,&
                 unitary,1,1,DESCZ(1:10),work_t,work_t_size,rwork_t,&
                 rwork_t_size,iwork_t,iwork_t_size,infoconv)
    work_t_size = INT(ABS(work_t(1)))
    iwork_t_size = INT(ABS(iwork_t(1)))
    rwork_t_size = INT(ABS(rwork_t(1)))
    DEALLOCATE(work_t,iwork_t,rwork_t,matr_lin,unitary,evals)

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
    INTEGER                           :: infoconv
    ALLOCATE(work_t(work_t_size),rwork_t(rwork_t_size),&
             iwork_t(iwork_t_size),evals(nlin))
    CALL PZHEEVD('V','L',nlin,matr_in,1,1,DESCA(1:10),evals,evecs,1,1,DESCZ(1:10),&
                  work_t,work_t_size,rwork_t,rwork_t_size,iwork_t,iwork_t_size,infoconv)
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
    ALLOCATE(tmatr(nstloc_x,nstloc_y),eigen_h(nlin))
    CALL eigenvecs(imatr,tmatr,eigen_h)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO it = 1,nstloc_x
      DO jt = 1,nstloc_y
        iy = globalindex_y(jt)-npmin(my_iso)+1
        tmatr(it,jt) = tmatr(it,jt)*sqrt(eigen_h(iy))
      ENDDO
    ENDDO 
    CALL PZGEMM('N','C',nlin,nlin,nlin,cmplxone,tmatr,1,1,DESCA(1:10),&
                tmatr,1,1,DESCZ(1:10),cmplxzero,smatr,1,1,DESCC(1:10))
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
