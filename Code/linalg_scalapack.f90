MODULE LINALG
  USE Params,   ONLY: db,cmplxzero,cmplxone
  USE Levels
  USE Parallel, ONLY: nb,mb,contxt,contxt1d,nstloc_x,nstloc_y,globalindex_x,&
                      globalindex_y,nb_psi,npsi_loc,npmin_loc,psiloc_x,psiloc_y
  USE Grids,    ONLY: wxyz
!
  IMPLICIT NONE
  INTEGER                 :: nlin(2)
  INTEGER                 :: desca(2,10),descz(2,10),descc(2,10),desc_to(2,10),&
                             desc_psi1d(2,10),desc_psi2d(2,10),&
                             work_t_size(2),iwork_t_size(2),rwork_t_size(2)
!
  REAL(db)   ,ALLOCATABLE :: rwork_t(:),evals(:)  
  COMPLEX(db),ALLOCATABLE :: work_t(:),matr_lin(:,:),unitary(:,:)
  INTEGER    ,ALLOCATABLE :: iwork_t(:)
!
  CONTAINS
  !************************************************************
  SUBROUTINE init_linalg
    INTEGER                 :: iq,infoconv
    DO iq=1,2
      nlin(iq)=npsi(iq)-npmin(iq)+1
      CALL DESCINIT(DESCA(iq,1:10),npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,&
                    NB,MB,0,0,CONTXT,nstloc_x(iq),infoconv)
      CALL DESCINIT(DESCZ(iq,1:10),npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,&
                    NB,MB,0,0,CONTXT,nstloc_x(iq),infoconv)
      CALL DESCINIT(DESCC(iq,1:10),npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,&
                    NB,MB,0,0,CONTXT,nstloc_x(iq),infoconv)
      CALL DESCINIT(DESC_psi1D(iq,1:10),nx*ny*nz*2,npsi(iq)-npmin(iq)+1,nx*ny*nz*2,&
                      nb_psi,0,0,CONTXT1D,nx*ny*nz*2,infoconv)
      CALL DESCINIT(DESC_TO(iq,1:10),npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,&
                    npsi(iq)-npmin(iq)+1,nb_psi,0,0,CONTXT1D,npsi(iq)-npmin(iq)+1,infoconv)
      CALL DESCINIT(DESC_psi2D(iq,1:10),nx*ny*nz*2,npsi(iq)-npmin(iq)+1,&
                    NB,MB,0,0,CONTXT,psiloc_x(iq),infoconv)
      work_t_size(iq)  = -1
      iwork_t_size(iq) = -1
      rwork_t_size(iq) = -1
      ALLOCATE(work_t(1),iwork_t(1),rwork_t(1),matr_lin(nstloc_x(iq),nstloc_y(iq)),&
               unitary(nstloc_x(iq),nstloc_y(iq)),evals(nlin(iq)))
      CALL PZHEEVD('V','L',nlin(iq),matr_lin,1,1,DESCA(iq,1:10),evals,&
                   unitary,1,1,DESCZ(iq,1:10),work_t,work_t_size(iq),rwork_t,&
                   rwork_t_size(iq),iwork_t,iwork_t_size(iq),infoconv)
      work_t_size(iq) = INT(ABS(work_t(1)))
      iwork_t_size(iq) = INT(ABS(iwork_t(1)))
      rwork_t_size(iq) = INT(ABS(rwork_t(1)))
      DEALLOCATE(work_t,iwork_t,rwork_t,matr_lin,unitary,evals)
    END DO
  END SUBROUTINE init_linalg
  !************************************************************
  SUBROUTINE wf_1dto2d(psi_1d,psi_2d,iq)
    COMPLEX(db), INTENT(IN)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: psi_2d(:,:)
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_1d,1,1,desc_psi1d(iq,1:10),psi_2d,&
                  1,1,desc_psi2d(iq,1:10),contxt)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE wf_2dto1d(psi_2d,psi_1d,iq)
    COMPLEX(db), INTENT(OUT)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(IN)   :: psi_2d(:,:)
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_2d,1,1,desc_psi2d(iq,1:10),psi_1d,&
                  1,1,desc_psi1d(iq,1:10),contxt)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE calc_matrix(psi_1,psi_2,matrix,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: psi_1(:,:,:,:,:),psi_2(:,:,:,:,:)
    COMPLEX(db), INTENT(OUT):: matrix(:,:)
    INTEGER                 :: dim_x,dim_y
    COMPLEX(db),ALLOCATABLE :: psi_1_2d(:,:),psi_2_2d(:,:)
    ALLOCATE(psi_1_2d(psiloc_x(iq),psiloc_y(iq)),psi_2_2d(psiloc_x(iq),psiloc_y(iq)))    
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_1,1,1,desc_psi1d(iq,1:10),psi_1_2d,&
                  1,1,desc_psi2d(iq,1:10),contxt)
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_2,1,1,desc_psi1d(iq,1:10),psi_2_2d,&
                  1,1,desc_psi2d(iq,1:10),contxt)
    CALL PZGEMM('C','N',npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,nx*ny*nz*2,cmplxone,psi_1_2d,1,1,&
           desc_psi2d(iq,1:10),psi_2_2d,1,1,desc_psi2d(iq,1:10),cmplxzero,matrix,1,1,desca(iq,1:10))
    matrix=matrix*wxyz
    DEALLOCATE(psi_1_2d,psi_2_2d)
  END SUBROUTINE calc_matrix
  !************************************************************
  SUBROUTINE eigenvecs(matr_in,evecs,evals_out,iq)
    INTEGER,     INTENT(IN)           :: iq
    COMPLEX(db), INTENT(IN)           :: matr_in(:,:)
    COMPLEX(db), INTENT(OUT)          :: evecs(:,:)
    REAL(db),    INTENT(OUT),OPTIONAL :: evals_out(:)
    INTEGER                           :: infoconv
    ALLOCATE(work_t(work_t_size(iq)),rwork_t(rwork_t_size(iq)),&
             iwork_t(iwork_t_size(iq)),evals(nlin(iq)))
    CALL PZHEEVD('V','L',nlin(iq),matr_in,1,1,DESCA(iq,1:10),evals,evecs,1,1,DESCZ(iq,1:10),&
                  work_t,work_t_size(iq),rwork_t,rwork_t_size(iq),iwork_t,iwork_t_size(iq),infoconv)
    IF (PRESENT(evals_out))evals_out=evals
    DEALLOCATE(evals,work_t,iwork_t,rwork_t)
  END SUBROUTINE eigenvecs
  !************************************************************
  SUBROUTINE loewdin(imatr,smatr,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: imatr(:,:)
    COMPLEX(db),INTENT(OUT) :: smatr(:,:)
    COMPLEX(db),ALLOCATABLE :: tmatr(:,:)
    REAL(db)   ,ALLOCATABLE :: eigen_h(:)
    INTEGER                 :: i,it,jt,iy
    ALLOCATE(tmatr(nstloc_x(iq),nstloc_y(iq)),eigen_h(nlin(iq)))
    CALL eigenvecs(imatr,tmatr,eigen_h,iq)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO it = 1,nstloc_x(iq)
      DO jt = 1,nstloc_y(iq)
        iy = globalindex_y(jt,iq)-npmin(iq)+1
        tmatr(it,jt) = tmatr(it,jt)*sqrt(eigen_h(iy))
      ENDDO
    ENDDO 
    CALL PZGEMM('N','C',nlin(iq),nlin(iq),nlin(iq),cmplxone,tmatr,1,1,DESCA(iq,1:10),&
                tmatr,1,1,DESCZ(iq,1:10),cmplxzero,smatr,1,1,DESCC(iq,1:10))
    DEALLOCATE(tmatr,eigen_h)
  END SUBROUTINE
  !************************************************************
  SUBROUTINE comb_orthodiag(unitary_1,unitary_2,unitary,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: unitary_1(:,:),unitary_2(:,:)
    COMPLEX(db),INTENT(OUT) :: unitary(:,:)  
    CALL PZGEMM('T','T',nlin(iq),nlin(iq),nlin(iq),cmplxone,unitary_1,1,1,DESCA(iq,1:10),unitary_2,&
                    1,1,DESCZ(iq,1:10),cmplxzero,unitary,1,1,DESCC(iq,1:10))
  END SUBROUTINE
  !************************************************************
  SUBROUTINE recombine(psi_in,matrix,psi_out,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: psi_in(:,:,:,:,:),matrix(:,:)
    COMPLEX(db),INTENT(OUT) :: psi_out(:,:,:,:,:)
    COMPLEX(db),ALLOCATABLE :: psi_in_2d(:,:),psi_out_2d(:,:)
    ALLOCATE(psi_in_2d(psiloc_x(iq),psiloc_y(iq)),psi_out_2d(psiloc_x(iq),psiloc_y(iq)))    
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_in,1,1,desc_psi1d(iq,1:10),psi_in_2d,&
                  1,1,desc_psi2d(iq,1:10),contxt)
    CALL PZGEMM('N','T',nx*ny*nz*2,npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,cmplxone,psi_in_2d,1,1,desc_psi2d(iq,1:10),&
           matrix,1,1,desca(iq,1:10),cmplxzero,psi_out_2d,1,1,desc_psi2d(iq,1:10))
    CALL PZGEMR2D(nx*ny*nz*2,npsi(iq)-npmin(iq)+1,psi_out_2d,1,1,desc_psi2d(iq,1:10),psi_out,&
                  1,1,desc_psi1d(iq,1:10),contxt)
  END SUBROUTINE
END MODULE LINALG
