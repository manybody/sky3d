MODULE LINALG
  USE Params, ONLY: db,cmplxzero,cmplxone
  USE Levels
  USE Grids, ONLY: wxyz
!
  IMPLICIT NONE
  INTEGER :: nlin(2)  
!
  CONTAINS
  !************************************************************
  SUBROUTINE init_linalg
    INTEGER :: iq
    DO iq=1,2
      nlin(iq)=npsi(iq)-npmin(iq)+1
    END DO
  END SUBROUTINE init_linalg
  !************************************************************
  SUBROUTINE calc_matrix(psi_1,psi_2,matrix,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: psi_1(:,:,:,:,:),psi_2(:,:,:,:,:)
    COMPLEX(db), INTENT(OUT):: matrix(:,:)
    CALL ZGEMM('C','N',npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,nx*ny*nz*2,cmplxone,psi_1,nx*ny*nz*2,&
           psi_2,nx*ny*nz*2,cmplxzero,matrix,npsi(iq)-npmin(iq)+1)
           matrix=matrix*wxyz
  END SUBROUTINE calc_matrix
  !************************************************************
  SUBROUTINE eigenvecs(matr_in,evecs,evals_out,iq)
    INTEGER,     INTENT(IN)           :: iq
    COMPLEX(db), INTENT(IN)           :: matr_in(:,:)
    COMPLEX(db), INTENT(OUT)          :: evecs(:,:)
    REAL(db),    INTENT(OUT),OPTIONAL :: evals_out(:)
    INTEGER                           :: infoconv
    REAL(db),   ALLOCATABLE           :: evals(:)
    COMPLEX(db),ALLOCATABLE           :: cwork(:)
    REAL(db)   ,ALLOCATABLE           :: rwork(:)
    INTEGER    ,ALLOCATABLE           :: iwork(:)
    evecs=matr_in
    ALLOCATE(cwork(2*nlin(iq)*nlin(iq)),rwork(2*nlin(iq)*nlin(iq)+5*nlin(iq)+1),&
             iwork(5*nlin(iq)+3),evals(nlin(iq)))
    CALL zheevd('V','L',nlin(iq),evecs,nlin(iq),evals,cwork,nlin(iq)*nlin(iq)*2,&
                rwork,2*nlin(iq)*nlin(iq)+5*nlin(iq)+1,iwork,5*nlin(iq)+3,infoconv)
    IF (PRESENT(evals_out))evals_out=evals
    DEALLOCATE(cwork,rwork,iwork,evals)
  END SUBROUTINE eigenvecs
  !************************************************************
  SUBROUTINE loewdin(imatr,smatr,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: imatr(:,:)
    COMPLEX(db),INTENT(OUT) :: smatr(:,:)
    COMPLEX(db),ALLOCATABLE :: tmatr(:,:)
    REAL(db)   ,ALLOCATABLE :: eigen_h(:)
    INTEGER                 :: i
    ALLOCATE(tmatr(nlin(iq),nlin(iq)),eigen_h(nlin(iq)))
    CALL eigenvecs(imatr,tmatr,eigen_h,iq)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO i=1,nlin(iq)
      tmatr(:,i)=tmatr(:,i)*sqrt(eigen_h(i))
    END DO
    CALL zgemm('N','C',nlin(iq),nlin(iq),nlin(iq),cmplxone,tmatr,nlin(iq),&
               tmatr,nlin(iq),cmplxzero,smatr,nlin(iq))  
  END SUBROUTINE
  !************************************************************
  SUBROUTINE comb_orthodiag(unitary_1,unitary_2,unitary,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: unitary_1(:,:),unitary_2(:,:)
    COMPLEX(db),INTENT(OUT) :: unitary(:,:)  
    CALL zgemm('T','T',nlin(iq),nlin(iq),nlin(iq),cmplxone,&
               unitary_1,nlin(iq),unitary_2,nlin(iq),cmplxzero,unitary,nlin(iq))
  END SUBROUTINE
  !************************************************************
  SUBROUTINE recombine(psi_in,matrix,psi_out,iq)
    INTEGER,    INTENT(IN)  :: iq
    COMPLEX(db),INTENT(IN)  :: psi_in(:,:,:,:,:),matrix(:,:)
    COMPLEX(db),INTENT(OUT) :: psi_out(:,:,:,:,:)     
    CALL ZGEMM('N','T',nx*ny*nz*2,npsi(iq)-npmin(iq)+1,npsi(iq)-npmin(iq)+1,cmplxone,psi_in,nx*ny*nz*2,&
           matrix,npsi(iq)-npmin(iq)+1,cmplxzero,psi_out,nx*ny*nz*2)
  END SUBROUTINE
END MODULE LINALG
