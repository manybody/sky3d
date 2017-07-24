MODULE LINALG
  USE Params, ONLY: db,cmplxzero,cmplxone
  USE Levels
  USE Grids, ONLY: wxyz
  USE Parallel, ONLY : globalindex_diag_x,globalindex_diag_y,globalindex_x,globalindex_y,&
                       psiloc_x,psiloc_y,nstloc_diag_x,nstloc_diag_y,nstloc_x,nstloc_y
!
  IMPLICIT NONE
  INTEGER :: nlin           !<number of wave functions for local isospin.
!
  CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: init_linalg
!> @brief
!!This subroutine initializes BLACS descriptors for ScaLAPACK routines and 
!!initializes diagonalization routines.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_linalg(is)
    INTEGER,INTENT(IN) :: is
    INTEGER            :: i,noffset
    IF(is==0) RETURN
    nlin=npsi(is)-npmin(is)+1
    noffset=npmin(is)-1
    psiloc_x=nlin
    psiloc_y=nlin
    nstloc_x=nlin
    nstloc_y=nlin
    nstloc_diag_x=nlin
    nstloc_diag_y=nlin
    DO i=npmin(is),npsi(is)
      globalindex_x(i-noffset)=i
      globalindex_y(i-noffset)=i
      globalindex_diag_x(i-noffset)=i
      globalindex_diag_y(i-noffset)=i
    END DO
  END SUBROUTINE init_linalg
!---------------------------------------------------------------------------  
! DESCRIPTION: wf_1dto2d
!> @brief
!!This subroutine transfers wave functions from 1d distribution to 2d distribution.
!--------------------------------------------------------------------------- 
  SUBROUTINE wf_1dto2d(psi_1d,psi_2d)
    COMPLEX(db), INTENT(IN)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(OUT) :: psi_2d(:,:)
    STOP 'Should not be called in sequential version'
  END SUBROUTINE
!---------------------------------------------------------------------------  
! DESCRIPTION: wf_2dto1d
!> @brief
!!This subroutine transfers wave functions from 2d distribution to 1d distribution.
!--------------------------------------------------------------------------- 
  SUBROUTINE wf_2dto1d(psi_2d,psi_1d)
    COMPLEX(db), INTENT(OUT)  :: psi_1d(:,:,:,:,:)
    COMPLEX(db), INTENT(IN)   :: psi_2d(:,:)
    STOP 'Should not be called in sequential version'
  END SUBROUTINE
!---------------------------------------------------------------------------  
! DESCRIPTION: matrix_split
!> @brief
!!This subroutine transfers matrices to isospin subgroups
!--------------------------------------------------------------------------- 
  SUBROUTINE matrix_split(rhomatr_lin,hmatr_lin,rhomatr_lin_d,hmatr_lin_d)
    COMPLEX(db), INTENT(OUT)  :: rhomatr_lin_d(:,:),hmatr_lin_d(:,:)
    COMPLEX(db), INTENT(IN)   :: rhomatr_lin(:,:),hmatr_lin(:,:)
    STOP 'Should not be called in sequential version'
  END SUBROUTINE matrix_split
!---------------------------------------------------------------------------  
! DESCRIPTION: matrix_gather
!> @brief
!!This subroutine transfers matrices from isospin subgroups to full isospin groups.
!--------------------------------------------------------------------------- 
  SUBROUTINE matrix_gather(unitary_rho,unitary_h,unitary_rho_d,unitary_h_d)
    COMPLEX(db), INTENT(IN)  :: unitary_rho_d(:,:),unitary_h_d(:,:)
    COMPLEX(db), INTENT(OUT) :: unitary_rho(:,:),unitary_h(:,:)
    STOP 'Should not be called in sequential version'
  END SUBROUTINE matrix_gather
!---------------------------------------------------------------------------  
! DESCRIPTION: calc_matrix
!> @brief
!!This subroutine calculates overlap or Hamiltonian matrices.
!--------------------------------------------------------------------------- 
  SUBROUTINE calc_matrix(psi_1,psi_2,matrix)
    COMPLEX(db),INTENT(IN)  :: psi_1(:,:),psi_2(:,:)
    COMPLEX(db), INTENT(OUT):: matrix(:,:)
    CALL ZGEMM('C','N',nlin,nlin,nx*ny*nz*2,cmplxone,psi_1,nx*ny*nz*2,&
           psi_2,nx*ny*nz*2,cmplxzero,matrix,nlin)
    matrix=matrix*wxyz
  END SUBROUTINE calc_matrix
!---------------------------------------------------------------------------  
! DESCRIPTION: eigenvecs
!> @brief
!!This subroutine calculates eigenvectors and optional eigen values.
!--------------------------------------------------------------------------- 
  SUBROUTINE eigenvecs(matr_in,evecs,evals_out)
    COMPLEX(db), INTENT(IN)           :: matr_in(:,:)
    COMPLEX(db), INTENT(OUT)          :: evecs(:,:)
    REAL(db),    INTENT(OUT),OPTIONAL :: evals_out(:)
    INTEGER                           :: infoconv
    REAL(db),   ALLOCATABLE           :: evals(:)
    COMPLEX(db),ALLOCATABLE           :: cwork(:)
    REAL(db)   ,ALLOCATABLE           :: rwork(:)
    INTEGER    ,ALLOCATABLE           :: iwork(:)
    evecs=matr_in
    ALLOCATE(cwork(2*nlin*nlin),rwork(2*nlin*nlin+5*nlin+1),&
             iwork(5*nlin+3),evals(nlin))
    CALL zheevd('V','L',nlin,evecs,nlin,evals,cwork,nlin*nlin*2,&
                rwork,2*nlin*nlin+5*nlin+1,iwork,5*nlin+3,infoconv)
    IF (PRESENT(evals_out))evals_out=evals
    DEALLOCATE(cwork,rwork,iwork,evals)
  END SUBROUTINE eigenvecs
!---------------------------------------------------------------------------  
! DESCRIPTION: loewdin
!> @brief
!!This subroutine calculates the Loewdin matrix.
!--------------------------------------------------------------------------- 
  SUBROUTINE loewdin(imatr,smatr)
    COMPLEX(db),INTENT(IN)  :: imatr(:,:)
    COMPLEX(db),INTENT(OUT) :: smatr(:,:)
    COMPLEX(db),ALLOCATABLE :: tmatr(:,:)
    REAL(db)   ,ALLOCATABLE :: eigen_h(:)
    INTEGER                 :: i
    ALLOCATE(tmatr(nlin,nlin),eigen_h(nlin))
    CALL eigenvecs(imatr,tmatr,eigen_h)
    eigen_h=1.0d0/sqrt(eigen_h)
    DO i=1,nlin
      tmatr(:,i)=tmatr(:,i)*sqrt(eigen_h(i))
    END DO
    CALL zgemm('N','C',nlin,nlin,nlin,cmplxone,tmatr,nlin,&
               tmatr,nlin,cmplxzero,smatr,nlin)  
  END SUBROUTINE
!---------------------------------------------------------------------------  
! DESCRIPTION: comb_orthodiag
!> @brief
!!This subroutine combines Loewdin and diagonalization matrices.
!--------------------------------------------------------------------------- 
  SUBROUTINE comb_orthodiag(unitary_1,unitary_2,unitary)
    COMPLEX(db),INTENT(IN)  :: unitary_1(:,:),unitary_2(:,:)
    COMPLEX(db),INTENT(OUT) :: unitary(:,:)  
    CALL zgemm('T','T',nlin,nlin,nlin,cmplxone,&
               unitary_1,nlin,unitary_2,nlin,cmplxzero,unitary,nlin)
  END SUBROUTINE
!---------------------------------------------------------------------------  
! DESCRIPTION: recombine
!> @brief
!!This subroutine performs matrix vector multiplication for .
!--------------------------------------------------------------------------- 
  SUBROUTINE recombine(psi_in,matrix,psi_out)
    COMPLEX(db),INTENT(IN)  :: psi_in(:,:),matrix(:,:)
    COMPLEX(db),INTENT(OUT) :: psi_out(:,:)  
    CALL ZGEMM('N','T',nx*ny*nz*2,nlin,nlin,cmplxone,psi_in,nx*ny*nz*2,&
           matrix,nlin,cmplxzero,psi_out,nx*ny*nz*2)

  END SUBROUTINE
END MODULE LINALG
