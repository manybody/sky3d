!------------------------------------------------------------------------------
! MODULE: Trivial
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains some basic calculations with wave functions and
!!densities, which are similar enough to be grouped together.
!------------------------------------------------------------------------------
MODULE Trivial
  USE Params, ONLY: db
  USE Grids, ONLY: nx,ny,nz,wxyz
  IMPLICIT NONE  
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: cmulx
!> @brief
!!The routine does a matrix multiplication of a real matrix with a complex 
!!wave function along the x-direction. 
!!\c ifadd allows accumulating results but is not used in the present code.
!>
!> @details
!!The operation carried out here includes a loop over the spin index.
!>
!> @param[in] xmat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the x-direction.
!> @param[in] pinn
!> COMPLEX(db), array, takes the input wave function.
!> @param[out] pout
!> COMPLEX(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if nonzero, the result of the multiplication is added to the 
!! input wave function.
!--------------------------------------------------------------------------- 

  PURE SUBROUTINE cmulx(xmat,pinn,pout,ifadd)  
    REAL(db) :: xmat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: xmat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz
    IF(ifadd==0) pout=0.0D0
    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(xmat(ix,:)*pinn(:,iy,iz,is))
    END FORALL
  END SUBROUTINE cmulx
!---------------------------------------------------------------------------  
! DESCRIPTION: cmuly
!> @brief
!!The routine does a matrix multiplication of a real matrix with a complex 
!!wave function along the y-direction. 
!!\c ifadd allows accumulating results but is not used in the present code.
!>
!> @details
!!The operation carried out here includes a loop over the spin index.
!>
!> @param[in] ymat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the y-direction.
!> @param[in] pinn
!> COMPLEX(db), array, takes the input wave function.
!> @param[out] pout
!> COMPLEX(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if nonzero, the result of the multiplication is added to the 
!! input wave function.
!--------------------------------------------------------------------------- 
  PURE SUBROUTINE cmuly(ymat,pinn,pout,ifadd)  
    REAL(db) :: ymat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: ymat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz
    IF(ifadd==0) pout=0.0D0
    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(ymat(iy,:)*pinn(ix,:,iz,is))
    END FORALL
  END SUBROUTINE cmuly
!---------------------------------------------------------------------------  
! DESCRIPTION: cmulz
!> @brief
!!The routine does a matrix multiplication of a real matrix with a complex 
!!wave function along the z-direction. 
!!\c ifadd allows accumulating results but is not used in the present code.
!>
!> @details
!!The operation carried out here includes a loop over the spin index.
!>
!> @param[in] zmat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the z-direction.
!> @param[in] pinn
!> COMPLEX(db), array, takes the input wave function.
!> @param[out] pout
!> COMPLEX(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if nonzero, the result of the multiplication is added to the 
!! input wave function.
!--------------------------------------------------------------------------- 
  PURE SUBROUTINE cmulz(zmat,pinn,pout,ifadd)  
    REAL(db) :: zmat(:,:)  
    COMPLEX(db) :: pinn(:,:,:,:),pout(:,:,:,:)  
    INTEGER :: ifadd
    INTENT(IN) :: zmat,pinn,ifadd
    INTENT(INOUT) :: pout
    INTEGER :: is,ix,iy,iz,izz
    IF(ifadd==0) pout=0.0D0
    !
    !    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
    !       pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+SUM(zmat(iz,:)*pinn(ix,iy,:,is))
    !    END FORALL
    DO is=1,2
       DO iz=1,nz
          DO izz=1,nz
             DO iy=1,ny
                DO ix=1,nx
                   pout(ix,iy,iz,is)=pout(ix,iy,iz,is)+ &
                        zmat(iz,izz)*pinn(ix,iy,izz,is)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE cmulz
!---------------------------------------------------------------------------  
! DESCRIPTION: rpsnorm
!> @brief
!!This function with one argument calculates the norm of a wave function 
!>
!> @details
!!i.e.,
!!\f[ |\psi|=\sum_{s=\pm{1\over2}}\int |\psi(\vec r,s)|^2\,\D^3 r\rightarrow
!!{\tt wxyz}\sum_{s=1}^2 \sum_{i=1}^{nx}\sum_{j=1}^{ny}\sum_{k=1}^{nz}
!!|{\tt psi(i,j,k,s)}|^2, \f]
!!with \c wxyz the volume element.
!>
!> @param[in] ps
!> COMPLEX(db), array, takes the wave function.
!--------------------------------------------------------------------------- 
  PURE FUNCTION rpsnorm(ps) RESULT(r)
    COMPLEX(db),INTENT(IN) :: ps(:,:,:,:)
    REAL(db) :: r
    r=wxyz*SUM(REAL(CONJG(ps)*ps))
  END FUNCTION rpsnorm
!---------------------------------------------------------------------------  
! DESCRIPTION: overlap
!> @brief
!!This calculates the overlap of two wave functions \f$ psi_L \f$ and \f$ psi_R \f$
!>
!> @details
!!This is defined as
!!\f{eqnarray}{
!!  \langle\psi_L|\psi_R\rangle&=&\sum_{s=\pm{1\over2}}\int \psi_L^*(\vec
!!  r)\,\psi_R(\vec r)\,\D^3 r \\ 
!!  &\rightarrow&
!!  {\tt wxyz}\sum_{s=1}^2
!!  \sum_{i=1}^{nx}\sum_{j=1}^{ny}\sum_{k=1}^{nz} 
!!  {\tt CONJG(pl(i,j,k,s))pr(i,j,k,s)}.
!!\f}
!>
!> @param[in] pl
!> COMPEX(db), array, takes \f$ psi_L \f$.
!> @param[in] pr
!> COMPEX(db), array, takes \f$ psi_R \f$.
!--------------------------------------------------------------------------- 
  PURE FUNCTION overlap(pl,pr)  RESULT(c)
    COMPLEX(db) :: c,pl(:,:,:,:),pr(:,:,:,:)
    INTENT(IN) :: pl,pr
    c=wxyz*SUM(CONJG(pl)*pr)
  END FUNCTION overlap
!---------------------------------------------------------------------------  
! DESCRIPTION: rmulx
!> @brief
!!The routine does a matrix multiplication of a real matrix with a real 
!!field along the x-direction. 
!!\c ifadd allows accumulating or substracting results.
!>
!> @param[in] xmat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the x-direction.
!> @param[in] finn
!> REAL(db), array, takes the input wave function.
!> @param[out] fout
!> REAL(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if positive, the result of the multiplication is added to the 
!! input wave function, if negative it is substracted.
!--------------------------------------------------------------------------- 
  PURE SUBROUTINE rmulx(xmat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: xmat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: xmat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(xmat(ix,:)*finn(:,iy,iz))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(xmat(ix,:)*finn(:,iy,iz))
       END FORALL
    ENDIF
  END SUBROUTINE rmulx
!---------------------------------------------------------------------------  
! DESCRIPTION: rmuly
!> @brief
!!The routine does a matrix multiplication of a real matrix with a real 
!!field along the y-direction. 
!!\c ifadd allows accumulating or substracting results.
!>
!> @param[in] ymat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the y-direction.
!> @param[in] finn
!> REAL(db), array, takes the input wave function.
!> @param[out] fout
!> REAL(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if positive, the result of the multiplication is added to the 
!! input wave function, if negative it is substracted.
!--------------------------------------------------------------------------- 
  PURE SUBROUTINE rmuly(ymat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: ymat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: ymat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(ymat(iy,:)*finn(ix,:,iz))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(ymat(iy,:)*finn(ix,:,iz))
       END FORALL
    ENDIF
  END SUBROUTINE rmuly
!---------------------------------------------------------------------------  
! DESCRIPTION: rmulz
!> @brief
!!The routine does a matrix multiplication of a real matrix with a real 
!!field along the z-direction. 
!!\c ifadd allows accumulating or substracting results.
!>
!> @param[in] zmat
!> REAL(db), takes the matrix dimensioned 
!! as a square matrix with the number of points in the z-direction.
!> @param[in] finn
!> REAL(db), array, takes the input wave function.
!> @param[out] fout
!> REAL(db), array, returns the multiplied wave function.
!> @param[in] ifadd
!> INTEGER, if positive, the result of the multiplication is added to the 
!! input wave function, if negative it is substracted.
!---------------------------------------------------------------------------   
  PURE SUBROUTINE rmulz(zmat,finn,fout,ifadd)  
    INTEGER :: ifadd
    REAL(db) :: zmat(:,:),finn(:,:,:),fout(:,:,:)
    INTENT(IN) :: zmat,finn,ifadd
    INTENT(INOUT) :: fout
    INTEGER :: ix,iy,iz
    IF(ifadd==0) fout=0.D0  
    IF(ifadd>=0) THEN  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)+SUM(zmat(iz,:)*finn(ix,iy,:))
       END FORALL
    ELSE  
       FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
          fout(ix,iy,iz)=fout(ix,iy,iz)-SUM(zmat(iz,:)*finn(ix,iy,:))
       END FORALL
    ENDIF
    RETURN  
  END SUBROUTINE rmulz
  !***********************************************************************
END MODULE Trivial
