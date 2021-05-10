! This is a collection of subroutines for analysis during
! TDHF computation. Each subroutein can be activated by
! moving it to 'user.f90' and inserting a CALL at the
! appropriate place.



!***************************************************************************
! This routine computes the expectation value of the spin-dipole 
! and writes the result to the file 'spindipole.res'.
! The routine can be placed immediately in 'util.f90'.
!
  SUBROUTINE spindipole()
    USE Grids, ONLY: wxyz,x,y,z
    USE Params, ONLY: time

    REAL(db) :: spindip(3,3)
    INTEGER :: ix,iy,iz
    LOGICAL,SAVE :: tinit=.TRUE.

    IF(tinit) THEN
      OPEN(unit=scratch,file='spindipole.res')
      WRITE(scratch,'(a)') &
       '# time    r_1*sigma_1 r_2*sigma_1  r_3*sigma_1  r_1*sigma_2 ....'
      CLOSE(unit=scratch)
      tinit=.FALSE.
    END IF

    spindip = 0D0
    DO iz=1,nz  
      DO iy=1,ny  
        DO ix=1,nx
          spindip(1,:) = spindip(1,:) + x(ix)*sdens(ix,iy,iz,:,1)  
          spindip(2,:) = spindip(2,:) + y(iy)*sdens(ix,iy,iz,:,1)  
          spindip(3,:) = spindip(3,:) + z(iz)*sdens(ix,iy,iz,:,1)  
        END DO
      END DO
    END DO
    spindip = wxyz*spindip

    OPEN(unit=scratch,file='spindipole.res', POSITION='APPEND')  
    WRITE(scratch,'(f10.2,9(1pg13.5))') time,spindip
    CLOSE(unit=scratch)

  END SUBROUTINE spindipole


! This routine tests the ortho-normality of the wavefunctions
! in the array 'psi'. 
! The subroutine was part of a test version of 'levels.f90'
! and uses variables available in 'levels.f90'. Some new 'USE' 
! statements may be necessary before usign it in 'user.f90'.
!
  SUBROUTINE test_orthonorm(deviat_norm)
    USE Trivial, ONLY: rpsnorm, overlap
    REAL(db),INTENT(OUT) :: deviat_norm(2)
    COMPLEX(db) :: oij
    REAL(db) :: acc
    INTEGER :: nst1,nst2,j,iq
    DO iq=1,2
       acc = 0D0
       DO nst1=npmin(iq),npsi(iq)
          DO nst2=nst1,npsi(iq)
             oij=overlap(psi(:,:,:,:,nst1),psi(:,:,:,:,nst2))
             IF(nst1==nst2) THEN
                acc = acc + (REAL(oij)-1D0)**2
             ELSE
                acc = acc + 2D0*(REAL(oij)**2+AIMAG(oij)**2)
             END IF
          END DO
       END DO
       deviat_norm(iq) = SQRT(acc/(npsi(iq)-npmin(iq)+1)**2)
    END DO
  END SUBROUTINE test_orthonorm


! This routine computes the wavefunctions in momentum space
! and prints their absolute values at maxmimum momentum. This
! serves to check whether dynamics hits the upper bounds of
! momentum space. Results are written on unit 923 (unnamed).
! The subroutine was part of a test version of 'levels.f90'
! and uses variables available in 'levels.f90'. Some new 'USE' 
! statements may be necessary before usign it in 'user.f90'.
!
  SUBROUTINE pmomsmax(psin,timact)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    REAL(db), INTENT(IN) :: timact
    COMPLEX(db), ALLOCATABLE:: ps1(:,:,:,:),ps2(:,:,:,:)
#ifdef CUDA
    COMPLEX(db), DEVICE :: psin_d(:,:,:,:),ps1_d(:,:,:,:),ps2_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iy
    ALLOCATE(ps1(nx,ny,nz,2),ps2(nx,ny,nz,2))
    kfac=(PI+PI)/(dy*ny)
#ifdef CUDA
    ! Copy to device, perform FFT, copy back to host
    psin_d = psin
    ps1_d = ps1
    ps2_d = ps2
    CALL cufftExecZ2Z(xforward,psin_d,ps1_d,CUFFT_FORWARD)
    CALL cufftExecZ2Z(yforward,ps1_d,ps2_d,CUFFT_FORWARD)
    CALL cufftExecZ2Z(zforward,ps2_d,ps1_d,CUFFT_FORWARD)
    psin = psin_d
    ps1 = ps1_d
    ps2 = ps2_d
#else
    CALL dfftw_execute_dft(xforward,psin,ps1)
    CALL dfftw_execute_dft(yforward,ps1,ps2)
    CALL dfftw_execute_dft(zforward,ps2,ps1)
#endif
    WRITE(921,'(f8.3,4(1pg13.5))') timact,SUM(ABS(ps1(nx/2,ny/2,nz/2,:))),&
         SUM(ABS(ps1(nx/2,1,1,:))),SUM(ABS(ps1(1,ny/2,1,:))),&
         SUM(ABS(ps1(1,1,nz/2,:)))
    CALL flush(921)
    DEALLOCATE(ps1,ps2)
  END SUBROUTINE pmomsmax


! This routine computes the wavefunctions in momentum space and 
! prints their absolute values along the bounds in momentum space
! in k_x, k_y, k_z direction. It serves to check on a broader 
! basis whether and where dynamics hits the upper bounds of 
! momentum space. Results are written on unit 922 (unnamed).
! The subroutine was part of a test version of 'levels.f90'
! and uses variables available in 'levels.f90'. Some new 'USE' 
! statements may be necessary before usign it in 'user.f90'.
!
  SUBROUTINE pmomscut(psin)  
    COMPLEX(db), INTENT(IN) :: psin(:,:,:,:)
    COMPLEX(db), ALLOCATABLE:: ps1(:,:,:,:),ps2(:,:,:,:)
#ifdef CUDA
    COMPLEX(db), DEVICE :: psin_d(:,:,:,:),ps1_d(:,:,:,:),ps2_d(:,:,:,:)
#endif
    REAL(db) :: kfac
    INTEGER :: iy,ix,iz
    LOGICAL, SAVE :: tfirst=.true.
    !    IF(.NOT.tfirst) RETURN
    ALLOCATE(ps1(nx,ny,nz,2),ps2(nx,ny,nz,2))
    kfac=(PI+PI)/(dy*ny)
#ifdef CUDA
    ! Copy to device, perform FFT, copy back to host
    psin_d = psin
    ps1_d = ps1
    ps2_d = ps2
    CALL cufftExecZ2Z(xforward,psin,ps1,CUFFT_FORWARD)
    CALL cufftExecZ2Z(yforward,ps1,ps2,CUFFT_FORWARD)
    CALL cufftExecZ2Z(zforward,ps2,ps1,CUFFT_FORWARD)
    psin = psin_d
    ps1 = ps1_d
    ps2 = ps2_d
#else
    CALL dfftw_execute_dft(xforward,psin,ps1)
    CALL dfftw_execute_dft(yforward,ps1,ps2)
    CALL dfftw_execute_dft(zforward,ps2,ps1)
#endif
    ! along x
    WRITE(922,'(a)') '# along x' 
    DO ix=nx/2+1,nx
       WRITE(922,'(f8.3,12(1pg13.5))') -(nx-ix+1)*kfac,&
            (ABS(ps1(ix,ny/2,nz/2,:))),&
            (ABS(ps1(ix,ny/2,1,:))),&
            (ABS(ps1(ix,ny/2+1,1,:))),&
            (ABS(ps1(ix,1,nz/2,:))),&
            (ABS(ps1(ix,1,nz/2+1,:))),&
            (ABS(ps1(ix,1,1,:)))
    ENDDO
    DO ix=1,nx/2
       WRITE(922,'(f8.3,12(1pg13.5))') (ix-1)*kfac,&
            (ABS(ps1(ix,ny/2,nz/2,:))),&
            (ABS(ps1(ix,ny/2,1,:))),&
            (ABS(ps1(ix,ny/2+1,1,:))),&
            (ABS(ps1(ix,1,nz/2,:))),&
            (ABS(ps1(ix,1,nz/2+1,:))),&
            (ABS(ps1(ix,1,1,:)))
    ENDDO
    WRITE(922,'(1x/1x)')
    ! along y
    WRITE(922,'(a)') '# along y' 
    DO iy=ny/2+1,ny
       WRITE(922,'(f8.3,12(1pg13.5))') -(ny-iy+1)*kfac,&
            (ABS(ps1(nx/2,iy,nz/2,:))),&
            (ABS(ps1(nx/2,iy,1,:))),&
            (ABS(ps1(nx/2,iy,nz/2+1,:))),&
            (ABS(ps1(1,iy,nz/2,:))),&
            (ABS(ps1(nx/2+1,iy,nz/2,:))),&
            (ABS(ps1(1,iy,1,:)))
    ENDDO
    DO iy=1,ny/2
       WRITE(922,'(f8.3,12(1pg13.5))') (iy-1)*kfac,&
            (ABS(ps1(nx/2,iy,nz/2,:))),&
            (ABS(ps1(nx/2,iy,1,:))),&
            (ABS(ps1(nx/2,iy,nz/2+1,:))),&
            (ABS(ps1(1,iy,nz/2,:))),&
            (ABS(ps1(nx/2+1,iy,nz/2,:))),&
            (ABS(ps1(1,iy,1,:)))
    ENDDO
    WRITE(922,'(1x/1x)')
    ! along z
    WRITE(922,'(a)') '# along z' 
    DO iz=nz/2+1,nz
       WRITE(922,'(f8.3,12(1pg13.5))') -(nz-iz+1)*kfac,&
            (ABS(ps1(nx/2,ny/2,iz,:))),&
            (ABS(ps1(1,ny/2,iz,:))),&
            (ABS(ps1(nx/2+1,ny/2,iz,:))),&
            (ABS(ps1(nx/2,1,iz,:))),&
            (ABS(ps1(nx/2,ny/2+1,iz,:))),&
            (ABS(ps1(1,1,iz,:)))
    ENDDO
    DO iz=1,nz/2
       WRITE(922,'(f8.3,12(1pg13.5))') (iz-1)*kfac,&
            (ABS(ps1(nx/2,ny/2,iz,:))),&
            (ABS(ps1(1,ny/2,iz,:))),&
            (ABS(ps1(nx/2+1,ny/2,iz,:))),&
            (ABS(ps1(nx/2,1,iz,:))),&
            (ABS(ps1(nx/2,ny/2+1,iz,:))),&
            (ABS(ps1(1,1,iz,:)))
    ENDDO
    !
    WRITE(922,'(1x/1x)')
    CALL flush(922)
    DEALLOCATE(ps1,ps2)
    !    tfirst=.false.
  END SUBROUTINE pmomscut
