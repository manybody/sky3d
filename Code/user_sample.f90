MODULE User
  USE Params
  USE Grids
  USE Levels
  IMPLICIT NONE
CONTAINS
  SUBROUTINE init_user
    INTEGER :: ic,ix,iy,iz
    REAL(db) :: center(3),d,r
    NAMELIST /user/ d,r
    READ(5,user)
    center=(/-d,0.D0,d/)
    WRITE(*,'(A,F5.2,A,F5.4)') '3 alpha cluster initialization with distance ' &
         ,d,' and radii ',r
    wocc=1.D0
    DO ic=1,3
      FORALL(ix=1:nx,iy=1:ny,iz=1:nz)
        psi(ix,iy,iz,1,ic)= &
             EXP(-((x(ix)-center(ic))**2+y(iy)**2+z(iz)**2)/(2.D0*r**2))
      END FORALL
      psi(:,:,:,2,ic)=0.D0
      psi(:,:,:,1,ic+3)=0.D0
      psi(:,:,:,2,ic+3)=psi(:,:,:,1,ic)
      psi(:,:,:,:,ic+6)=psi(:,:,:,:,ic)
      psi(:,:,:,:,ic+9)=psi(:,:,:,:,ic+3)
    END DO
  END SUBROUTINE init_user
END MODULE User
