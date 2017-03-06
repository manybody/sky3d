MODULE User
  USE Params
  USE Grids
  USE Levels
  USE Densities, ONLY: localization
  IMPLICIT NONE
CONTAINS
  SUBROUTINE init_user
  IMPLICIT NONE
    CALL localize()    
  END SUBROUTINE init_user
  SUBROUTINE localize()
    !
    IMPLICIT NONE  
    !
    !***********************************************************************
    !                                                                      *
    !       localize:                                                      *
    !        computes the localization criterion from                      *
    !             Becke and Egdecomb, JCP {\bf 92} (1990) 5397.            *
    !        The (nabla rho)/2 was accumulated on 'nablarho' in the        *
    !        'subroutine densit'. It has to be computed exactly in the     *
    !        same manner as the current to guarantee compensation.         *
    !        The localization is composed first as                         *
    !          C = tau*rho - 1/4*(nabla rho)**2 -curr**2  ,                *
    !        then the 'C' is rescaled into the final criterion             *
    !            1/(1+C/tau_TF)                                            *
    !        where 'tau_TF' is the kinetic energy density in               *
    !        Thomas-Fermi approximation.                                   *
    !                                                                      *
    !***********************************************************************
    !
    INTEGER      :: ixyz,iq,i,is
    REAL(db)     :: rhos(nx,ny,nz,2,2),taus(nx,ny,nz,2,2),&
                    drhos_sq(nx,ny,nz,2,2),currs_sq(nx,ny,nz,2,2),&
                    drhos(nx,ny,nz,2,2,3),currs(nx,ny,nz,2,2,3),time
    COMPLEX(db)  :: ps1(nx,ny,nz,2) 
    CHARACTER(10) :: filename

  DO is=1,2
    DO iq=1,2
      rhos(:,:,:,is,iq)=0.0d0
      taus(:,:,:,is,iq)=0.0d0
      drhos_sq(:,:,:,is,iq)=0.0d0
      drhos(:,:,:,is,iq,:)=0.0d0
      currs(:,:,:,is,iq,:)=0.0d0
      DO i=npmin(iq),npsi(iq)
        rhos(:,:,:,is,iq)=rhos(:,:,:,is,iq)+wocc(i)*(psi(:,:,:,is,i)*CONJG(psi(:,:,:,is,i)))
        CALL cdervx(psi(:,:,:,:,i),ps1)
        taus(:,:,:,is,iq)=taus(:,:,:,is,iq)+wocc(i)*ps1(:,:,:,is)*CONJG(ps1(:,:,:,is))
        drhos(:,:,:,is,iq,1)=drhos(:,:,:,is,iq,1)+wocc(i)*REAL(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
        currs(:,:,:,is,iq,1)=currs(:,:,:,is,iq,1)+wocc(i)*AIMAG(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
        CALL cdervy(psi(:,:,:,:,i),ps1)
        taus(:,:,:,is,iq)=taus(:,:,:,is,iq)+wocc(i)*ps1(:,:,:,is)*CONJG(ps1(:,:,:,is))
        drhos(:,:,:,is,iq,2)=drhos(:,:,:,is,iq,2)+wocc(i)*REAL(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
        currs(:,:,:,is,iq,2)=currs(:,:,:,is,iq,2)+wocc(i)*AIMAG(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
        CALL cdervz(psi(:,:,:,:,i),ps1)
        taus(:,:,:,is,iq)=taus(:,:,:,is,iq)+wocc(i)*ps1(:,:,:,is)*CONJG(ps1(:,:,:,is))
        drhos(:,:,:,is,iq,3)=drhos(:,:,:,is,iq,3)+wocc(i)*REAL(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
        currs(:,:,:,is,iq,3)=currs(:,:,:,is,iq,3)+wocc(i)*AIMAG(CONJG(psi(:,:,:,is,i))*ps1(:,:,:,is))
      END DO
      drhos_sq(:,:,:,is,iq)=drhos(:,:,:,is,iq,1)**2+drhos(:,:,:,is,iq,2)**2+drhos(:,:,:,is,iq,3)**2
      currs_sq(:,:,:,is,iq)=currs(:,:,:,is,iq,1)**2+currs(:,:,:,is,iq,2)**2+currs(:,:,:,is,iq,3)**2
    END DO
  END DO
  rhos(:,:,:,1,:)=rhos(:,:,:,1,:)+rhos(:,:,:,2,:)
  taus(:,:,:,1,:)=taus(:,:,:,1,:)+taus(:,:,:,2,:)
  drhos_sq(:,:,:,1,:)=drhos_sq(:,:,:,1,:)+drhos_sq(:,:,:,2,:)
  currs_sq(:,:,:,1,:)=currs_sq(:,:,:,1,:)+currs_sq(:,:,:,2,:)
  localization=1.0d0/(1.0d0+((rhos(:,:,:,1,:)*taus(:,:,:,1,:)-drhos_sq(:,:,:,1,:)-currs_sq(:,:,:,1,:))/&
               (3.0d0/5.0d0*(6.0d0*pi**2.0d0)**(2.0d0/3.0d0)*rhos(:,:,:,1,:)*2.0d0*(rhos(:,:,:,1,:)/2.0d0)**(5.0d0/3.0d0)))**2)
    !
    !  preliminary print along axes
    ! 
    OPEN(33,file='localization.res',status='unknown')
    WRITE(33,'(/a)')  '# localization'
    WRITE(33,'(/a)')  '#   x    neutrons      protons  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (x(ixyz),(localization(ixyz,NY/2,NZ/2,iq),iq=1,2),ixyz=1,NX)
    WRITE(33,*)
    WRITE(33,'(/a)')  '#   y    neutrons      protons  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (y(ixyz),(localization(NX/2,ixyz,NZ/2,iq),iq=1,2),ixyz=1,NY)
    WRITE(33,*)
    WRITE(33,'(/a)')  '#   z    neutrons      protons  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (z(ixyz),(localization(NX/2,NY/2,ixyz,iq),iq=1,2),ixyz=1,NZ)
    CLOSE(33)
    RETURN

  END SUBROUTINE localize
END MODULE User
