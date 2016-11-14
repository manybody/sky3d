MODULE User
  USE Params
  USE Grids
  USE Levels
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
                    drhos(nx,ny,nz,2,2,3),currs(nx,ny,nz,2,2,3),&
                    localization(nx,ny,nz,2,2),time
    COMPLEX(db)  :: ps1(nx,ny,nz,2) 
    CHARACTER(10) :: filename
    !
    !***********************************************************************
    !                                                                      *
    !       compose pre-localization                                       *
    !       use (nabla rho)/2 already accumulated in 'subr. densit'        *
    !       use first component of 'nablarho' as workspace                 *
    !                                                                      *
    !***********************************************************************
    !
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
  localization=1.0d0/(1.0d0+((rhos*taus-drhos_sq-currs_sq)/&
               (3.0d0/5.0d0*(6.0d0*pi**2.0d0)**(2.0d0/3.0d0)*rhos**(8.0d0/3.0d0)))**2)
    !
    !  preliminary print along axes
    ! 
    OPEN(33,file='localization.res',status='unknown')
    WRITE(33,'(/a)')  '# localization'
    WRITE(33,'(/a)')  '#   x    neutrons_up neutrons_down protons_up  protons_down  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (x(ixyz),((localization(ixyz,NY/2,NZ/2,is,iq),is=1,2),iq=1,2),ixyz=1,NX)
    WRITE(33,*)
    WRITE(33,'(/a)')  '#   y    neutrons_up neutrons_down protons_up  protons_down  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (y(ixyz),((localization(NX/2,ixyz,NZ/2,is,iq),is=1,2),iq=1,2),ixyz=1,NY)
    WRITE(33,*)
    WRITE(33,'(/a)')  '#   z    neutrons_up neutrons_down protons_up  protons_down  '
    WRITE(33,'(1x,f6.2,4g13.5)') &
         (z(ixyz),((localization(NX/2,NY/2,ixyz,is,iq),is=1,2),iq=1,2),ixyz=1,NZ)
    CLOSE(33)
    
    WRITE(filename,'(I6.6,A4)') iter,'.loc'
    time=0.D0
    OPEN(UNIT=scratch,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE')
    WRITE(scratch) iter,time,nx,ny,nz
    WRITE(scratch) dx,dy,dz,wxyz,x,y,z
    CALL write_density('Rho_n_up',rhos(:,:,:,1,1))
    CALL write_density('Rho_n_down',rhos(:,:,:,2,1))
    CALL write_density('Rho_p_up',rhos(:,:,:,1,2))
    CALL write_density('Rho_p_down',rhos(:,:,:,2,2))
    CALL write_density('Loc_n_up',localization(:,:,:,1,1))
    CALL write_density('Loc_n_down',localization(:,:,:,2,1))
    CALL write_density('Loc_p_up',localization(:,:,:,1,2))
    CALL write_density('Loc_p_down',localization(:,:,:,2,2))
    CLOSE(UNIT=scratch)

    RETURN
    !
  END SUBROUTINE localize
!
  SUBROUTINE write_density(name,values)
    CHARACTER(*),INTENT(IN) :: name
    REAL(db),INTENT(IN) :: values(nx,ny,nz)
    CHARACTER(10) :: stored_name
    stored_name=name
    WRITE(scratch) stored_name,.FALSE.,.FALSE.
    WRITE(scratch) values
  END SUBROUTINE write_density
  !
END MODULE User
