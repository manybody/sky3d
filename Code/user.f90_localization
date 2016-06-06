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
    USE Densities, ONLY: rho, tau, current
  USE Trivial, ONLY: cmulx,cmuly,cmulz
!    USE Trivial, ONLY: rmulx,rmuly,rmulz
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
    INTEGER :: ixyz,iq,nst,n
    REAL(db) :: tf_fac                  ! Thomas-Fermi factor
    REAL(db),ALLOCATABLE :: nablarho(:,:,:,:,:)
    COMPLEX(db),ALLOCATABLE :: ps1(:,:,:,:)  
    !
    !***********************************************************************
    !                                                                      *
    !       compose pre-localization                                       *
    !       use (nabla rho)/2 already accumulated in 'subr. densit'        *
    !       use first component of 'nablarho' as workspace                 *
    !                                                                      *
    !***********************************************************************
    !
    ALLOCATE(nablarho(nx,ny,nz,3,2))
    ALLOCATE(ps1(nx,ny,nz,2)  )
    nablarho=0D0
    DO nst=1,nstmax
      iq=isospin(nst)
      IF(TFFT) THEN
        CALL cdervx(psi(:,:,:,:,nst),ps1)  
      ELSE
        CALL cmulx(der1x,psi(:,:,:,:,nst),ps1,0)  
      ENDIF
      nablarho(:,:,:,1,iq)=nablarho(:,:,:,1,iq)+wocc(nst)* &
         REAl(ps1(:,:,:,1)*CONJG(psi(:,:,:,1,nst))+ &
              ps1(:,:,:,2)*CONJG(psi(:,:,:,2,nst))   )
      IF(TFFT) THEN
        CALL cdervy(psi(:,:,:,:,nst),ps1)  
      ELSE
        CALL cmuly(der1y,psi(:,:,:,:,nst),ps1,0)  
      ENDIF
      nablarho(:,:,:,2,iq)=nablarho(:,:,:,2,iq)+wocc(nst)* &
         REAl(ps1(:,:,:,1)*CONJG(psi(:,:,:,1,nst))+ &
              ps1(:,:,:,2)*CONJG(psi(:,:,:,2,nst))   )
      IF(TFFT) THEN
        CALL cdervz(psi(:,:,:,:,nst),ps1)  
      ELSE
        CALL cmulz(der1z,psi(:,:,:,:,nst),ps1,0)  
      ENDIF
      nablarho(:,:,:,3,iq)=nablarho(:,:,:,3,iq)+wocc(nst)* &
         REAl(ps1(:,:,:,1)*CONJG(psi(:,:,:,1,nst))+ &
              ps1(:,:,:,2)*CONJG(psi(:,:,:,2,nst))   )
    ENDDO
    DEALLOCATE(ps1)
!    DO iq=1,2
!       CALL rmulx(der1x,rho(:,:,:,iq),nablarho(:,:,:,1,iq),0)
!       CALL rmuly(der1y,rho(:,:,:,iq),nablarho(:,:,:,2,iq),0)
!       CALL rmulz(der1z,rho(:,:,:,iq),nablarho(:,:,:,3,iq),0)
!    ENDDO
!    nablarho=nablarho/2D0
!    WRITE(33,'(/a)')  '#   x  rho_n  nablarho_n  '
!    WRITE(33,'(1x,f6.2,4g13.5)') &
!         (x(ixyz),rho(ixyz,NY/2,NZ/2,2),(nablarho(ixyz,NY/2,NZ/2,n,2),n=1,3),ixyz=1,NX)

    DO iq = 1,2  
       nablarho(:,:,:,1,iq) = tau(:,:,:,iq) - &
            ( (nablarho(:,:,:,1,iq)**2+nablarho(:,:,:,2,iq)**2  &
            +nablarho(:,:,:,3,iq)**2)  + &
            (current(:,:,:,1,iq)**2+current(:,:,:,2,iq)**2  &
            +current(:,:,:,3,iq)**2)      &
            )/rho(:,:,:,iq)
    ENDDO
    !
    !***********************************************************************
    !                                                                      *
    !       rescale to final localization criterion                        *
    !                                                                      *
    !***********************************************************************
    !
    tf_fac = 0.6D0*(3.0D0*pi**2)**0.666666666667D0
    nablarho(:,:,:,1,:) =    1.0D0/(1.0D0   &
         +nablarho(:,:,:,1,:)/(tf_fac*rho(:,:,:,:)**1.666666666667D0))
    !
    !  preliminary print along axes
    ! 
    OPEN(33,file='localization.res',status='unknown')
    WRITE(33,'(/a)')  '# localization'
    WRITE(33,'(/a)')  '#   x    neutrons  protons  '
    WRITE(33,'(1x,f6.2,2g13.5)') &
         (x(ixyz),(nablarho(ixyz,NY/2,NZ/2,1,iq),iq=1,2),ixyz=1,NX)
    WRITE(33,'(/a)')  '#   y    neutrons  protons  '
    WRITE(33,'(1x,f6.2,2g13.5)') &
         (y(ixyz),(nablarho(NX/2,ixyz,NZ/2,1,iq),iq=1,2),ixyz=1,NY)
    WRITE(33,'(/a)')  '#   z    neutrons  protons  '
    WRITE(33,'(1x,f6.2,2g13.5)') &
         (z(ixyz),(nablarho(NX/2,NY/2,ixyz,1,iq),iq=1,2),ixyz=1,NZ)
    CLOSE(33)
    !
    DEALLOCATE(nablarho)
    RETURN
    !
  END SUBROUTINE localize

  !
END MODULE User
