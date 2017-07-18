MODULE Moment
  USE Params
  USE Grids, ONLY: nx,ny,nz,x,y,z,wxyz
  USE Levels, ONLY: nprot,nneut
  IMPLICIT NONE
  PRIVATE
  REAL(db) :: pnr(2),pnrtot,cm(3,2),cmtot(3),pcm(3,2)
  REAL(db) :: rms(2),rmstot,q20(2),q20tot,q22(2),q22tot,q20T1, &
       x2m(3,2),x2mtot(3),beta20tot,beta22tot,beta,gamma
  PUBLIC :: pnr,pnrtot,cm,cmtot,pcm,rmstot,beta,gamma, &
       moments,moment_print,moment_shortprint
CONTAINS
  !***********************************************************
  SUBROUTINE moments
    USE Densities, ONLY: rho,current
    INTEGER :: ix,iy,iz,iq
    REAL(db) :: xx(3),x2(3),vol,radius,r0rms
    REAL(db) :: qmat(3,3,2),qmtot(3,3)
    pnr=0.D0
    cm=0.D0
    pcm=0.D0
    DO iq=1,2  
       DO iz=1,nz  
          xx(3)=z(iz)  
          DO iy=1,ny  
             xx(2)=y(iy)  
             DO ix=1,nx  
                xx(1)=x(ix)
                pnr(iq)=pnr(iq)+wxyz*rho(ix,iy,iz,iq)
                cm(:,iq)=cm(:,iq)+wxyz*xx*rho(ix,iy,iz,iq)
                IF(tdynamic) pcm(:,iq)=pcm(:,iq)+wxyz*current(ix,iy,iz,:,iq)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    pnrtot=pnr(1)+pnr(2)  
    cmtot=(cm(:,1)+cm(:,2))/pnrtot  
    DO iq=1,2
       cm(:,iq)=cm(:,iq)/pnr(iq)
    ENDDO
    !***********************************
    rms=0.D0
    qmat=0.D0
    x2m=0.D0
    DO iq=1,2  
       DO iz=1,nz  
          xx(3)=z(iz)-cm(3,iq)  
          x2(3)=xx(3)**2  
          DO iy=1,ny  
             xx(2)=y(iy)-cm(2,iq)  
             x2(2)=xx(2)**2
             DO ix=1,nx  
                xx(1)=x(ix)-cm(1,iq)  
                x2(1)=xx(1)**2  
                vol=wxyz*rho(ix,iy,iz,iq)  
                rms(iq)=vol*SUM(x2)+rms(iq)
                qmat(1,1,iq)=qmat(1,1,iq)+vol*(x2(1)+x2(1)-x2(2)-x2(3))
                qmat(1,2,iq)=qmat(1,2,iq)+3.D0*vol*xx(1)*xx(2)
                qmat(1,3,iq)=qmat(1,3,iq)+3.D0*vol*xx(1)*xx(3)
                qmat(2,2,iq)=qmat(2,2,iq)+vol*(x2(2)+x2(2)-x2(1)-x2(3))
                qmat(2,3,iq)=qmat(2,3,iq)+3.D0*vol*xx(2)*xx(3)
                qmat(3,3,iq)=qmat(3,3,iq)+vol*(x2(3)+x2(3)-x2(1)-x2(2))
                x2m(:,iq)=vol*x2(:)+x2m(:,iq)  
             ENDDO
          ENDDO
       ENDDO
       qmat(2,1,iq)=qmat(1,2,iq)
       qmat(3,1,iq)=qmat(1,3,iq)
       qmat(3,2,iq)=qmat(2,3,iq)
    ENDDO
    rmstot=SQRT((rms(1)+rms(2))/pnrtot)
    rms=SQRT(rms/pnr)  
    x2mtot=(x2m(:,1)+x2m(:,2))/pnrtot  
    DO iq=1,2
       x2m(:,iq)=x2m(:,iq)/pnr(iq)
    ENDDO
    qmtot=qmat(:,:,1)+qmat(:,:,2) 
    IF(printnow.AND.wflag) WRITE(*,'(/A)') 'Cartesian quadrupole tensor,&
         &  principal values, and axes:'
    CALL q2diag(qmat(:,:,1),q20(1),q22(1),'Neutrons ')
    CALL q2diag(qmat(:,:,2),q20(2),q22(2),'Protons  ')
    CALL q2diag(qmtot,q20tot,q22tot,'Total    ')
    q20T1 = -q20(1)/nneut+q20(2)/nprot
    r0rms = r0*SQRT(0.6)           ! convert from box to r.m.s. radius
    radius=rmstot!r0rms*pnrtot**(1.D0/3.D0)
    beta20tot=q20tot*(4.0D0*PI/(5.0D0*radius**2*pnrtot))
    beta22tot=q22tot*(4.0D0*PI/(5.0D0*radius**2*pnrtot))
    beta=SQRT(beta20tot**2+2.0*beta22tot**2)
    gamma=ABS(ATAN2(SQRT(2.0)*beta22tot,beta20tot)*180.0D0/PI)
    IF(gamma>120.D0) THEN
       gamma=gamma-120.D0
    ELSEIF(gamma>60.D0) THEN
       gamma=120.D0-gamma
    ENDIF
    IF(printnow.AND.wflag) WRITE(*,'(4(A,F8.4)/)') &
         ' Beta20: ',beta20tot,' Beta22: ',beta22tot,' Beta: ',beta, &
         ' Gamma: ',gamma
  END SUBROUTINE moments
  !***********************************************************
  SUBROUTINE moment_shortprint
    REAL(db),SAVE :: q2old(4,3),asig(4)=(/1D0,1D0,1D0,1D0/)
    INTEGER, SAVE :: countcall=0
    INTEGER :: i
    countcall = 1+countcall
    q2old(:,1) = q2old(:,2)
    q2old(:,2) = q2old(:,3)
    q2old(1:2,3) = q20
    q2old(3,3) = q20tot
    q2old(4,3) = q20T1
    IF(countcall>2) THEN
       DO i=1,4
          IF(q2old(i,1)>q2old(i,2) .AND. q2old(i,3)>q2old(i,2)) &
               asig(i) = -asig(i)
       END DO
    END IF
    OPEN(unit=scratch,file=monopolesfile,POSITION='APPEND')  
    WRITE(scratch,'(4F10.2,E14.5)') time,rms,rmstot,rms(1)-rms(2)
    CLOSE(unit=scratch)
    OPEN(unit=scratch,file=quadrupolesfile,POSITION='APPEND')  
    WRITE(scratch,'(4F10.2,E14.5)') time,q20,q20tot,x2m
    CLOSE(unit=scratch)
  END SUBROUTINE moment_shortprint
  !***********************************************************
  SUBROUTINE moment_print
    INTEGER :: iq
    CHARACTER(11),PARAMETER :: Name(2)=(/ '  Neutron: ','   Proton: '/)
    Write(*,'(A)') '              Part.Num.   rms-radius   q20         &
         &<x**2>      <y**2>      <z**2>        <x>            &
         &<y>            <z>'    
    WRITE(*,'(a,2f12.4,1p,4e12.4,3e15.7)') '    Total: ',pnrtot,rmstot, &
         q20tot,x2mtot,cmtot
    DO iq=1,2
       WRITE(*,'(a,2f12.4,1p,4e12.4,3e15.7)') name(iq),pnr(iq),rms(iq),q20(iq), &
            x2m(:,iq),cm(:,iq)
    ENDDO
  END SUBROUTINE moment_print
  !***********************************************************
  SUBROUTINE q2diag(q_mat,q20x,q22x,title)
    REAL(db),INTENT(INOUT) :: q_mat(3,3)
    REAL(db),INTENT(OUT) :: q20x,q22x
    CHARACTER(LEN=*),INTENT(IN) :: title
    REAL(db) :: q_eig(3),q_vec(3,3),fv1(20)
    INTEGER :: info,i, j,k
    IF(printnow.AND.wflag) WRITE(*,'(3(f12.5,1x))') ((q_mat(j,k),k=1,3),j=1,3)    
    CALL DSYEV('V','U',3,q_mat,3,q_eig,fv1,20,info)
    q_vec=q_mat
    IF(info/=0) THEN
       q_mat=0.0d0
       q_eig=0.0d0
    END IF
    IF(printnow.AND.wflag) THEN
       WRITE(*,'(1X,A,3(F10.2,''('',3F8.4,'')''))') &
            title,(REAL(q_eig(i)),REAL(q_vec(:,i)),i=3,1,-1)
    ENDIF
    q20x=SQRT(5.D0/(16.D0*pi))*REAL(q_eig(3))
    q22x=SQRT(5.D0/(96.D0*pi))*(REAL(q_eig(2))-REAL(q_eig(1)))
  END SUBROUTINE q2diag
END MODULE Moment
