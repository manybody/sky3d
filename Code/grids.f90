MODULE Grids
  USE Params, ONLY: db,pi,wflag
  IMPLICIT NONE
  SAVE
  INTEGER :: nx,ny,nz
  LOGICAL :: periodic
  REAL(db) :: dx,dy,dz
  REAL(db) :: wxyz
  REAL(db),POINTER ::  x(:),y(:),z(:)
  REAL(db),POINTER,DIMENSION(:,:) ::  der1x,der2x,cdmpx, &
       der1y,der2y,cdmpy,der1z,der2z,cdmpz
  PRIVATE :: init_coord, sder, sder2, setdmc, gauss
CONTAINS
  !***************************************************
  ! Initialization of all the grid quantities,
  !***************************************************
  SUBROUTINE init_grid
    NAMELIST /Grid/ nx,ny,nz,dx,dy,dz,periodic
    dx=0.D0
    dy=0.D0
    dz=0.D0
    READ(5,Grid)
    IF(MOD(nx,2)/=0.OR.MOD(ny,2)/=0.OR.MOD(nz,2)/=0) THEN
      IF(wflag) WRITE(*,'(A,3I4)') 'Dimensions must be even: ',nx,ny,nz
      STOP
    END IF
    IF(wflag) THEN
      WRITE(*,*) '***** Grid Definition *****'
      IF(periodic) THEN
        WRITE(*,*) 'Grid is periodic'
      ELSE
        WRITE(*,*) 'Grid is not periodic'
      END IF
    END IF
    IF(dx*dy*dz<=0.D0) THEN
      IF(dx<=0.D0) STOP 'Grid spacing given as zero'
      dy=dx
      dz=dx
    END IF
    CALL init_coord('x',nx,dx,x,der1x,der2x,cdmpx)
    CALL init_coord('y',ny,dy,y,der1y,der2y,cdmpy)
    CALL init_coord('z',nz,dz,z,der1z,der2z,cdmpz)
    wxyz=dx*dy*dz
  END SUBROUTINE init_grid
  !***************************************************
  ! Initialization for one coordinate direction
  !***************************************************
  SUBROUTINE init_coord(name,nv,dv,v,der1v,der2v,cdmpv)
    CHARACTER(*) :: name
    INTEGER :: nv
    REAL(db) :: dv
    INTENT(IN) :: name,nv,dv
    REAL(db),POINTER :: v(:),der1v(:,:), &
         der2v(:,:),cdmpv(:,:)
    INTEGER :: i
    ALLOCATE(v(nv),der1v(nv,nv),der2v(nv,nv),cdmpv(nv,nv))
    v=(/ ((i-1)*dv-0.5D0*FLOAT(nv-1)*dv,i=1,nv) /)
    IF(wflag) THEN
       WRITE(*,'(1X,A,I3,A,F8.4,2(A,F8.4))') name // ' direction: ',nv, &
            ' points, spacing:',dv,' ranging from ',v(1),' to ',v(nv)
    ENDIF
    CALL sder(der1v,nv,dv)
    CALL sder2(der2v,nv,dv)
  END SUBROUTINE init_coord
  !***************************************************
  ! Computation of first derivative matrix
  !***************************************************
  PURE SUBROUTINE sder(der,nmax,d)
    INTEGER :: nmax
    REAL(db) :: d,der(:,:)
    INTENT(IN) :: nmax,d
    INTENT(OUT) :: der
    REAL(db) :: afac,sum
    INTEGER :: i,k,j,icn
    icn=(nmax+1)/2
    afac=pi/icn
    DO k=1,nmax
       DO i=1,nmax
          sum=0.0D0
          DO j=1,icn-1
             sum=sum-j*SIN(j*afac*(k-i))
          ENDDO
          sum=sum-0.5D0*icn*SIN(icn*afac*(k-i))
          der(i,k)=-afac*sum/(icn*d)
       ENDDO
    ENDDO
  END SUBROUTINE sder
  !***************************************************
  ! Computation of second-derivative matrices
  !***************************************************
  PURE SUBROUTINE sder2(der,nmax,d)
    INTEGER :: nmax
    REAL(db) :: d,der(1:nmax,1:nmax)
    INTENT(IN) :: nmax,d
    INTENT(OUT) :: der
    REAL(db) :: afac,sum
    INTEGER :: i,k,j,icn
    icn=(nmax+1)/2
    afac=pi/icn
    DO k=1,nmax
       DO i=1,nmax
          sum=0.0D0
          DO j=1,icn-1
             sum=sum+j**2*COS(j*afac*(k-i))
          ENDDO
          sum=sum+0.5D0*icn**2*COS(icn*afac*(k-i))
          der(i,k)=-(afac*afac)*sum/(icn*d*d)
       ENDDO
    ENDDO
  END SUBROUTINE sder2
  !***************************************************
  ! Calculation of damping matrices
  !***************************************************
  SUBROUTINE setup_damping(e0dmp)
    REAL(db),INTENT(IN) :: e0dmp
    CALL setdmc(der2x,nx,cdmpx,e0dmp)
    CALL setdmc(der2y,ny,cdmpy,e0dmp)
    CALL setdmc(der2z,nz,cdmpz,e0dmp)
  END SUBROUTINE setup_damping
  !***************************************************
  ! Calculation of damping matrix, one direction
  !***************************************************
  PURE SUBROUTINE setdmc(der2,nmax,cdmp,e0dmp)
    USE Forces, ONLY: h2ma
    REAL(db),INTENT(IN) :: der2(:,:)
    REAL(db),INTENT(OUT) :: cdmp(:,:)
    INTEGER,INTENT(IN) :: nmax
    REAL(db),INTENT(IN) :: e0dmp
    REAL(db) :: unit(nmax,nmax)
    INTEGER :: i
    cdmp=0.D0
    IF(e0dmp<=0.D0) RETURN
    unit=-h2ma*der2/e0dmp  
    FORALL(i=1:nmax)
       unit(i,i)=1.0D0+unit(i,i)  
       cdmp(i,i)=1.0D0  
    END FORALL
    CALL gauss(unit,cdmp,nmax)
  END SUBROUTINE setdmc
  !***************************************************
  PURE SUBROUTINE gauss(a,b,n)
    INTEGER,INTENT(IN) :: n
    REAL(db),INTENT(IN) :: a(n,n)
    REAL(db),INTENT(INOUT) :: b(n,n)
    REAL(db) :: c(n,2*n), d(2*n)
    INTEGER :: k,ks,ip(1),i,j
    ! copy matrix a into left half of m
    c(:,1:n)=a
    ! fill right-hand side
    c(:,n+1:2*n)=b
    ! start triangle decomposition
    DO k=1,n
       ! look for pivot in rows k..n
       ip=MAXLOC(ABS(c(k:n,k)))
       ! ks=index of pivot row
       ks=k-1+ip(1)
       ! interchange with row k if k/=ks
       ! using d as intermediate storage
       IF(k/=ks) THEN
          d(k:2*n)=c(k,k:2*n)
          c(k,k:2*n)=c(ks,k:2*n) 
          c(ks,k:2*n)=d(k:2*n)
       END IF
       c(k,k+1:2*n)=c(k,k+1:2*n)/c(k,k)
       FORALL(i=k+1:n,j=k+1:2*n) c(i,j)=&
            c(i,j)-c(i,k)*c(k,j)
    END DO
    b=c(:,n+1:2*n)
    DO k=n-1,1,-1
       DO i=1,n
          b(1:k,i)=b(1:k,i)-c(:k,k+1)*b(k+1,i)
       END DO
    END DO
  END SUBROUTINE gauss
  !***************************************************
END MODULE Grids
