!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module deals with the definition of the spatial grid and associated operations.
!------------------------------------------------------------------------------
MODULE Grids
  USE Params, ONLY: db,pi,wflag,tabc_myid,tabc_nprocs,tfft
  IMPLICIT NONE
  SAVE
  INTEGER  :: nx            !< Number of points in x-direction. Must be even to preserve
  !!reflection symmetry.
  INTEGER  :: ny            !< Number of points in y-direction. Must be even to preserve
  !!reflection symmetry.
  INTEGER  :: nz            !< Number of points in z-direction. Must be even to preserve
  !!reflection symmetry.
  INTEGER  :: tabc_x=0      !< Number of TABC points in x-direction.
  INTEGER  :: tabc_y=0      !< Number of TABC points in y-direction.
  INTEGER  :: tabc_z=0      !< Number of TABC points in z-direction.
  REAL(db) :: bangx         !< bloch phase in x-direction.
  REAL(db) :: bangy         !< bloch phase in y-direction.
  REAL(db) :: bangz         !< bloch phase in z-direction.
  LOGICAL  :: periodic      !< logical variable indicating whether the situation is triply 
  !!periodic in three-dimensional space.
  REAL(db) :: dx            !< The spacing between grid points (in fm) in x-directions.
  REAL(db) :: dy            !< The spacing between grid points (in fm) in y-directions.
  REAL(db) :: dz            !< The spacing between grid points (in fm) in z-directions.
  REAL(db) :: wxyz          !< the volume element <tt> wxyz=dx*dy*dz</tt>.
  REAL(db),POINTER :: x(:)  !<array containing the actual x-coordinate values in fm, dimensioned
  !!as \c x(nx) and allocated dynamically, thus allowing dynamic dimensioning through input 
  !!values of \c nx.
  REAL(db),POINTER :: y(:)  !<array containing the actual y-coordinate values in fm, dimensioned
  !!as \c x(ny) and allocated dynamically, thus allowing dynamic dimensioning through input 
  !!values of \c ny.
  REAL(db),POINTER :: z(:)  !<array containing the actual z-coordinate values in fm, dimensioned
  !!as \c x(nz) and allocated dynamically, thus allowing dynamic dimensioning through input 
  !!values of \c nz.
  REAL(db),POINTER,DIMENSION(:,:) ::  &
              der1x, &      !<matrices describing the first spatial derivatives in the x-direction, 
  !!dynamically allocated with dimensions <tt> (nx,nx) </tt> and calculated in subroutine \c sder.
              der2x, &      !<matrices describing the second spatial derivatives in the x-direction, 
  !!dynamically allocated with dimensions <tt> (nx,nx) </tt> and calculated in subroutine \c sder2.
              cdmpx, &      !<matrices describing the damping operation in the x-direction, 
  !!dynamically allocated with dimensions <tt> (nx,nx) </tt> and calculated in subroutine \c setdmc.
              der1y, &      !<matrices describing the first spatial derivatives in the y-direction, 
  !!dynamically allocated with dimensions <tt> (ny,ny) </tt> and calculated in subroutine \c sder.
              der2y, &      !<matrices describing the second spatial derivatives in the y-direction, 
  !!dynamically allocated with dimensions <tt> (ny,ny) </tt> and calculated in subroutine \c sder2.
              cdmpy, &      !<matrices describing the damping operation in the y-direction, 
  !!dynamically allocated with dimensions <tt> (ny,ny) </tt> and calculated in subroutine \c setdmc.
              der1z, &      !<matrices describing the first spatial derivatives in the z-direction, 
  !!dynamically allocated with dimensions <tt> (nz,nz) </tt> and calculated in subroutine \c sder.
              der2z, &      !<matrices describing the second spatial derivatives in the z-direction, 
  !!dynamically allocated with dimensions <tt> (nz,nz) </tt> and calculated in subroutine \c sder2.
              cdmpz         !<matrices describing the damping operation in the z-direction, 
  !!dynamically allocated with dimensions <tt> (nz,nz) </tt> and calculated in subroutine \c setdmc.
  PRIVATE :: init_coord, sder, sder2, setdmc, gauss
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: init_grid
!> @brief
!!This subroutine is called during the initialization for both static
!!and dynamic calculations. It reads the grid dimension and spacing
!!information using namelist \c Grid, allocates the necessary arrays,
!!and calculates the coordinate values and the derivative and damping
!!matrices.
!>
!> @details
!!This is done by calling the subroutine \c init_coord once for each
!!direction. It does everything needed except the calculation of
!!the volume element. 
!--------------------------------------------------------------------------- 
  SUBROUTINE init_grid
    NAMELIST /Grid/ nx,ny,nz,dx,dy,dz,periodic,bangx,bangy,bangz,&
                    tabc_x,tabc_y,tabc_z
    dx=0.D0
    dy=0.D0
    dz=0.D0
    bangx=0.0d0
    bangy=0.0d0
    bangz=0.0d0
    READ(5,Grid)
    IF(wflag) THEN
       WRITE(*,*)
       WRITE(*,*) '***** Grid Definition *****'
       IF(periodic) THEN
          WRITE(*,*) 'Grid is periodic'
       ELSE
          WRITE(*,*) 'Grid is not periodic'
       END IF
    END IF
    IF(.NOT.TFFT.AND.(abs(bangx)>0.00001.OR.abs(bangy)>0.00001.OR.abs(bangz)>0.00001)) &
      STOP 'Bloch boundaries cannot be used without TFFT'
    IF(tabc_nprocs>1) CALL tabc_init_blochboundary
    IF(tabc_nprocs==1.AND.(tabc_x/=0.OR.tabc_y/=0.OR.tabc_z/=0)) &
      STOP 'No TABC possible with tabc_nprocs=1!!!' 
    bangx=bangx*PI
    bangy=bangy*PI
    bangz=bangz*PI
    IF(MOD(nx,2)/=0.OR.MOD(ny,2)/=0.OR.MOD(nz,2)/=0) THEN
       IF(wflag) WRITE(*,'(A,3I4)') 'Dimensions must be even: ',nx,ny,nz
       STOP
    END IF
    IF(wflag) THEN
      WRITE(*,*) 'Bloch twist x-direction: ',bangx
       WRITE(*,*) 'Bloch twist y-direction: ',bangy
       WRITE(*,*) 'Bloch twist z-direction: ',bangz
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
!---------------------------------------------------------------------------  
! DESCRIPTION: Routinename
!> @brief
!!In this subroutine the defining information for a grid direction
!!(generically called \c v, which can be replaced by \c x, \c y, or \c z) 
!!in the form of the number of points \c nv and the
!!spacing \c dv is used to generate the associated data.  The arrays
!!of coordinate values, derivative and damping matrices are allocated.
!!Since all the quantities that are later used in the code are passed as
!!arguments, this subroutine can handle all three directions in a
!!unified way. To print the information intelligibly, it is also passed
!!the name of the coordinate as \c name.
!>
!> @details
!!It is assumed that the coordinate zero is in the center of the grid,
!!i.e., since the dimension is even the number of points to each side
!!of zero is equal andthe origin is in the center of a cell. The special
!!position of the origin is used in static calculations, e.g., for the
!!parity determination. In other situations, the position of the center
!!of mass is more important, this is defined in module \c Moment.
!!
!!If a difference location of the origin in the grid is desired, it can
!!be done by changing the statement generating the values of \c v.
!!
!!Finally the derivative matrices and damping matrix are computed using
!!\c sder, \c sder2, and \c setdmc. 
!>
!> @param[in] name
!> CHARACTER, array, takes x, y, or z as sirection.
!> @param[in] nv
!> INTEGER, takes the number of grid points.
!> @param[in] dv
!> REAL(db), takes the grid spacing.
!> @param[out] v
!> REAL(db), array, returns the coordinates.
!> @param[out] der1v
!> REAL(db), array, returns matrix for first derivative.
!> @param[out] der2v
!> REAL(db), array, returns matrix for second derivative.
!> @param[out] cdmpv
!> REAL(db), array, returns matrix for damping.
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: sder
!> @brief
!!This subroutine calculates the matrix for the first derivative. 
!> @param[out] der
!> REAL(db), array, returns derivative matrix.
!> @param[in] nmax
!> INTEGER, takes the dimension.
!> @param[in] d
!> REAL(db), takes the grid spacing. 
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: sder2
!> @brief
!!This subroutine calculates the matrix for the second derivative. 
!> @param[out] der
!> REAL(db), array, returns derivative matrix.
!> @param[in] nmax
!> INTEGER, takes the dimension.
!> @param[in] d
!> REAL(db), takes the grid spacing. 
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: setup_damping
!> @brief
!!This sets up the damping matrices by calls to \c setdmc for each
!!coordinate direction. 
!>
!> @details
!!The reason for not including this in \c init_grid 
!!is that it used only in the static calculation and
!!requires the damping parameter \c e0dmp, which is in the static
!!module. It has to be passed as a parameter because 
!!circular dependence of modules would result otherwise.
!>
!> @param[in] e0dmp
!>REAL(db), takes the damping parameter.
!--------------------------------------------------------------------------- 

  SUBROUTINE setup_damping(e0dmp)
    REAL(db),INTENT(IN) :: e0dmp
    CALL setdmc(der2x,nx,cdmpx,e0dmp)
    CALL setdmc(der2y,ny,cdmpy,e0dmp)
    CALL setdmc(der2z,nz,cdmpz,e0dmp)
  END SUBROUTINE setup_damping
!---------------------------------------------------------------------------  
! DESCRIPTION: setdmc
!> @brief
!!This subroutine calculates the matrices corresponding to the
!!one-dimensional operators
!!\f[ \frac{1}{1+\hat t/E_0} {\rm~with~} \hat
!!t=-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}. \f] 
!>
!> @details
!!Here \f$ E_0 \f$ is
!!the damping parameter called \c e0dmp in the code. Since the
!!kinetic energy operator \f$ \hat t \f$ contains the parameter \f$ \hbar^2/2m \f$,
!!which is force-dependent, this subroutine depends on module \c Forces.
!!
!!The calculation proceeds simply by constructing the unit matrix,
!!adding the operator to it to form the denominator, and then
!!calculating the inverse matrix using subroutine \c gauss.
!>
!> @param[in] der2
!> REAL(db), array, takes second derivative matrix.
!> @param[in] nmax
!> INTEGER, takes the number of grid points
!> @param[out] cdmp
!> REAL(db), array, returns the damping matrix. 
!> @param[in] e0dmp
!> REAL(db), takes the damping parameter.
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: gauss
!> @brief
!!This is a Fortran 95 implementation of the standard Gauss algorithm
!!with pivoting. It is simplified for the special case of computing
!!\f$ B=B^{-1}A \f$ with both matrices dimensioned <tt> (n,n) </tt>.
!> @param[in] a
!> REAL(db), array, takes the matrix a.
!> @param[in,out] b
!> REAL(db), array, takes and returns the matrix b.
!> @param[in] n
!> INTEGER, takes the number of grid points. 
!--------------------------------------------------------------------------- 
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
!---------------------------------------------------------------------------  
! DESCRIPTION: tabc_init_blochboundary
!> @brief
!!This routine calculates the bloch phases for TABC calculations
!--------------------------------------------------------------------------- 
  SUBROUTINE tabc_init_blochboundary
    INTEGER :: xbloch,ybloch,zbloch,nbloch
    nbloch=MAX(1,abs(tabc_x))*MAX(1,abs(tabc_y))*MAX(1,abs(tabc_z))
    WRITE(*,'(X,I4,A,I3,A,I3,A,I3,A)') nbloch,' sets of bloch twists (x:',abs(tabc_x),&
    ', y:',abs(tabc_y),', z:',abs(tabc_z),')' 
    IF (nbloch/=tabc_nprocs) &
      STOP 'number of processes not adequate for this setup of TABC'
    IF(tabc_x/=0) xbloch=MOD(tabc_myid,abs(tabc_x))
    IF(tabc_y/=0) ybloch=MOD(tabc_myid/MAX(1,abs(tabc_x)),abs(tabc_y))
    IF(tabc_z/=0) zbloch=MOD(tabc_myid/MAX(1,abs(tabc_x))/MAX(1,abs(tabc_y)),abs(tabc_z))
!
    IF (tabc_x<0) THEN
      bangx=(REAL(xbloch)+0.5d0)/REAL(abs(tabc_x))
    ELSE IF(tabc_x>0) THEN
      bangx=-1.0+(REAL(xbloch)+0.5d0)*2.0d0/REAL(abs(tabc_x))
    END IF
!
    IF (tabc_y<0) THEN
      bangy=(REAL(ybloch)+0.5d0)/REAL(abs(tabc_y))
    ELSE IF(tabc_y>0) THEN
      bangy=-1.0+(REAL(ybloch)+0.5d0)*2.0d0/REAL(abs(tabc_y))
    END IF
!
    IF (tabc_z<0) THEN
      bangz=(REAL(zbloch)+0.5d0)/REAL(abs(tabc_z))
    ELSE IF(tabc_z>0) THEN
      bangz=-1.0+(REAL(zbloch)+0.5d0)*2.0d0/REAL(abs(tabc_z))
    END IF
!
    WRITE(*,'(X,A,I4,A,I3,A,I3,A,I3)') 'TABC-ranks: ',tabc_myid,' x=',xbloch,' y=',ybloch,' z=',zbloch
    WRITE(*,'(X,A,F8.2,A,F8.2,A,F8.2)')  'local values: x:', bangx,' y:',bangy,' z:',bangz
  END SUBROUTINE tabc_init_blochboundary
END MODULE Grids


