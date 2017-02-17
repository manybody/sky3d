!------------------------------------------------------------------------------
! MODULE: Constraint
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module adds the fields and subroutines which are required to
!!run a static iteration with multipole constraints.
!>
!>@details
!!The wanted deformations are read in from namelist 'constraint'.
!!The code runs without constraint if this nmelist is no found in
!!the input.
!------------------------------------------------------------------------------
Module Constraint
  USE Params, ONLY: db,pi,iter
  USE Densities
  USE Grids, ONLY: nx,ny,nz,x,y,z,wxyz
  USE Levels, ONLY: nstmax,psi,nstloc,wocc,isospin,mass_number,schmid
  USE Parallel, ONLY: globalindex
  IMPLICIT NONE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: constr_field   !<this is the array of local
  !!constraining fields
  REAL(db),ALLOCATABLE,DIMENSION(:) :: lambda_crank !<vector of Lagrange parameters
  REAL(db),ALLOCATABLE,DIMENSION(:) :: dlambda_crank !<vector of correction on Lagrange parameters
  REAL(db),ALLOCATABLE,DIMENSION(:) :: goal_crank !<this is the vector of wanted expectation values
  REAL(db),ALLOCATABLE,DIMENSION(:) :: actual_crank !<this is the vector of actual expectation values
  REAL(db),ALLOCATABLE,DIMENSION(:) :: old_crank !<this is the vector of previous expectation values
  REAL(db),ALLOCATABLE,DIMENSION(:) :: actual_crank2 !<this is the vector of actual variances
  REAL(db),ALLOCATABLE,DIMENSION(:) :: qcorr 
  REAL(db) :: alpha20_wanted=-1D99,alpha22_wanted=-1D99 !<wanted expectation values
  REAL(db) :: c0constr=0.8D0,d0constr=0.1D-4,qepsconstr=0.3D0 !<parameters for Q-corrective step 
  REAL(db) :: dampgamma=1D0,damprad=6D0 !<parameters for damping multipole constraints
  REAL(db) :: corrlambda,actual_numb
  LOGICAL :: tconstraint=.FALSE.
  INTEGER :: numconstraint=0,iconstr
  INTEGER,PRIVATE:: is,ix,iy,iz,nst
  LOGICAL,PRIVATE,PARAMETER :: tprintcons=.TRUE.
CONTAINS
  !***********************************************************************
  SUBROUTINE init_constraint() 
    REAL(db), PARAMETER :: prefac20=SQRT(pi/5D0),prefac22=SQRT(1.2D0*pi)
    REAL(db), PARAMETER :: r0rms=0.93D0  ! corresponds to box R0=1.2
    REAL(db) :: fac20,fac22
    NAMELIST /constraint/alpha20_wanted,alpha22_wanted,&
           c0constr,d0constr,qepsconstr,damprad,dampgamma

    READ(5,constraint,ERR=99,END=98)

    IF(alpha20_wanted>-1D99) numconstraint=1+numconstraint
    IF(alpha22_wanted>-1D99) numconstraint=1+numconstraint
    IF(numconstraint==0) GO TO 98
    ALLOCATE(constr_field(numconstraint,nx,ny,nz,2))
    ALLOCATE(lambda_crank(numconstraint),dlambda_crank(numconstraint),&
             goal_crank(numconstraint))
    ALLOCATE(actual_crank(numconstraint),actual_crank2(numconstraint),&
             old_crank(numconstraint),qcorr(numconstraint))
    lambda_crank=-0.2D0
    OPEN(803,file='constraint_iter')
    WRITE(803,*) &
    'Constraint: iter iconstr old-value new-value variance  corrlambda lambda'

!    defining constrains, here preliminarily unscaled Cartesian quadrupoles
    iconstr=0
    IF(alpha20_wanted>-1D99) THEN       ! isoscalar axial quadrupole Q_20
        iconstr=1+iconstr
        fac20=prefac20/(r0rms**2*mass_number**1.666666667D0)
        goal_crank(iconstr)=alpha20_wanted
        DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx
                constr_field(iconstr,ix,iy,iz,1) = &
                   fac20*(2D0*z(iz)**2-x(ix)**2-y(iy)**2)/ &
                   (1D0+EXP((SQRT(z(iz)**2+y(iy)**2+x(ix)**2))-damprad)/dampgamma)
                constr_field(iconstr,ix,iy,iz,2) = constr_field(iconstr,ix,iy,iz,1) 
             ENDDO
          ENDDO
       ENDDO
!       IF(tprintcons) WRITE(*,*) 'fac20,field:',&
!           fac20,SUM(constr_field(iconstr,:,:,:,:))
    END IF
    IF(alpha22_wanted>-1D99) THEN       ! isoscalar triaxial quadrupole Q_22
        iconstr=1+iconstr
        fac22=prefac22/(r0rms**2*mass_number**1.666666667D0)
        goal_crank(iconstr)=alpha22_wanted
        DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx
                constr_field(iconstr,ix,iy,iz,1) = &
                   fac22*(x(ix)**2-y(iy)**2)/ &
                   (1D0+EXP((SQRT(z(iz)**2+y(iy)**2+x(ix)**2))-damprad)/dampgamma)
                constr_field(iconstr,ix,iy,iz,2) = constr_field(iconstr,ix,iy,iz,1) 
             ENDDO
          ENDDO
       ENDDO
       IF(tprintcons) WRITE(*,*) 'fac22,field:',&
           fac22,SUM(constr_field(iconstr,:,:,:,:))
    END IF
    IF(iconstr.NE.numconstraint) &
      STOP 'INIT_CONSTRAINT: inconsistent nr. of constraints'

    tconstraint = .TRUE.

    WRITE(6,'(a,2(1pg13.5))') 'With constraints: alpha_20,alpha_22=',&
           alpha20_wanted,alpha22_wanted
    RETURN

98  CONTINUE       ! override constraint
      WRITE(6,*)  'Without constraints'
    RETURN

99 STOP 'error in reading NAMELIST constraint'

  END SUBROUTINE init_constraint

  SUBROUTINE add_constraint(potential) 
    REAL(db),INTENT(IN OUT) :: potential(nx,ny,nz,2)

    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
       potential(ix,iy,iz,is)=potential(ix,iy,iz,is)- &
         SUM(lambda_crank(:)*constr_field(:,ix,iy,iz,is))
    END FORALL

!    IF(tprintcons) WRITE(6,*) 'constraint added to potential'

  END SUBROUTINE add_constraint

  SUBROUTINE tune_constraint(e0act,x0act)
    REAL(db), INTENT(IN) :: e0act,x0act
    REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: corrfield 
    ! compute Q-moments
    ! simplify variance by Q**2 moment, ignore non-diagonal elements
    actual_numb = wxyz*SUM(rho(:,:,:,:))
    DO iconstr=1,numconstraint
      actual_crank(iconstr) = wxyz*SUM(rho(:,:,:,:)*constr_field(iconstr,:,:,:,:))
      actual_crank2(iconstr) = wxyz*SUM(rho(:,:,:,:)*constr_field(iconstr,:,:,:,:)**2)
      actual_crank2(iconstr)= ABS(actual_crank2(iconstr)-actual_crank(iconstr)**2/actual_numb)
      qcorr(iconstr) = c0constr*(actual_crank(iconstr)-goal_crank(iconstr))/ &
                                (2D0*actual_crank2(iconstr)+d0constr)
    END DO    

    
!   correct waverfunctions towards wanted deformation
    ALLOCATE(corrfield(nx,ny,nz,2))
    FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)
      corrfield(ix,iy,iz,is) = EXP(-SUM(qcorr(:)*constr_field(:,ix,iy,iz,is)))
    END FORALL
    DO nst=1,nstmax
       psi(:,:,:,:,nst) = corrfield(:,:,:,:)*psi(:,:,:,:,nst)
    ENDDO
    DEALLOCATE(corrfield)
    CALL schmid
    rho=0.0D0
    tau=0.0D0
    current=0.0D0
    sdens=0.0D0
    sodens=0.0D0
    DO nst=1,nstloc
       CALL add_density(isospin(globalindex(nst)),wocc(globalindex(nst)),&
                        psi(:,:,:,:,nst),rho,tau,current,sdens,sodens)  
    ENDDO

!   update Lagrange parameters
    DO iconstr=1,numconstraint
      corrlambda = -qepsconstr*(MAX(e0act,1D0)/x0act)* &
               (actual_crank(iconstr)-old_crank(iconstr))/ &
               (2D0*actual_crank2(iconstr)+d0constr)
!      corrlambda=-corrlambda
!      corrlambda = corrlambda/lambda_crank(iconstr)
!      corrlambda = SIGN(MIN(ABS(corrlambda),2D0*ABS(qepsconstr*lambda_crank(iconstr))),corrlambda)
      lambda_crank(iconstr) = lambda_crank(iconstr)+corrlambda
      WRITE(803,'(2i4,6(1pg13.5))') iter,iconstr, &
        (actual_crank(iconstr)-old_crank(iconstr)), &
        actual_crank(iconstr), &
        actual_crank2(iconstr),corrlambda, &
        lambda_crank(iconstr),qcorr(iconstr)
    ENDDO

  END SUBROUTINE tune_constraint


  SUBROUTINE before_constraint(dens) 
    REAL(db),INTENT(IN OUT) :: dens(nx,ny,nz,2)
!   actual_numb = wxyz*SUM(rho(:,:,:,:))
    DO iconstr=1,numconstraint
       old_crank(iconstr) = wxyz*&
          SUM(rho(:,:,:,:)*constr_field(iconstr,:,:,:,:))
    END DO    
 
    RETURN

  END SUBROUTINE before_constraint


END Module Constraint
