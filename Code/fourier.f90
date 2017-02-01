!------------------------------------------------------------------------------
! MODULE: Fourier
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module initializes the \c FFTW3 package for doing Fourier
!!transforms \cite Fri05a. It needs the file \c fftw3.f from the
!!\c FFTW3 package to define the symbolic parameters for the
!!initialization calls.
!>
!>@details
!!Note that \c FFTW is used only for the complex transforms of the
!!wave functions and the Coulomb solver; the Fourier derivatives of the
!!real fields are handled by explicit matrix multiplication (see module
!!\c Trivial, routines \c rmulx, \c rmuly, and \c rmulz).
!!The definition of the plans for the Coulomb solver is contained in
!!subroutine \c coulinit in module \c Coulomb.
!------------------------------------------------------------------------------
MODULE Fourier
  USE params, ONLY: db,wflag
  USE Grids, ONLY: nx,ny,nz
  USE ISO_C_BINDING
  IMPLICIT NONE
  !>@name plans for full three-dimensional forward and backward transforms for both spin components.
  !>@{
  INTEGER(C_LONG),SAVE :: pforward,pbackward
  !>@}
  !>@name one-dimensional forward and backward transform in the x-direction,
  !!but for all values of \f$ y \f$, \f$ z \f$, and spin \f$ s \f$.
  !>@{
  INTEGER(C_LONG),SAVE :: xforward,xbackward
  !>@}
  !>@name one-dimensional forward and backward transform in the y-direction,
  !!but for all values of \f$ x \f$, \f$ z \f$, and spin \f$ s \f$.
  !>@{
  INTEGER(C_LONG),SAVE :: yforward,ybackward
  !>@}
  !>@name one-dimensional forward and backward transform in the z-direction,
  !!but for all values of \f$ x \f$, \f$ y \f$, and spin \f$ s \f$.
  !>@{
  INTEGER(C_LONG),SAVE :: zforward,zbackward
  !>@}
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: init_fft
!> @brief
!!In this subroutine the \c FFTW system is initialized
!>
!> @details
!!and a 2-wave
!!function array \c p is allocated temporarily for performing the
!!tests to make the plans efficient (we do not use dynamic allocation,
!!since having this array on the stack may produce different results
!!than on the heap - a point which has not been tested, though).
!!
!!This is followed by the calls to set up the plans. These have
!!different variations, since the way in which the index or indices over
!!which the transform is done are intertwined with the unaffected
!!indices is quite different for the different directions. Understanding
!!this Section requires a thorough familiarity with the FFTW
!!documentation and Fortran indexing.
!!
!!<b> An important point to consider when modifying the code:</b> the
!!setting-up of a plan must agree with its later use with respect to
!!whether the input and output arrays are the same (in-place transform)
!!or different. This is why \c p sometimes has a last index of 2 in
!!these calls: in the code all transforms are in-place except for \c xforward,
!!\c yforward, and \c zforward. It was found that
!!very strange things happen if this rule is not obeyed.
!--------------------------------------------------------------------------- 
  SUBROUTINE init_fft
    INCLUDE 'fftw3.f'
    COMPLEX(db),ALLOCATABLE :: p(:,:,:,:,:)
    INTEGER,SAVE :: FFTW_planflag             
! set option for FFTW setup here
!    FFTW_planflag=FFTW_ESTIMATE
!    FFTW_planflag=FFTW_MEASURE
    FFTW_planflag=FFTW_PATIENT
!    FFTW_planflag=FFTW_EXHAUSTIVE
!    FFTW_planflag=FFTW_MEASURE+FFTW_UNALIGNED
    ALLOCATE(p(nx,ny,nz,2,2))
    CALL dfftw_plan_dft_3d(pforward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_FORWARD, FFTW_planflag)
    CALL dfftw_plan_dft_3d(pbackward,nx,ny,nz,p(:,:,:,1,1),p(:,:,:,1,1), &
         FFTW_BACKWARD, FFTW_planflag)
    CALL dfftw_plan_many_dft(xforward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,2),0,1,nx, &
         FFTW_FORWARD, FFTW_planflag)
    CALL dfftw_plan_many_dft(xbackward,1,(/nx/),2*ny*nz, &
         p(:,:,:,:,1),0,1,nx,p(:,:,:,:,1),0,1,nx, &
         FFTW_BACKWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(yforward,1,(/ny/),nx, &
       p(:,:,:,:,1),0,nx,1,p(:,:,:,:,2),0,nx,1, &
       FFTW_FORWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(ybackward,1,(/ny/),nx, &
       p(:,:,:,:,1),0,nx,1,p(:,:,:,:,1),0,nx,1, &
       FFTW_BACKWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(zforward,1,(/nz/),nx*ny, &
       p(:,:,:,:,1),0,nx*ny,1,p(:,:,:,:,2),0,nx*ny,1, &
       FFTW_FORWARD, FFTW_planflag)
  CALL dfftw_plan_many_dft(zbackward,1,(/nz/),nx*ny, &
       p(:,:,:,:,1),0,nx*ny,1,p(:,:,:,:,1),0,nx*ny,1, &
       FFTW_BACKWARD, FFTW_planflag)
    IF(wflag) WRITE(*,*) '***** FFTW3 plans established *****'
    DEALLOCATE(p)
  END SUBROUTINE init_fft
END MODULE Fourier
