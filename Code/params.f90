!------------------------------------------------------------------------------
! MODULE: Params
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains some general parameters critical to controlling
!!the code and some mathematical and physical constants used everywhere.
!------------------------------------------------------------------------------
MODULE Params
  IMPLICIT NONE
  !**********************************************************************
  !     data type definition                                            *
  !**********************************************************************
  !>@name General paramters
  !>@{
  INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100) !< this is a constant 
  !!determining the precision of real numbers in the code in a portable manner. 
  !!In practice it will usually be <tt> REAL(8) </tt>.
  !
  !**********************************************************************
  !     useful constants                                                *
  !**********************************************************************
  REAL(db),PARAMETER :: pi=3.14159265358979D0  !< the constant \f$ \pi \f$.
  REAL(db),PARAMETER :: hbc=197.32164D0        !< the constant \f$ \hbar c \f$ in units of MeV*fm.
  REAL(db),PARAMETER :: e2=1.43989D0           !< the electron charge squared. It is calculated as
  !!\f$ \alpha\hbar c \f$ and has units of MeV*fm.
  !>@}
  !
  !**********************************************************************
  !  names of files and units to be used                                *
  !**********************************************************************
  !>@name File names and units
  !>These variables allow the user to change the names of some of the
  !>files used by the code. Except for \c wffile, output is produced
  !>for an iteration or a time step at selected intervals.
  !>@{
  CHARACTER(LEN=80) :: wffile='none'                     !<file to contain 
  !!the static single-particle
  !!wave functions plus some additional data. It is also used for
  !!restarting an interrupted calculation and is rewritten regularly
  !!(see variables \c trestart and \c mrest). It can be turned off
  !!completely using the name <tt>'NONE'</tt>.
  CHARACTER(LEN=80) :: converfile='conver.res'           !<contains convergence 
  !!information for the static calculation.
  CHARACTER(LEN=80) :: monopolesfile='monopoles.res'     !<contains moment 
  !!values of monopole type.
  CHARACTER(LEN=80) :: dipolesfile='dipoles.res'         !<contains moment values of dipole type.
  CHARACTER(LEN=80) :: momentafile='momenta.res'         !<contains components of the total
  !!momentum.
  CHARACTER(LEN=80) :: energiesfile='energies.res'       !<energy data for time-dependent mode.
  CHARACTER(LEN=80) :: quadrupolesfile='quadrupoles.res' !<contains moment values of dipole type.
  CHARACTER(LEN=80) :: spinfile='spin.res'               !<time-dependent total, orbital, and spin
  !!angular-momentum data as three-dimensional vectors. 
  CHARACTER(LEN=80) :: extfieldfile='extfield.res'       !<contains the time-dependence of the
  !!expectation value of the external field. Present only if an
  !!external field for boost or time-dependent excitation is used.
  INTEGER,PARAMETER :: scratch=11                        !<unit number used for temporary storage.
  INTEGER,PARAMETER :: scratch2=12                       !<unit number used for temporary storage.
  !>@}
  !**********************************************************************
  !     basic parameters controlling the job                            *
  !**********************************************************************
  !>@name Switches
  !>@{
  LOGICAL :: tcoul=.TRUE.     !<indicates whether the Coulomb field should be included or not.
  LOGICAL :: tstatic          !<\c true for a static job. Not input
  !!directly but from the input variable \c imode, which is 1 for the static.
  LOGICAL :: tdynamic         !<\c true for a dynamic job. Not input
  !!directly but from the input variable \c imode, which is 1 for the static.
  LOGICAL :: tfft=.TRUE.      !<if \c true, the derivatives of the wave
  !!functions, but not of the densities, are done directly through FFT.
  !!Otherwise matrix multiplication is used, but with the matrix also
  !!obtained from FFT.
  LOGICAL :: trestart=.FALSE. !< if \c true, restarts the calculation
  !!from the \c wffile.
  LOGICAL :: dconstr=.FALSE. !< if \c true, performs the density constraint
  !!during the time evolution.
  !>@}
  !**********************************************************************
  ! parameters controlling printout frequency etc.
  !**********************************************************************
  !>@name Output control
  !>@{
  INTEGER :: mprint=100 !<control for printer output. If \c mprint
  !!is greater than zero, more detailed output is produced every \c
  !!mprint iterations or time steps on standard output.
  INTEGER :: mplot=0    !<if \c mplot is greater than zero, a
  !!printer plot is produced and the densities are dumped onto <tt>
  !!*.tdd</tt> files every \c  mplot time steps or iterations.
  INTEGER :: mrest=0    !<if greater than zero, a \c wffile is produced
  !!every \c mrest iterations or time steps.
  INTEGER :: mconstr=0  !<if greater than zero and \c dconstr is \c true,
  !!the density constraint is performed every \c mconstr time steps.
  !>@}
  !>@name Globally used variables
  !>@{
  INTEGER :: iter     !<number of the current iteration.
  !!Used in the static and dynamic mode.
  INTEGER :: diter     !<number of the current time step.
  !!Used in the dynamic mode to save time step.
  REAL(db) :: time    !<the simulation time of the current step in
  !!fm/c; only meaningful in a dynamic calculation.
  LOGICAL :: wflag    !<indicates whether printing is allowed. This
  !!is necessary for the parallel job to have only one processor print
  !!and concerns both the standard output and the <tt> *.res</tt> files.
  LOGICAL :: printnow !<this variable is set to true if conditions
  !!for printing are met, such as a certain interval in iteration number.
  !>@}
  !>@name Field output control
  !>@{
  INTEGER,PARAMETER :: nselect=10           !< parameter limiting how many fields can be
  !!selected for binary output, it is just the length of the character string \c writeselect.
  CHARACTER(LEN=nselect) :: writeselect='r' !< it is used to determine which fields
  !!should be output under the control of \c mplot.
  LOGICAL :: write_isospin=.FALSE.          !<if this is <tt>.FALSE.</tt>, the proton
  !!and neutron contributions of a field are added up before output. Otherwise both are written.
  !>@}
  !>@name Fragment number parameters
  !>@{
  INTEGER,PARAMETER ::  mnof=4 !<number of fragments for the initialization maximum allowed. 
  INTEGER :: nof               !<actual number of fragments for the initialization.               
  REAL(db) :: r0=1.2D0         !<nuclear radius parameter. The nuclear radius
  !!\f$ R=r_0A^{1/3} \f$ is used to compute the \f$ \beta \f$ and \f$ \gamma \f$ deformation parameters 
  !!in subroutine \c moments. Units: fm.
  !>@}
  !>@name Density constraint parameters
  !>@{
  REAL(db) :: c0=1.90D0   !<c0 parameter for density constraint
  REAL(db) :: d0=5e-5     !<d0 parameter for density constraint exchange term
  !>@}
END MODULE Params
