!---------------------------------------------------------------------------  
! MODULE: User
!> @brief
!!This is included as a place to insert arbitrary user initialization of
!!the wave functions. 
!>
!> @details
!!It really does the same job as subroutine \c harmosc, which is a 
!!relatively complicated example. For this reason
!!a sample user initialization is also provided in the file \c
!!user_sample.f90. It produces Gaussian wave functions for three
!!alpha-particle-like nuclei separated by a distance \c d and with
!!radii \c r.
!!
!!The only routine that has to be defined is \c user_init
!!\c user_init, but for more complicated initializations there can 
!!be any number of additional procedures accompanying it in \c user.f90.
!!
!!The setup assumes that the relevant particle numbers \c nneut and
!!\c nprot and numbers of states \c npmin, \c npsi, and \c nstmax are 
!!set correctly using the static input. In this case we
!!assume <tt>nprot=6</tt>, <tt>nneut=6</tt>, <tt>npmin=1,7</tt>, <tt>npsi=6,12</tt>,
!!and <tt>nstmax=12</tt>.  Note that the occupation numbers <tt>wocc</tt>
!!still have to be set explicitly, in this case they are all unity.
!!
!!The routine then reads the parameters for the setup from namelist \c user. 
!!This namelist can be used to read anything desired. If there
!!is no user initialization, it is simply omitted from the input file.
!!
!!Now there is a loop over center positions with index \c ic. For
!!each of them, the appropriate Gaussian is calculated and put into wave
!!function # \c ic in the spin-up component; the spin-down component
!!is set to zero. Then the Gaussian is copied into index position <tt>ic+3</tt>, 
!!spin-down component, and finally the complete wave functions
!!are copied to the proton indices by adding 6.
!!
!!There is no need to orthonormalize the wave functions, since \c schmid 
!!is called before the static iterations are started. User
!!initialization for the dynamic case does not appear useful; it could
!!be easily done without modifying the code by running one static
!!iteration and using the wave-function file generated at the beginning
!!to initialize the dynamic calculation. 
!---------------------------------------------------------------------------
MODULE User
  USE Params
  USE Grids
  USE Levels
  IMPLICIT NONE
CONTAINS
  SUBROUTINE init_user
    STOP 'No user-supplied initialization'
  END SUBROUTINE init_user
END MODULE User
