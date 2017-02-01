!------------------------------------------------------------------------------
! MODULE: Coulomb
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module offers the subroutine \c poisson, which calculates the Coulomb 
!!potential from the proton density. The two boundary conditions allowed are 
!!distinguished by the value of the variable periodic. If it is true, the mesh 
!!is assumed to be periodic in all three directions and the Poisson equation 
!!is solved using the \f$\frac{1}{k^2}\f$ Green’s function in momentum space and 
!!assuming a homogeneous negative background annulling the total charge (jellium 
!!approximation). Otherwise the boundary condition is for an isolated charge 
!!distribution. In this case the potential is calculated by the method of 
!!doubling the grid in all directions with zero density, folding with the 
!!1/r-potential in momentum space and then restricting to the physical region.
!>
!>@details 
!!The Coulomb solver follows the ideas of \cite Eas79. We summarize
!!here the way it is realized in the code without detailed proofs.
!!
!!Two grids are used. We will call them here \f$\mathcal{G}_1\f$ and
!!\f$\mathcal{G}_2\f$. The grid \f$\mathcal{G}_1\f$ is our standard working grid.
!!We repeat here the essentials
!!\f{eqnarray*}{
!!\mathcal{G}_1:\;&
!!  \nu_x=1,...,{\tt nx}\;,\;
!!  \nu_y=1,...,{\tt ny}\;,\;
!!  \nu_z=1,...,{\tt nz}\;;\;
!!\\&
!!  \mbox{coordinate spacing: } {\tt dx}, {\tt dy}, {\tt dz};
!!\\&
!!  \mbox{momentum spacing: } 
!!  {\tt dk}_x
!!  =
!!  \frac{2\pi}{{\tt nx}\!\cdot\!{\tt dx}},
!!  {\tt dk}_y
!!  =
!!  \frac{2\pi}{{\tt ny}\!\cdot\!{\tt dy}},
!!  {\tt dk}_z
!!  =
!!  \frac{2\pi}{{\tt nz}\!\cdot\!{\tt dz}};
!!\f}
!!The second grid amounts to a doubled box in each direction, thus
!!\f{eqnarray*}{
!!\mathcal{G}_2:\;&
!!  \nu_x=1,...,2{\tt nx}\;,\;
!!  \nu_y=1,...,2{\tt ny}\;,\;
!!  \nu_z=1,...,2{\tt nz}\;;\;
!!\\&
!!  x_{\nu_x}=(\nu_x\!-\!{\tt nx}){\tt dx}\,,\,
!!  y_{\nu_y}=(\nu_y\!-\!{\tt ny}){\tt dy}\,,\,
!!  z_{\nu_z}=(\nu_z\!-\!{\tt nz}){\tt dz}\,;
!!\\&
!!  \mbox{momentum spacing: } 
!!  {\tt dk}_x
!!  =
!!  \frac{2\pi}{{2\tt nx}\!\cdot\!{\tt dx}},
!!  {\tt dk}_y
!!  =
!!  \frac{2\pi}{{2\tt ny}\!\cdot\!{\tt dy}},
!!  {\tt dk}_z
!!  =
!!  \frac{2\pi}{{2\tt nz}\!\cdot\!{\tt dz}};
!!\f}
!!The momenta \f$k_{n_i}\f$ in \f$\mathcal{G}_2\f$ are arranged around \f$k=0\f$
!!similar as in \f$\mathcal{G}_1\f$, but with half spacing \f${\tt dk}_i\f$.
!!For simplicity, we will use compact notation
!!\f$\mathbf{r}_{\bnu}=\left(x_{\nu_x},y_{\nu_y},z_{\nu_z}\right)\f$
!!and similarly for \f$\mathbf{k}_\mathbf{n}\f$.
!!
!!First, the Coulomb solver has to be initialized by defining the
!!appropriate Greens function which is done on \f$\mathcal{G}_2\f$.
!!The two necessary steps are:
!!  -# Prepare the Greens function in r-space as
!!     \f[ G(\mathbf{r}_{\bnu})=\frac{e^2}{|\mathbf{r}_{\bnu}|} \f]
!!     Special consideration has to be given to the singularity of the \f$ 1/r \f$
!!     potential. In the discretized version, the value at \f$ \vec r=0 \f$ should
!!     be replaced by the average of \f$ 1/r \f$ over a cuboid of size \f$\Delta
!!     x\times\Delta y\times\Delta z \f$. For equal spacing in all three
!!     directions this would lead to a value of \f$ 2.38/\Delta x \f$. Practical
!!     experimentation, however, showed that the underlying assumption of a
!!     point charge can be improved by some smeared-out density for nuclear
!!     applications; a value of \f$ 2.84/\Delta x \f$ was found to be optimal.  
!!     We use the expression 
!!     \f[ G(0) = \frac{2.84\sqrt{3}}{\sqrt{{\tt dx}^2+{\tt dy}^2+{\tt dz}^2}} \quad.\f]
!!  -# Prepare the Greens function in \f$ k \f$-space by
!!     3D fast Fourier transformation (FFT) on the double grid \f$\mathcal{G}_2\f$.
!!     \f$\tilde{G}(\mathbf{k}_{\bmu})={\tt FFT}\{G(\mathbf{r}_{\bnu})\}\f$.
!!     The array \f$ \tilde{G}(\mathbf{k}_{\bmu}) \f$ is stored 
!!     for further continued use.
!!.
!!Once properly prepared, the practical steps for computing
!!the Coulomb field \f$ U_\mathrm{Coul}(\mathbf{r}_{\bnu}) \f$ for a 
!!density \f$\rho(\mathbf{r}_{\bnu})\f$ given on \f$\mathcal{G}_1\f$ are
!!  -# Extend \f$ \rho_{\bnu} \f$ from  \f$ \mathcal{G}_1 \f$ to  \f$ \mathcal{G}_2 \f$ 
!!     by zero padding:
!!     \f[ \rho_2(\mathbf{r}_{\bnu}) = \left\{\begin{array}{lcl} \rho(\mathbf{r}_{\bnu}) &\mbox{ if }& 1\leq \nu_i \leq {\tt n}i \mbox{ for } i\in\{x,y,z\}\\ 0 &\mbox{ else }\end{array}\right. \f]
!!  -# Fourier transform the density in  \f$ \mathcal{G}_2 \f$:
!!     \f[ \rho_2(\mathbf{r}_{\bnu})\ \longrightarrow\ \tilde{\rho}_2(\mathbf{k}_{\bmu})={\tt FFT}\{\rho_2(\mathbf{r}_{\bnu})\} \f]
!!  -# Compute solution in \f$ k \f$-space by multiplication with the
!!     formerly prepared Greens function
!!     \f[ \tilde{U}_2(\mathbf{k}_{\bmu}) =\tilde{G}_2(\mathbf{k}_{\bmu})\tilde{\rho}_2(\mathbf{k}_{\bmu}) \f]
!!  -# Compute solution in \f$ r \f$-space by Fourier back transformation
!!     in \f$\mathcal{G}_2\f$:
!!     \f[ {U}_2(\mathbf{r}_{\bnu})={\tt FFT}^{-1}\{\tilde{U}_2(\mathbf{k}_{\bmu})\} \f]
!!  -# Map to standard grid  \f$ \mathcal{G}_1 \f$:
!!     \f[ {U}_\mathrm{Coul}(\mathbf{r}_{\bnu})=U_2(\mathbf{r}_{\bnu})\quad\mbox{for}\quad\bnu\in\mathcal{G}_1 \f]
!!.
!------------------------------------------------------------------------------
MODULE Coulomb
  USE Params, ONLY: db,pi,e2
  USE Grids, ONLY: nx,ny,nz,dx,dy,dz,wxyz,periodic
  USE Densities, ONLY: rho
  USE ISO_C_BINDING
  IMPLICIT NONE
  !>@name Dimensions of the grid on which the Fourier transform is calculated. 
  !>For the periodic case they are identical to the regular dimensions, for the 
  !>isolated case they are doubled.
  !>@{
  INTEGER,PRIVATE :: nx2,ny2,nz2 
  !>@}
  !>@name Plans for FFTW complex forward and reverse transforms 
  !>with array dimensions depending on the boundary condition.
  !>@{
  INTEGER(C_LONG),PRIVATE,SAVE :: coulplan1,coulplan2
  !>@}  
  REAL(db),ALLOCATABLE,SAVE :: wcoul(:,:,:)  !< the Coulomb potential as a three-dimensional array. Units: MeV.
  COMPLEX(db),PRIVATE,ALLOCATABLE,SAVE :: q(:,:,:) !<array for the complex Green’s function (isolated) or array
  !!of 1/r values. Its dimension also depends on the boundary condition.
  PUBLIC :: poisson,coulinit,wcoul
  PRIVATE :: initiq
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: poisson
!> @brief
!!This subroutine solves the Poisson equation by the Fourier method. 
!>
!> @details
!!For the periodic case, this means to just Fourier transform, multiply by
!!\f$1/k^2\f$, and the transform back. Note that the coefficient for
!!momentum zero must be zero also: this means that the total charge in
!!the box vanishes, corresponding to the jellium approximation.
!!In the non-periodic case we use use the trick of doubling the box size
!!in every direction and then folding with \f$1/r\f$ my multiplying the
!!Fourier transforms. The result in the physical box is then the correct
!!solution for the boundary condition of zero potential at infinity.
!!The key to success is not to use simply \f$1/k^2\f$ for the Fourier
!!transform of \f$1/r\f$, but to compute it on the actual grid the same way
!!as the densities and potentials are transformed
!---------------------------------------------------------------------------  
  SUBROUTINE poisson
    COMPLEX(db),ALLOCATABLE :: rho2(:,:,:)
    ALLOCATE(rho2(nx2,ny2,nz2))
    ! put proton density into array of same or double size, zeroing rest
    IF(.NOT.periodic) rho2=(0.D0,0.D0)
    rho2(1:nx,1:ny,1:nz)=rho(:,:,:,2)
    ! transform into momentum space
    CALL dfftw_execute_dft(coulplan1,rho2,rho2)
    ! add charge factor and geometric factors
    ! note that in the periodic case q has only a real part
    IF(periodic) THEN
       rho2=4.D0*pi*e2*REAL(q)*rho2
    ELSE
       rho2=e2*wxyz*q*rho2
    END IF
    ! transform back to coordinate space and return in wcoul
    CALL dfftw_execute_dft(coulplan2,rho2,rho2)
    wcoul=REAL(rho2(1:nx,1:ny,1:nz))/(nx2*ny2*nz2)
    DEALLOCATE(rho2)
  END SUBROUTINE poisson
!---------------------------------------------------------------------------  
! DESCRIPTION: coulinint
!> @brief
!!This subroutine does the necessary initialization. 
!!
!>
!> @details
!!It calculates the dimension for the Fourier transform depending on 
!!the boundary condition, allocates the necessary arrays and sets up 
!!the FFTW plans.
!!Then it composes the array \c q. For the periodic case this
!!contains the values of the Green's function \f$ 1/k^2 \f$, with zero at the
!!origin, while for the isolated case it first calculates the inverse
!!shortest distance $1/r$ from the origin with indices (1,1,1), replaces
!!the value at the origin and then
!!Fourier-transforms this to use it for folding the density with the
!!coordinate-space Green's function.   
!>
!--------------------------------------------------------------------------- 
  SUBROUTINE coulinit
    INCLUDE 'fftw3.f'
    REAL(db),ALLOCATABLE :: iqx(:),iqy(:),iqz(:)
    INTEGER :: i,j,k
    IF(ALLOCATED(q)) RETURN ! has been initialized already
    ! dimensions will be doubled for isolated distribution
    IF(periodic) THEN
       nx2=nx; ny2=ny; nz2=nz
    ELSE
       nx2=2*nx; ny2=2*ny; nz2=2*nz
    END IF
    ! allocated helper arrays
    ALLOCATE (wcoul(nx,ny,nz),q(nx2,ny2,nz2),iqx(nx2),iqy(ny2),iqz(nz2))
    ! set up FFTW plans
    CALL dfftw_plan_dft_3d(coulplan1,nx2,ny2,nz2,q,q, &
         FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
    CALL dfftw_plan_dft_3d(coulplan2,nx2,ny2,nz2,q,q, &
         FFTW_BACKWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
    ! calculate coordinate contributions to 1/k^2 or 1/r^2, respectively
    CALL initiq(nx2,dx,iqx)
    CALL initiq(ny2,dy,iqy)
    CALL initiq(nz2,dz,iqz)
    ! combine proper values, using dummy value at (1,1,1) to avoid div0
    FORALL(k=1:nz2,j=1:ny2,i=1:nx2) q(i,j,k)=iqx(i)+iqy(j)+iqz(k)
    q(1,1,1)=1.D0
    IF(periodic) THEN
       q=1.D0/REAL(q)
       q(1,1,1)=(0.D0,0.D0)
    ELSE
       q=1.D0/SQRT(REAL(q))
       q(1,1,1)=2.84D0/(dx*dy*dz)**(1.D0/3.D0)
       CALL dfftw_execute_dft(coulplan1,q,q)
    END IF
    DEALLOCATE(iqx,iqy,iqz)
  END SUBROUTINE coulinit
!---------------------------------------------------------------------------  
! DESCRIPTION: initiq
!> @brief
!!This subroutine calculates the contributions of each Cartesian index
!!to \f$r^2\f$ and \f$k^{-2}\f$. 
!>
!> @details
!!This depends on the boundary condition.  For
!!the periodic case, the values of the momenta are given by
!!\f[ k_i=\frac{2\pi}{n\Delta x}\bigl(0,\;1,\;2,\;\ldots \frac{n}2-1,\;
!!-\frac{n}2,\;-\frac{n}2+1,\;\ldots -1\bigr).\f]
!!For an isolated distribution the shortest distances to the point with
!!index 1 are calculated (periodicity used):
!!\f[ d_i=d\bigl(0,\;1,\;2,\;\ldots \frac{n}2-1,\;
!!-\frac{n}2,-\frac{n}2+1,\;\;\ldots 1\bigr).\f]
!!The input is the dimension along the coordinate direction given, the
!!output is the one-dimensional array \c iq containing these values
!!squared.
!>
!> @param[in] n
!> INTEGER, takes the number of grid points.
!> @param[in] d
!> REAL(db), takes the grid spacing.
!> @param[out] iq
!> REAL(db), array, returns the values for \f$r^2\f$ and \f$k^{-2}\f$.
!--------------------------------------------------------------------------- 
  SUBROUTINE initiq(n,d,iq)
    INTEGER,INTENT(IN) :: n
    REAL(db),INTENT(IN) :: d
    REAL(db),INTENT(OUT) :: iq(:)
    INTEGER :: i,ii
    DO i=1,n
       IF(i<=n/2) THEN
          ii=i-1
       ELSE
          ii=i-n-1
       ENDIF
       IF(periodic) THEN
          iq(i)=(2.D0*pi*ii/(n*d))**2
       ELSE
          iq(i)=(d*ii)**2
       END IF
    ENDDO
  END SUBROUTINE initiq
END MODULE Coulomb
