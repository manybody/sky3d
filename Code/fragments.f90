!------------------------------------------------------------------------------
! MODULE: Fragments
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module is concerned with setting up the initial condition from
!!data files containing fragment wave functions, usually obtained in a
!!previous static calculation. 
!>
!>@details 
!!A typical application is the initialization of a heavy-ion reaction.
!!Beyond that, any number of fragments (limited by the parameter
!!variable \c mnof) can be put into the grid at prescribed positions
!!with given initial velocities. A condition is, however, that the grid
!!they are defined in is smaller than the new grid. If they are put
!!close to the boundary, they may have density appearing on the other
!!side because of periodicity; this may be acceptable if the boundary
!!condition is periodic. A calculation may be restarted on a larger
!!grid, but this may be accompanied by some loss of accuracy, since the
!!wave function outside the original regions is set to a constant small value.
!!
!!For the case of parallel \c MPI calculations, there is logic to read wave functions
!!distributed over a series of files as described in connection with
!!subroutine \c write_wavefunctions. Since only the dynamic case
!!can run under \c MPI at present, the only application of this is to
!!restart a dynamic calculation. The number of processors may be
!!different in the restart.

!------------------------------------------------------------------------------
MODULE Fragments
  USE Params
  USE Grids 
  USE Forces, ONLY: f,nucleon_mass
  USE Levels
  IMPLICIT NONE 
  SAVE
  PRIVATE
  LOGICAL :: fix_boost                 !<if this is set to true, the boost
  !!(initial kinetic energy) values are used unchanged from the input;
  !!otherwise they are calculated from the initial kinetic energy of
  !!relative motion \c ecm and the impact parameter \c b. This
  !!implies a two-body initial configuration. The initial motion is
  !!assumed to be in the \f$ (x,z) \f$-plane in this case.  Note that for two
  !!initial fragments \c fix_boost can also be set to true for special
  !!initial conditions.
  REAL(db) :: ecm                      !< The kinetic energy of relative motion in MeV.
  REAL(db) :: b                        !< The impact parameter in fm.
  CHARACTER(64) :: filename(mnof)      !<for each fragment this indicates the
  !!name of the file with the associated wave functions. One file can be
  !!used several times in the case of identical fragments, abut the code
  !!does not treat that as a special case.
  REAL(db) :: fcent(3,mnof)            !< for fragment no. \c i the initial
  !!three-dimensional position is determined by <tt> fcent(:,i) </tt>.
  REAL(db) :: fboost(3,mnof)           !< gives the initial motion <tt>
  !!fboost(:,i)</tt> of fragment \c i in 3 dimensions. This is not a
  !!velocity, but the total kinetic energy of the fragment <em> in that
  !!Cartesian direction</em> in MeV. The sign indicates the direction,
  !!positive or negative along the corresponding axis. The sum of
  !!absolute values thus is the total kinetic energy of the fragment. These
  !!values are used only if \c fix_boost is true.
  REAL(db) :: fcmtot(3,mnof)           !< Centers of mass of the fragments.
  REAL(db) :: fmass(mnof)              !< Masses of the fragments.
  REAL(db) :: fcharge(mnof)            !< Charges of the fragments.
  INTEGER  :: fnneut(mnof)             !< Neutron number of the fragments.
  INTEGER  :: fnprot(mnof)             !< Proton number of the fragments.
  INTEGER  :: fnstmax(mnof)            !< Last index number of the wave functions of the fragments.
  INTEGER  :: fnpmin(2,mnof)           !< First index numbers for protons and neutrons of the fragments.
  INTEGER  :: fnpsi(2,mnof)            !< Last index numbers for protons and neutrons of the fragments.
  INTEGER  :: fnumber(2,mnof)          !< gives the number of wave functions in the
  !!fragment for each isospin. For initialization of a dynamic
  !!calculation only wave functions with a nonzero occupation are taken
  !!into account.
  INTEGER  :: fnewnpmin(2,mnof)        !<are the starting indices for each fragment's wave 
  !!functions in the total set of wave functions combining all the fragments.
  INTEGER  :: fnewnpsi(2,mnof)         !<are the ending indices for each fragment's wave 
  !!functions in the total set of wave functions combining all the fragments.
  INTEGER  :: fnx(mnof)                !<indicate the grid size in x-direction on which
  !!the fragment wave functions are defined.
  INTEGER  :: fny(mnof)                !<indicate the grid size in y-direction on which
  !!the fragment wave functions are defined.
  INTEGER  :: fnz(mnof)                !<indicate the grid size in z-direction on which
  !!the fragment wave functions are defined.
  PUBLIC   :: getin_fragments, read_fragments
  PRIVATE  :: locate,phases
CONTAINS 
!---------------------------------------------------------------------------  
! DESCRIPTION: getin_fragments
!> @brief
!!This subroutine reads the input variables, makes some consistency
!!checks, and then prepares for the initialization.
!>
!> @details
!!A restart is set up by placing one fragment taken from \c wffile at
!!the origin with zero velocity.
!!
!!The main loop is over fragments. The fragment files are opened a first
!!time to obtain information about the fragments, which is stored in
!!fragment-specific arrays. The checks include agreement of the forces
!!and grid spacings as well as that the fragment grid is not larger than
!!the new grid.
!!
!!The code at this point makes a distinction between static and dynamic
!!modes: for a static calculation it is assumed that the data may be
!!needed for a restart. In this case all wave functions are read in even
!!if not occupied even fractionally. For the dynamic case only those
!!with non-zero occupation are input as determined from \c fwocc.
!!This yields a reduced value for \c fnumber.
!!
!!<b> Note that at present it is assumed that the static wave functions
!!  are ordered in ascending energy, so that the occupied ones will
!!  start at index one and all empty states will be at the uppermost
!!  index positions.</b>
!!
!!Following this, the index positions in the new wave function array are
!!calculated in \c fnewnpmin and \c fnewnpsi by adding the
!!numbers of wave functions of each fragment-specific successively. Note
!!that at this stage the proton indices are still counted from one.
!!
!!After this input loop, the indices \c npmin and \c npsi as well
!!as the total number \c nstmax for the combined system are
!!calculated and finally the proton indices where the fragment wave
!!functions should be inserted are shifted to behind the neutrons,
!!i.e., starting at \c npmin(2).
!!
!!At the end the two-body initialization \c twobody_init is called
!!for the dynamic case with two fragments and <tt> fix_boost=.FALSE.</tt>.
!!This calculates the \c fboost values from \c ecm and \c b. 
!--------------------------------------------------------------------------- 
  SUBROUTINE getin_fragments
    INTEGER :: i
    REAL(db) :: fdx,fdy,fdz,fwxyz
    CHARACTER(8) :: forcename
    NAMELIST /fragments/ filename,fcent,fboost,ecm,b,fix_boost
    IF(trestart) THEN
       filename(1)=wffile
       fcent(:,1)=0.D0
       fboost(:,1)=0.D0
       fix_boost=.TRUE.
       IF(wflag) WRITE(*,*) 'Reading restart as one fragment from file ',wffile
    ELSE
       fix_boost=.FALSE.
       READ(5,fragments)
    END IF
    IF(nof/=2.AND..NOT.fix_boost) THEN  
       IF(wflag) WRITE(*,"(//,a,//)") "Non-fixed boost only for nof=2"  
       STOP  
    ENDIF
    IF(wflag) THEN
       IF(fix_boost) THEN
          WRITE(*,*) '***** Fixed boost values used *****'
       ELSE
          WRITE(*,*) '***** Boost values computed from twobody kinematics *****'
       END IF
    END IF
    fnewnpmin(:,1)=1
    DO i=1,nof
       OPEN(UNIT=scratch,FILE=filename(i),STATUS='old',FORM=&
            'unformatted')
       READ(scratch) iter,time,forcename,fnstmax(i),fnneut(i),fnprot(i), &
            fnumber(:,i),fnpsi(:,i),fcharge(i),fmass(i),fcmtot(:,i)
       IF(trestart) fcmtot(:,i)=0.D0
       IF(forcename/=f%name) THEN
          IF(wflag) WRITE(*,*) 'Forces do not agree: ',forcename, &
               ' : ',f%name,' Fragment:',filename(i)
          STOP
       END IF
       READ(scratch) fnx(i),fny(i),fnz(i),fdx,fdy,fdz,fwxyz
       CLOSE(UNIT=scratch)
       IF(fnx(i)>nx.OR.fny(i)>ny.OR.fnz(i)>nz) &
            CALL errormsg('Fragment dimensioning larger than present', &
            filename(i))
       IF(MAX(ABS(fdx-dx),ABS(fdy-dy),ABS(fdz-dz)) &
            >0.01D0) CALL errormsg('Grid spacing different for', &
            filename(i))
       IF(wflag) THEN
          WRITE(*,"(A,I2,2A)") " ***** Data for fragment #",i,' from file ',&
               filename(i)
          WRITE(*,"(A,1P,3E12.4)") " Location(x,y,z):  ",fcent(:,i)
          WRITE(*,"(2(A,F9.4))") " Mass number: ",fmass(i), &
               ' Charge number: ',fcharge(i)   
          IF(fix_boost) &
               WRITE(*,"(A,1P,3E12.4)") " Boost(x,y,z):     ",fboost(:,i)
       END IF
       ! for static restart read also unoccupied wave functions
       IF(tstatic.AND.nof==1) THEN
          fnumber(1,i)=fnpsi(1,i)
          fnumber(2,i)=fnpsi(2,i)-fnpsi(1,i)
       ENDIF
       ! get isospin starting points in fragment file
       fnpmin(1,i)=1
       fnpmin(2,i)=fnpsi(1,i)+1
       IF(i<nof) &
            fnewnpmin(:,i+1)=fnewnpmin(:,i)+fnumber(:,i)
       fnewnpsi(:,i)=fnewnpmin(:,i)+fnumber(:,i)-1
    END DO
    ! determine total wave function numbers
    nneut=SUM(fnneut)
    nprot=SUM(fnprot)
    npsi(1)=SUM(fnumber(1,:))
    npsi(2)=SUM(fnumber(2,:))
    nstmax=npsi(1)+npsi(2)
    npmin(1)=1
    npmin(2)=npsi(1)+1
    npsi(2)=nstmax
    fnewnpmin(2,:)=fnewnpmin(2,:)+npsi(1)
    fnewnpsi(2,:)=fnewnpsi(2,:)+npmin(2)-1
    charge_number=nprot
    mass_number=nneut+nprot
    IF(.NOT.fix_boost.AND.nof==2.AND.tdynamic) CALL twobody_init
  END SUBROUTINE getin_fragments
!---------------------------------------------------------------------------  
! DESCRIPTION: read_fragments
!> @brief
!!This subroutine prints summary information about which fragment wave
!!functions occupy which range of indices. Then it does a loop over
!!fragments to read in their wave functions using 
!!\c read_one_fragment, followed by applying the boost to them using
!!\c boost_fragment.
!--------------------------------------------------------------------------- 
  SUBROUTINE read_fragments
    INTEGER :: iff
    IF(wflag) &
         WRITE(*,*) '***** Input of fragment wave functions *****'
    DO iff=1,nof  
       CALL read_one_fragment(iff)
       CALL boost_fragment(iff)  
    ENDDO
    IF(wflag) &
         WRITE(*,*) '***** All fragments loaded and boosted *****'
  END SUBROUTINE read_fragments
!---------------------------------------------------------------------------  
! DESCRIPTION: read_one_fragment
!> @brief
!!This subroutine reopens a fragment file
!!for the fragment indexed by \c iff and reads the wave functions,
!!inserting them at the correct index into the new wave function array
!!while also moving them to the desired center-of-mass position.
!>
!> @details
!!The principal loop is over isospin. The index limits in the fragment
!!file for that isospin are put into \c il and \c iu, and the new
!!indices into \c newil and \c newiu. In the next step the grid
!!coordinates and the single-particle quantities are read. The latter
!!are copied into the new arrays and the isospin is also recorded. Then
!!this information is printed.
!!
!!Once this loop is concluded, the spatial shift is prepared. The shift
!!is calculated from the difference between the desired position \c fcent 
!!with respect to the origin of the new coordinates given by
!!\c x etc., and the fragment center-of-mass \c fcmtot with
!!respect to the origin in fragment coordinates \c fx etc.
!!Subroutine \c phases is used to calculate essentially the shift
!!phase factor \f$ \exp(-\I\vec k\cdot\Delta\vec r) \f$, which is a product of
!!phases in each coordinate direction \c akx, \c aky, and \c akz. 
!!These have an index corresponding to \f$ \vec k \f$ in the Fourier transform.
!!
!!Now the index arrays \c fnode and \c flocalindex are read,
!!which indicate where the wave functions are to be found in case of an
!!MPI job. The logical variable \c multifile records whether this is
!!the case by testing whether any of the wave functions was stored on
!!other than node 0.
!!
!!The following loop runs over the new index positions.
!!If the wave function is not on the current node, nothing is done
!!except for ignoring the input record.
!!
!!For the case of multiple files for one fragment (which is recognized
!!by not all the indices \c fnode being zero, a short subroutine \c locate 
!!is used to position input at the correct location.
!!
!!Otherwise the wave function is read from the fragment file using the
!!fragment grid dimensions into a variable \c ps1 defined with the
!!full new dimension, then it is Fourier transformed, multiplied by the
!!phase factor, and transformed back. Finally it is inserted into the
!!total wave function array \c psi, where any zeroes are replaced by
!!a small number, presently set to \f$ 10^{-20} \f$. 
!>
!> @param[in] iff
!> INTEGER, takes the umber of the fragment which is read in.
!--------------------------------------------------------------------------- 
  SUBROUTINE read_one_fragment(iff)
    USE Parallel, ONLY: node,mpi_myproc,localindex
    USE Fourier
    INTEGER,INTENT(IN) :: iff
    LOGICAL :: multifile
    INTEGER :: ipn
    COMPLEX(db) :: ps1(nx,ny,nz,2),akx(nx),aky(ny),akz(nz)
    INTEGER :: iq,is,nst,oldnst,newnst,ix,iy,iz,iold,inew
    REAL(db) :: cmi(3)
    REAL(db) :: fx(fnx(iff)),fy(fny(iff)),fz(fnz(iff))
    REAL(db),DIMENSION(fnstmax(iff)) :: fwocc,fsp_energy,fsp_parity, &
         fsp_norm,fsp_kinetic,fsp_efluct1
    INTEGER,DIMENSION(fnstmax(iff)) :: fnode,flocalindex
    OPEN(UNIT=scratch,FILE=filename(iff),STATUS='old',FORM=&
         'unformatted')
    READ(scratch) 
    READ(scratch) 
    READ(scratch) fx,fy,fz
    ! read the s.p.quantities, reinserting them into the correct positions
    READ(scratch) fwocc,fsp_energy,fsp_parity,fsp_norm,fsp_kinetic,fsp_efluct1
    ps1=0.0D0
    IF(wflag) THEN
       WRITE(*,'(A,I3)') ' ***** S.p. levels for fragment #',iff
       WRITE(*,'(A,3I5,A,2I5)') '   Iso   #new  #old  Occup.     E_sp        &
            &Parity       Norm       E_kin      E_fluct' 
    END IF
    DO iq=1,2
       DO inew=fnewnpmin(iq,iff),fnewnpsi(iq,iff)
          iold=fnpmin(iq,iff)+inew-fnewnpmin(iq,iff)
          wocc(inew)=fwocc(iold)
          sp_energy(inew)=fsp_energy(iold)
          sp_parity(inew)=fsp_parity(iold)
          sp_norm(inew)=fsp_norm(iold)
          sp_kinetic(inew)=fsp_kinetic(iold)
          sp_efluct1(inew)=fsp_efluct1(iold)
          isospin(inew)=iq
          IF(wflag) THEN
             WRITE(*,'(3I6,6G12.5)') &
                  iq,inew,iold,wocc(inew),sp_energy(inew),sp_parity(inew), &
                  sp_norm(inew),sp_kinetic(inew),sp_efluct1(inew)
          ENDIF
       END DO
    ENDDO
    cmi(1)=(fcent(1,iff)-fcmtot(1,iff)-x(1)+fx(1))/(nx*dx)
    cmi(2)=(fcent(2,iff)-fcmtot(2,iff)-y(1)+fy(1))/(ny*dy)
    cmi(3)=(fcent(3,iff)-fcmtot(3,iff)-z(1)+fz(1))/(nz*dz)
    IF(nof==2.AND.wflag) WRITE(*,'(A,I3,3F10.5)') &
         ' Translation for fragment ',iff,nx*cmi(1),ny*cmi(2),nz*cmi(3)
    CALL phases(nx,akx,cmi(1))
    CALL phases(ny,aky,cmi(2))
    CALL phases(nz,akz,cmi(3))
    ! Code for multiple files: find out which file each w.f. is on and at which index
    READ(scratch) fnode,flocalindex
    multifile=ANY(fnode/=0)
    ! start reading wave functions
    DO iq=1,2
       ! skip unoccupied neutron states in one-file case
       IF(iq==2.AND..NOT.multifile) THEN
          DO nst=fnumber(1,iff)+1,fnpsi(1,iff)
             READ(scratch)
          END DO
       END IF
       ! Read occupied states
       DO nst=1,fnumber(iq,iff)
          newnst=fnewnpmin(iq,iff)+nst-1
          oldnst=fnpmin(iq,iff)+nst-1 ! includes empty states
          IF(node(newnst)==mpi_myproc) THEN           
             ipn=localindex(newnst)
             ! open correct file and position for multifile case
             IF(multifile) CALL locate(fnode(oldnst),flocalindex(oldnst))
             ps1=0.D0
             READ(scratch) ps1(1:fnx(iff),1:fny(iff),1:fnz(iff),:)
             DO is=1,2
                CALL dfftw_execute_dft(pforward,ps1(:,:,:,is),ps1(:,:,:,is))
                FORALL(ix=1:nx,iy=1:ny,iz=1:nz) 
                   ps1(ix,iy,iz,is)=ps1(ix,iy,iz,is)*akx(ix)*aky(iy)*akz(iz) &
                        /DBLE(nx*ny*nz)
                END FORALL
                CALL dfftw_execute_dft(pbackward,ps1(:,:,:,is),ps1(:,:,:,is))
                WHERE(ABS(ps1)==0.D0) 
                   psi(:,:,:,:,ipn)=1.0D-20
                ELSEWHERE
                   psi(:,:,:,:,ipn)=ps1
                ENDWHERE
             ENDDO
          ELSE
             IF(.NOT.multifile) READ(scratch)
          ENDIF
       ENDDO
    ENDDO
    CLOSE(unit=scratch)
  END SUBROUTINE read_one_fragment
!---------------------------------------------------------------------------  
! DESCRIPTION: locate
!> @brief
!!This has the task to position the file for reading wave functions at
!!the correct place. It has as arguments the number of the file
!!(corresponding to the node number in the previous calculation) and the
!!index of the wave function in this file. The variable \c presentfile
!!keeps track of which file is currently open.
!>
!> @details
!!If we are starting to read a new file (which is always true initially,
!!as no wave functions are written into the header file \c wffile
!!itself), it closes the present file, composes the file name for the
!!new one and opens it.  Then the proper number of records are ignored
!!to position at the correct one. If the file has not changed, it simply
!!returns and reading continues sequentially. The logic of course
!!assumes that files are stored sequentially in each partial file, but
!!the division into partial files need not be the same as in the
!!previous run.
!>
!> @param[in] fileno
!> INTEGER, takes the number of the file.
!> @param[in] pos
!> INTEGER, takes the index of the wave function in this file.
!--------------------------------------------------------------------------- 
  ! position scratch file in correct location
  SUBROUTINE locate(fileno,pos)
    INTEGER,INTENT(IN) :: fileno,pos
    INTEGER,SAVE :: presentfile=-1
    INTEGER :: i
    CHARACTER(120) :: rsfp
    IF(fileno/=presentfile) THEN
       CLOSE(scratch)
       WRITE(rsfp,'(I3.3,''.'',A)') fileno,wffile
       OPEN(scratch,FORM='unformatted',FILE=rsfp,STATUS='OLD')
       DO i=1,pos-1
          READ(scratch)
       END DO
    END IF
  END SUBROUTINE locate
!---------------------------------------------------------------------------  
! DESCRIPTION: phases
!> @brief
!!This subroutine calculates the phase factors for a one-dimensional
!!translation \f$ \Delta x \f$ as
!!\f[ \exp\left(-\frac{2\pi\I k_j \Delta x}{L}\right). \f]
!>
!> @details
!!The shift was calculated in \c read_one_fragment, the
!!argument \c c includes a denominator \f$ L={\tt nx*dx} \f$, the total
!!length of the grid. The momentum \f$ k_j \f$ is determined in the usual way
!!for the finite Fourier transform.
!>
!> @param[in] n
!> INTEGER, takes the number of grid points.
!> @param[out] a
!> COMPLEX(db), array, returns the phases.
!> @param[in] c
!> REAL(db), takes the shift  
!--------------------------------------------------------------------------- 
  PURE SUBROUTINE phases(n,a,c)
    INTEGER :: n,i,k
    COMPLEX(db),INTENT(OUT) :: a(n)
    REAL(db) :: c
    INTENT(IN) :: n,c
    DO i=1,n
       k=i-1
       IF(i>(n+1)/2) k=k-n
       a(i)=EXP(CMPLX(0.D0,-2.D0*pi*k*c,db))
    END DO
  END SUBROUTINE phases
!---------------------------------------------------------------------------  
! DESCRIPTION: twobody_init
!> @brief
!!The purpose of this subroutine is to calculate the boost values from
!!\c ecm and \c b in the two-body case. The calculation can be
!!followed with an elementary understanding of the two-body system
!!kinematics, so we just give a brief overview.
!>
!> @details
!!The reduced mass \c xmu, the relative velocity \c vrel and the
!!angular momentum \c xli are calculated first. Then the components
!!of the vector linking the two centers of mass \c dix and \c diz
!!are divided by the length of this vector \c roft to get the
!!direction cosines.
!!
!!The Coulomb energy is calculated assuming two point charges at
!!distance \c roft and is subtracted from \c ecm to yield the
!!kinetic energy remaining at this distance, from which the relative 
!!velocity \c vrel_d is calculated. Since the total center
!!of mass is assumed to be at rest, the velocities of the fragments \c v1
!!and \c v2 can then be simply obtained. The instantaneous impact parameter
!!\c b_d  results from the angular momentum and the relative
!!velocity, and the angle by which the two fragments need to miss each other in
!!order to realize this impact parameter is computed in \c sint and \c cost.
!!
!!This now makes it possible to calculate the boost velocity components,
!!which are converted into kinetic energies to conform to the definition
!!of \c fboost as given in the input. Note that these energies are
!!signed to indicate the direction of motion.
!!
!!Both velocities and energies are printed. 
!--------------------------------------------------------------------------- 
  SUBROUTINE twobody_init
    REAL(db) :: vrel,vrel_d,b_d,v1,v2,ec,dix,diz,totmass,sint,cost, &
         xli,xmu,roft
    INTEGER :: i
    IF(ABS(fcent(2,1)-fcent(2,2))>1.D-5) &
         STOP 'Two fragments must be in the x-z plane'
    fboost=0.D0
    totmass=SUM(fmass(1:2))
    xmu=fmass(1)*fmass(2)/totmass*nucleon_mass
    vrel=SQRT(2.D0*ecm/xmu)
    xli=xmu*vrel*b/hbc
    dix=fcent(1,1)-fcent(1,2)
    diz=fcent(3,1)-fcent(3,2)
    roft=SQRT(dix**2+diz**2)
    dix=dix/roft; diz=diz/roft
    ec=e2*fcharge(1)*fcharge(2)/roft
    IF(ec>ecm) STOP 'Not enough energy to reach this distance'
    vrel_d=SQRT(2.D0*(ecm-ec)/xmu)
    v1=fmass(2)/totmass*vrel_d
    v2=fmass(1)/totmass*vrel_d
    b_d=xli*hbc/(xmu*vrel_d)
    sint=b_d/roft
    cost=SQRT(1.D0-sint**2)
    fboost(1,1)=-v1*(dix*cost-diz*sint)
    fboost(3,1)=-v1*(dix*sint+diz*cost)
    fboost(1,2)=v2*(dix*cost-diz*sint)
    fboost(3,2)=v2*(dix*sint+diz*cost)
    IF(wflag) THEN
       WRITE(*,*) '***** Two-body initialization *****'
       WRITE(*,"(2(A,F12.4),A)") " c. m. Energy:",ecm, &
            ' MeV. Impact parameter b: ',b,' fm'
       WRITE(*,"(A,F12.4)") " xli/hbar    :",xli  
       WRITE(*,*) 'Computed boost velocities in units of c'
       WRITE(*,'(A,I2,3G15.6)') (' Fragment #',i,fboost(:,i),i=1,2)
    END IF
    !convert into kinetic energy for that direction
    FORALL(i=1:2)
       fboost(:,i)=SIGN(0.5D0*fmass(i)*fboost(:,i)**2*nucleon_mass, &
            fboost(:,i))
    END FORALL
    IF(wflag) THEN
       WRITE(*,*) 'Total translational kinetic energies in units of MeV'
       WRITE(*,'(A,I2,3G15.6)') (' Fragment #',i,fboost(:,i),i=1,2)
    END IF
  END SUBROUTINE twobody_init
!---------------------------------------------------------------------------  
! DESCRIPTION: boost_fragment
!> @brief
!!This subroutine multiplies the configuration-space wave functions of
!!fragment no. \c iff by plane-wave factors to give them
!!translational motion. 
!>
!> @details
!!This is done by calculating the wave number
!!vector \c akf from the kinetic energy. There is one subtle point:
!!since the boost is applied to <em> single-particle</em> wave functions,
!!the momentum \f$ k \f$ should be the correct one for <em> one nucleon </em>. Thus
!!the total kinetic energy of the fragment is
!!\f[ T=A\frac{\hbar^2}{2m}k^2, \f]
!!where \f$ A \f$ is the mass number of the fragment and \f$ m \f$ the nucleon mass.
!!Solving for \f$ k \f$ and taking the sign into account yields the expression
!!in the code.
!!
!!The rest is then straightforward.
!>
!> @param[in] iff
!> INTEGER, takes the number of the fragment to be boosted. 
!--------------------------------------------------------------------------- 
  SUBROUTINE boost_fragment(iff)  
    USE Parallel, ONLY: node,localindex,mpi_myproc
    INTEGER,INTENT(IN) :: iff
    REAL(db) :: sb(3),tb(3),akf(3)
    INTEGER :: iq,nst,is,ix,iy,iz
    sb=0.0D0  
    tb=fboost(:,iff)  
    WHERE(tb/=0.0D0) 
       sb=tb/ABS(tb)  
    END WHERE
    DO iq=1,2  
       akf=sb*SQRT(1.0D0/f%h2m(iq)*ABS(tb)/fmass(iff))  
       DO nst=fnewnpmin(iq,iff),fnewnpsi(iq,iff)
          IF(node(nst)==mpi_myproc) THEN
             FORALL(is=1:2,ix=1:nx,iy=1:ny,iz=1:nz)             
                psi(ix,iy,iz,is,localindex(nst))=&
                     psi(ix,iy,iz,is,localindex(nst)) * EXP( &
                     CMPLX(0.D0,akf(1)*x(ix)+akf(2)*y(iy)+akf(3)*z(iz),db))
             END FORALL
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE boost_fragment
!---------------------------------------------------------------------------  
! DESCRIPTION: errormsg
!> @brief
!!Writes out error messages.
!!
!> @param[in] msg1
!> CHARACTER, array, takes the first message.
!> @param[in] msg2
!> CHARACTER, array, takes the second message.
!--------------------------------------------------------------------------- 
  SUBROUTINE errormsg(msg1,msg2)
    CHARACTER(*),INTENT(IN) :: msg1,msg2
    IF(wflag) WRITE(*,*) msg1,msg2
    STOP
  END SUBROUTINE errormsg
END MODULE Fragments
