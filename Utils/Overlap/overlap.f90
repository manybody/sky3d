PROGRAM Analyze
  IMPLICIT NONE
  INTEGER :: nstmax,nneut,nprot,nx,ny,nz
  INTEGER :: iter,number(2),npsi(2)
  INTEGER :: nstmax2,nneut2,nprot2,nx2,ny2,nz2
  REAL(8) :: time,charge_number,mass_number,cm(3),cm2(3),dx,dy,dz, &
       wxyz,wxyz2
  CHARACTER(8) :: forcename,forcename2
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: psi1,psi2
  COMPLEX(8) :: norm,on,op
  COMPLEX(8),ALLOCATABLE :: ovn(:,:),ovp(:,:)
  REAL(8) :: fn,fp,dnmax,dnmin,dpmax,dpmin,onmax,opmax
  INTEGER :: i,j,m
  CHARACTER(40) :: file1,file2
  DO
     WRITE(*,*) 'Enter two file names on separate lines -'
     READ(*,'(A80)',END=1) file1
     READ(*,'(A80)',END=1) file2
! open files and check for compatibility
     OPEN(UNIT=11,FILE=file1,STATUS ='old',FORM='unformatted')
     OPEN(UNIT=12,FILE=file2,STATUS ='old',FORM='unformatted')
     READ(11) iter,time,forcename,nstmax,nneut,nprot,number,npsi, &
          charge_number,mass_number,cm  
     READ(12) iter,time,forcename2,nstmax2,nneut2,nprot2,number,npsi, &
          charge_number,mass_number,cm2
     IF(nstmax/=nstmax2.OR.nneut/=nneut2.OR.nprot/=nprot2) THEN
        WRITE(*,'(6(A,I4))') '# states ',nstmax,'/=',nstmax2, &
             ' or neutron number ',nneut,'/=',nneut2, &
             ' or proton number ',nprot,'/=',nprot2
        STOP
     END IF
     IF(forcename/=forcename2) THEN
        WRITE(*,'(4A)') 'Forces different :',forcename,' : ',forcename2
        STOP
     END IF
     READ(11) nx,ny,nz,dx,dy,dz,wxyz
     READ(12) nx2,ny2,nz2,dx,dy,dz,wxyz2
     IF(nx/=nx2.OR.ny/=ny2.OR.nz/=nz2.OR.ABS(wxyz-wxyz2)>0.01D0) THEN
        WRITE(*,'(2(A,3I4,F8.3))') 'Grid dimensions do not agree: ', &
             nx,ny,nz,wxyz,' : ',nx2,ny2,nz2,wxyz2
        STOP
     END IF
     ALLOCATE(psi1(nx,ny,nz,2,nstmax),psi2(nx,ny,nz,2,nstmax), &
          ovn(nneut,nneut),ovp(nprot,nprot))
     ! ignore unnecessary data
     DO i=1,3
        READ(11)
        READ(12)
     END DO
! read wave function values
     DO i=1,nstmax
        READ(11) psi1(:,:,:,:,i)
        READ(12) psi2(:,:,:,:,i)
     END DO
     CLOSE(11)
     CLOSE(12)
     WRITE(*,'(5A,F10.4)') 'Files: ',TRIM(file1),',',TRIM(file2), &
          ' Distance of c.m.: ',SQRT(SUM((cm-cm2)**2))
     WRITE(*,'(A,6F8.4)') 'CM: ',cm,cm2
! normalize wave functions
     DO m=1,nstmax
        norm=SCALPROD(psi1(:,:,:,:,m),psi1(:,:,:,:,m))
        psi1(:,:,:,:,m)=psi1(:,:,:,:,m)/SQRT(norm)
        norm=SCALPROD(psi2(:,:,:,:,m),psi2(:,:,:,:,m))
        psi2(:,:,:,:,m)=psi2(:,:,:,:,m)/SQRT(norm)
     ENDDO
! calculate neutron overlap matrix
     dnmin=2.D0
     dnmax=0.D0
     onmax=0.D0
     DO j=1,nneut
        DO i=1,nneut
           ovn(i,j)=SCALPROD(psi1(:,:,:,:,j),psi2(:,:,:,:,i))
           IF(i==j) THEN
              dnmax=MAX(dnmax,ABS(ovn(i,j)))
              dnmin=MIN(dnmin,ABS(ovn(i,j)))
           ELSE
              onmax=MAX(onmax,ABS(ovn(i,j)))
           END IF
        ENDDO
     ENDDO
     WRITE(*,*) 'Diagonal max, min, off-diagonal max: '
     WRITE(*,'(A/3G15.6)') 'Neutrons: ',dnmax,dnmin,onmax
! calculate proton overlap matrix
     dpmin=2.D0
     dpmax=0.D0
     opmax=0.D0
     DO j=1,nprot
        DO i=1,nprot
           ovp(i,j)=SCALPROD(psi1(:,:,:,:,npsi(1)+j),psi2(:,:,:,:,npsi(1)+i))
           IF(i==j) THEN
              dpmax=MAX(dpmax,ABS(ovp(i,j)))
              dpmin=MIN(dpmin,ABS(ovp(i,j)))
           ELSE
              opmax=MAX(opmax,ABS(ovp(i,j)))
           END IF
        ENDDO
     ENDDO
     WRITE(*,'(A/3G15.6)') ' Protons:',dpmax,dpmin,opmax
! now  calculate antisymmetrized overlap for protons and neutrons
     on=determinant(ovn,nneut)
     op=determinant(ovp,nprot)
     WRITE(*,'(A/3G15.6)') 'Neutron, proton, and total overlap:', &
          ABS(on),ABS(op),ABS(on)*ABS(op)
     DEALLOCATE(psi1,psi2,ovn,ovp)
  END DO
1 STOP
CONTAINS
  COMPLEX(8) FUNCTION scalprod(a,b)
    IMPLICIT NONE
    COMPLEX(8),INTENT(IN),DIMENSION(:,:,:,:) :: a,b
    scalprod=SUM(CONJG(a)*b)*wxyz
  END FUNCTION scalprod
  COMPLEX(8) function determinant(a,n)
    INTEGER,INTENT(IN) :: n
    COMPLEX(8) a(n,n),det(2),work(n)
    integer ipvt(n),info
    call zgefa(a,n,n,ipvt,info)
    if(info.ne.0) then
       write(*,*) ' Determinant error ',info
       stop
    endif
    call zgedi(a,n,n,ipvt,det,work,10)
    determinant=det(1)*10.D0**real(det(2))
  end function determinant
END PROGRAM ANALYZE
