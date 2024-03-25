!------------------------------------------------------------------------------
! MODULE: Spherical Harmonics
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the definition of the spherical harmonics function \f$Y_{LM}\f$
!>
!>@details 
!! It contains the definition \f$Y_{LM}\f$ in terms of Associate Legendre Polinomial (\f$P_{LM}\f$)
Module Spherical_Harmonics
   USE Params, ONLY: db,pi
CONTAINS
! DESCRIPTION: Fact(n)
!> @brief
!!Function Fact Calculates the Factorial.
   REAL(db) function Fact(n)result(Fct) 
      
      IMPLICIT NONE
      integer:: n,i
      ! real(db)::Fct
      Fct=1.0d0
      do i =1,n
         Fct = Fct*i
      end do
   end function
! DESCRIPTION: Fact2(n)
!> @brief
!!Function Fact2 Calculates the double Factorial.
   REAL(db) function Fact2(n)result(Fct2) 
      
      IMPLICIT NONE
      integer :: n,i
      ! real(db)::Fct2
      Fct2=1.0d0
      do while (n>0)
         Fct2 = Fct2*n
         n=n-2
      end do
   end function

! DESCRIPTION: Plm
!> @brief
!!Function Plm Calculates the Associate Legendre Polinomial
   real(db) function Plm(x,l,m)result(P_lm) 
      IMPLICIT NONE
      integer, intent(in) :: l,m
      integer :: em_i,em_i_lower,em_i_lower2,el,em
      real(db):: P00,P10,P11,x
      real(db), ALLOCATABLE :: Ps(:,:)

      if (m>l)then
         P_lm= 0
         end if
      if (l==0)then
         P_lm= 1
      end if
      P00=1
      P10=x
      P11 = -sqrt(1-x**2)
      
      ALLOCATE(Ps(l+1,(2*l+1)))
      Ps=0.0D0
      Ps(1,1)=P00
      Ps(2,1) = P11
      Ps(2,2) = P10
      Ps(2,3) = -Ps(2,1)/2
      do el = 2,l
         Ps(el+1,el-el+1) = ((-1)**el)*Fact2(2*el-1)*(1-x**2)**(el/2) 
         do em = el-1,-el,-1
            em_i = el-em+1
            if (em>=0)then
                  em_i_lower = el-em
                  em_i_lower2 = el-1-em
                  Ps(el+1,em_i) = (x*(2*el-1)*Ps(el,em_i_lower)-(el+em-1)*Ps(el-1,em_i_lower2))/(el-em)
            else if (em<0)then
                  Ps(el+1, em_i) = ((-1)**abs(m))*(Fact(l-m)/Fact(l+m))*Ps(el+1,el+1+em)
            end if
         end do
      end do
      P_lm = Ps(l+1,l+1-m)
   end function
! DESCRIPTION: Y_lm
!> @brief
!!Function Y_lm Calculates the spherical harmonics (Real part)
   real(db) function Y_lm(l,m,x,y,z)result(Ylm) 
      IMPLICIT NONE
      integer, intent(in) :: l,m
      real(db), intent(in) :: x,y,z
      real(db) :: r,cos_theta,mphi,const,ylm_im,ylm_re
      ! PI = 4.0d0*atan(1.0d0)

      r=sqrt(x**2+y**2+z**2)
      cos_theta = z/r
      if (y .ge. 0)then
         mphi = abs(m)*acos(x/(sqrt(x**2+y**2)))
      else
         mphi = -1*abs(m)*acos(x/(sqrt(x**2+y**2)))
      end if
      if (m .lt. 0) then
         const = sqrt(((2*l+1)*(Fact(l-abs(m))))/(4*PI*Fact(l+abs(m))))
         ylm_im = const*Plm(cos_theta,l,abs(m))*sin(mphi)
         Ylm = ylm_im!*((-1)**m)
         ! write(*,*) 'negative m',l,m,abs(m),mphi*180/pi,acos(cos_theta)*180/pi,x,y,z,sin(mphi)
      else
         const = sqrt(((2*l+1)*(Fact(l-m)))/(4*PI*Fact(l+m)))
         ylm_re = const*Plm(cos_theta,l,m)*cos(mphi)
         Ylm = ylm_re
         ! write(*,*) 'positive m',l,m,abs(m),mphi*180/pi,acos(cos_theta)*180/pi,x,y,z,cos(mphi)
      end if
      

   end function

End Module Spherical_Harmonics