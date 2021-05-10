! 
!     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
!
! NVIDIA CORPORATION and its licensors retain all intellectual property
! and proprietary rights in and to this software, related documentation
! and any modifications thereto.
!
!
!    These example codes are a portion of the code samples from the companion
!    website to the book "CUDA Fortran for Scientists and Engineers":
!
! http://store.elsevier.com/product.jsp?isbn=9780124169708
!

!
! Define the INTERFACE to the NVIDIA CUFFT routines
!
module cufft_m
#ifdef CUDA
  integer, public :: CUFFT_FORWARD = -1
  integer, public :: CUFFT_INVERSE =  1
  integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
  integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
  integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
  integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
  integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
  integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex


  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy') 
       use iso_c_binding
       integer(C_LONG) :: plan
     end subroutine cufftDestroy
  end interface cufftDestroy
  
  interface cufftSetStream
     subroutine cufftSetStream(plan, stream) bind(C,name='cufftSetStream') 
       use iso_c_binding
       use cudafor
       integer(C_LONG) :: plan
       integer(kind=cuda_stream_kind),value:: stream
     end subroutine cufftSetStream
  end interface cufftSetStream

  interface cufftExecC2C
     subroutine cufftExecC2C(plan, idata, odata, direction) &
          bind(C,name='cufftExecC2C') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       integer,value:: direction
       !pgi$ ignore_tr idata,odata
       complex(singlePrecision),device:: idata(*),odata(*)
     end subroutine cufftExecC2C
  end interface cufftExecC2C

  interface cufftExecZ2Z
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          bind(C,name='cufftExecZ2Z') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       integer,value:: direction
       !pgi$ ignore_tr idata,odata
       complex(doublePrecision),device:: idata(*),odata(*)
     end subroutine cufftExecZ2Z
  end interface cufftExecZ2Z

  interface cufftExecR2C
     subroutine cufftExecR2C(plan, idata, odata) &
          bind(C,name='cufftExecR2C') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       integer,value:: direction
       !pgi$ ignore_tr idata,odata
       real(singlePrecision),device:: idata(*)
       complex(singlePrecision),device:: odata(*)
     end subroutine cufftExecR2C
  end interface cufftExecR2C

  interface cufftExecD2Z
     subroutine cufftExecD2Z(plan, idata, odata) &
          bind(C,name='cufftExecD2Z') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       integer,value:: direction
       !pgi$ ignore_tr idata,odata
       real(doublePrecision),device:: idata(*)
       complex(doublePrecision),device:: odata(*)
     end subroutine cufftExecD2Z
  end interface cufftExecD2Z


  interface cufftExecR2Cinplace
     subroutine cufftExecR2Cinplace(plan, idata, odata) &
          bind(C,name='cufftExecR2C') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       integer,value:: direction
       !pgi$ ignore_tr idata,odata
       real(singlePrecision),device:: idata(*)
       real(singlePrecision),device:: odata(*)
     end subroutine cufftExecR2Cinplace
  end interface cufftExecR2Cinplace

  interface cufftExecD2Zinplace
     subroutine cufftExecD2Zinplace(plan, idata, odata) &
          bind(C,name='cufftExecD2Z') 
       use iso_c_binding
       use precision_m
       integer(C_LONG) :: plan
       !pgi$ ignore_tr idata,odata
       real(doublePrecision),device:: idata(*)
       real(doublePrecision),device:: odata(*)
     end subroutine cufftExecD2Zinplace
  end interface cufftExecD2Zinplace

  interface cufftPlan1d
     subroutine cufftPlan1d(plan, nx, type, batch) &
          bind(C,name='cufftPlan1d') 
       use iso_c_binding
       integer(C_LONG) :: plan
       integer,value:: nx, batch,type
     end subroutine cufftPlan1d
  end interface cufftPlan1d
  
  interface cufftPlanMany
     subroutine cufftPlanMany(plan, rank, n, inembed, istride, idist, & 
          onembed, ostride, odist,  &
          type, batch) bind(C,name='cufftPlanMany')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr n, inembed, onembed       
       integer(C_LONG) :: plan
       integer :: n, inembed, onembed
       integer, value:: rank, istride, ostride, idist, odist, type, batch
     end subroutine cufftPlanMany
  end interface cufftPlanMany
  
  interface cufftPlan2d
     subroutine cufftPlan2d(plan, nx, ny, type) &
          bind(C,name='cufftPlan2d')
       use iso_c_binding
       integer(C_LONG) :: plan
       integer,value:: nx, ny, type
     end subroutine cufftPlan2d
  end interface cufftPlan2d

  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type) &
          bind(C,name='cufftPlan3d')
       use iso_c_binding
       integer(C_LONG) :: plan
       integer,value:: nx, ny, nz, type
     end subroutine cufftPlan3d
  end interface cufftPlan3d
#endif
end module cufft_m
