! 
!     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
!
! NVIDIA CORPORATION and its licensors retain all intellectual property
! and proprietary rights in and to this software, related documentation
! and any modifications thereto.  Any use, reproduction, disclosure or
! distribution of this software and related documentation without an express
! license agreement from NVIDIA CORPORATION is strictly prohibited.
! 

! This module provides a simple facility for changing between single
! and double precision

module precision_m
#ifdef CUDA
  integer, parameter, public :: singlePrecision = kind(0.0) ! Single precision
  integer, parameter, public :: doublePrecision = kind(0.0d0) ! Double precision
  
  ! Comment out one of the lines below
  !integer, parameter, public :: fp_kind = singlePrecision
  integer, parameter, public :: fp_kind = doublePrecision
#endif
end module precision_m
