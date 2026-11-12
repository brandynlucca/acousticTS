! Adapted from the param module distributed with prolate_swf.
! Original code Copyright (c) 2021 Arnie Lee Van Buren, licensed under MIT.
! See inst/COPYRIGHTS for the complete upstream license notice.
module param_quad
#ifdef __GFORTRAN__
    use iso_c_binding, only: c_float128
#endif
    implicit none
    public

#ifdef __GFORTRAN__
    ! Preserve binary128 precision and identify the interoperable callback kind.
    integer, parameter :: knd = c_float128
#else
    integer, parameter :: knd = selected_real_kind(33)  ! quad precision (128-bit)
#endif
    logical, parameter :: debug = .false.
    logical, parameter :: warn = .false.
    logical, parameter :: output = .false.
    logical, parameter :: suffix = .false.
end module param_quad
