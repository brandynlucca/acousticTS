! Adapted from the param module distributed with prolate_swf.
! Original code Copyright (c) 2021 Arnie Lee Van Buren, licensed under MIT.
! See inst/COPYRIGHTS for the complete upstream license notice.
module param
    integer, parameter :: knd = selected_real_kind(15)
    logical, parameter :: debug = .false.
    logical, parameter :: warn = .false.
    logical, parameter :: output = .false.
    logical, parameter :: suffix = .false.
end module param
