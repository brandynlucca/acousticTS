module param_quad
    integer, parameter :: knd = selected_real_kind(33)  ! quad precision (128-bit)
    logical, parameter :: debug = .false.
    logical, parameter :: warn = .false.
    logical, parameter :: output = .false.
    logical, parameter :: suffix = .false.
end module param_quad