#ifdef USE_QUAD
module prolate_swf_quad
 use param_quad
#else
module prolate_swf
 use param
#endif

 ! ---------------------------------------------------------------------------
 ! Lightweight caches for repeated Legendre- and quadrature-related work.
 ! The PSMS backend frequently asks for nearby spheroidal orders at the same
 ! size parameter, so caching these auxiliary tables avoids rebuilding the same
 ! support arrays on every call into profcn.
 ! Modified by: Brandyn M. Lucca; March 2026
 ! ---------------------------------------------------------------------------
 integer, parameter :: pleg_cache_slots = 48
 integer, parameter :: qleg_cache_slots = 48

 type :: pleg_cache_entry
   logical :: valid = .false.
   integer :: m = -1
   integer :: iopd = -1
   integer :: narg = 0
   integer :: maxt = 0
   integer :: maxp = 0
   integer :: lim = 0
   integer :: ndec = 0
   integer :: nex = 0
   real(knd), allocatable :: barg(:)
   real(knd), allocatable :: pr(:,:), pdr(:,:)
   real(knd), allocatable :: pdnorm(:), pnorm(:)
   integer, allocatable :: ipdnorm(:), ipnorm(:)
   real(knd), allocatable :: alpha(:), beta(:), gamma(:)
   real(knd), allocatable :: coefa(:), coefb(:), coefc(:), coefd(:), coefe(:)
 end type

 type :: qleg_cache_entry
   logical :: valid = .false.
   integer :: m = -1
   integer :: lnum = 0
   integer :: limq = 0
   integer :: maxq = 0
   integer :: ndec = 0
   integer :: iqdml = 0
   integer :: iqml = 0
   integer :: itermpq = 0
   real(knd) :: x1 = 0.0_knd
   real(knd) :: qdml = 0.0_knd
   real(knd) :: qml = 0.0_knd
   real(knd) :: termpq = 0.0_knd
   real(knd), allocatable :: qdr(:), qr(:), qdl(:), ql(:)
   integer, allocatable :: iqdl(:), iql(:)
 end type

 type :: gauss_cache_entry
   logical :: valid = .false.
   integer :: ndec = 0
   integer :: n = 0
   real(knd), allocatable :: x(:), w(:)
 end type

 type(pleg_cache_entry), save :: pleg_cache(pleg_cache_slots)
 type(qleg_cache_entry), save :: qleg_cache(qleg_cache_slots)
 type(gauss_cache_entry), save :: gauss_cache
 integer, save :: pleg_cache_next = 1
 integer, save :: qleg_cache_next = 1

 contains

  logical function cache_real_equal(a, b)
    real(knd), intent(in) :: a, b
    real(knd) :: scale, tol

    ! Use a precision-scaled tolerance so cache lookup is robust to the small
    ! roundoff differences introduced by repeated calls from C++.
    scale = max(1.0_knd, max(abs(a), abs(b)))
    tol = 64.0_knd * epsilon(1.0_knd) * scale
    cache_real_equal = abs(a - b) <= tol
  end function

  logical function cache_real_vector_equal(a, b, n)
    real(knd), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    integer :: i

    cache_real_vector_equal = .false.
    if(size(a) < n .or. size(b) < n) return
    do i = 1, n
      if(.not. cache_real_equal(a(i), b(i))) return
    end do
    cache_real_vector_equal = .true.
  end function

  subroutine clear_pleg_cache_slot(idx)
    integer, intent(in) :: idx

    ! Reset and deallocate one cache slot before it is reused.
    pleg_cache(idx)%valid = .false.
    pleg_cache(idx)%m = -1
    pleg_cache(idx)%iopd = -1
    pleg_cache(idx)%narg = 0
    pleg_cache(idx)%maxt = 0
    pleg_cache(idx)%maxp = 0
    pleg_cache(idx)%lim = 0
    pleg_cache(idx)%ndec = 0
    pleg_cache(idx)%nex = 0
    if(allocated(pleg_cache(idx)%barg)) deallocate(pleg_cache(idx)%barg)
    if(allocated(pleg_cache(idx)%pr)) deallocate(pleg_cache(idx)%pr)
    if(allocated(pleg_cache(idx)%pdr)) deallocate(pleg_cache(idx)%pdr)
    if(allocated(pleg_cache(idx)%pdnorm)) deallocate(pleg_cache(idx)%pdnorm)
    if(allocated(pleg_cache(idx)%pnorm)) deallocate(pleg_cache(idx)%pnorm)
    if(allocated(pleg_cache(idx)%ipdnorm)) deallocate(pleg_cache(idx)%ipdnorm)
    if(allocated(pleg_cache(idx)%ipnorm)) deallocate(pleg_cache(idx)%ipnorm)
    if(allocated(pleg_cache(idx)%alpha)) deallocate(pleg_cache(idx)%alpha)
    if(allocated(pleg_cache(idx)%beta)) deallocate(pleg_cache(idx)%beta)
    if(allocated(pleg_cache(idx)%gamma)) deallocate(pleg_cache(idx)%gamma)
    if(allocated(pleg_cache(idx)%coefa)) deallocate(pleg_cache(idx)%coefa)
    if(allocated(pleg_cache(idx)%coefb)) deallocate(pleg_cache(idx)%coefb)
    if(allocated(pleg_cache(idx)%coefc)) deallocate(pleg_cache(idx)%coefc)
    if(allocated(pleg_cache(idx)%coefd)) deallocate(pleg_cache(idx)%coefd)
    if(allocated(pleg_cache(idx)%coefe)) deallocate(pleg_cache(idx)%coefe)
  end subroutine

  integer function find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    integer, intent(in) :: m, iopd, ndec, nex, narg, maxt, lim, maxp
    real(knd), intent(in) :: barg(:)
    integer :: i

    ! Only reuse a slot when the order, options, precision metadata, argument
    ! vector, and retained dimensions are all compatible.
    find_pleg_cache = 0
    do i = 1, pleg_cache_slots
      if(.not. pleg_cache(i)%valid) cycle
      if(pleg_cache(i)%m /= m) cycle
      if(pleg_cache(i)%iopd /= iopd) cycle
      if(pleg_cache(i)%ndec /= ndec) cycle
      if(pleg_cache(i)%nex /= nex) cycle
      if(pleg_cache(i)%narg /= narg) cycle
      if(pleg_cache(i)%maxt /= maxt) cycle
      if(pleg_cache(i)%lim < lim) cycle
      if(pleg_cache(i)%maxp < maxp) cycle
      if(.not. cache_real_vector_equal(pleg_cache(i)%barg, barg, narg)) cycle
      find_pleg_cache = i
      return
    end do
  end function

  subroutine load_pleg_cache(idx, maxt, narg, lim, maxp, barg, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: idx, maxt, narg, lim, maxp
    real(knd), intent(out) :: barg(maxt), pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(out) :: ipdnorm(maxt), ipnorm(maxt)

    ! Copy cached arrays back into the caller's workspaces.
    barg(1:narg) = pleg_cache(idx)%barg(1:narg)
    pr(1:maxt, 1:lim) = pleg_cache(idx)%pr(1:maxt, 1:lim)
    pdr(1:maxt, 1:lim) = pleg_cache(idx)%pdr(1:maxt, 1:lim)
    pdnorm(1:maxt) = pleg_cache(idx)%pdnorm(1:maxt)
    ipdnorm(1:maxt) = pleg_cache(idx)%ipdnorm(1:maxt)
    pnorm(1:maxt) = pleg_cache(idx)%pnorm(1:maxt)
    ipnorm(1:maxt) = pleg_cache(idx)%ipnorm(1:maxt)
    alpha(1:lim) = pleg_cache(idx)%alpha(1:lim)
    beta(1:lim) = pleg_cache(idx)%beta(1:lim)
    gamma(1:lim) = pleg_cache(idx)%gamma(1:lim)
    coefa(1:lim) = pleg_cache(idx)%coefa(1:lim)
    coefb(1:lim) = pleg_cache(idx)%coefb(1:lim)
    coefc(1:lim) = pleg_cache(idx)%coefc(1:lim)
    coefd(1:lim) = pleg_cache(idx)%coefd(1:lim)
    coefe(1:lim) = pleg_cache(idx)%coefe(1:lim)
  end subroutine

  subroutine store_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: m, iopd, ndec, nex, narg, maxt, lim, maxp
    real(knd), intent(in) :: barg(maxt), pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(in) :: ipdnorm(maxt), ipnorm(maxt)
    integer :: idx

    ! Reuse a matching slot when possible; otherwise take a free or round-robin
    ! replacement slot.
    idx = find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    if(idx == 0) then
      idx = 1
      do while(idx <= pleg_cache_slots)
        if(.not. pleg_cache(idx)%valid) exit
        idx = idx + 1
      end do
      if(idx > pleg_cache_slots) then
        idx = pleg_cache_next
        pleg_cache_next = pleg_cache_next + 1
        if(pleg_cache_next > pleg_cache_slots) pleg_cache_next = 1
      end if
    end if

    call clear_pleg_cache_slot(idx)
    allocate(pleg_cache(idx)%barg(narg))
    allocate(pleg_cache(idx)%pr(maxt, maxp))
    allocate(pleg_cache(idx)%pdr(maxt, maxp))
    allocate(pleg_cache(idx)%pdnorm(maxt))
    allocate(pleg_cache(idx)%pnorm(maxt))
    allocate(pleg_cache(idx)%ipdnorm(maxt))
    allocate(pleg_cache(idx)%ipnorm(maxt))
    allocate(pleg_cache(idx)%alpha(maxp))
    allocate(pleg_cache(idx)%beta(maxp))
    allocate(pleg_cache(idx)%gamma(maxp))
    allocate(pleg_cache(idx)%coefa(maxp))
    allocate(pleg_cache(idx)%coefb(maxp))
    allocate(pleg_cache(idx)%coefc(maxp))
    allocate(pleg_cache(idx)%coefd(maxp))
    allocate(pleg_cache(idx)%coefe(maxp))

    pleg_cache(idx)%valid = .true.
    pleg_cache(idx)%m = m
    pleg_cache(idx)%iopd = iopd
    pleg_cache(idx)%narg = narg
    pleg_cache(idx)%maxt = maxt
    pleg_cache(idx)%maxp = maxp
    pleg_cache(idx)%lim = lim
    pleg_cache(idx)%ndec = ndec
    pleg_cache(idx)%nex = nex
    pleg_cache(idx)%barg(1:narg) = barg(1:narg)
    pleg_cache(idx)%pr = 0.0_knd
    pleg_cache(idx)%pdr = 0.0_knd
    pleg_cache(idx)%alpha = 0.0_knd
    pleg_cache(idx)%beta = 0.0_knd
    pleg_cache(idx)%gamma = 0.0_knd
    pleg_cache(idx)%coefa = 0.0_knd
    pleg_cache(idx)%coefb = 0.0_knd
    pleg_cache(idx)%coefc = 0.0_knd
    pleg_cache(idx)%coefd = 0.0_knd
    pleg_cache(idx)%coefe = 0.0_knd
    pleg_cache(idx)%pr(1:maxt, 1:lim) = pr(1:maxt, 1:lim)
    pleg_cache(idx)%pdr(1:maxt, 1:lim) = pdr(1:maxt, 1:lim)
    pleg_cache(idx)%pdnorm(1:maxt) = pdnorm(1:maxt)
    pleg_cache(idx)%pnorm(1:maxt) = pnorm(1:maxt)
    pleg_cache(idx)%ipdnorm(1:maxt) = ipdnorm(1:maxt)
    pleg_cache(idx)%ipnorm(1:maxt) = ipnorm(1:maxt)
    pleg_cache(idx)%alpha(1:lim) = alpha(1:lim)
    pleg_cache(idx)%beta(1:lim) = beta(1:lim)
    pleg_cache(idx)%gamma(1:lim) = gamma(1:lim)
    pleg_cache(idx)%coefa(1:lim) = coefa(1:lim)
    pleg_cache(idx)%coefb(1:lim) = coefb(1:lim)
    pleg_cache(idx)%coefc(1:lim) = coefc(1:lim)
    pleg_cache(idx)%coefd(1:lim) = coefd(1:lim)
    pleg_cache(idx)%coefe(1:lim) = coefe(1:lim)
  end subroutine

  subroutine pleg_cached(m, lim, maxp, limcsav, iopd, ndec, nex, barg, narg, maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: m, lim, maxp, iopd, ndec, nex, narg, maxt
    integer, intent(inout) :: limcsav
    real(knd), intent(inout) :: barg(maxt)
    real(knd), intent(out) :: pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(out) :: ipdnorm(maxt), ipnorm(maxt)
    integer :: idx

    ! Cached wrapper around pleg: load when available, otherwise compute once
    ! and store for later calls at the same settings.
    idx = find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    if(idx /= 0) then
      call load_pleg_cache(idx, maxt, narg, lim, maxp, barg, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
      limcsav = max(limcsav, lim)
      return
    end if

    call pleg(m, lim, maxp, limcsav, iopd, ndec, nex, barg, narg, maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    call store_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
  end subroutine

  subroutine clear_qleg_cache_slot(idx)
    integer, intent(in) :: idx

    qleg_cache(idx)%valid = .false.
    qleg_cache(idx)%m = -1
    qleg_cache(idx)%lnum = 0
    qleg_cache(idx)%limq = 0
    qleg_cache(idx)%maxq = 0
    qleg_cache(idx)%ndec = 0
    qleg_cache(idx)%iqdml = 0
    qleg_cache(idx)%iqml = 0
    qleg_cache(idx)%itermpq = 0
    qleg_cache(idx)%x1 = 0.0_knd
    qleg_cache(idx)%qdml = 0.0_knd
    qleg_cache(idx)%qml = 0.0_knd
    qleg_cache(idx)%termpq = 0.0_knd
    if(allocated(qleg_cache(idx)%qdr)) deallocate(qleg_cache(idx)%qdr)
    if(allocated(qleg_cache(idx)%qr)) deallocate(qleg_cache(idx)%qr)
    if(allocated(qleg_cache(idx)%qdl)) deallocate(qleg_cache(idx)%qdl)
    if(allocated(qleg_cache(idx)%ql)) deallocate(qleg_cache(idx)%ql)
    if(allocated(qleg_cache(idx)%iqdl)) deallocate(qleg_cache(idx)%iqdl)
    if(allocated(qleg_cache(idx)%iql)) deallocate(qleg_cache(idx)%iql)
  end subroutine

  integer function find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    integer, intent(in) :: m, lnum, limq, maxq, ndec
    real(knd), intent(in) :: x1
    integer :: i

    find_qleg_cache = 0
    do i = 1, qleg_cache_slots
      if(.not. qleg_cache(i)%valid) cycle
      if(qleg_cache(i)%m /= m) cycle
      if(qleg_cache(i)%ndec /= ndec) cycle
      if(.not. cache_real_equal(qleg_cache(i)%x1, x1)) cycle
      if(qleg_cache(i)%lnum < lnum) cycle
      if(qleg_cache(i)%limq < limq) cycle
      if(qleg_cache(i)%maxq < maxq) cycle
      find_qleg_cache = i
      return
    end do
  end function

  subroutine store_qleg_cache(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    integer, intent(in) :: m, lnum, limq, maxq, ndec, iqdml, iqml, itermpq
    real(knd), intent(in) :: x1, qdr(maxq), qdml, qdl(lnum), qr(maxq), qml, ql(lnum), termpq
    integer, intent(in) :: iqdl(lnum), iql(lnum)
    integer :: idx

    idx = find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    if(idx == 0) then
      idx = 1
      do while(idx <= qleg_cache_slots)
        if(.not. qleg_cache(idx)%valid) exit
        idx = idx + 1
      end do
      if(idx > qleg_cache_slots) then
        idx = qleg_cache_next
        qleg_cache_next = qleg_cache_next + 1
        if(qleg_cache_next > qleg_cache_slots) qleg_cache_next = 1
      end if
    end if

    call clear_qleg_cache_slot(idx)
    allocate(qleg_cache(idx)%qdr(maxq))
    allocate(qleg_cache(idx)%qr(maxq))
    allocate(qleg_cache(idx)%qdl(lnum))
    allocate(qleg_cache(idx)%ql(lnum))
    allocate(qleg_cache(idx)%iqdl(lnum))
    allocate(qleg_cache(idx)%iql(lnum))

    qleg_cache(idx)%valid = .true.
    qleg_cache(idx)%m = m
    qleg_cache(idx)%lnum = lnum
    qleg_cache(idx)%limq = limq
    qleg_cache(idx)%maxq = maxq
    qleg_cache(idx)%ndec = ndec
    qleg_cache(idx)%x1 = x1
    qleg_cache(idx)%qdml = qdml
    qleg_cache(idx)%iqdml = iqdml
    qleg_cache(idx)%qml = qml
    qleg_cache(idx)%iqml = iqml
    qleg_cache(idx)%termpq = termpq
    qleg_cache(idx)%itermpq = itermpq
    qleg_cache(idx)%qdr(1:maxq) = qdr(1:maxq)
    qleg_cache(idx)%qr(1:maxq) = qr(1:maxq)
    qleg_cache(idx)%qdl(1:lnum) = qdl(1:lnum)
    qleg_cache(idx)%ql(1:lnum) = ql(1:lnum)
    qleg_cache(idx)%iqdl(1:lnum) = iqdl(1:lnum)
    qleg_cache(idx)%iql(1:lnum) = iql(1:lnum)
  end subroutine

  subroutine qleg_cached(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    integer, intent(in) :: m, lnum, limq, maxq, ndec
    real(knd), intent(in) :: x1
    real(knd), intent(out) :: qdr(maxq), qdml, qdl(lnum), qr(maxq), qml, ql(lnum), termpq
    integer, intent(out) :: iqdml, iqdl(lnum), iqml, iql(lnum), itermpq
    integer :: idx

    idx = find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    if(idx /= 0) then
      qdr(1:maxq) = qleg_cache(idx)%qdr(1:maxq)
      qr(1:maxq) = qleg_cache(idx)%qr(1:maxq)
      qdl(1:lnum) = qleg_cache(idx)%qdl(1:lnum)
      ql(1:lnum) = qleg_cache(idx)%ql(1:lnum)
      iqdl(1:lnum) = qleg_cache(idx)%iqdl(1:lnum)
      iql(1:lnum) = qleg_cache(idx)%iql(1:lnum)
      qdml = qleg_cache(idx)%qdml
      iqdml = qleg_cache(idx)%iqdml
      qml = qleg_cache(idx)%qml
      iqml = qleg_cache(idx)%iqml
      termpq = qleg_cache(idx)%termpq
      itermpq = qleg_cache(idx)%itermpq
      return
    end if

    call qleg(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    call store_qleg_cache(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
  end subroutine

  subroutine gauss_cached(ndec, n, x, w)
    integer, intent(in) :: ndec, n
    real(knd), intent(out) :: x(n), w(n)

    if(gauss_cache%valid) then
      if(gauss_cache%ndec == ndec .and. gauss_cache%n == n) then
        x(1:n) = gauss_cache%x(1:n)
        w(1:n) = gauss_cache%w(1:n)
        return
      end if
      if(allocated(gauss_cache%x)) deallocate(gauss_cache%x)
      if(allocated(gauss_cache%w)) deallocate(gauss_cache%w)
    end if

    call gauss(ndec, n, x, w)
    allocate(gauss_cache%x(n))
    allocate(gauss_cache%w(n))
    gauss_cache%valid = .true.
    gauss_cache%ndec = ndec
    gauss_cache%n = n
    gauss_cache%x(1:n) = x(1:n)
    gauss_cache%w(1:n) = w(1:n)
  end subroutine

  subroutine profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
           r1c, ir1e, r1dc, ir1de, r2c, &
           ir2e, r2dc, ir2de, naccr, &
           s1c, is1e, s1dc, is1de, naccs)

!      version March 2026
!
!  Subroutine version of the fortran program profcn originally developed
!  about 2000 by arnie lee van buren and jeffrey boisvert. Updated
!  several times since then. For more information see the GitHub
!  repository: GitHub.com/MathieuandSpheroidalWaveFunctions/Prolate_swf.
!  Especially see the readme file, example input and output files and two
!  journal articles describing the methods used in profcn. These subroutines 
!  were modified by Brandyn M. Lucca to improve computational efficiency and 
!  compatibility with the psms.cpp C++ interface used by the acousticTS 
!  R-package. The scaffolding for the actual calculations and numerical 
!  methods otherwise remain mostly the same as version 1.15 (Dec 2023) of 
!  this file.
!
!  purpose:     To calculate the first and second kind prolate
!               radial functions r1 and r2 and their first
!               derivatives r1d and r2d for a given order m,
!               a range of degrees l beginning at m, and for a
!               specific size parameter c and shape parameter x.
!               [Note that profcn inputs x1 = x - 1]
!               To calculate the first kind prolate angular
!               functions and their first derivatives with
!               respect to eta for a range of values for l
!               and eta for specified values of c and m.
!
!  Profcn can be run in either double precision or quadruple precision
!  arithmetic. The choice is set in the module param provided in the github
!  repository. If this is not available, then create param as follows:
!    module param
!    integer, parameter :: knd = selected_real_kind(8)
!    logical, parameter :: debug = .true.
!    logical, parameter :: warn = .true.
!    logical, parameter :: output = .true.
!    end module param
!  Set the value of knd in the parenthesis to either 8 for double
!  precision or 16 for quadruple precision arithmetic. The logicals
!  in param are described in the readme file and below in the discussion
!  of the output files.
!
!  Profcn provides accurate results over very wide parameter ranges when
!  using double precision. It provides higher accuracy using quadruple
!  precision but run times are considerable greater.
!
!    Input and output parameters
!
!          c      : desired value of the size parameter (= kd/2, where
!                   k = wavenumber and d = interfocal length) (either
!                   real*8 or real*16)
!          m      : desired value for the order m (integer)
!          lnum   : number of values desired for the degree l equal
!                   to m, m + 1, m + 2, ..., m + lnum - 1 (integer)
!          ioprad : (integer)
!                 : =0 if radial functions are not computed
!                 : =1 if radial functions of only the first kind
!                      and their first derivatives are computed
!                 : =2 if radial functions of both kinds and
!                      their first derivatives are computed
!          x1     : value of the radial coordinate x minus 1.0. This
!                   choice is made to avoid subtraction errors in
!                   calculating quantities containing x - 1 when x
!                   is close to unity. Note that when x1 = 0.0,
!                   ioprad must be set = 1, since r2 and r2d are
!                   infinite in value in this case.
!                   (a nominal value can be entered for x1 if ioprad
!                   = 0)
!          iopang : (integer)
!                 : =0 if angular functions are not computed
!                 : =1 if angular functions of the first kind
!                      are computed
!                 : =2 if angular functions of the first kind and
!                      their first derivatives are computed
!          iopnorm: (integer)
!                 : =0 if not scaled. The angular functions have
!                      the same norm as the corresponding associated
!                      legendre function [i.e., we use the Meixner and
!                      Schafke normalization scheme.] This norm
!                      becomes very large as m becomes large. The
!                      angular functions are computed below as
!                      a characteristic and an exponent to avoid
!                      overflow.
!                 : =1 if angular functions of the first kind
!                      (and their first derivatives if computed)
!                      are scaled by the square root of the
!                      normalization of the corresponding
!                      associated Legendre function. The resulting
!                      scaled angular functions have unity norm.
!                      This is very useful since it removes the
!                      need to calculate a normalization factor
!                      when using the angular function values given
!                      here. It also eliminates any chance for
!                      overflow when the characteristics and exponents
!                      are combined to form the angular functions.
!          narg   : number of values of the angular coordinate eta for
!                   which angular functions are calculated (integer)
!          arg:     vector containing the values of eta for which
!                   angular functions are desired (real*8 or real*16)
!          r1c   :  either real*8 or real*16 vectors of length lnum
!          r1dc     containing the characteristics for the radial
!                   radial functions of the first kind r1 and their
!                   first derivatives
!          ir1e   : integer vectors of length lnum containing the
!          ir1de    exponents corresponding to r1c and r1dc
!          r2c    : real*8 or real*16 vectors of length lnum containing
!          r2dc     the characteristics for the radial functions of the
!                   second kind r2 and their first derivatives
!          ir2e   : integer vectors of length lnum containing the
!          ir2de    exponents corresponding to r2c and r2dc
!          naccr  : integer vector of length lnum containing the estimated
!                   accuracy of the radial functions
!          s1c,   : two-dimensional arrays s1c(lnum,narg) and
!          s1dc     s1dc(lnum,narg) that contain narg calculated
!                   characteristics for the angular functions and
!                   their first derivatives for each of the lnum
!                   values of l (real*8 or real*16)
!                   For example, s1(10,1) is the characteristic
!                   of the angular function for l = m +10 -1 and
!                   for the first value of eta given by arg(1)
!          is1e,  : integer arrays is1e(lnum,narg) and is1de(lnum,narg)
!          is1de    containing the exponents corresponding to s1c and
!                   s1dc
!          naccs  : two-dimensional array naccs(lnum,narg) containing
!                   narg estimated accuracy values for the angular functions
!                   for each of the lnum values of l
!
!  Profcn offers several several output files: Fort.20 and fort.30
!  list the calculated radial and angular functions. Fort.40 and
!  fort.50 are diagnostic files. Fort.60 provides warning whenever the
!  estimated accuracy falls below a specified minimum, currently set
!  equal to 6. Writing to these files is controlled by logicals specified
!  above in the module param. False suppresses the file; true enables it.
!  Debug controls fort.30 and fort.40, warn controls fort.60 and output
!  controls fort.20 and fort.30. Information about these files as well
!  as a discussion about accuracy, expansion d coefficients and eigenvalues
!  is given in the readme file.

    real(knd), intent (in) :: c, x1, arg(narg)
    integer, intent (in)  :: m, lnum, ioprad, iopang, iopnorm, narg
    real(knd), intent (out) :: r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                  s1c(lnum, narg), s1dc(lnum, narg)
    integer, intent (out)  :: ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                  is1e(lnum, narg), is1de(lnum, narg), naccr(lnum), &
                  naccs(lnum, narg)
    real(knd) :: qr1_tmp(lnum, 1), qr1d_tmp(lnum, 1), qr2_tmp(lnum, 1), qr2d_tmp(lnum, 1), &
                 s1_tmp(lnum, narg, 1), s1d_tmp(lnum, narg, 1)
    integer :: ir1_tmp(lnum, 1), ir1d_tmp(lnum, 1), ir2_tmp(lnum, 1), ir2d_tmp(lnum, 1), &
               nar_tmp(lnum, 1), is1_tmp(lnum, narg, 1), is1d_tmp(lnum, narg, 1), &
               naccs_tmp(lnum, narg, 1)

!       Here is where the user sets kindd, the value for kind that
!       corresponds to double precision data on the users computer.
!       Similarly, this is where kindq, the value of kind for quadruple
!       precision data, is set. These values are set below to 8 and 16,
!       respectively. They should be changed to the kind values for double
!       precision and quadruple precision if those values are different than
!       these.

5    kindd = 8
    kindq = 16

!       set the minimum desired accuray minacc to 8 for real*8
!       arithmetic and to 15 for real*16 arithmetic. These can be
!       changed if desired. See the readme file.

    if(knd == kindd) minacc = 8
    if(knd == kindq) minacc = 8

!       ndec: the maximum number of decimal digits available in real(knd)
!             arithmetic.
!       nex:  the maximum exponent available in real(knd) arithmetic.

    ndec = precision(c)
    nex = range(c) - 1

!       open input and output files
!    if (output) then
!     open(20, file='fort.20')
!      open(30, file='fort.30')
!    end if
!    if (debug) then
!      open(40, file='fort.40')
!      open(50, file='fort.50')
!    end if
!    if (warn) then
!      open(60, file='fort.60')
!    end if
!
!       set array dimensions
    mnum = 1
    mmin = m
    minc = 0
    maxm = mmin

    maxint = lnum + 3 * ndec + int(c) + 5
    maxj = maxint + maxm
    maxp = maxint
    maxn = maxp + maxm
    maxpdr = 4 * ndec + 5
    neta = 993
    ngau = 200
    if(ioprad == 2) then
      lnump = max(lnum + maxm, 64)
      if(x1 >= 0.00065e0_knd) maxn = 2 * (lnump * (-18.5e0_knd - 20.e0_knd * &
              log10(x1)) + 5 * ndec + 4 * maxm + c + 5000) + maxm + 5

      if(x1 > 0.08e0_knd) maxn = 2 * (lnump * (0.5e0_knd - 3.0e0_knd * &
              log10(x1)) + 5 * ndec + 4 * maxm + c + 1000) + maxm + 5

      if(x1 > 1.0e0_knd) maxn = 2 * (lnump * 0.5e0_knd + 5 * ndec &
              + 4 * maxm + c + 500) + maxm + 5
      maxp = max(maxn, maxp)
      if(x1 < 1.0e-3_knd) ngau = 200 - 50 * int(log10(x1) - 1.0e-30_knd)
      if(x1 < 1.0e-10_knd) ngau = 250 - 50 * int(log10(x1) - 1.0e-30_knd)
      if(x1 < 1.0e-11_knd) ngau =1200
      if(x1 < 1.0e-12_knd) ngau =2400
      if(x1 <= 0.5e0_knd) maxpdr = maxpdr + int(2.e0_knd * c + 100.0e0_knd * x1) + 400
    end if
    maxq = maxint + maxm + maxm
    maxdr = maxpdr / 2 + 1
    maxp = max(maxp, maxpdr)
    maxd = maxp / 2 + 1
    maxlp = lnum + maxm + 5
    maxmp = maxm +5
    maxt = 1
    jnenmax = 10
    if(iopang /= 0) maxt = narg

     call main (mmin, minc, mnum, lnum, c, ioprad, iopang, iopnorm, minacc, &
          x1, ngau, arg, narg, neta, maxd, maxdr, maxint, maxj, maxlp, &
          maxm, maxmp, maxn, maxp, maxpdr, maxq, maxt, jnenmax, &
          kindd, kindq, ndec, nex, qr1_tmp, ir1_tmp, qr1d_tmp, ir1d_tmp, &
          qr2_tmp, ir2_tmp, qr2d_tmp, ir2d_tmp, nar_tmp, s1_tmp, is1_tmp, &
          s1d_tmp, is1d_tmp, naccs_tmp)

    r1c = qr1_tmp(:, 1)
    ir1e = ir1_tmp(:, 1)
    r1dc = qr1d_tmp(:, 1)
    ir1de = ir1d_tmp(:, 1)
    r2c = qr2_tmp(:, 1)
    ir2e = ir2_tmp(:, 1)
    r2dc = qr2d_tmp(:, 1)
    ir2de = ir2d_tmp(:, 1)
    naccr = nar_tmp(:, 1)
    s1c = s1_tmp(:, :, 1)
    is1e = is1_tmp(:, :, 1)
    s1dc = s1d_tmp(:, :, 1)
    is1de = is1d_tmp(:, :, 1)
    naccs = naccs_tmp(:, :, 1)

    end subroutine

  subroutine profcn_batch(c, m_start, m_count, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
           r1c, ir1e, r1dc, ir1de, r2c, &
           ir2e, r2dc, ir2de, naccr, &
           s1c, is1e, s1dc, is1de, naccs)

!  purpose:     Batched wrapper around profcn for consecutive orders
!               m = m_start, ..., m_start + m_count - 1.
!               The underlying mathematics is unchanged; the benefit is that
!               setup inside main can be amortized across several nearby orders.

    real(knd), intent (in) :: c, x1, arg(narg)
    integer, intent (in)  :: m_start, m_count, lnum, ioprad, iopang, iopnorm, narg
    real(knd), intent (out) :: r1c(lnum, m_count), r1dc(lnum, m_count), &
                  r2c(lnum, m_count), r2dc(lnum, m_count), &
                  s1c(lnum, narg, m_count), s1dc(lnum, narg, m_count)
    integer, intent (out)  :: ir1e(lnum, m_count), ir1de(lnum, m_count), &
                  ir2e(lnum, m_count), ir2de(lnum, m_count), naccr(lnum, m_count), &
                  is1e(lnum, narg, m_count), is1de(lnum, narg, m_count), &
                  naccs(lnum, narg, m_count)
5    kindd = 8
    kindq = 16

    if(knd == kindd) minacc = 8
    if(knd == kindq) minacc = 8

    ndec = precision(c)
    nex = range(c) - 1

    ! Request a consecutive block of orders from main.
    mnum = m_count
    mmin = m_start
    minc = 1
    maxm = mmin + minc * (mnum - 1)

    maxint = lnum + 3 * ndec + int(c) + 5
    maxj = maxint + maxm
    maxp = maxint
    maxn = maxp + maxm
    maxpdr = 4 * ndec + 5
    neta = 993
    ngau = 200
    if(ioprad == 2) then
      lnump = max(lnum + maxm, 64)
      if(x1 >= 0.00065e0_knd) maxn = 2 * (lnump * (-18.5e0_knd - 20.e0_knd * &
              log10(x1)) + 5 * ndec + 4 * maxm + c + 5000) + maxm + 5

      if(x1 > 0.08e0_knd) maxn = 2 * (lnump * (0.5e0_knd - 3.0e0_knd * &
              log10(x1)) + 5 * ndec + 4 * maxm + c + 1000) + maxm + 5

      if(x1 > 1.0e0_knd) maxn = 2 * (lnump * 0.5e0_knd + 5 * ndec &
              + 4 * maxm + c + 500) + maxm + 5
      maxp = max(maxn, maxp)
      if(x1 < 1.0e-3_knd) ngau = 200 - 50 * int(log10(x1) - 1.0e-30_knd)
      if(x1 < 1.0e-10_knd) ngau = 250 - 50 * int(log10(x1) - 1.0e-30_knd)
      if(x1 < 1.0e-11_knd) ngau =1200
      if(x1 < 1.0e-12_knd) ngau =2400
      if(x1 <= 0.5e0_knd) maxpdr = maxpdr + int(2.e0_knd * c + 100.0e0_knd * x1) + 400
    end if
    maxq = maxint + maxm + maxm
    maxdr = maxpdr / 2 + 1
    maxp = max(maxp, maxpdr)
    maxd = maxp / 2 + 1
    maxlp = lnum + maxm + 5
    maxmp = maxm +5
    maxt = 1
    jnenmax = 10
    if(iopang /= 0) maxt = narg

     call main (mmin, minc, mnum, lnum, c, ioprad, iopang, iopnorm, minacc, &
          x1, ngau, arg, narg, neta, maxd, maxdr, maxint, maxj, maxlp, &
          maxm, maxmp, maxn, maxp, maxpdr, maxq, maxt, jnenmax, &
          kindd, kindq, ndec, nex, r1c, ir1e, r1dc, ir1de, r2c, ir2e, &
          r2dc, ir2de, naccr, s1c, is1e, s1dc, is1de, naccs)

    end subroutine

    subroutine main (mmin, minc, mnum, lnum, c, ioprad, iopang, iopnorm, &
             minacc, x1, ngau, barg, narg, neta, maxd, maxdr, &
             maxint, maxj, maxlp, maxm, maxmp, maxn, maxp, maxpdr, &
             maxq, maxt, jnenmax, kindd, kindq, ndec, nex, qr1, ir1, &
             qr1d, ir1d, qr2, ir2, qr2d, ir2d, nar, s1, is1, s1d, is1d, &
             nas)

!  purpose:     To coordinate the calculation of both the prolate
!               spheroidal radial and angular functions and their
!               first derivatives using various algorithms.
!
!  parameters:
!
!     input:    mmin   : minimum desired value of m
!               minc   : increment in m used to compute other values
!               mnum   : number of values of m that are desired
!               lnum   : desired number of values of l = m, m + 1, ...,
!                        m + lnum - 1
!               c      : size parameter
!               ioprad : equal to 0 if no radial functions are desired;
!                        equal to 1 if only radial functions of the
!                          first kind and their first derivatives are
!                          desired;
!                        equal to 2 if radial functions of both kinds
!                          and their first derivatives are desired
!               iopang : equal to 0 if no angular functions are desired;
!                        equal to 1 if only angular functions of the
!                          first kind are desired;
!                        equal to 2 if angular functions of the first
!                          kind and their first derivatives are desired
!               iopnorm: equal to 0 when the angular functions have
!                        the same norm as the corresponding associated
!                        Legendre functions;
!                        equal to 1 when the angular functions are
!                        scaled by the square root of the normalization
!                        of the corresponding Legendre function, giving
!                        them unity norm
!               minacc : desired minimum accuracy for the radial
!                        functions
!               x1     : radial coordinate x minus 1
!               ngau   : order of the Gaussian quadrature to be used in
!                        computing integrals in subroutine pint for use
!                        in subroutine r2int where the integal method
!                        is used to calculate r2 and r2d
!               barg   : array containing the values of eta for which
!                        angular functions are to be computed (named
!                        arg in subroutine calling statement)
!               narg   : number of desired eta values; dimension of
!                        barg.
!               neta   : number of values available for eta in the
!                        variable eta method for calculating r2 and r2d
!                        (subroutine r2eta); set equal to 993 above
!               maxd   : dimension of enr array containing ratios of
!                        the expansion d coefficients
!               maxdr  : dimension of drhor array containing special d
!                        coefficient ratios used in subroutine r2leg
!                        when computing the sum of Legendre functions of
!                        the first kind that appear in the Legendre
!                        function expansion for r2 and r2d
!               maxint : maximum number of terms available for computing
!                        r2 and r2d in the subroutine r2int; dimension
!                        of the arrays of integrals computed in
!                        subroutine pint
!               maxj   : equal to the dimension of the array of ratios
!                        of spherical Bessel functions of the first kind
!                        and the array of ratios of the first derivative
!                        of this Bessel function to the corresponding
!                        Bessel function
!               maxlp  : maximum value desired for l
!               maxm   : maximum value desired for m
!               maxmp  : maxm + 5; dimension of the integer array norme
!                        used in scaling of the Neumann functions in
!                        the integrands in subroutine pint
!               maxn   : dimension of the arrays of Neumann function
!                        ratios used in computing r2 and r2d
!               maxp   : dimension of arrays of Legendre functions of
!                        the first kind used in computing angular
!                        functions, in computing integrands in
!                        subroutine pint and in computing r2 and r2d in
!                        subroutine r2eta
!               maxpdr : dimension of the arrays of ratios of both
!                        Legendre functions of the first kind and their
!                        first derivatives used in the sum of these
!                        functions that contribute to r2 and r2d in
!                        subroutine r2leg
!               maxq   : dimension of arrays of ratios of Legendre
!                        functions of the second kind and ratios of
!                        their first derivatives used in their sum in
!                        subroutine r2leg
!               maxt   : equal to narg if angular functions are
!                        computed where it is the maximum value of the
!                        first index in the arrays of Legendre functions
!                        used in subroutine s1leg;
!                        otherwise equal to 1 to specify the
!                        first index for the Legendre functions used
!                        in the variable eta method for computing r2
!                        and r2d in subroutine r2eta
!               jneumax: number of arrays of ratios of Legendre and
!                        Neumann functions stored as eta is varied in
!                        subroutine r2eta; set equal to 10 so that the
!                        previous 10 sets of ratios are available
!                        to use without recalculating them when one of
!                        thes previous values for eta is used again for
!                        a later value of l
!               kindd  : kind value for double precision real data
!               kindq  : kind value for quadruple precision real data
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent for real(knd)
!
!     output:   qr1    : array of lnum values for the characteristic
!                        of the radial function of the first kind
!                        (r1c in call to subroutine profcn)
!               ir1    : array of exponents corresponding to qr1
!                        (ir1e in call to subroutine profcn)
!               qr1d   : array of lnum values for the characteristic
!                        of the derivative of the radial function of
!                        the first kind (r1dc in call to subroutine
!                        profcn)
!               ir1d   : array of exponents corresponding to qr1d
!                        (ir1de in call to subroutine profcn)
!               qr2    : array of lnum values for the characteristic
!                        of the radial function of the second kind
!                        (r2c in call to subroutine profcn)
!               ir2    : array of exponents corresponding to qr2
!                        (ir2e in call to subroutine profcn)
!               qr2d   : array of lnum values for the characteristic
!                        of the derivative of the radial function of
!                        the second kind (r2dc in call to subroutine
!                        profcn)
!               ir2d   : array of exponents corresponding to qr1d
!                        (ir2de in call to subroutine profcn)
!               nar    : vector of lnum values for the estimated radial
!                        function accuracy
!               s1     : two dimensional array of the characteristics
!                        of the angular functions of the first kind
!                        s1(i,j) is the characteristic for the jth
!                        value of eta and the degree = m + i -1
!                        (s1c in call to subroutine profcn)
!               is1    : array of exponents corresponding to s1 (is1e)
!               s1d    : two dimensional array of the characteristics
!                        of the first derivatives of the angular
!                        functions of the first kind; s1d(i,j) is the
!                        characteristic for the jth value of eta and
!                        the degree = m + i -1 (s1dc in call to
!                        subroutine profcn)
!               is1d   : array of exponents corresponding to s1d
!                       (is1de in call to subroutine profcn)
!               nas    : two dimensional array nas(lnum,narg) of the
!                        estimated accuracy values for the narg angular
!                        function values for each of the lnum values of l
!
!    use param
!
!  real(knd) scalars
    real(knd) aj1, aj2, ang, apcoef, apcoefn, c, c2, c4, coefn, &
         coefme, coefmo, coefr1e, coefr1o, dec, dfnorm, dmfnorm, &
         dmsnorm, dmlf, dmlmf, dmlms, dmlms1, dneg, d01, dold, dnew, &
         eigval, eigvalp, eig1, eig2, eig3, eig4, eig5, etaval, factor, &
         pcoefe, pcoefet, pcoefn, pcoefo, pdcoefe, pdcoefet, pdcoefo, &
         pi, qdml, qml, rl, rm, rm2, r1c, r1dc, r2c, r2dc, r2ec, r2dec, &
         r2ic, r2dic, r2lc, r2dlc, r2nc, r2dnc, sgn, termpq, x, xb, &
         xbninp, x1, wm, wronc, wront, wronca, wroncb
    character chr
!
!  integer and real(knd) arrays with dimension lnum
    dimension iqdl(lnum), iql(lnum), ifajo(lnum)
    real(knd) qdl(lnum), ql(lnum), fajo(lnum), eig(lnum)
    dimension ir1(lnum, mnum), ir1d(lnum, mnum), ir2(lnum, mnum), ir2d(lnum, mnum), &
             nar(lnum, mnum)
    real(knd) qr1(lnum, mnum), qr1d(lnum, mnum), qr2(lnum, mnum), qr2d(lnum, mnum)
!
!  real(knd) and integer arrays with dimensions lnum and narg
    real(knd) s1(lnum, narg, mnum), s1d(lnum, narg, mnum)
    dimension is1(lnum, narg, mnum), is1d(lnum, narg, mnum), nas(lnum, narg, mnum)
!
!  real(knd) arrays with dimension maxd
    real(knd) enr(maxd), bliste(maxd), gliste(maxd), &
         blisto(maxd), glisto(maxd)
!
!  real(knd) arrays with dimension maxdr
    real(knd) drhor(maxdr)
!
!  real(knd) arrays with dimension maxint
    real(knd) pint1(maxint), pint2(maxint), pint3(maxint), &
         pint4(maxint), rpint1(maxint), rpint2(maxint)
!
!  real(knd) array with dimension maxj
    real(knd) sbesf(maxj), sbesdf(maxj)
!
!  integer and real(knd) arrays with dimension maxlp
    dimension ibese(maxlp), ineue(maxlp), ineuee(maxlp), &
         ipnormint(maxlp), ineuesv(jnenmax, maxlp)
    real(knd) pnormint(maxlp), sbesdr(maxlp), sbesn(maxlp), &
         sneun(maxlp), sneune(maxlp), sneudr(maxlp), &
         sneudre(maxlp), sneunsv(jnenmax, maxlp), &
         sneudrsv(jnenmax, maxlp)
!
!  integers and real(knd) arrays with dimension maxmp
    real(knd) enrneg(maxmp)
    dimension norme(maxmp)
!
!  real(knd) arrays with dimension maxn
    real(knd) sneuf(maxn), sneudf(maxn), sneufe(maxn), sneudfe(maxn), &
         sneufsv(jnenmax, maxn), sneudfsv(jnenmax, maxn)
!
!  real(knd) arrays with dimension given by maxp
    real(knd) alpha(maxp), beta(maxp), coefa(maxp), coefb(maxp), &
         coefc(maxp), coefd(maxp), coefe(maxp), gamma(maxp), &
         pdr(maxt, maxp), pdrat(maxt, maxp), pdratt(maxp), &
         pr(maxt, maxp), prat(maxt, maxp), pratb(maxp), pratt(maxp), &
         prat1(maxp), pratbsv(jnenmax, maxp), &
         prattsv(jnenmax, maxp), pdratsv(jnenmax, maxp)
!
!  real(knd) arrays with dimension maxpdr
    real(knd) prx(maxpdr), pdrx(maxpdr)
!
!  real(knd) arrays with dimension maxq
    real(knd) qr(maxq), qdr(maxq)
!
!  real(knd) and integer arrays with dimension maxt
    real(knd) barg(maxt), etainp(maxt), pdnorm(maxt), pdnorma(maxt), &
         pnorm(maxt), pnorma(maxt), pdtempe(maxt), pdtempo(maxt), &
         ptempe(maxt), ptempo(maxt), s1c(maxt), s1dc(maxt), &
         xin(maxt), xlninp(maxt)
    dimension ipdnorm(maxt), ipdnorma(maxt), ipnorm(maxt), &
         ipnorma(maxt), ipdtempe(maxt), ipdtempo(maxt), &
         iptempe(maxt), iptempo(maxt), is1e(maxt), is1de(maxt), &
         naccs(maxt)
!
!  real(knd) arrays with dimension neta
    real(knd) eta(neta), wmeta2(neta), xbn(neta), xln(neta)
!
!  real(knd) arrays with dimension ngau
    real(knd) wr(ngau), xr(ngau)
!
!  miscellaneous integer arrays
    dimension nees(100), naccsav(100), neeb(jnenmax), limpsv(jnenmax), &
         limnsv(jnenmax), jelimsv(jnenmax)
!
    dec = 10.0e0_knd ** (-ndec - 1)
    if(ioprad /= 0) x = x1 + 1.0e0_knd
    jtest = ndec - minacc - 2
    pi = acos(-1.0e0_knd)
    c2 = c * c
    c4 = c2 * c2
    nbp = int(2.0e0_knd * c / 3.14e0_knd)
    if(knd == kindd) legtest = 8
    if(knd == kindq) legtest = min(15, minacc)
!
!  begin loops
     igau = 0
if (debug) then
     if(knd == kindd .and. ioprad /= 0) write(40, 25) x, c
25    format(1x,'x = ',e23.14,/,1x,'c = ',e23.14)
     if(knd == kindq .and. ioprad /= 0) write(40, 30) x, c
30    format(1x,'x = ',e39.30,/,1x,'c = ',e39.30)
end if
     nc = int(log10(c))
     if(nc < 0) nc = 0
     if(ioprad == 2) wront = 1.0e0_knd / (c * x1 * (x1 + 2.0e0_knd))
     ibflag1 = 0
      do 900 mi = 1, mnum
      m = mmin + minc * (mi - 1)
      m2 = m + m
if (debug) then
       if(knd == kindd .and. iopang /= 0) write(50, 35) c, m
 35     format(1x,'c = ',e23.14,'; m = ',i5)
       if(knd == kindq .and. iopang /= 0) write(50, 40) c, m
 40     format(1x,'c = ',e39.30,'; m = ',i5)
       if(ioprad /= 0) write(40, 50) m
50      format(1x,'m = ',i5)
end if
if (output) then
       if(knd == kindd .and. iopang /= 0) write(30, 60) c, m
 60     format(1x, e23.14, i5)
       if(knd == kindq .and. iopang /= 0) write(30, 70) c, m
 70     format(1x, e39.30, i5)
end if
      rm = m
      rm2 = m + m
      iopleg = 0
      iopneu = 0
      iopeta = 0
      iopint = 1
      jintm = 0
      iopd = 3
      limcsav = 0
      jjjflag = 0
      if(ioprad /= 2) go to 80
      if(x1 <= 0.4e0_knd .and. c <= 10.0e0_knd) iopleg = 1
      if(x1 > 0.4e0_knd .and. c <= 10.0e0_knd) iopneu = 1
      if(x1 <= 0.4e0_knd .and. c <= 15.0e0_knd .and. minacc <= 16) &
         iopleg = 1
      if(x1 > 0.4e0_knd .and. c <= 20.0e0_knd .and. minacc <= 16) &
         iopneu = 1
      if(iopleg == 1 .or. iopneu == 1) iopint = 0
      ioppsum = 1
      iopqnsum = 1
       if(knd == kindd) then
       neest = 897
       if(x1 > 0.01e0_knd) neest = 769
       if(x1 > 0.03e0_knd) neest = 705
       if(x1 > 0.04e0_knd) neest = 641
       if(x1 > 0.05e0_knd) neest = 577
       if(x1 > 0.06e0_knd) neest = 513
       if(x1 > 0.07e0_knd) neest = 449
       if(x1 > 0.08e0_knd) neest = 385
       if(x1 > 0.09e0_knd) neest = 1
       end if
       if(knd == kindq) then
       neest = 905
       if(x1 > 0.01e0_knd) neest = 1
       end if
      nee = neest
      jnen = 0
      incnee = 64
      if(knd == kindd .and. x1 < 0.2e0_knd) incnee = 32
      msearch = 0
80     continue
      if(iopang == 0) go to 90
      limps1 = lnum + 3 * ndec + int(c)
      if((limps1 + 3) > maxp) limps1 = maxp - 3
      iopd = 0
      if(iopang == 2) iopd = 1
      call pleg_cached(m, limps1, maxp, limcsav, iopd, ndec, nex, barg, narg, &
           maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, &
           beta, gamma, coefa, coefb, coefc, coefd, coefe)
      limcsav = limps1
      iopd = 3
90     if(ioprad == 0 .or. mi /= 1 .or. x1 == 0.0e0_knd) go to 100
      limj = lnum + 3 * ndec + int(c) + maxm
      xb = sqrt(x1 * (x1 + 2.0e0_knd))
      call sphbes(c, xb, limj, maxj, maxlp, sbesf, sbesdf, sbesn, ibese, &
            sbesdr)
100     eig1 = 0.0e0_knd
      eig2 = 0.0e0_knd
      eig3 = 0.0e0_knd
      eig4 = 0.0e0_knd
      eig5 = 0.0e0_knd
      iflag = 0
      ibflag2 = 0
      legflag = 0
      jflagleg = 0
      legstart = m
      nflag = 0
      lowacc = ndec
      lowtest = minacc
      nacctest = minacc
      naccintp = minacc + 1
      naccleg = 0
      naccneu = 0
      nacclegp = 0
      naccneup = 0
      nacceta = 0
      naccr = minacc + 1
      ietacount = 0
      incnflag = 0
      iplflag = 0
      factor = 1.0e0_knd
      ijnet = 0
      naccrp = 0
      iflagnee = 0
      istartr2 = 1
      jeta = 0
      iflagq = 0
      iflagp = 0
      jbes = 3 * ndec + int(c)
110     continue
if (output) then
      if(knd == kindd .and. ioprad /= 0) write(20, 115) x, c, m
115     format(1x, e23.14, e23.14, i5)
      if(knd == kindq .and. ioprad /= 0) write(20, 120) x, c, m
120     format(1x, e39.30, e39.30, i5)
end if
       do 850 li = 1, lnum
       l = m + (li - 1)
if (output) then
       if(iopang /= 0) write(30, 140) l
140      format(1x, i6)
end if
if (debug) then
       if(iopang /= 0) write(50, 150) l
150      format(1x,'l = ',i6)
end if
       ix = l - m - 2 * ((l - m) / 2)
       iopnee = 0
        if(iflagnee == 1) then
        incnee = 8
         if(knd == kindd) then
         if(x1 >= 0.05e0_knd) incnee = 16
         if(x1 >= 0.2e0_knd) incnee = 32
         end if
         if(knd == kindq) then
         if(x1 >= 0.05e0_knd) incnee = 16
         if(x1 >= 0.1e0_knd) incnee = 32
         end if
        iflagnee = 2
        end if
       naccetas = minacc
       if(li == 1) naccrsav = minacc
       if(li > 1) naccrsav = naccr
       naccr = -1
       nacce = 0
       limdrad = 3 * ndec + int(c)
       if(ioprad /= 0 .and. li /= 1) limdrad = jbes + jbes + 20+ &
                         int(sqrt(c))
       if(iopint /= 0 .and. li /= 1 .and. jintm > jbes) &
         limdrad = jintm + jintm + 20 + int(sqrt(c))
       limdang = 3 * ndec + int(c)
       if(iopang /= 0 .and. li /= 1) limdang = jang + jang + 20+ &
                         int(sqrt(c))
       if(iopang == 0) limd = limdrad
       if(ioprad == 0) limd = limdang
       if(iopang /= 0 .and. ioprad /= 0) limd = max(limdang, limdrad)
       if(li == 1) limmf = limdang
       if(li > 1) limmf = jmf + jmf + 20 + int(sqrt(c))
       limd = max(limd, limmf)
       if(ioprad /= 2) go to 155
       if(iopleg == 1) limdleg = l - m + 3 * ndec + int(c)
       if(iopleg == 2) limdleg = jleg + jleg + 20 + int(sqrt(c))
       if(iopleg /= 0) limd = max(limd, limdleg)
       limdneu = limd
       lplus = max(l, lnum + maxm)
       if(x1 >= 0.00065e0_knd) limdneu = 2 * ((lplus) * (-18.5e0_knd- &
                    20.0e0_knd * log10(x1))+ &
                    5 * ndec + 4 * m + c + 01000)
       if(x1 > 0.08e0_knd) limdneu = 2 * ((lplus) * (0.5e0_knd- &
                    3.0e0_knd * log10(x1))+ &
                    5 * ndec + 4 * m + c + 01000)
       if(x1 > 1.0e0_knd) limdneu = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                     4 * m + c + 00500)
       if(iopneu == 2 .and. naccneu > 0) &
           limdneu = jneu + jneu + 20 + int(sqrt(c))
       if(iopneu /= 0) limd = max(limd, limdneu)
       limdeta = limd
       if(x1 >= 0.00065e0_knd) limdeta = 2 * ((lplus) * (-18.5e0_knd- &
                      20.0e0_knd * log10(x1)) + 5 * ndec+ &
                      4 * m + c + 05000)
       if(x1 > 0.08e0_knd) limdeta = 2 * ((lplus) * (0.5e0_knd- &
                      3.0e0_knd * log10(x1)) + 5 * ndec+ &
                      4 * m + c + 01000)
       if(x1 > 1.0e0_knd) limdeta = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                     4 * m + c + 00500)
       if(iopeta == 3 .and. naccrsav > minacc) &
               limdeta = jeta + jeta + 500 + c / 10
       if(iopeta == 3 .and. naccrsav <= minacc) &
               limdeta = jeta + jeta + 500 + c
       if(iopeta /= 0) limd = max(limd, limdeta)
155      continue
       if(limd > maxp) limd = maxp
       call geteig (l, m, c, eig2, eig3, eig4, eig5, eigval)
       eig1 = eig2
       eig2 = eig3
       eig3 = eig4
       eig4 = eig5
!
!  use Bouwkamp procedure to obtain accurate eigenvalues
       if(l == m) ienre = (3 * ndec + int(c)) / 2
       if(l == m) jlowe = 1
       if(l == m) limdle = 2
       if(l == m + 1) ienro = (3 * ndec + int(c)) / 2
       if(l == m + 1) jlowo = 1
       if(l == m + 1) limdlo = 3
!
!  compute the coeficients in the Bouwkamp method
       if(ix == 1) go to 160
!
!  beta coefficients (bliste) for l-m even
       if(limdle > limd) go to 163
       j = jlowe
        do 158 i = limdle, limd, 2
        i2 = i + i
        bliste(j) = c4 * real(i, knd) * real((i - 1), knd)* &
             real((m2 + i), knd) * real((m2 + i - 1), knd)/ &
             (real((m2 + i2 - 1), knd) * real((m2 + i2 - 1), knd)* &
             real((m2 + i2 - 3), knd) * real((m2 + i2 + 1), knd))
        j = j + 1
158       continue
!
!  gamma coeficients (gliste) for l-m even
       j = jlowe
        do 159 i = limdle - 1, limd + 1, 2
        i2 = i + i
        gliste(j) = real((m + i - 1), knd) * real((m + i), knd) + 0.5e0_knd* &
             c2 * ((1.0e0_knd - real((m2 * m2 - 1), knd)/ &
             (real((m2 + i2 - 3), knd) * real((m2 + i2 + 1), knd))))
        j = j + 1
159       continue
       go to 163
160      continue
!
!  beta coefficients (blisto) for l-m odd
       if(limdlo > limd) go to 163
       j = jlowo
        do 161 i = limdlo, limd, 2
        i2 = i + i
        blisto(j) = c4 * real(i, knd) * real((i - 1), knd)* &
             real((m2 + i), knd) * real((m2 + i - 1), knd)/ &
             (real((m2 + i2 - 1), knd) * real((m2 + i2 - 1), knd)* &
             real((m2 + i2 - 3), knd) * real((m2 + i2 + 1), knd))
        j = j + 1
161       continue
!
!  gamma coeficient (glisto) for l-m odd
       j = jlowo
        do 162 i = limdlo - 1, limd + 1, 2
        i2 = i + i
        glisto(j) = real((m + i - 1), knd) * real((m + i), knd) + 0.5e0_knd* &
             c2 * (1.0e0_knd - real((m2 * m2 - 1), knd)/ &
             (real((m2 + i2 - 3), knd) * real((m2 + i2 + 1), knd)))
       j = j + 1
162      continue
163      continue
       if(ix == 0) call conver(l, m, c, limd, bliste, gliste, eig1, &
               eig3, eig4, ndec, maxd, eigval, eig5, enr, ienre)
       if(ix == 1) call conver(l, m, c, limd, blisto, glisto, eig1, &
               eig3, eig4, ndec, maxd, eigval, eig5, enr, ienro)
       eig(li) = eigval
if (debug) then
       if(knd == kindd .and. ioprad /= 0) write(40, 165) l, eigval
165      format(1x,'l =',i6, 5x,'eigenvalue =',e23.14)
       if(knd == kindq .and. ioprad /= 0) write(40, 170) l, eigval
170      format(1x,'l =',i6, 5x,'eigenvalue =',e39.30)
end if
       if(ix == 1) go to 175
       limdle = limd + 2
       if(2 * (limd / 2) /= limd) limdle = limd + 1
       jlowe = limd / 2 + 1
       go to 176
175      limdlo = limd + 1
       if(2 * (limd / 2) /= limd) limdlo = limd + 2
       jlowo = (limd - 1) / 2 + 1
176      call dnorm (l, m, c, ndec, nex, limd, maxd, enr, sgn, d01, id01, &
             dmfnorm, idmfe, dmlmf, idmlmfe, dmsnorm, idmse, &
             dmlms, idmlmse, jmf, jsub)
       if(li /= 1 .and. eigval <= eigvalp) go to 900
       if(l == m .and. jsub > jtest .and. iopneu /= 0) iopneu = 0
       if(l == m .and. jsub > jtest .and. iopleg /= 0) iopleg = 0
       eigvalp = eigval
!
!  determine prolate radial functions of the first kind
       qr1(li, mi) = 0.0e0_knd
       qr1d(li, mi) = 0.0e0_knd
       ir1(li, mi) = 0
       ir1d(li, mi) = 0
       qr2(li, mi) = 0.0e0_knd
       qr2d(li, mi) = 0.0e0_knd
       ir2(li, mi) = 0
       ir2d(li, mi) = 0
       if(ioprad == 0) go to 720
if (debug) then
       write(40, 178)
178      format(4x,'r1 and r1d calculation')
end if
!  calculation of r1 and r1d for x = 1
   if(x1 == 0.0e0_knd .and. m == 0) then
!   calculation of dfnorm
!   forward summation of series
    mml = ix - 1
    lm2 = l/2
    if(lm2 == 0) limfl = 1.5 * ndec + int(0.5e0_knd*c)
    dold = 1.0e0_knd
    dfnorm = dold
     do j = lm2 + 1, limfl
     jj = j + j + ix
     dnew = -dold * enr(j) * real((jj + mml), knd) / real(jj - ix, knd)
     dfnorm = dfnorm + dnew
     if(abs(dnew / dfnorm) < dec) exit
     dold = dnew
     jmax=j
     end do
if (debug) then
    write(40, 179) j, limfl
179   format(15x,'Flammer norm. series converged in ',i6,' terms; ', &
        i6,' available.')
end if
    limfl = jmax + 10
! backward summation of series
    if(lm2 >= 1) then
    dold = 1.0e0_knd
     do j = lm2, 1,-1
     jj = j + j + ix
     dnew = -dold * (jj - ix) / (real((jj + mml), knd) &
        *enr(j))
     dfnorm = dfnorm + dnew
     if(abs(dnew / dfnorm) < dec) exit
     dold = dnew
     end do
    end if
    iterm = int(log10(abs(dfnorm)))
    dfnorm = dfnorm * (10.0e0_knd ** (-iterm))
    idfe = iterm
    dmlf = 1.0e0_knd / dfnorm
    idmlfe = -idfe
    if(l == 0) coefr1e = 1.0e0_knd
    if(l == 1) coefr1o = 1.0e0_knd
    rl = real(l, knd)
    if(ix == 0) then
     if(l > 0) coefr1e = coefr1e * rl / (rl - 1.0e0_knd)
     r1c = coefr1e * d01 * dmlf
     iterm = int(log10(abs(r1c)))
     r1c = r1c * (10.0e0_knd ** (-iterm))
     ir1e = id01 + idmlfe + iterm
     r1dc = c * c * coefr1e * d01 * dmlf * ((enr(1) / 15.0e0_knd)-&
       (1.0e0_knd / 3.0e0_knd))
     iterm = int(log10(abs(r1dc)))
     r1dc = r1dc * (10.0e0_knd ** (-iterm))
     ir1de = id01 + idmlfe + iterm
    end if
    if(ix == 1) then
     if(l > 1) coefr1o = coefr1o * (rl - 1.0e0_knd) / rl
     r1c = c * coefr1o * d01 * dmlf / 3.0e0_knd
     iterm = int(log10(abs(r1c)))
     r1c = r1c * (10.0e0_knd ** (-iterm))
     ir1e = id01 + idmlfe + iterm
     r1dc = c * c * c * coefr1o * d01 * dmlf * ((enr(1) / 35.0e0_knd)- &
       (1.0e0_knd / 15.0e0_knd) + (1.0e0_knd / (3.0e0_knd * c * c)))
     iterm = int(log10(abs(r1dc)))
     r1dc = r1dc * (10.0e0_knd ** (-iterm))
     ir1de = id01 + idmlfe + iterm
    end if
    if(abs(r1c) < 1.0e0_knd) then
     r1c = r1c * 10.0e0_knd
     ir1e = ir1e - 1
    end if
    if(abs(r1dc) < 1.0e0_knd) then
     r1dc = r1dc * 10.0e0_knd
     ir1de = ir1de - 1
    end if
   go to 680
   end if
    if(x1 == 0.0e0_knd .and. m /= 0) then
     r1c = 0.0e0_knd
     ir1e = 0
     r1dc = 0.0e0_knd
     ir1de = 0
     go to 680
    end if
! calculation of r1 and r1d for x /= 1
       if(li == 1) limr1 = 3 * ndec + int(c)
       if(li /= 1) limr1 = jbes + jbes + 20 + int(sqrt(c))
       call r1bes(l, m, c, x1, limr1, ndec, maxd, enr, maxj, maxlp, &
             nex, iflag, sbesf, sbesdf, sbesn, ibese, sbesdr, &
             d01, id01, r1c, ir1e, r1dc, ir1de, dfnorm, jbes, &
             factor)
       iterm = int(log10(abs(dfnorm)))
       dfnorm = dfnorm * (10.0e0_knd ** (-iterm))
       idfe = iterm
       dmlf = 1.0e0_knd / dfnorm
       idmlfe = -idfe
       if(ioprad == 2) ir2est = int(log10(wront)) - ir1de + 1
if (debug) then
       if(knd == kindd) write(40, 180) r1c, ir1e, r1dc, ir1de
       if(knd == kindq) write(40, 185) r1c, ir1e, r1dc, ir1de
180      format(10x,'r1 = ', f17.14, i6, 5x,'r1d = ',f17.14, i6)
185      format(10x,'r1 = ', f33.30, i6, 5x,'r1d = ',f33.30, i6)
end if
       if(ioprad /= 2) go to 680
!
!  determine prolate radial functions of the second kind
!
if (debug) then
       write(40, 187)
187      format(4x,'r2 and r2d calculation')
end if
!  calculation using integration technique
       if(iopint == 0) go to 230
       if(iopint == 2) go to 190
       limint = lnum + 3 * ndec + int(c)
       if(igau == 0) call gauss_cached(ndec, ngau, xr, wr)
       igau = 1
       ngqs = 10
       if(c > 2000.0e0_knd) ngqs = ngqs * (c / 2000.0e0_knd)* &
                     (c / 2000.0e0_knd)
       call pint(c, m, lnum, x1, limint, maxint, maxlp, maxmp, ndec, &
            wr, xr, ngau, ngqs, rpint1, rpint2, pint1, &
            pint2, pint3, pint4, norme, pnormint, ipnormint, &
            coefme, coefmo)
190      continue
       if(iopint == 1) limint = 3 * ndec + int(c)
       if(iopint == 2) limint = jintm + jintm + 20 + int(sqrt(c))
       call r2int(l, m, c, x, limint, ndec, nex, maxd, enr, d01, id01, &
             maxint, maxmp, maxlp, rpint1, rpint2, pint1, pint2, &
             pint3, pint4, norme, pnormint, ipnormint, coefme, &
             coefmo, r2ic, ir2ie, r2dic, ir2die, jint, coefn, &
             icoefn)
       iopint = 2
       if(jint > jintm) jintm = jint
       wronca = r1c * r2dic * 10.0e0_knd ** (ir1e+ir2die)
       wroncb = r2ic * r1dc * 10.0e0_knd ** (ir2ie+ir1de)
       wronc = wronca - wroncb
       naccint = -int(log10(abs((wronc - wront) / wront) + dec))
       if(naccint < 0) naccint = 0
       nsubw = -int(log10(abs(wronc / wronca) + dec))
       if(nsubw < 0) nsubw = 0
       if(naccint > 1) naccint = naccint + nsubw
       if(naccint > ndec - 1) naccint = ndec - 1
if (debug) then
        if(nsubw > 0) then
        write(40, 200) nsubw
        end if
200      format(15x,'sub. error in forming wronskian = ',i3, &
           ' digits.')
end if
        if(naccint >= naccr) then
        naccr = naccint
        r2c = r2ic
        ir2e = ir2ie
        r2dc = r2dic
        ir2de = ir2die
        nacce = 0
        end if
       if(naccint >= minacc .and. iopneu /= 0) iopneu = 4
       if(naccint >= minacc .and. (iopeta == 1 .or. iopeta == 2 .or. &
          iopeta == 4)) iopeta = 4
       istartr2 = 1
        if(naccint >= minacc .and. naccintp >= minacc) then
        istartr2 = 0
        iopneu = 0
        iopeta = 0
        end if
       if(naccint >= minacc .and. ndec - jsub <= naccint .and. &
           iopleg /= 0) iopleg = 0
       if(naccint < minacc .and. x1 <= 0.1e0_knd .and. &
         iopleg == 0 .and. l >= legstart .and. &
         jsub <= ndec - naccint .and. jsub <= ndec - naccrp) iopleg = 1
       if(naccint == 0 .and. naccintp == 0) iopint = 0
if (debug) then
       if(knd == kindd) write(40, 210) naccint, r2ic, ir2ie, r2dic, &
                       ir2die
       if(knd == kindq) write(40, 220) naccint, r2ic, ir2ie, r2dic, &
                       ir2die
210      format(15x,'accuracy in decimal digits = ',i2,/,10x, &
           'r2 = ',f17.14, i6, 5x,'r2d = ',f17.14, i6)
220      format(15x,'accuracy in decimal digits = ',i2,/,10x, &
           'r2 = ',f33.30, i6, 5x,'r2d = ',f33.30, i6)
end if
230      continue
!
!  calculation using Legendre expansion and joining factor
       if(iopleg == 0) go to 360
       if(jflagleg == 1) go to 310
       jflagleg = 1
       limdr = c + 2 * ndec + 50.0e0_knd * x1 + 200
       if(limdr > maxdr - 2) limdr = maxdr - 2
       if(ioppsum == 0) go to 250
       xin(1) = x
       limpleg = limdr + limdr
       iopd = 3
       call pleg_cached(m, limpleg, maxp, limcsav, iopd, ndec, nex, xin, 1, maxt, &
            prat, pdrat, pdnorma, ipdnorma, pnorma, ipnorma, &
            alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
       limcsav = max(limcsav, limpleg)
        do jj = 1, limpleg
        prx(jj) = prat(1, jj)
        pdrx(jj) = pdrat(1, jj)
        end do
250      limq = lnum + 3 * ndec + int(c)
       call qleg_cached(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, &
            iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
       fajo(1) = c / (rm2 - 1.0e0_knd)
       ifajo(1) = 0
       if(m == 0) go to 280
        do im = 1, m
        fajo(1) = fajo(1) * (im + im) / c
        if(abs(fajo(1)) < 1.0e+10_knd) go to 260
        fajo(1) = fajo(1) * (1.0e-10_knd)
        ifajo(1) = ifajo(1) + 10
260       continue
        if(abs(fajo(1)) > 1.0e-10_knd) go to 270
        fajo(1) = fajo(1) * (1.0e+10_knd)
        ifajo(1) = ifajo(1) - 10
270       continue
        end do
280      continue
       fajo(2) = -c * fajo(1) / (rm2 - 3.0e0_knd)
       ifajo(2) = ifajo(1)
        do jl = 3, lnum - 1, 2
        fajo(jl) = fajo(jl - 2) * real((jl + m + m - 1), knd) / (jl - 2)
        ifajo(jl) = ifajo(jl - 2)
        if(abs(fajo(jl)) < 1.0e10_knd) go to 290
        fajo(jl) = fajo(jl) * 1.0e-10_knd
        ifajo(jl) = ifajo(jl) + 10
290       fajo(jl + 1) = fajo(jl - 1) * real((jl + m + m - 1), knd) / (jl)
        ifajo(jl + 1) = ifajo(jl - 1)
        if(abs(fajo(jl + 1)) < 1.0e10_knd) go to 300
        fajo(jl + 1) = fajo(jl + 1) * 1.0e-10_knd
        ifajo(jl + 1) = ifajo(jl + 1) + 10
300       end do
       if(2 * (lnum / 2) == lnum .or. lnum == 2) go to 310
       fajo(lnum) = fajo(lnum - 2) * real((lnum + m + m - 1), knd) / (lnum - 2)
       ifajo(lnum) = ifajo(lnum - 2)
310      continue
!
       limleg = l - m + 3 * ndec + int(c)
       limdr = c + ndec + 50.0e0_knd * x1 + 200
       if(iopleg == 2) limleg = jleg + jleg + 20 + int(sqrt(c))
       if(iopleg == 2) limdr = jlegp + 10 + int(0.5e0_knd * sqrt(c))
       if(limdr > maxdr) limdr = maxdr
       nsdneg = 0
       dneg = 1.0e0_knd
       idneg = 0
       call dalt(l, m, c, limdr, maxdr, maxmp, ndec, nex, ioppsum, eigval, &
            enrneg, drhor, dneg, idneg, nsdneg, nsdrho)
       call r2leg(l, m, c, x1, lnum, limleg, limdr, ndec, nex, &
             maxd, maxmp, maxpdr, maxdr, maxq, enr, enrneg, drhor, &
             nsdrho, d01, id01, dneg, idneg, nsdneg, dfnorm, idfe, &
             dmfnorm, idmfe, prx, pdrx, qdr, qdml, iqdml, qdl, iqdl, &
             qr, qml, iqml, ql, iql, fajo, ifajo, jsub, termpq, &
             itermpq, ioppsum, iopqnsum, r1c, ir1e, r1dc, &
             ir1de, wront, minacc, r2lc, ir2le, r2dlc, ir2dle, &
             jleg, jlegp, naccleg, nsubw, jflagl, iflagq, iflagp)
if (debug) then
        if(nsubw > 0) then
        write(40, 200) nsubw
        end if
end if
       if(nacclegp - naccleg > 8) naccleg = nacclegp
       if(naccleg <= naccr) go to 320
       naccr = naccleg
       r2c = r2lc
       ir2e = ir2le
       r2dc = r2dlc
       ir2de = ir2dle
       nacce = jflagl
320      continue
       if(naccleg >= naccrsav) then
       iopleg = 2
       else
       if(iopleg == 1 .and. l /= m) iopleg = 0
       legstart = l + naccrsav - naccleg
       end if
       if(naccleg >= minacc .and. nacclegp >= minacc) then
       iopleg = 2
       iopneu = 0
       iopeta = 0
       end if
       if(iopint /= 0 .and. naccleg > max(naccint, naccintp) &
         .and. nacclegp > max(naccint, naccintp) &
         .and. jflagl == 0) iopint = 0
       nacclegp = naccleg
if (debug) then
       if(knd == kindd) write(40, 210) naccleg, r2lc, ir2le, r2dlc, &
                       ir2dle
       if(knd == kindq) write(40, 220) naccleg, r2lc, ir2le, r2dlc, &
                       ir2dle
end if
360      continue
!
!  calculation using conventional Neumann expansion (eta=1)
       if(iopneu == 0) go to 420
       if(iopneu == 2) go to 380
       if(ibflag1 == 1) go to 370
       ibflag1 = 1
       lnump = max(lnum + maxm, 64)
       limn1 = 2 * (lnump * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
          5 * ndec + 4 * m + c + 01000) + maxm
       if(x1 > 0.08e0_knd) limn1 = 2 * (lnump * (0.5e0_knd - 3.0e0_knd* &
                    log10(x1)) + 5 * ndec + 4 * m + c + 01000)+ &
                     maxm
       if(x1 > 1.0e0_knd) limn1 = 2 * (lnump * 0.5e0_knd + 5 * ndec + 4 * m + c+ &
                    00500) + maxm
       if(limn1 > maxn) limn1 = maxn
       call sphneu(c, x, limn1, maxn, maxlp, sneuf, sneun, ineue, sneudf, &
             sneudr)
370      if(ibflag2 == 1) go to 380
       ibflag2 = 1
       lp = max(lnum + m, 64)
       limp1 = 2 * (lp * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
          5 * ndec + 4 * m + c + 01000)
       if(x1 > 0.08e0_knd) limp1 = 2 * (lp * (0.5e0_knd - 3.0e0_knd* &
                    log10(x1)) + 5 * ndec + 4 * m + c + 01000)
       if(x1 > 1.0e0_knd) limp1 = 2 * (lp * 0.5e0_knd + 5 * ndec + 4 * m + c+ &
                    00500)
       if(limp1 > maxp) limp1 = maxp
       prat1(1) = 1.0e0_knd
       prat1(2) = rm2 + 1.0e0_knd
        do jp = 3, limp1
        aj1 = jp - 1
        aj2 = jp - 2
        prat1(jp) = (rm2 + aj1) * (rm2 + aj2) / (aj1 * aj2)
        end do
       pcoefn = x1 * (x1 + 2.0e0_knd) / (x * x)
       apcoefn = (rm / 2.0e0_knd) * log10(pcoefn)
       ipcoefn = int(apcoefn)
       pcoefn = 10.0e0_knd ** (apcoefn - ipcoefn)
380      continue
       lplus = max(l, lnum + maxm)
       limneu = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
           5 * ndec + 4 * m + c + 01000)
       if(x1 > 0.08e0_knd) limneu = 2 * ((lplus) * (0.5e0_knd- &
              3.0e0_knd * log10(x1)) + 5 * ndec + 4 * m + c + 01000)
       if(x1 > 1.0e0_knd) limneu = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                     4 * m + c + 00500)
       if(iopneu == 2 .and. naccneu > 0) limneu = jneu + jneu + 20+ &
                   int(sqrt(c)) + int(1.0e0_knd / x1)
       if(limneu > limp1 - 2) limneu = limp1 - 2
       call r2neu(l, m, c, x1, limneu, ndec, nex, maxd, maxlp, maxn, maxp, &
             minacc, enr, sneuf, sneun, ineue, sneudf, sneudr, &
             prat1, pcoefn, ipcoefn, dmfnorm, idmfe, r1dc, ir1de, &
             r2nc, ir2ne, r2dnc, ir2dne, jneu)
       wronca = r1c * r2dnc * 10.0e0_knd ** (ir1e+ir2dne)
       wroncb = r2nc * r1dc * 10.0e0_knd ** (ir2ne+ir1de)
       wronc = wronca - wroncb
       naccneu = -int(log10(abs((wronc - wront) / wront) + dec))
       if(naccneu < 0) naccneu = 0
       nsubw = -int(log10(abs(wronc / wronca) + dec))
       if(nsubw < 0) nsubw = 0
       if(naccneu > 1) naccneu = naccneu + nsubw
       if(naccneu > ndec - 1) naccneu = ndec - 1
if (debug) then
        if(nsubw > 0) then
        write(40, 200) nsubw
        end if
end if
       if(naccneup - naccneu > 8) naccneu = naccneup
       naccneup = naccneu
       if(naccneu <= naccr) go to 390
       naccr = naccneu
       iopneu = 2
       r2c = r2nc
       ir2e = ir2ne
       r2dc = r2dnc
       ir2de = ir2dne
       nacce = 0
390      continue
       if(naccneu > minacc) then
       nflag = 1
       iopneu = 2
       iopeta = 0
       iopint = 0
       end if
       if(naccneu == minacc) then
       if(iopeta /= 0) iopeta = 4
       end if
       if(iopeta == 0 .and. naccr < minacc .and. nflag == 0) then
       iopeta = 1
       end if
if (debug) then
       if(knd == kindd) write(40, 210) naccneu, r2nc, ir2ne, r2dnc, &
                       ir2dne
       if(knd == kindq) write(40, 220) naccneu, r2nc, ir2ne, r2dnc, &
                       ir2dne
end if
420      continue
!
!  calculation using the variable eta expansion
       if(iopeta == 0 .or. iopeta == 4) go to 670
        do 430 inn = 1, 100
        nees(inn) = 0
        naccsav(inn) = 0
430       continue
       inen = 0
       neemark = nee
       naccetamax = 0
       neemax = nee
       naccnmax = 0
       nacctemp = 0
       kounte = 0
       netatry = 1
       naccdp = 0
       if(iopeta > 1) go to 440
       kounter = 0
        if(ijnet == 0) then
         do jnet = 1, neta
         ang = (neta + 1 - jnet) * pi * 0.5e0_knd / (neta + 1)
         eta(jnet) = cos(ang)
         wmeta2(jnet) = 2.0e0_knd * (1.0e0_knd + eta(jnet))* &
                (sin(0.5e0_knd * ang) ** 2)
         xbn(jnet) = sqrt(x1 * (x1 + 2.0e0_knd) + eta(jnet)* &
              eta(jnet))
         xln(jnet) = eta(jnet) * (x1 + 1.0e0_knd) / xbn(jnet)
         end do
        ijnet = 1
        end if
       iopeta = 2
440      if(iopeta == 3) go to 540
       etaval = eta(nee)
450      xbninp = xbn(nee)
       netainp = 1
       etainp(1) = eta(nee)
       xlninp(1) = xln(nee)
       lplus = max(l, lnum + maxm)
       limn = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
          5 * ndec + 10 * incnee + 4 * m + c + 05000) + m
       if(x1 > 0.08e0_knd) limn = 2 * ((lplus) * (0.5e0_knd- &
         3.0e0_knd * log10(x1)) + 10 * ndec + 10 * incnee+4 * m + c + 01000) + m
       if(x1 > 1.0e0_knd) limn = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                    10 * incnee + 4 * m + c + 00500) + m
       if(limn > maxn - 2) limn = maxn - 2
       limp = limn - m
       if(jnen == 0) go to 510
       jnenlim = jnen
       if(jnen > jnenmax) jnenlim = jnenmax
       limplim = limp
       limnlim = limn
        do 500 jn = 1, jnenlim
        if(nee /= neeb(jn)) go to 500
        if(limplim > limpsv(jn)) limplim = limpsv(jn)
        if(limnlim > limnsv(jn)) limnlim = limnsv(jn)
         do 460 je = 1, limplim
         if(je <= limpd) pratb(je) = pratbsv(jn, je)
         pratt(je) = prattsv(jn, je)
         pdratt(je) = pdratsv(jn, je)
460        continue
         do 470 je = 1, limnlim
         sneufe(je) = sneufsv(jn, je)
         sneudfe(je) = sneudfsv(jn, je)
470        continue
         jelim = maxlp
         if(maxlp > limn + 1) jelim = limn + 1
         if(jelim > jelimsv(jn)) jelim = jelimsv(jn)
         do 480 je = 1, jelim
         sneune(je) = sneunsv(jn, je)
         sneudre(je) = sneudrsv(jn, je)
         ineuee(je) = ineuesv(jn, je)
480        continue
if (debug) then
        write(40, 490) etaval
490       format(8x,'r2eta: reused expansion functions for eta =', f13.9,'.')
end if
        go to 530
500       continue
510      continue
       jnen = jnen + 1
       jnencur = jnen - (jnenmax * int((jnen - 1) / jnenmax))
       neeb(jnencur) = nee
       call sphneu(c, xbninp, limn, maxn, maxlp, sneufe, sneune, &
             ineuee, sneudfe, sneudre)
        do je = 1, limn
        sneufsv(jnencur, je) = sneufe(je)
        sneudfsv(jnencur, je) = sneudfe(je)
        limnsv(jnencur) = limn
        end do
       jelim = maxlp
       if(maxlp > limn + 1) jelim = limn + 1
        do 520 je = 1, jelim
        sneunsv(jnencur, je) = sneune(je)
        sneudrsv(jnencur, je) = sneudre(je)
        ineuesv(jnencur, je) = ineuee(je)
520       continue
       jelimsv(jnencur) = jelim
       iopd = 3
       call pleg_cached(m, limp, maxp, limcsav, iopd, ndec, nex, xlninp, &
            netainp, maxt, prat, pdrat, pdnorma, ipdnorma, pnorma, &
            ipnorma, alpha, beta, gamma, coefa, coefb, coefc, &
            coefd, coefe)
       limcsav = max(limcsav, limp)
        do je = 1, limp
        pratt(je) = prat(1, je)
        pdratt(je) = pdrat(1, je)
        prattsv(jnencur, je) = pratt(je)
        pdratsv(jnencur, je) = pdratt(je)
        limpsv(jnencur) = limp
        end do
       limpd = 2 * (lnum + int(c) + ndec)
       if(limpd > limp) limpd = limp
       iopd = 2
       call pleg_cached(m, limpd, maxp, limcsav, iopd, ndec, nex, etainp, &
            netainp, maxt, prat, pdrat, pdnorma, ipdnorma, pnorma, &
            ipnorma, alpha, beta, gamma, coefa, coefb, coefc, &
            coefd, coefe)
       iopd = 3
        do je = 1, limpd
        pratb(je) = prat(1, je)
        pratbsv(jnencur, je) = pratb(je)
        end do
       pratb(limpd + 1) = 0.0e0_knd
       pratb(limpd + 2) = 0.0e0_knd
       pratbsv(jnencur, limpd + 1) = 0.0e0_knd
       pratbsv(jnencur, limpd + 2) = 0.0e0_knd
530      continue
       pcoefe = ((x1 * (x1 + 2.0e0_knd)) / (x1 * (x1 + 2.0e0_knd)+ &
           eta(nee) ** 2))
       apcoef = (rm / 2.0e0_knd) * log10(pcoefe)
       ipcoefe = int(apcoef)
       pcoefe = 10.0e0_knd ** (apcoef - ipcoefe)
       pcoefo = pcoefe * pratt(2) / pratb(2)
       ipcoefo = ipcoefe
       pdcoefe = pcoefe
       if(m /= 0) pdcoefe = -pcoefe * rm * xln(nee) * xbn(nee) * xbn(nee)/ &
                 (x1 * (x1 + 2.0e0_knd) * wmeta2(nee))
       ipdcoefe = ipcoefe
       pdcoefo = pdcoefe * pdratt(2) / pratb(2)
       ipdcoefo = ipdcoefe
       if(li < 3) go to 540
        do jl = 3, li + ix, 2
        pcoefe = pcoefe * pratt(jl) / pratb(jl)
        iterm = log10(abs(pcoefe))
        pcoefe = pcoefe * 10.0e0_knd ** (-iterm)
        ipcoefe = ipcoefe + iterm
        pdcoefe = pdcoefe * pdratt(jl) / pratb(jl)
        iterm = log10(abs(pdcoefe))
        pdcoefe = pdcoefe * 10.0e0_knd ** (-iterm)
        ipdcoefe = ipdcoefe + iterm
        end do
       continue
       if(li < 4) go to 540
        do jl = 4, li + 1 - ix, 2
        pcoefo = pcoefo * pratt(jl) / pratb(jl)
        iterm = log10(abs(pcoefo))
        pcoefo = pcoefo * 10.0e0_knd ** (-iterm)
        ipcoefo = ipcoefo + iterm
        pdcoefo = pdcoefo * pdratt(jl) / pratb(jl)
        iterm = log10(abs(pdcoefo))
        pdcoefo = pdcoefo * 10.0e0_knd ** (-iterm)
        ipdcoefo = ipdcoefo + iterm
        end do
540      continue
       if(ix == 0) go to 550
       pcoefet = pcoefo
       ipcoefet = ipcoefo
       pcoefo = pcoefo * pratt(li + 2) / pratb(li + 2)
       iterm = int(log10(abs(pcoefo)))
       pcoefo = pcoefo * 10.0e0_knd ** (-iterm)
       ipcoefo = ipcoefo + iterm
       pdcoefet = pdcoefo
       ipdcoefet = ipdcoefo
       pdcoefo = pdcoefo * pdratt(li + 2) / pratb(li + 2)
       iterm = int(log10(abs(pdcoefo)))
       pdcoefo = pdcoefo * 10.0e0_knd ** (-iterm)
       ipdcoefo = ipdcoefo + iterm
       go to 560
550      pcoefet = pcoefe
       ipcoefet = ipcoefe
       pcoefe = pcoefe * pratt(li + 2) / pratb(li + 2)
       iterm = int(log10(abs(pcoefe)))
       pcoefe = pcoefe * 10.0e0_knd ** (-iterm)
       ipcoefe = ipcoefe + iterm
       pdcoefet = pdcoefe
       ipdcoefet = ipdcoefe
       pdcoefe = pdcoefe * pdratt(li + 2) / pratb(li + 2)
       iterm = int(log10(abs(pdcoefe)))
       pdcoefe = pdcoefe * 10.0e0_knd ** (-iterm)
       ipdcoefe = ipdcoefe + iterm
560      continue
       lplus = max(l, lnum + maxm)
       limeta = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
           5 * ndec + 4 * m + c + 05000)
       if(x1 > 0.08e0_knd) limeta = 2 * ((lplus) * (0.50e0_knd- &
          3.0e0_knd * log10(x1)) + 5 * ndec + 4 * m + c + 01000)
       if(x1 > 1.0e0_knd) limeta = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                     4 * m + c + 00500)
       if(iopeta == 3 .and. naccrsav > minacc) &
               limeta = jeta + jeta + 500 + c
       if(iopeta == 3 .and. naccrsav <= minacc) &
               limeta = jeta + jeta + 500 + c
       if(iopeta == 2) limeta = max(limeta, jeta + jeta + 500 + int(c))
       if(limeta > limp - 2) limeta = limp - 2
       if(limeta > limd) limeta = limd
       wm = wmeta2(nee)
       call r2eta(l, m, c, x1, etaval, nee, incnee, limeta, ndec, nex, &
             maxd, maxlp, maxn, maxp, minacc, wm, enr, sneufe, &
             sneune, ineuee, sneudfe, sneudre, pdratt, &
             pratb, pratt, pcoefet, ipcoefet, pdcoefet, &
             ipdcoefet, r1c, ir1e, r1dc, ir1de, naccetamax, naccr, &
             r2ec, ir2ee, r2dec, ir2dee, nacceta, nacciop, jeta, &
             iopnee, neemark, naccd, naccn, naccnmax, naccns)
       netatry = netatry + 1
       naccetas = nacceta
       if(naccetas == 0 .and. naccmax == 0 .and. iopnee == 0 &
          .and. naccnmax < 3) neemax = nee
        if(naccetas == naccmax .and. naccetas > 0) then
        kounte = kounte + 1
        if(kounte == 1) nee1 = nee
        if(kounte == 2) nee2 = nee
         if(kounte > 2) then
         neemax = nee1
         nee1 = nee2
         nee2 = nee
         end if
        end if
        if(naccetas > naccmax) then
        naccmax = naccetas
        kounte = 0
        neemax = nee
        end if
        if(naccetas < naccmax) then
        kounte = 0
        end if
if (debug) then
       if(nacciop == 0) write(40, 590)
590      format(15x,'r2eta accuracy is calculated using the', &
           ' wronskian.')
       if(nacciop == 1) write(40, 600)
600      format(15x,'r2eta accuracy = estimated numerator', &
           ' accuracy minus sub. error in forming Wronskian.')
end if
       iopeta = 3
       if(naccetas > nacctemp) nacctemp = naccetas
        if(naccetas > naccr .or. (naccetas == naccr .and. &
           nacciop == 0)) then
        naccr = naccetas
        r2c = r2ec
        ir2e = ir2ee
        r2dc = r2dec
        ir2de = ir2dee
        nacce = nacciop
        end if
       if(iopint /= 0 .and. naccetas > minacc .and. naccint < &
         6 .and. naccintp < 6) iopint = 0
if (debug) then

       if(knd == kindd) write(40, 620) naccetas, etaval, nee, r2ec, &
                     ir2ee, r2dec, ir2dee
       if(knd == kindq) write(40, 630) naccetas, etaval, nee, r2ec, &
                     ir2ee, r2dec, ir2dee
620      format(15x,'accuracy in decimal digits = ',i2, 5x, &
          'eta = ',f17.14,'; nee = ',i4/,10x, &
           'r2 = ',f17.14, i5, 5x,'r2d = ',f17.14, i5)
630      format(15x,'accuracy in decimal digits = ',i2, 5x, &
          'eta = ',f17.14,'; nee = ',i4/,10x, &
           'r2 = ',f33.30, i5, 5x,'r2d = ',f33.30, i5)
end if
       if(naccetas > naccetamax .or. (naccetas == 0 .and. naccetamax &
         == 0 .and. iopnee == 0)) neemax = nee
       if(naccetas > naccetamax) naccetamax = naccetas
       if(naccetas >= minacc) ietacount = ietacount + 1
       if(ietacount >= 5) incnflag = 1
        if(naccetas >= minacc .and. iplflag == 0) then
        nee = nee - incnee
        iopeta = 2
        if(nee < 1) nee = 1
        iplflag = 1
        end if
        if(naccetas >= minacc) go to 660
       iopeta = 2
       if(iplflag == 1 .and. incnflag == 1 .and. netatry == 2) &
           iopnee = 0
       ietacount = 0
       insflag = 0
       if(naccns == ndec .and. x <= 1.11e0_knd) insflag = 1
       if((naccd >= naccdp .or. naccd >= minacc .or. naccd >= naccrp &
          .or. naccd >= naccr) .and. (naccd + naccdp > 0) .and. &
         nee /= neta .and. insflag == 0) iopnee = 0
       naccdp = naccd
       if(iopnee == 0) go to 650
       if(iopnee == 2) nee = neemax - incnee
       if(iopnee == 1) nee = max(neemark, neemax) - incnee
       if(nee < 1) nee = 1
       if(iopnee == 2 .and. naccetas < lowtest - 1) go to 640
       incnee = 8
        if(knd == kindd) then
        if(x1 >= 0.01e0_knd) incnee = 16
        end if
        if(knd == kindq) then
        if(x1 >= 0.05e0_knd) incnee = 16
        if(x1 >= 0.1e0_knd) incnee = 32
        end if
640      if(nacctemp >= lowtest - 1) msearch = 1
       if(msearch == 0) iopeta = 0
       if(msearch == 0) lowtest = lowtest - 1
       go to 665
650      if(nee == neta) go to 660
       nee = nee + incnee
       if(nee > neta) nee = neta
       if(msearch /= 0) kounter = 0
       go to 440
660      continue
       if(naccetas < minacc .and. nee == neta) &
           nee = nee - incnee
       if(naccetas < minacc) iopeta = 2
       if(nee /= neta) msearch = 1
       if(naccetas >= minacc) kounter = kounter + 1
       if(kounter >= (2 * incnee) .and. msearch /= 0) &
           incnee = 2 * incnee
       if(incnee > 64) incnee = 64
        if(knd == kindd) then
        if(x1 <= 0.2e0_knd .and. incnee > 32) incnee = 32
        if(x1 <= 0.15e0_knd .and. incnee > 16) incnee = 16
        end if
        if(knd == kindq) then
        if(x1 <= 0.1e0_knd .and. incnee > 32) incnee = 32
        end if
       if(iopint /= 0 .and. naccetas < lowacc) iopeta = 0
       if(iopeta == 0) nacctest = naccetas
       if(naccetas < minacc) iplflag = 0
       naccetap = naccetamax
665      if(naccetas < minacc .and. iflagnee == 0) iflagnee = 1
       if(nee < neest) nee = neest
670      if(naccr > 0) go to 680
       naccr = 0
       r2c = 0.0e0_knd
       ir2e = 0
       r2dc = 0.0e0_knd
       ir2de = 0
680      continue
if (output) then
       if(nacce /= 1) chr = 'w'
       if(nacce == 1) chr = 'e'
       if(ioprad == 2) write(20, 690)l, r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, chr
       if(ioprad == 1) write(20, 710) l, r1c, ir1e, r1dc, ir1de
690      format(1x, i6, 2x, 4(f17.14, 1x, i6, 2x), i2, a)
710      format(1x, i6, 2x, 2(f17.14, 1x, i6, 2x))
end if
       qr1(li, mi) = r1c
       ir1(li, mi) = ir1e
       qr1d(li, mi) = r1dc
       ir1d(li, mi) = ir1de
        if(ioprad == 2) then
        qr2(li, mi) = r2c
        ir2(li, mi) = ir2e
        qr2d(li, mi) = r2dc
        ir2d(li, mi) = ir2de
        nar(li, mi) = naccr
        end if
       if(ioprad /= 2) go to 720
       if(lowacc > naccr) lowacc = naccr
       if(iopint /= 0) naccintp = naccint
        if(istartr2 == 1) then
        if(ndec - jsub > naccr .and. x1 <= 0.4e0_knd .and. &
          iopleg == 0 .and. l >= legstart) iopleg = 1
        if(ndec - jsub > naccr .and. x1 >= 0.00065e0_knd .and. &
          iopneu == 0 .and. iopleg /= 2) iopneu = 1
        if(iopeta == 0 .and. x1 >= 0.00065e0_knd .and. iopneu == 0 &
          .and. iopleg /= 2) iopeta = 1
        end if
       if(iopeta == 4) iopeta = 1
       if((naccr < legtest .or. naccrp < legtest) .and. &
         x1 <= 0.01e0_knd .and. l > legstart .and. iopleg == 0) &
         iopleg = 1
       if(iopneu == 0 .and. naccr < minacc .and. naccrp < minacc &
         .and. ndec - jsub > naccr .and. ndec - jsub > naccrp .and. &
         x1 >= 0.00065e0_knd) iopneu = 1
       naccrp = naccr
       if (warn) then
        if(ioprad == 2 .and. naccr < 6) then
        write(60,*) ' est. acc. = ',naccr, ' digits for x = ', &
             x,' c = ', c,' m = ',m,' l = ',l
        end if
        end if
720      if(iopang == 0) go to 850
!
!  determine first kind prolate angular function
       if(l == m) lims1 = 3 * ndec + int(c)
       if(l /= m) lims1 = jang + jang + 20 + int(sqrt(c))
       if(lims1 > maxp) lims1 = maxp
       call s1leg(l, m, c, iopang, iopnorm, barg, narg, lims1, ndec, maxt, &
             maxd, maxp, enr, pr, pdr, pdnorm, ipdnorm, pnorm, &
             ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo, &
             ptempe, iptempe, ptempo, iptempo, dmlms, idmlmse, &
             s1c, is1e, s1dc, is1de, dmlms1, idmlms1e, naccs, jang)
        do 810 jarg = 1, narg
        s1(li, jarg, mi) = s1c(jarg)
        s1d(li, jarg, mi) = s1dc(jarg)
        is1(li, jarg, mi) = is1e(jarg)
        is1d(li, jarg, mi) = is1de(jarg)
        nas(li, jarg, mi) = naccs(jarg)
if (debug) then
        if(knd == kindd) write(50, 740) barg(jarg), naccs(jarg)
        if(knd == kindq) write(50, 745) barg(jarg), naccs(jarg)
740       format(1x,'eta = ',f17.14,'  accuracy = ',i5, ' digits.')
745       format(1x,'eta = ',f33.30,'  accuracy = ',i5, ' digits.')
end if
if (output) then
        if(iopang == 1) write(30, 750) barg(jarg), s1c(jarg), is1e(jarg), naccs(jarg)
        if(iopang == 2) write(30, 760) barg(jarg), s1c(jarg), is1e(jarg), s1dc(jarg), is1de(jarg), naccs(jarg)
750       format(1x, f19.14, 2x, f17.14, 2x, i5, 2x,', ',i2)
760       format(1x, f19.14, 2x, f17.14, 2x, i5, 2x, f17.14, 2x, i5, 2x, i2)
end if
if (debug) then
        if(knd == kindd .and. iopang == 1) write(50, 770) s1c(jarg), is1e(jarg)
        if(knd == kindd .and. iopang == 2) write(50, 780) s1c(jarg), is1e(jarg), s1dc(jarg), is1de(jarg)
        if(knd == kindq .and. iopang == 1) write(50, 790) s1c(jarg), is1e(jarg)
        if(knd == kindq .and. iopang == 2) write(50, 800)s1c(jarg), is1e(jarg), s1dc(jarg), is1de(jarg)
770       format(12x,'s1 = ',f17.14, 2x, i5)
780       format(12x,'s1 = ',f17.14, 2x, i5, 5x,'s1d = ',f17.14, 2x, i5)
790       format(12x,'s1 = ',f33.30, 2x, i5)
800       format(12x,'s1 = ',f33.30, 2x, i5,/,10x,'s1d = ',f33.30, 2x, i5)
end if
810       continue
850      continue
900     continue
     return
     end subroutine
!
!
    subroutine s1leg (l, m, c, iopang, iopnorm, barg, narg, lims1, ndec, &
             maxt, maxd, maxp, enr, pr, pdr, pdnorm, ipdnorm, &
             pnorm, ipnorm, pdtempe, ipdtempe, pdtempo, &
             ipdtempo, ptempe, iptempe, ptempo, iptempo, &
             dmlms, idmlmse, s1c, is1e, s1dc, is1de, dmlms1, &
             idmlms1e, naccs, jang)
!
!  purpose:     To calculate the prolate angular functions of the first
!               kind and their first derivatives with respect to eta.
!
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               iopang  : index = 1 when angular functions of the
!                         first kind are calculated; = 2 when the
!                         first derivatives with respect to eta are
!                         also calculated
!               iopnorm : = 1 when the angular functions (and
!                         first derivatives) are scaled by the
!                         square root of the normalization of the
!                         corresponding Legendre function, giving them
!                         unity norm; iopnorm = 0 otherwise
!               barg    : array of eta values for which angular
!                         functions are desired
!               narg    : number of eta values
!               lims1   : approximately twice the maximum number
!                         of terms available to be taken in the
!                         sums
!               ndec    : number of decimal digits for real(knd)
!               maxt    : dimension of barg, pdnorm, ipdnorm, pnorm,
!                         ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo,
!                         ptempe, iptempe, ptempo, iptempo, s1c, is1e,
!                         s1dc, is1de, and naccs arrays. first dimension
!                         of the doubly dimensioned arrays pr and pdr
!               maxd    : dimension of enr array
!               maxp    : second dimension of pr and pdr arrays
!               enr     : array of d coefficient ratios
!               pr      : array of ratios of successive first kind
!                         associated Legendre functions of the same
!                         parity
!               pdr     : array of ratios of successive derivatives of
!                         first kind associated Legendre functions of
!                         the same parity
!               pdnorm  : array of characteristics of the first
!                         derivatives of associated Legendre functions
!                         of the first kind of order m and degree m
!               ipdnorm : array of exponents corresponding to pdnorm
!               pnorm   : array of characteristics of the associated
!                         Legendre functions of the first kind of order
!                         m and degree m
!               ipnorm  : array of exponents corresponding to pnorm
!               pdtempe : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is even
!                         or odd, to the first derivative of the
!                         function of order m and degree m
!               ipdtempe: array of exponents corresponding to pdtempe
!               pdtempo : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is odd
!                         or even, to the first derivtive of the
!                         function of order m and degree m
!               ipdtempo: array of exponents corresponding to pdtempo
!               ptempe  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is even or odd, to the function
!                         of order m and degree m
!               iptempe : array of exponents corresponding to ptempe
!               ptempo  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is odd or even, to the
!                         function of order m and degree m
!               iptempo : array of exponents corresponding to ptempo
!               dmlms   : characteristic of d coefficient with subscript
!                         l - m using the Meixner and Schafke
!                         normalization
!               idmlmse : exponent associated with dmlms
!
!
!     output:   s1c    : array of characteristics of prolate
!                        angular functions of the first kind
!               is1e   : array of exponents of prolate angular
!                        functions of the first kind
!               s1dc   : array of characteristics of derivative with
!                        respect to eta of prolate angular functions
!                        of the first kind
!               is1de  : array of exponents of derivatives with respect
!                        to eta of prolate angular functions of first
!                        kind
!               dmlms1 : characteristic of the d coefficient with
!                        subscript l - m for unit normalization of
!                        the angular functions
!               idmlms1e: exponent associated with dmlms1
!               naccs  : array of integer estimates of the number of
!                        accurate decimal digits in the values obtained
!                        for s1 and s1d
!               jang   : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) adec, aj, c, dcon, dec, dnew, dmlms, dmlms1, dnewd, dold, &
         doldd, factor, fterm, rm2, rm2m1, rm2m3, rm2p1, s1, s1d
    real(knd) barg(maxt), enr(maxd), pdr(maxt, maxp), pdnorm(maxt), &
         pnorm(maxt), pr(maxt, maxp), pdtemp(maxt), ptemp(maxt), &
         pdtempe(maxt), ptempe(maxt), pdtempo(maxt), ptempo(maxt), &
         s1c(maxt), s1dc(maxt)
!
!  integer arrays
    integer ipdnorm(maxt), ipnorm(maxt), ipdtemp(maxt), iptemp(maxt), &
        ipdtempe(maxt), iptempe(maxt), ipdtempo(maxt), &
        iptempo(maxt), is1de(maxt), is1e(maxt), naccs(maxt)
!
    dec = 10.0e0_knd ** (-ndec - 1)
    dcon = dec
    adec = 1000.0e0_knd * dec
    rm2 = m + m
    rm2m1 = m + m - 1
    rm2p1 = m + m + 1
    rm2m3 = m + m - 3
    if(l > (m + 1)) go to 30
     do 20 k = 1, narg
     if(pnorm(k) == 0.0e0_knd) go to 20
     if(l == (m + 1)) go to 10
     ptempe(k) = pr(k, 1)
     iptempe(k) = 0
     ptemp(k) = ptempe(k)
     iptemp(k) = 0
     if(iopang /= 2) go to 20
     pdtempe(k) = pdr(k, 1)
     ipdtempe(k) = 0
     pdtemp(k) = pdtempe(k)
     ipdtemp(k) = 0
     go to 20
10    ptempo(k) = pr(k, 2)
     iptempo(k) = 0
     ptemp(k) = ptempo(k)
     iptemp(k) = 0
     if(iopang /= 2) go to 20
     pdtempo(k) = pdr(k, 2)
     ipdtempo(k) = 0
     pdtemp(k) = pdtempo(k)
     ipdtemp(k) = 0
20    continue
30   continue
    lm2 = (l - m) / 2
    ix = l - m - 2 * lm2
    ixx = ix - 1
    ixx2 = ixx + 2
    if(l < (m + 2)) go to 110
     do 100 k = 1, narg
     if(pnorm(k) == 0.0e0_knd) go to 100
     if(ix /= 0) go to 60
     ptempe(k) = ptempe(k) * pr(k, l - m + 1)
     if(abs(ptempe(k)) < 1.0e+10_knd) go to 40
     ptempe(k) = ptempe(k) * (1.0e-10_knd)
     iptempe(k) = iptempe(k) + 10
40    ptemp(k) = ptempe(k)
     iptemp(k) = iptempe(k)
     if(iopang /= 2) go to 100
     if(abs(barg(k)) < adec) go to 100
     pdtempe(k) = pdtempe(k) * pdr(k, l - m + 1)
     if(abs(pdtempe(k)) < 1.0e+10_knd) go to 50
     pdtempe(k) = pdtempe(k) * (1.0e-10_knd)
     ipdtempe(k) = ipdtempe(k) + 10
50    pdtemp(k) = pdtempe(k)
     ipdtemp(k) = ipdtempe(k)
     go to 100
60    if(abs(barg(k)) < adec) go to 80
     ptempo(k) = ptempo(k) * pr(k, l - m + 1)
     if(abs(ptempo(k)) < 1.0e+10_knd) go to 70
     ptempo(k) = ptempo(k) * (1.0e-10_knd)
     iptempo(k) = iptempo(k) + 10
70    ptemp(k) = ptempo(k)
     iptemp(k) = iptempo(k)
     if(iopang /= 2) go to 100
80    pdtempo(k) = pdtempo(k) * pdr(k, l - m + 1)
     if(abs(pdtempo(k)) < 1.0e+10_knd) go to 90
     pdtempo(k) = pdtempo(k) * (1.0e-10_knd)
     ipdtempo(k) = ipdtempo(k) + 10
90    pdtemp(k) = pdtempo(k)
     ipdtemp(k) = ipdtempo(k)
100    continue
110   continue
    lim = lims1 / 2 - ix
    jlow = lm2 + 1
    jang = 0
!
!  compute the associated Legendre function normalization factor
    factor = 1.0e0_knd
    ifactor = 0
    if(iopnorm == 0) go to 210
    if(m == 0) go to 170
     do 160 j = 1, m
     aj = j
     factor = factor * (aj + aj) * (aj + aj - 1.0e0_knd)
     if(factor < 1.0e100_knd) go to 160
     factor = factor * 1.0e-100_knd
     ifactor = ifactor + 100
160    continue
170   if(l == m) go to 190
     do 180 j = 1, l - m
     aj = j
     factor = factor * (rm2 + aj) / (aj)
     if(factor < 1.0e100_knd) go to 180
     factor = factor * 1.0e-100_knd
     ifactor = ifactor + 100
180    continue
190   factor = factor * 2.0e0_knd / (l + l + 1.0e0_knd)
    factor = sqrt(factor)
    ifactor = ifactor / 2
    iterm = int(log10(factor))
    factor = factor * (10.0e0_knd ** (-iterm))
    ifactor = ifactor + iterm
    dmlms1 = dmlms / factor
    idmlms1e = idmlmse - ifactor
if (debug) then
    write(50, 200)
200   format(5x,'s1 is normalized to have unit norm.')
end if
210 continue
if (debug) then
    if(iopnorm == 0) write(50, 215)
215   format(5x,'s1 has the same normalization as the', &
        ' corresponding Legendre function.')
end if
!
!  compute the angular function s1
     do 380 k = 1, narg
     if(pnorm(k) == 0.0e0_knd) go to 220
     if((ix == 1) .and. (abs(barg(k)) < adec)) go to 220
     if(((abs(abs(barg(k)) - 1.0e0_knd)) < adec) &
        .and. (m /= 0)) go to 220
     go to 230
220    s1c(k) = 0.0e0_knd
     is1e(k) = 0
     naccs(k) = ndec
     go to 300
230    dold = 1.0e0_knd
     s1 = dold
     fterm = 1.0e0_knd
      do 240 j = jlow, lim
      dnew = dold * enr(j) * pr(k, j + j + ixx2)
      s1 = s1 + dnew
      if(abs(dnew) > fterm) fterm = abs(dnew)
      if(abs(dnew / s1) < dcon) go to 250
      dold = dnew
240     continue
250     if(j > jang) jang = j
if (debug) then
     write(50, 260) barg(k), j
260    format(8x,'s1 calculation for eta = ',f17.14,' converged in ', &
         i6,' terms.')
end if
     if(lm2 < 1) go to 280
     dold = 1.0e0_knd
     j = lm2
      do 270 jj = 1, lm2
      dnew = dold / (pr(k, j + j + ixx2) * enr(j))
      s1 = s1 + dnew
      if(abs(dnew) > fterm) fterm = abs(dnew)
      if(abs(dnew / s1) < dcon) go to 280
      dold = dnew
      j = j - 1
270     continue
280    s1c(k) = s1 * dmlms * ptemp(k) * pnorm(k) / factor
     if(s1c(k) /= 0.0e0_knd) iterm = int(log10(abs(s1c(k))))
     if(s1c(k) == 0.0e0_knd) iterm = 0
     s1c(k) = s1c(k) * (10.0e0_knd ** (-iterm))
     is1e(k) = iptemp(k) + ipnorm(k) + iterm + idmlmse - ifactor
     if(abs(s1c(k)) >= 1.0e0_knd) go to 290
     s1c(k) = s1c(k) * 10.0e0_knd
     is1e(k) = is1e(k) - 1
     if(s1c(k) == 0.0e0_knd) naccs(k) = 0
290    if(s1c(k) /= 0.0e0_knd) naccs(k) = ndec - 2- &
                 log10(abs((fterm) / (s1)))
     if(naccs(k) > 0) go to 300
     naccs(k) = 0
     s1c(k) = 0.0e0_knd
     is1e(k) = 0
     s1dc(k) = 0.0e0_knd
     is1de(k) = 0
     go to 380
!
!       compute the first derivative of the anguar function when
!       iopang equals 2
300    if(iopang /= 2) go to 380
     if(pnorm(k) == 0.0e0_knd) go to 310
     if((ix == 0) .and. (abs(barg(k)) < adec)) go to 310
     if(((abs(abs(barg(k)) - 1.0e0_knd)) < adec) .and. (m /= 0) &
        .and. (m /= 2)) go to 310
     go to 320
310    s1dc(k) = 0.0e0_knd
     is1de(k) = 0
     go to 370
320    doldd = 1.0e0_knd
     s1d = doldd
     if(l == 0) s1d = 0.0e0_knd
      do 330 j = jlow, lim
      dnewd = doldd * enr(j) * pdr(k, j + j + ixx2)
      s1d = s1d + dnewd
      if(abs(dnewd / s1d) < dcon) go to 340
      doldd = dnewd
330     continue
340    if(lm2 < 1) go to 360
     doldd = 1.0e0_knd
     j = lm2
     ja = lm2
     if(m == 0 .and. ix == 0) ja = lm2 - 1
     if(ja == 0) go to 360
      do 350 jj = 1, ja
      dnewd = doldd / (pdr(k, j + j + ixx2) * enr(j))
      s1d = s1d + dnewd
      if(abs(dnewd / s1d) < dcon) go to 360
      doldd = dnewd
      j = j - 1
350     continue
360    s1dc(k) = s1d * dmlms * pdtemp(k) * pdnorm(k) / factor
     if(s1dc(k) /= 0.0e0_knd) iterm = int(log10(abs(s1dc(k))))
     if(s1dc(k) == 0.0e0_knd) iterm = 0
     s1dc(k) = s1dc(k) * 10.0e0_knd ** (-iterm)
     is1de(k) = ipdtemp(k) + ipdnorm(k) + iterm + idmlmse - ifactor
     if(abs(s1dc(k)) >= 1.0e0_knd) go to 370
     s1dc(k) = s1dc(k) * 10.0e0_knd
     is1de(k) = is1de(k) - 1
370    continue
380    continue
    return
    end subroutine
!
!
    subroutine r1bes (l, m, c, x1, limr1, ndec, maxd, enr, maxj, maxlp, &
             nex, iflag, sbesf, sbesdf, sbesn, ibese, sbesdr, &
             d01, id01, r1c, ir1e, r1dc, ir1de, dfnorm, jbes, &
             factor)
!
!  purpose:     To calculate the prolate radial function of the
!               first kind and its first derivative with respect
!               to x, using an expansion of spherical Bessel
!               functions of the first kind with argument
!               c*sqrt(x*x-1).
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x1     : x-1
!               limr1  : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits for real(knd)
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               maxj   : dimension of sbesf and sbesdf arrays
!               maxlp  : maximum  l value desired; dimension
!                        of the sbesdr, sbesn, and ibese arrays
!               nex    : maximum exponent available in real(knd)
!                        arithmetic
!               iflag  : integer = 1 if forward series not needed;
!                        =0 if the forward series is computed
!               sbesf  : array of ratios of successive first kind
!                        spherical Bessel functions of the same parity
!               sbesdf : array of ratios of successive derivatives of
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents corresponding to sbesn
!               sbesdr : value of ratio of first derivative of
!                        spherical Bessel function to the corresponding
!                        Bessel function
!               d01    : characteristic of ratio of first d coefficient
!                        to the d coefficient for n = l-m
!               id01   : exponent associated with d01
!
!     output  : r1c    : characteristic of prolate radial function
!                        of the first kind
!               ir1e   : exponent of prolate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               dfnorm : Flammer normalization factor of the
!                        d coefficients. equal to the reciprocal of
!                        the value of the d coefficient d(n = l - m)
!                        using this normalization for the angular
!                        functions
!               jbes   : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!               factor : coefficient used in alternative expression for
!                        contribution to r1d from the first order
!                        Bessel function and its first derivative.
!                        Used to avoid subtraction error that occurs
!                        for m = 0 when x1 is very small, especially
!                        for values of l. Value of the factor for a
!                        given l is used to obtain the value for the
!                        next l.
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, con, dec, dfnorm, dnew, dnewd, dold, doldd, d01, factor, &
         r1bot, r1c, r1d, r1dc, r1dstore, r1temp, r1top, r1dtop, &
         r1topd, r1dtopd, rj1, rj2, term, termd, teste, testeo, x, x1
    real(knd) enr(maxd), sbesdf(maxj), sbesdr(maxlp), sbesf(maxj), &
         sbesn(maxlp)
!
!  integer array
    dimension ibese(maxlp)
!
!  convergence ratio dec is set according to the requested accuracy
    dec = 10.0e0_knd ** (-ndec - 1)
    lm2 = (l - m) / 2
!
!  ix=0 for l-m even, ix=1 for l-m odd
    ix = l - m - 2 * lm2
    lim = limr1 / 2 - ix
    x = x1 + 1.0e0_knd
    con = x / sqrt(x1 * (x1 + 2.0e0_knd))
    nfac = nex / 3
    teste = 10.0e0_knd ** nfac
    testeo = 1.0e0_knd / teste
    ir1tope = 0
    mml = m + m - 1 + ix
    iflagd = 0
    if(x1 < 0.1e0_knd .and. ix == 1 .and. m == 0) iflagd = 1
    if(iflagd == 1 .and. l /= 1) factor = factor * real(l, knd) / (l - 1)
!
!
!  compute radial function of the first kind r1 and its first
!  derivative r1d
!
!  forward summation of numerator series for both r1 and r1d
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
    r1top = dold
    r1dtop = doldd
if (debug) then
    if(iflag == 1) write(40, 10) lm2
10   format(8x,'r1bes: numerator forward series not used; backward ', &
        'series has 'i5,' terms.')
end if
    jtop = lm2
    if(iflag == 1) go to 50
     do 20 j = lm2 + 1, lim
     jj = j + j + ix
     rj1 = jj - ix
     rj2 = jj + mml
     dnew = dold * enr(j) * sbesf(jj + m) * rj2 / rj1
     dnewd = doldd * enr(j) * sbesdf(jj + m) * rj2 / rj1
     r1top = r1top + dnew
     r1dtop = r1dtop + dnewd
     if((abs(dnew / r1top) + abs(dnewd / r1dtop)) < dec) go to 30
     dold = dnew
     doldd = dnewd
20    continue
30   continue
    jtop = min(lim, j)
if (debug) then
    write(40, 40) jtop, lim
40   format(8x,'r1bes: numerator series converged in ',i6,' terms; ', &
        i6,' available.' )
end if
    if(iflagd == 0 .or. l /= 1) go to 45
    r1topd = r1top - 1.0e0_knd
    r1dtopd = r1dtop - 1.0e0_knd
45   continue
50   continue
!
!  backward summation of numerator series for r1 and r1d
    if (lm2 < 1) go to 80
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
     do 70 j = lm2, 1,-1
     jj = j + j + ix
     rj1 = jj - ix
     rj2 = jj + mml
     dnew = dold * rj1 / (sbesf(jj + m) * rj2 * enr(j))
     dnewd = doldd * rj1 / (sbesdf(jj + m) * rj2 * enr(j))
     r1top = r1top + dnew
     r1dtop = r1dtop + dnewd
     if((abs(dnew / r1top) + abs(dnewd / r1dtop)) < dec) go to 75
     if(abs(r1top) < teste) go to 60
     r1top = r1top * testeo
     dnew = dnew * testeo
     ir1tope = ir1tope + nfac
     r1dtop = r1dtop * testeo
     dnewd = dnewd * testeo
60    dold = dnew
     doldd = dnewd
70    continue
     if(jj /= 3) iflagd = 0
     if(iflagd == 0) go to 80
     r1topd = r1top - dnew
     r1dtopd = r1dtop - dnewd
     go to 80
75    continue
     iflagd = 0
80    continue
!
!  forward summation of denominator series for r1 and r1d
!  the denominator series is the Flammer normalization constant dfnorm
    dold = 1.0e0_knd
    r1bot = dold
     do 90 j = lm2 + 1, lim
     jj = j + j + ix
     dnew = -dold * enr(j) * real((jj + mml), knd) / real(jj - ix, knd)
     r1bot = r1bot + dnew
     if(abs(dnew / r1bot) < dec) go to 100
     dold = dnew
90    continue
100   continue
    jbot = min(j, lim)
    jbes = max(jtop, jbot)
if (debug) then
    write(40, 110) jbot, lim
110   format(15x,'denominator series converged in ',i6,' terms; ', &
        i6,' available.')
end if
!
!  backward summation of denominator series for both r1 and r1d
    if(lm2 < 1) go to 130
    dold = 1.0e0_knd
     do 120 j = lm2, 1,-1
     jj = j + j + ix
     dnew = -dold * (jj - ix) / (real((jj + mml), knd) &
        *enr(j))
     r1bot = r1bot + dnew
     if(abs(dnew / r1bot) < dec) go to 130
     dold = dnew
120    continue
130   dfnorm = r1bot
!
!  compute r1 and r1d
    r1temp = r1top * sbesn(l + 1) / r1bot
    if(ix == 1) r1temp = r1temp * con
    iterm = int(log10(abs(r1temp)))
    ir1e = ir1tope + ibese(l + 1) + iterm
    r1c = r1temp * (10.0e0_knd ** (-iterm))
    if(abs(r1c) >= 1.0e0_knd) go to 140
    r1c = r1c * 10.0e0_knd
    ir1e = ir1e - 1
140   if(iflagd == 1) r1temp = r1temp * r1topd / r1top
    if(iflagd == 1) r1dtop = r1dtopd
    r1d = r1dtop * sbesn(l + 1) * c * con * sbesdr(l + 1) / r1bot
    r1dstore = r1d * con
    ndsub = 0
    if(ix == 1) r1d = r1d * con - r1temp / (x1 * x * (x1 + 2.0e0_knd))
    if(ix == 1) ndsub = -int(log10(abs(r1d / r1dstore)))
    if(ndsub < 0) ndsub = 0
    if(iflagd == 0) go to 150
    term = x1 * (x1 + 2.0e0_knd) * sbesdr(2) * sbesn(2)
    termd = term - sbesn(3) * (10.0e0_knd ** (ibese(3) - ibese(2)))
    ndsub1 = -log10(abs(termd / term))
    if(ndsub1 < 0) ndsub1 = 0
    termd = termd * (c * d01 / (x1 * (x1 + 2.0e0_knd) * r1bot * factor))* &
        (10.0e0_knd ** (id01 + ibese(2) - ibese(l + 1) - ir1tope))
    r1d = r1d + termd
    ndsub2 = -log10(abs(r1d / termd))
    if(ndsub2 < 0) ndsub2 = 0
    ndsub = ndsub + ndsub1 + ndsub2
150   continue
if (debug) then
    if(ix == 1 .and. ndsub > 0) write(40, 160) ndsub
160   format(24x,'subtraction error in forming r1d =',i3,' digits.')
end if
    iterm = log10(abs(r1d))
    ir1de = ir1tope + ibese(l + 1) + iterm
    r1dc = r1d * (10.0e0_knd ** (-iterm))
    if(abs(r1dc) >= 1.0e0_knd) go to 170
    r1dc = r1dc * 10.0e0_knd
    ir1de = ir1de - 1
170   continue
    mfac = ir1e - ibese(l + 1)
    if(mfac > (ndec + 5)) iflag = 1
    if(mfac <= (ndec + 5)) iflag = 0
    return
    end subroutine
!
!
    subroutine r2int (l, m, c, x, limint, ndec, nex, maxd, enr, d01, id01, &
             maxint, maxmp, maxlp, rpint1, rpint2, &
             pint1, pint2, pint3, pint4, norme, pnorm, ipnorm, &
             coefme, coefmo, r2c, ir2e, r2dc, ir2de, jint, &
             coefn, icoefn)
!
!
!  purpose:     To calculate values of the radial function of the
!               second kind and its first derivative using an integral
!               representation of the radial functions in terms of the
!               angular function of the first kind together with a
!               Neumann function kernal. The angular function is
!               expanded in a series of associated Legendre functions.
!               Gaussian quadrature is used (in subroutine pint) to
!               evaluate the resulting integrals involving associated
!               Legendre functions times the Neumann function kernel.
!               This subroutine performs the summation of the
!               integrals times d coefficients to obtain r2 and r2d.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x      : x
!               limint : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent for real(knd)
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               d01    : characteristic of the first d coefficient,
!                        either d0 or d1, depending on whether l-m
!                        is even or odd
!               id01   : exponent (base 10) of the first d coefficient
!               maxint : dimension of pint and rpint arrays
!               maxmp  : dimension of norme array
!               maxlp  : maximum  l value desired; dimension
!                        of the pnorm and ipnorm arrays
!               rpint1 : arrays of ratios of successive integrals of
!                        either the first or the third kind, depending
!                        on whether l-m is even or odd
!               rpint2 : array of ratios of successive integrals of
!                        either the second or the fourth kind,
!                        depending on whether l-m is even or odd
!               pint1  : array of scaled values for the integrals of
!                        the first kind
!               pint2  : array of scaled values for the integrals of
!                        the second kind
!               pint3  : array of scaled values for the integrals of
!                        the third kind
!               pint4  : array of scaled values for the integrals of
!                        the fourth kind
!               norme  : array of exponents used to scale the Neumann
!                        function of order m involved in the integrals
!               pnorm  : array of characteristics of the scaling factors
!                        used for the associated Legendre functions in
!                        the integrals to avoid overflow
!               ipnorm : array of exponents (base 10) corresponding to
!                        pnorm
!               coefme : coefficient used to multiply the resulting
!                        sum to obtain r2 when l-m is even
!               coefmo : coefficient used to multiply the resulting
!                        sum to obtain r2 when l-m is odd
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the second kind
!               ir2e   : exponent of prolate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the second kind
!               jint   : maximum value of the index j in the forward
!                        sum for r2 and r2d, i.e., the highest enr(j)
!                        used
!               coefn  : characteristic of coefficient that is only
!                        calculated once (for l = m) and is then
!                        used for all values of l
!               icoefn : exponent for coefn
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, coefa, coefl, coefme, coefmo, coefn, dec, dcon, dnew, &
         dnewd, dold, doldd, d01, ri, rm, rm2, r2c, r2dc, r2dpos, &
         r2dtemp, r2pos, r2temp, x
    real(knd) enr(maxd), pnorm(maxlp), pint1(maxint), pint2(maxint), &
         pint3(maxint), pint4(maxint), rpint1(maxint), &
         rpint2(maxint)
!
!  integer arrays
    dimension norme(maxmp), ipnorm(maxlp)
!
    rm = m
    rm2 = rm + rm
    lm2 = (l - m) / 2
    ix = l - m - 2 * lm2
    ixx = ix - 1
    ixx2 = ixx + 2
    lim = limint / 2 - ix
!
!  compute the leading coefficient
    if(l > m) go to 20
    icoefn = norme(m + 1)
    coefn = 0.5e0_knd
    if(m == 0) go to 20
     do 10 i = 1, m
     ri = i
     coefn = coefn / (ri + ri)
     iterm = int(log10(abs(coefn)))
     coefn = coefn * 10.0e0_knd ** (-iterm)
     icoefn = icoefn + iterm
10    continue
20   continue
    if(ix == 0) coefa = (rm2 + 1.0e0_knd) * coefn
    if(ix == 1) coefa = (rm2 + 3.0e0_knd) * coefn
    if((ix == 0) .and. (2 * (lm2 / 2) /= lm2)) coefa = -coefa
    if((ix == 1) .and. (2 * ((l - m - 1) / 4) /= (l - m - 1) / 2)) coefa = -coefa
    coefl = coefa / d01
    icoefl = -id01 + icoefn
    dec = 10.0e0_knd ** (-ndec - 1)
    dcon = dec
    jlow = lm2 + 1
!
!  compute the integrals involving the angular functions by summing
!  d coefficients times corresponding integrals of Legendre
!  functions
!
!  forward summation of series for r2 and r2d
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
    r2dtemp = doldd
    r2temp = dold
    r2pos = dold
    r2dpos = doldd
     do 30 j = jlow, lim
     dnew = dold * enr(j) * rpint1(j + j + ixx2)
     dnewd = doldd * enr(j) * rpint2(j + j + ixx2)
     r2temp = r2temp + dnew
     r2dtemp = r2dtemp + dnewd
     if(dnew > 0.0e0_knd) r2pos = r2pos + dnew
     if(dnewd > 0.0e0_knd) r2dpos = r2dpos + dnewd
     if((abs(dnew / r2temp) + abs(dnewd / r2dtemp)) < dcon) go to 40
     dold = dnew
     doldd = dnewd
30    continue
!
!  backward summation of series for r2 and r2d
40   jint = j
    if(jint > lim) jint = lim
    if(lm2 < 1) go to 60
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
    j = lm2
     do 50 jj = 1, lm2
     dnew = dold / (rpint1(j + j + ixx2) * enr(j))
     dnewd = doldd / (rpint2(j + j + ixx2) * enr(j))
     r2temp = r2temp + dnew
     r2dtemp = r2dtemp + dnewd
     if(dnew > 00.0e0_knd) r2pos = r2pos + dnew
     if(dnewd > 00.0e0_knd) r2dpos = r2dpos + dnewd
     if((abs(dnew / r2temp) + abs(dnewd / r2dtemp)) < dcon) go to 60
     dold = dnew
     doldd = dnewd
     j = j - 1
50    continue
60   continue
    isub = int(log10(abs(r2pos / r2temp) + dec))
    isubd = int(log10(abs(r2dpos / r2dtemp) + dec))
    if(isubd < 0) isubd = 0
    r2temp = r2temp * coefl * pnorm(l - m + 1)
    if(ix == 0) r2temp = r2temp * pint1(l - m + 1)
    if(ix == 1) r2temp = r2temp * pint3(l - m + 1) * x
    iterm = int(log10(abs(r2temp)))
    ir2e = iterm + ipnorm(l - m + 1) + icoefl
    r2c = r2temp * 10.0e0_knd ** (-iterm)
    if(abs(r2c) >= 1.0e0_knd) go to 70
    r2c = r2c * 10.0e0_knd
    ir2e = ir2e - 1
70   r2dtemp = -r2dtemp * coefl * pnorm(l - m + 1) * c * x
    if(ix == 0) r2dtemp = r2dtemp * pint2(l - m + 1) + r2temp * coefme
    if(ix == 1) r2dtemp = r2dtemp * pint4(l - m + 1) * x + r2temp * coefmo
    if(ix == 0) jsub = int(log10(abs(r2temp * coefme / r2dtemp) + dec))
    if(ix == 1) jsub = int(log10(abs(r2temp * coefmo / r2dtemp) + dec))
    if(jsub < 0) jsub = 0
    isub = max(isub, isubd + jsub)
if (debug) then
    write(40, 80) jint, lim, isub, isubd + jsub
80   format(8x,'r2int: converged in ',i6,' terms; 'i6, &
        ' available; ',i3,' and ',i3,' digits of sub. error.')
end if
    jterm = int(log10(abs(r2dtemp)))
    ir2de = jterm + ipnorm(l - m + 1) + icoefl
    r2dc = r2dtemp * 10.0e0_knd ** (-jterm)
    if(abs(r2dc) >= 1.0e0_knd) go to 90
    r2dc = r2dc * 10.0e0_knd
    ir2de = ir2de - 1
90   continue
    return
    end subroutine
!
!
    subroutine r2leg (l, m, c, x1, lnum, limleg, limdr, ndec, nex, &
             maxd, maxmp, maxpdr, maxdr, maxq, enr, enrneg, drhor, &
             nsdrho, d01, id01, dneg, idneg, nsdneg, dfnorm, idfe, &
             dmfnorm, idmfe, prx, pdrx, qdr, qdml, iqdml, qdl, &
             iqdl, qr, qml, iqml, ql, iql, fajo, ifajo, jsub, &
             termpq, itermpq, ioppsum, iopqnsum, r1c, ir1e, &
             r1dc, ir1de, wront, minacc, r2c, ir2e, r2dc, ir2de, &
             jleg, jlegp, naccleg, nsubw, jflagl, iflagq, iflagp)
!
!  purpose:     To evaluate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x using the traditional expansion in associated
!               Legendre functions.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               x1      : x-1
!               lnum    : number of l values desired
!               limleg  : approximately twice the maximum number
!                         of terms available to be taken in qsum,
!                         (sum involving q's time d coefficients)
!               limdr   : maximum number of terms available to be
!                         taken in psum (sum involving p's time
!                         d rho coefficients)
!               ndec    : number of decimal digits for real(knd)
!               nex     : largest exponent available for real(knd)
!               maxd    : dimension of enr array
!               maxmp   : dimension of enrneg array
!               maxpdr  : dimension of prx and pdrx arrays
!               maxdr   : dimension of drhor array
!               maxq    : dimension of qr and qdr arrays
!               enr     : array of d coefficient ratios
!               enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               nsdrho  : maximum subtraction error in calculating
!                         drhor array
!               d01     : characteristic of the ratio of the first d
!                         coefficient with nonnegative subscript,
!                         either d0 or d1 depending on whether l-m is
!                         even or odd, to the d coefficient with
!                         subscript l - m
!               id01    : exponent (base 10) corresponding to d01
!               dneg    : characteristic of the ratio of the d
!                         coefficient with subscript -2m+ix to the
!                         d coefficient with subscript ix, where
!                         ix = 0 or 1 depending on whether l - m
!                         is even or odd
!               idneg   : exponent corresponding to dneg
!               nsdneg  : subtraction error in calculating dneg, also
!                         maximum subtraction error in the enrneg array
!               dfnorm  : characteristic of the Flammer normalization
!                         factor of the d coefficients. equal to the
!                         reciprocal of the value of the d coefficient
!                         d(n = l - m) using this normalization for
!                         the angular functions
!               idfe    : exponent associated with dfnorm
!               dmfnorm : characteristic of the Morse-Feshbach
!                         normalization factor of the d coefficients.
!                         equal to the reciprocal of the value of the
!                         d coefficient d(n = l - m) using this
!                         normalization for the angular functions
!               idmfe   : exponent associated with dmfnorm
!               prx     : ratios of successive Legendre functions of
!                         the first kind of the same parity
!               pdrx    : ratios of successive first derivatives of
!                         Legendre functions of the first kind of the
!                         same parity
!               qdr     : ratios of first derivatives of successive
!                         Legendre functions of the second kind
!               qdml    : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m-1
!               iqdml   : exponent corresponding to qdml
!               qdl     : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m
!               iqdl    : exponent corresponding to qdl
!               qr      : array of ratios of successive associated
!                         Legendre functions of the second kind
!               qml     : characteristic of the associated Legendre
!                         function of the second kind with order m
!                         and degree m-1
!               iqml    : exponent corresponding to qml
!               ql      : characteristic of the associated Legendre
!                         function of the second kind with order m and
!                         degree m
!               iql     : exponent corresponding to ql
!               fajo    : characteristic of the joining factor of the
!                         second kind
!               ifajo   : exponent corresponding to fajo
!               jsub    : subtraction error in forming fajo coming from
!                         dmfnorm
!               termpq  : characteristic of the relative size of the
!                         maximum terms in the positive degree q series
!                         and the p series used to calculate r2 and r2d
!               itermpq : exponent corresponding to termpq
!               ioppsum : integer flag = 0 if psum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if psum is computed
!               iopqnsum: integer flag = 0 if qnsum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if qnsum is computed
!               r1c     : charcteristic of corresponding radial function
!                         of the first kind
!               ir1e    : exponent of corresponding radial function of
!                         the first kind
!               r1dc    : charcteristic of corresponding first
!                         derivative of the radial function of the first
!                         kind
!               ir1de   : exponent of corresponding first derivative of
!                         the radial function of the first kind
!               wront   : theoretical Wronskian
!               minacc  : minimum desired accuracy in decimal digits
!
!     output:   r2c     : characteristic of prolate
!                         radial function of the second kind
!               ir2e    : exponent of prolate radial function of the
!                         second kind
!               r2dc    : characteristic of derivative with
!                         respect to x of prolate radial function
!                         of the second kind
!               ir2de   : exponent of derivative with respect to x of
!                         prolate radial function of second kind
!               jleg    : maximum number of terms taken in qsum
!               jlegp   : maximum number of terms taken in psum
!               naccleg : estimated accuracy of r2 and r2d
!               nsubw   : subtraction error in decimal digits that
!                         occurs in forming Wronskian for calculated
!                         radial functions and their first derivatives
!               jflagl  : equal to unity if Wronskian is used to
!                         calculate more accurate values for leading
!                         coefficient and thus more accurate values
!                         for r2 and r2d; equal to zero otherwise
!               iflagq  : set equal to zero at l = m. Remains equal
!                         to zero if set function values and accuracy
!                         to zero because terms in backward qsum
!                         become too large to allow accurate results.
!                         Set equal to unity when this first does not
!                         occur.
!               iflagp  : set equal to zero at l = m. Remains equal
!                         to zero if set function values and accuracy
!                         to zero because terms in psum become too large
!                         to allow accurate results. Set equal to unity
!                         when this first does not occur.
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, dconp, dconq, dconqn, dec, dec1, dfnorm, dmfnorm, dneg, &
         dnegjf, dnew, dnewd, dold, doldd, d01, fac, psum, psump, pdsum, &
         pdsump, qdml, qndsum, qndsump, qdsum, qdsump, qml, qnsum, &
         qnsump, qsum, qsump, r1c, r1dc, r2c, r2dc, rm, spsum, spsump, &
         spdsum, spdsump, ten, termpq, test, testd, teste, testeo, &
         testm, testdm, tm, wronca, wroncb, wronc, wront, x1
    real(knd) drhor(maxdr), enr(maxd), enrneg(maxmp), fajo(lnum), &
         prx(maxpdr), pdrx(maxpdr), qdl(lnum), qdr(maxq), ql(lnum), &
         qr(maxq)
!
!  integer arrays
    dimension ifajo(lnum), iqdl(lnum), iql(lnum)
!
    jflagl = 0
    ten = 10.0e0_knd
    nfac = nex / 3
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** (nfac)
    testeo = 1.0e0_knd / teste
    iscale = 0
    dec = ten ** (-ndec - 1)
    dec1 = ten ** (-ndec - 7)
    lm2 = (l - m) / 2
    ix = l - m - 2 * lm2
    imxp = m + m + ix
    ixx = 1 - ix
    lim1 = limleg / 2 - ix
    lim2 = limdr - 1
    if(ioppsum == 0) lim2 = 0
    rm = m
    tm = rm + rm
    dconq = dec
    dconqn = dec
    dconp = dec
    dnegjf = dneg * d01
    if(m == 0) dnegjf = d01
    idnegjf = idneg + id01
    if(m == 0) idnegjf = id01
    fajo(l - m + 1) = fajo(l - m + 1) * dmfnorm * dnegjf / dfnorm
    iterm = int(log10(abs(fajo(l - m + 1))))
    fajo(l - m + 1) = fajo(l - m + 1) * ten ** (-iterm)
    ifajo(l - m + 1) = ifajo(l - m + 1) + idnegjf + idmfe - idfe + iterm
    iterm = -int(log10(abs(c * x1 * (x1 + 2.0e0_knd) * r1dc)))
    itermq = int(log10(abs(fajo(l - m + 1) * termpq / ql(l - m + 1))))
    itestq = iterm + itermq - ir1de + itermpq + ifajo(l - m + 1) - iql(l - m + 1) + ndec + 3
    itermp = int(log10(abs(fajo(l - m + 1) / (dnegjf * termpq))))
    itestp = iterm + itermp - ir1de - idnegjf - itermpq + ifajo(l - m + 1) + ndec + 3
!
!  begin calculation of series for r2
!
!  calculate d*q sum over positive n using pyramid summation
!
!  backward summation
    qsum = 1.0e0_knd
    qdsum = 1.0e0_knd
    qsump = 1.0e0_knd
    qdsump = 1.0e0_knd
    if(lm2 < 1) go to 20
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
    j = lm2
     do 10 jj = 1, lm2
     dnew = dold / (qr(j + j + imxp) * qr(j + j + imxp - 1) * enr(j))
     qsum = qsum + dnew
     if(dnew > 0.0e0_knd) qsump = qsump + dnew
     dnewd = doldd / (qdr(j + j + imxp) * qdr(j + j + imxp - 1) * enr(j))
     qdsum = qdsum + dnewd
     if(dnewd > 0.0e0_knd) qdsump = qdsump + dnewd
      if(int(log10(abs(qsum))) + iscale > itestq .and. iflagq == 0) &
        then
      r2c = 0.0e0_knd
      r2dc = 0.0e0_knd
      ir2e = 0
      ir2de = 0
      nsub = ndec
      nsubd = ndec
      jleg = j
      jlegp = 0
      naccleg = 0
      go to 180
      end if
      if(abs(qsum) > teste) then
      dnew = dnew * testeo
      dnewd = dnewd * testeo
      qsum = qsum * testeo
      qdsum = qdsum * testeo
      qsump = qsump * testeo
      qdsump = qdsump * testeo
      iscale = iscale + nfac
      end if
     if((abs(dnew / qsum) + abs(dnewd / qdsum)) < dconq) go to 20
     dold = dnew
     doldd = dnewd
     j = j - 1
10    continue
20   continue
!
!  forward summation
     iflagq = 1
     if(iscale /= 0) then
     jleg = lm2
     go to 45
     end if
    jlow = lm2 + 1
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
     do 30 j = jlow, lim1
     dnew = dold * enr(j) * qr(j + j + imxp) * qr(j + j + imxp - 1)
     qsum = qsum + dnew
     if(dnew > 0.0e0_knd) qsump = qsump + dnew
     dnewd = doldd * enr(j) * qdr(j + j + imxp) * qdr(j + j + imxp - 1)
     qdsum = qdsum + dnewd
     if(dnewd > 0.0e0_knd) qdsump = qdsump + dnewd
     if((abs(dnew / qsum) + abs(dnewd / qdsum)) < dconq) go to 40
     dold = dnew
     doldd = dnewd
30    continue
40   continue
    jleg = j
    if(jleg > lim1) jleg = lim1
45   nsqsum = 0
    if(qsum * qsump /= 0.0e0_knd) nsqsum = int(log10(abs(qsump / qsum)))
    if(nsqsum > ndec) nsqsum = ndec
    if(nsqsum < 0) nsqsum = 0
    nsqdsum = 0
    if(qdsum * qdsump /= 0.0e0_knd) nsqdsum= &
                 int(log10(abs(qdsump / qdsum)))
    if(nsqdsum > ndec) nsqdsum = ndec
    if(nsqdsum < 0) nsqdsum = 0
    qsum = qsum * ql(l - m + 1) / (fajo(l - m + 1) * termpq)
    iterm = int(log10(abs(qsum)))
    qsum = qsum * (10.0e0_knd ** (-iterm))
    iqsum = iql(l - m + 1) - ifajo(l - m + 1) - itermpq + iscale + iterm
    qdsum = qdsum * qdl(l - m + 1) / (fajo(l - m + 1) * termpq)
    iterm = int(log10(abs(qdsum)))
    qdsum = qdsum * (10.0e0_knd ** (-iterm))
    iqdsum = iqdl(l - m + 1) - ifajo(l - m + 1) - itermpq + iscale + iterm
    qdsum = qdsum * (10.0e0_knd ** (iqdsum - iqsum))
!
!  calculate d*q sum over negative n
    if(m == 0) iopqnsum = 0
    qnsum = 0.0e0_knd
    qndsum = 0.0e0_knd
    qnsump = 0.0e0_knd
    qndsump = 0.0e0_knd
    iqnsum = 0
    iqndsum = 0
    nmterm = 0
    nsqnsum = 0
    nsqndsum = 0
    if(iopqnsum == 0 .or. m == 0) go to 90
    nmterm = m
    qnsum = enrneg(m)
    qndsum = enrneg(m)
    if(ix == 1) go to 50
    qnsum = qnsum * qr(m + m - 1)
    qndsum = qndsum * qdr(m + m - 1)
50   if(qnsum > 0.0e0_knd) qnsump = qnsum
    if(qndsum > 0.0e0_knd) qndsump = qndsum
    if(m == 1) go to 80
    dold = qnsum
    doldd = qndsum
     do j = 2, m
       dnew = dold * enrneg(m - j + 1) * qr(imxp - j - j + 1) * qr(imxp - j - j + 2)
       qnsum = qnsum + dnew
       if(dnew > 0.0e0_knd) qnsump = qnsump + dnew
       dnewd = doldd * enrneg(m - j + 1) * qdr(imxp - j - j + 1) * qdr(imxp - j - j + 2)
       qndsum = qndsum + dnewd
       if(dnewd > 0.0e0_knd) qndsump = qndsump + dnewd
       dold = dnew
       doldd = dnewd
     end do
70   nsqnsum = 0
    if(qnsum * qnsump /= 0.0e0_knd) nsqnsum= &
             int(log10(abs(qnsump / qnsum)))
    if(nsqnsum > ndec) nsqnsum = ndec
    if(nsqnsum < 0) nsqnsum = 0
    nsqnsum = max(nsqnsum, nsdneg)
    nsqndsum = 0
    if(qndsum * qndsump /= 0.0e0_knd) nsqndsum= &
             int(log10(abs(qndsump / qndsum)))
    if(nsqndsum > ndec) nsqndsum = ndec
    if(nsqndsum < 0) nsqndsum = 0
    nsqndsum = max(nsqndsum, nsdneg)
80   qnsum = qnsum * qml * d01 / (fajo(l - m + 1) * termpq)
    iterm = int(log10(abs(qnsum)))
    qnsum = qnsum * (10.0e0_knd ** (-iterm))
    iqnsum = iqml + id01 - ifajo(l - m + 1) - itermpq + iterm
    qnsum = qnsum * (10.0e0_knd ** (iqnsum - iqsum))
    qndsum = qndsum * qdml * d01 / (fajo(l - m + 1) * termpq)
    iterm = int(log10(abs(qndsum)))
    qndsum = qndsum * (10.0e0_knd ** (-iterm))
    iqndsum = iqdml + id01 - ifajo(l - m + 1) - itermpq + iterm
    qndsum = qndsum * (10.0e0_knd ** (iqndsum - iqsum))
90   continue
!
!       calculate d(rho|n)*p summation
    psum = 0.0e0_knd
    pdsum = 0.0e0_knd
    ipsum = 0
    ipdsum = 0
    jlegp = 0
    nspsum = 0
    nspdsum = 0
    if(ioppsum == 0) go to 160
    psum = prx(ixx + 1) * drhor(1)
    pdsum = pdrx(ixx + 1) * drhor(1)
    dold = psum
    doldd = pdsum
    if(m /= 0 .or. ix /= 1) go to 100
    pdsum = 0.0e0_knd
    doldd = drhor(1)
100   continue
    spsum = psum
    spdsum = pdsum
    psump = 0.0e0_knd
    if(psum > 0.0e0_knd) psump = psum
    pdsump = 0.0e0_knd
    if(pdsum > 0.0e0_knd) pdsump = pdsum
    spsump = psump
    spdsump = pdsump
    testm = 1.0e0_knd
    testdm = 1.0e0_knd
    jlegpf = 1
    jlegpd = 1
     do 130 j = 2, lim2
     dnew = dold * drhor(j) * prx(j + j - ix)
     psum = psum + dnew
     if(dnew > 0.0e0_knd) psump = psump + dnew
     dnewd = doldd * drhor(j) * pdrx(j + j - ix)
     pdsum = pdsum + dnewd
      if(int(log10(abs(psum))) > itestp .and. iflagp == 0) then
      r2c = 0.0e0_knd
      r2dc = 0.0e0_knd
      ir2e = 0
      ir2de = 0
      nsub = ndec
      nsubd = ndec
      naccleg = 0
      jlegp = 0
      go to 180
      end if
     if(dnewd > 0.0e0_knd) pdsump = pdsump + dnewd
     test = abs(dnew / psum)
     testd = abs(dnewd / pdsum)
     if(test > testm .or. test == 00.0e0_knd) go to 110
     testm = test
     spsum = psum
     spsump = psump
     jlegpf = j
110    if(testd > testdm .or. testd == 00.0e0_knd) go to 120
     testdm = testd
     spdsum = pdsum
     spdsump = pdsump
     jlegpd = j
120    if(test + testd < dconp) go to 150
     dold = dnew
     doldd = dnewd
130    continue
150   jlegp = max(jlegpf, jlegpd)
    iflagp = 1
    psum = spsum
    pdsum = spdsum
    psump = spsump
    pdsump = spdsump
    ntestm = -int(log10(testm))
    ntestdm = -int(log10(testdm))
    if(ntestm > ndec) ntestm = ndec
    if(ntestdm > ndec) ntestdm = ndec
    nspsum = 0
    if(psum * psump /= 0.0e0_knd) nspsum = int(log10(abs(psump / psum)))
    if(nspsum > ndec) nspsum = ndec
    if(nspsum < 0) nspsum = 0
    nspsum = max(ndec - ntestm, nspsum, nsdrho)
    nspdsum = 0
    if(pdsum * pdsump /= 0.0e0_knd) nspdsum= &
             int(log10(abs(pdsump / pdsum)))
    if(nspdsum > ndec) nspdsum = ndec
    if(nspdsum < 0) nspdsum = 0
    nspdsum = max(ndec - ntestdm, nspdsum, nsdrho)
    psum = psum * dnegjf * termpq / fajo(l - m + 1)
    iterm = 0
    if(psum /= 0.0e0_knd) iterm = int(log10(abs(psum)))
    psum = psum * (10.0e0_knd ** (-iterm))
    ipsum = idnegjf + itermpq - ifajo(l - m + 1) + iterm
    psum = psum * (10.0e0_knd ** (ipsum - iqsum))
    pdsum = pdsum * dnegjf * termpq / fajo(l - m + 1)
    if(m /= 0) pdsum = pdsum * rm * (x1 + 1.0e0_knd) / (x1 * (x1 + 2.0e0_knd))
    iterm = 0
    if(pdsum /= 0.0e0_knd) iterm = int(log10(abs(pdsum)))
    pdsum = pdsum * (10.0e0_knd ** (-iterm))
    ipdsum = idnegjf + itermpq - ifajo(l - m + 1) + iterm
    pdsum = pdsum * (10.0e0_knd ** (ipdsum - iqsum))
160   continue
    r2c = qsum + qnsum + psum
    r2dc = qdsum + qndsum + pdsum
    nqs = 0
    if(qsum / r2c /= 0.0e0_knd) nqs = int(log10(abs(qsum / r2c)))
    nqns = 0
    if(qnsum / r2c /= 0.0e0_knd) &
           nqns = int(log10(abs(qnsum / r2c)))
    nps = 0
    if(psum / r2c /= 0.0e0_knd) &
           nps = int(log10(abs(psum / r2c)))
    nsqsum = nsqsum + nqs
    if(nsqsum < 0) nsqsum = 0
    if(nsqsum > ndec) nsqsum = ndec
    if(qsum / r2c == 0.0e0_knd) nsqsum = 0
    nsqnsum = nsqnsum + nqns
    if(nsqnsum < 0) nsqnsum = 0
    if(nsqnsum > ndec) nsqnsum = ndec
    if(qnsum / r2c == 0.0e0_knd) nsqnsum = 0
    nspsum = nspsum + nps
    if(nspsum < 0) nspsum = 0
    if(nspsum > ndec) nspsum = ndec
    if(psum / r2c == 0.0e0_knd) nspsum = 0
    nsub = max(nsqsum, nsqnsum, nspsum)
    nqds = 0
    if(qdsum / r2dc /= 0.0e0_knd) nqds = int(log10(abs(qdsum / r2dc)))
    nqnds = 0
    if(qndsum / r2dc /= 0.0e0_knd) &
           nqnds = int(log10(abs(qndsum / r2dc)))
    if(qnsum == 0.0e0_knd .and. qndsum == 0.0e0_knd) iopqnsum = 0
    npds = 0
    if(pdsum / r2dc /= 0.0e0_knd) &
           npds = int(log10(abs(pdsum / r2dc)))
    if(psum == 0.0e0_knd .and. pdsum == 0.0e0_knd) ioppsum = 0
    nsqdsum = nsqdsum + nqds
    if(nsqdsum < 0) nsqdsum = 0
    if(nsqdsum > ndec) nsqdsum = ndec
    if(qdsum / r2dc == 0.0e0_knd) nsqdsum = 0
    nsqndsum = nsqndsum + nqnds
    if(nsqndsum < 0) nsqndsum = 0
    if(nsqndsum > ndec) nsqndsum = ndec
    if(qndsum / r2dc == 0.0e0_knd) nsqndsum = 0
    nspdsum = nspdsum + npds
    if(nspdsum < 0) nspdsum = 0
    if(nspdsum > ndec) nspdsum = ndec
    if(pdsum / r2dc == 0.0e0_knd) nspdsum = 0
    nsubd = max(nsqdsum, nsqndsum, nspdsum)
    if(qnsum == 0.0e0_knd .and. qndsum == 0.0e0_knd) iopqnsum = 0
    if(psum == 0.0e0_knd .and. pdsum == 0.0e0_knd) ioppsum = 0
    wronca = r1c * r2dc * 10.0e0_knd ** (ir1e+iqsum)
    wroncb = r2c * r1dc * 10.0e0_knd ** (iqsum + ir1de)
    wronc = wronca - wroncb
    naccleg = -int(log10(abs((wronc - wront) / wront) + dec))
    if(naccleg < 0) naccleg = 0
    if(naccleg > ndec - 1) naccleg = ndec - 1
    nsubw = -int(log10(abs(wronc / wronca) + dec))
    if(nsubw < 0) nsubw = 0
    if(naccleg > 1) naccleg = naccleg + nsubw
    itest = ndec - 2 - nsubw - max(nsub, nsubd)
    if(itest < 0) itest = 0
    if(naccleg < minacc .and. naccleg < itest .and. x1 <= 0.01e0_knd) &
       then
     fac = wront / wronc
     qsum = qsum * fac
     qdsum = qdsum * fac
     qnsum = qnsum * fac
     qndsum = qndsum * fac
     psum = psum * fac
     pdsum = pdsum * fac
     r2c = qsum + qnsum + psum
     r2dc = qdsum + qndsum + pdsum
     jflagl = 1
     naccleg = itest
     end if
    if(naccleg > 3 .and. nps < (-ndec - 1) .and. npds < (-ndec - 1)) &
       ioppsum = 0
    if(naccleg > 3 .and. nqns < (-ndec - 1) .and. nqnds < (-ndec - 1)) &
      iopqnsum = 0
    if(naccleg < 0) naccleg = 0
     if(jflagl == 0) then
     nsub = max(nsub, jsub, nsdneg)
     nsubd = max(nsubd, jsub, nsdneg)
     end if
    iterm = int(log10(abs(r2c)))
    r2c = r2c * (10.0e0_knd ** (-iterm))
    ir2e = iqsum + iterm
    if(abs(r2c) >= 1.0e0_knd) go to 170
    r2c = r2c * 10.0e0_knd
    ir2e = ir2e - 1
170   continue
    iterm = int(log10(abs(r2dc)))
    r2dc = r2dc * (10.0e0_knd ** (-iterm))
    ir2de = iqsum + iterm
    if(abs(r2dc) >= 1.0e0_knd) go to 180
    r2dc = r2dc * 10.0e0_knd
    ir2de = ir2de - 1
180   continue
if (debug) then
    if(ioppsum == 1 .and. iopqnsum == 1) write(40, 190) jleg, jlegp, &
                      m, lim1, lim2, m, nsub, nsubd
190   format(8x,'r2leg: qsum, psum and qnsum series converged in ',i6, &
       ',' i6,' and ',i4,' terms; ',i6,',' i6,' and ' i4, &
       ' terms avail.',/,15x, i2,' and ',i2,' digits of sub.', &
       ' error in r2 and r2d.')
    if(ioppsum == 1 .and. iopqnsum == 0) write(40, 200) jleg, jlegp, &
                      lim1, lim2, nsub, nsubd
200   format(8x,'r2leg: qsum and psum series converged in ',i6, &
       ' and ',i6,' terms; ',i6,' and ',i6,' terms avail.',/, &
       15x, i2,' and ',i2,' digits of sub. error in r2 and r2d;', &
       ' qnsum is negligible.')
    if(ioppsum == 0 .and. iopqnsum == 1) write(40, 210) jleg, m, &
                      lim1, m, nsub, nsubd
 210   format(8x,'r2leg: qsum and qnsum series converged in ',i6, &
       ' and ',i4,' terms; ',i6,' and ',i4,' terms avail.',/, &
        15x, i2,' and ',i2,' digits of sub. error in r2 and r2d;' &
        ' psum is negligible.')
     if(ioppsum == 0 .and. iopqnsum == 0) write(40, 220) jleg, lim1, &
                         nsub, nsubd
 220   format(8x,'r2leg: qsum series converged in ',i6,' terms with ', &
        i6,' terms avail.; 'i2,' and ',i2,' digits of',/,15x, &
        'sub. error in r2 and r2d; psum and qnsum are ', &
        'negligible.')
     if(jflagl == 1) write(40, 230)
 230   format(15x,'Wronskian used to improve accuracy of the', &
        ' joining factor including dmfnorm and dneg.')
end if
    return
    end subroutine
!
!
    subroutine r2neu (l, m, c, x1, limneu, ndec, nex, maxd, maxlp, maxn, &
             maxp, minacc, enr, sneuf, sneun, ineue, sneudf, &
             sneudr, prat1, pcoefn, ipcoefn, dmfnorm, idmfe, &
             r1dc, ir1de, r2c, ir2e, r2dc, ir2de, jneu)
!
!  purpose:     To calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using the traditional expansions in terms of
!               spherical Neumann functions.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x1     : x-1
!               limneu : maximum number of terms to be taken in the
!                        series summations for r2 and r2d
!               ndec   : number of decimal digits available in
!                        real(knd)
!               nex    : maximum exponent in real(knd) arithmetic
!               maxd   : dimension of enr array
!               maxlp  : maximum  l value desired; dimension
!                        of the sneun, sneudn, ineue, and ineude arrays
!               maxn   : dimension of sneuf and sneudf arrays
!               maxp   : dimension of prat1 array
!               minacc : number of decimal digits of desired accuracy
!                        of the resulting radial functions
!               enr    : array of ratios of successive d coefficients
!               sneuf  : array of ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : array of characteristics for Neumann functions
!               ineue  : array of exponents for Neumann functions
!               sneudf : array of ratios of successive first derivatives
!                        of spherical Neumann functions of same parity
!               sneudr : array of ratios of first derivatives of Neumann
!                        functions to the corresponding functions
!               prat1  : array of ratios of successive coefficients in
!                        r2 and r2d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r2 and r2d sums that contains Neumann function
!                        of order l
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm: characteristic of the Morse-Feshbach
!                        normalization factor of the d coefficients.
!                        equal to the reciprocal of the value of the d
!                        coefficient d(n = l - m) using this
!                        normalization for the angular functions
!               idmfe  : exponent associated with dmfnorm
!               r1dc   : charcteristic of corresponding first
!                        derivative of the radial function of the first
!                        kind
!               ir1de  : exponent of corresponding first derivative of
!                        the radial function of the first kind
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the second kind
!               ir2e   : exponent of prolate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the second kind
!               jneu   : index of term where best convergence is
!                        achieved for r2 or for r2d, whichever term is
!                        larger
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, dconb, dconf, dconi, dmfnorm, dnew, dnewd, dold, doldd, &
         pcoefn, rm, rm2, r1dc, r2c, r2dc, r2dcoef, r2dtemp, r2est, &
         r2temp, r2test, sr2temp, sr2dtemp, sumcoef, sump, sumdp, &
         ten, test, testd, testdm, testm, teste, testeo, tx, txd, x1
    real(knd) enr(maxd), sneudr(maxlp), sneun(maxlp), &
         prat1(maxp), sneuf(maxn), sneudf(maxn)
!
!  integer arrays
    dimension ineue(maxlp)
!
    rm = m
    ten = 10.0e0_knd
    nfac = nex / 2
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** (nfac)
    testeo = 1.0e0_knd / teste
    iscale = 0
    dconf = ten ** (-ndec - 1)
    dconi = ten ** (ndec + 2)
    sumcoef = (ten ** (-ir1de - ineue(l + 1) - ipcoefn + idmfe))/ &
        (c * x1 * (x1 + 2.0e0_knd) * r1dc * sneun(l + 1) * pcoefn)
    r2est = abs(sumcoef * dmfnorm)
    dconb = r2est / dconi
    r2test = r2est * dconi
    r2dcoef = rm / ((x1 + 1.0e0_knd) * (x1 + 2.0e0_knd) * x1)
    rm2 = rm * 2.0e0_knd
    lm2 = (l - m) / 2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
    ix = l - m - 2 * lm2
    lim = limneu / 2 - ix
!
!  compute radial function of the second kind
!
!  backward series
    r2temp = 1.0e0_knd
    r2dtemp = 1.0e0_knd
    sump = 1.0e0_knd
    sumdp = 1.0e0_knd
    if(r2est > dconi) go to 20
    if (lm2 < 1) go to 20
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
     do 10 j = lm2, 1,-1
     jj = j + j + ix
     dnew = -dold / (sneuf(jj + m) * prat1(jj + 1) * enr(j))
     dnewd = -doldd / (sneudf(jj + m) * prat1(jj + 1) * enr(j))
     r2temp = r2temp + dnew
     r2dtemp = r2dtemp + dnewd
     if(dnew > 0.0e0_knd) sump = sump + dnew
     if(dnewd > 0.0e0_knd) sumdp = sumdp + dnewd
     if(abs(dnew / r2temp) + abs(dnewd / r2dtemp) < dconb) go to 20
     dold = dnew
     doldd = dnewd
10   continue
20   continue
!
!  forward series
    dold = 1.0e0_knd
    doldd = 1.0e0_knd
    testm = 1.0e0_knd
    testdm = 1.0e0_knd
    sr2temp = r2temp
    sr2dtemp = r2dtemp
    tx = sump
    txd = sumdp
    js = lim
    jds = lim
     do 70 j = lm2 + 1, lim
     jj = j + j + ix
     dnew = -dold * enr(j) * sneuf(jj + m) * prat1(jj + 1)
     dnewd = -doldd * enr(j) * sneudf(jj + m) * prat1(jj + 1)
     r2temp = r2temp + dnew
     r2dtemp = r2dtemp + dnewd
     if(dnew > 0.0e0_knd) sump = sump + dnew
     if(dnewd > 0.0e0_knd) sumdp = sumdp + dnewd
     test = abs(dnew / r2temp)
     testd = abs(dnewd / r2dtemp)
     if(abs(r2temp) > r2test) go to 80
     if(test < testm) go to 30
     go to 40
30    testm = test
     sr2temp = r2temp
     tx = sump
     js = j
40    continue
     if(testd < testdm) go to 50
     go to 60
50    testdm = testd
     sr2dtemp = r2dtemp
     txd = sumdp
     jds = j
60    continue
     if(test + testd < dconf) go to 90
      if(abs(r2temp) > teste) then
      r2temp = r2temp * testeo
      r2dtemp = r2dtemp * testeo
      sr2temp = sr2temp * testeo
      sr2dtemp = sr2dtemp * testeo
      sump = sump * testeo
      sumdp = sumdp * testeo
      dnew = dnew * testeo
      dnewd = dnewd * testeo
      iscale = iscale + nfac
      r2test = r2test * testeo
      tx = tx * testeo
      txd = txd * testeo
      end if
     dold = dnew
     doldd = dnewd
70    continue
80   r2temp = sr2temp
    r2dtemp = sr2dtemp
    sump = tx
    sumdp = txd
90   continue
    jneu = max(js, jds)
    naccs1 = int(log10(abs(tx / r2temp) + dconf))
    if(naccs1 < 0) naccs1 = 0
    if(naccs1 > ndec) naccs1 = ndec
    naccs2 = int(log10(abs(txd / r2dtemp) + dconf))
    if(naccs2 < 0) naccs2 = 0
    if(naccs2 > ndec) naccs2 = ndec
    jtestm = -int(log10(testm + dconf))
    if(jtestm < 0) jtestm = 0
    if(jtestm > ndec) jtestm = ndec
    jtestdm = -int(log10(testdm + dconf))
    if(jtestdm < 0) jtestdm = 0
    if(jtestdm > ndec) jtestdm = ndec
if (debug) then
    write(40, 100) j, lim, js, jtestm, naccs1, jds, jtestdm, naccs2
100   format(8x,'r2neu: numerator converged in ',i6,' terms; ', &
        i6,' terms available.',/,15x,'best r2 at ',i6,' terms', &
        ' with convergence to ',i3,' digits and sub error', &
        ' of',i3,' digits.',/,15x,'best r2d at ',i6,' terms', &
        ' with convergence to ',i3,' digits and sub. error of', &
        i3,' digits.')
end if
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
    r2c = r2temp * sneun(l + 1) * pcoefn / dmfnorm
    iterm = int(log10(abs(r2c)))
    ir2e = ineue(l + 1) + ipcoefn - idmfe + iterm + iscale
    r2c = r2c * 10.0e0_knd ** (-iterm)
    if(abs(r2c) >= 1.0e0_knd) go to 110
    r2c = r2c * 10.0e0_knd
    ir2e = ir2e - 1
110   continue
    r2dc = r2dcoef * r2c + (c * r2dtemp * sneun(l + 1) * sneudr(l + 1) * pcoefn/ &
       dmfnorm) * 10.0e0_knd ** (ineue(l + 1) + ipcoefn - idmfe+iscale-ir2e)
    iterm = int(log10(abs(r2dc)))
    ir2de = ir2e + iterm
    r2dc = r2dc * 10.0e0_knd ** (-iterm)
    if(abs(r2dc) >= 1.0e0_knd) go to 120
    r2dc = r2dc * 10.0e0_knd
    ir2de = ir2de - 1
120   continue
    return
    end subroutine
!
!
    subroutine r2eta (l, m, c, x1, eta, nee, incnee, limeta, ndec, nex, maxd, &
             maxlp, maxn, maxp, minacc, wm, enr, sneuf, sneun, &
             ineue, sneudf, sneudr, pdratt, pratb, pratt, &
             pcoefn, ipcoefn, pdcoefn, ipdcoefn, r1c, ir1e, r1dc, &
             ir1de, naccmax, naccr, r2c, ir2e, r2dc, ir2de, &
             nacceta, nacciop, jeta, iopnee, neemark, naccd, &
             naccn, naccnmax, naccns)
!
!  purpose:     To calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using an expansion of spherical Neumann
!               functions.
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               x1      : x-1
!               eta     : value for eta used in calculation
!               nee     : index in the array of eta values in the main
!                         program that corresponds to the value of eta
!                         used in r2eta calculations
!               incnee  : increment in nee
!               limeta  : maximum number of terms available in the sums
!                         for r2 and r2d
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent in real(knd) arithmetic
!               maxd    : dimension of enr array
!               maxlp   : maximum  l value desired; dimension
!                         of the sneun, sneudr, and ineue arrays
!               maxn    : dimension of sneuf and sneudf arrays
!               maxp    : dimension of pdratt, pratb, and pratt arrays
!               minacc  : minimum number of accurate decimal digits
!                         that are requested
!               wm      : 1 - eta*eta = sin(nee)*sin(nee)
!               enr     : array of ratios of successive d coefficients
!               sneuf   : array of ratios of successive spherical
!                         Neumann functions of the same parity
!               sneun   : array of characteristics for Neumann functions
!               ineue   : array of exponents corresponding to sneun
!               sneudf  : array of ratios of successive first
!                         derivatives of spherical Neumann functions of
!                         the same parity
!               sneudr  : array of ratios of first derivatives of the
!                         spherical Neumann functions to the
!                         corresponding functions
!               pdratt  : array of ratios of successive first
!                         derivatives of the associated Legendre
!                         functions of the first kind of the same parity
!                         (used in numerator series)
!               pratb   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in denominator series)
!               pratt   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in numerator series)
!               pcoefn  : characteristic of the ratio of the numerator
!                         and denominator associated Legendre functions
!                         of the first kind of order m and degree l
!               ipcoefn : exponent corresponding to pcoefn
!               pdcoefn : characteristic of the ratio of the first
!                         derivative of the associated Legendre function
!                         of the first kind in the numerator and the
!                         associated Legendre function of the first kind
!                         in the denominator, both of order m and
!                         degree l
!               ipdcoefn: exponent corresponding to pdcoefn
!               r1c     : characteristic of the radial function of the
!                         first kind (calculated in r1bes)
!               irie    : exponent corresponding to r1c
!               r1dc    : characteristic of the first derivative with
!                         respect to x of the radial function of the
!                         first kind (calculated in r1bes)
!               ir1de   : exponent corresponding to r1dc
!               naccmax : maximum accuracy (in decimal digits) obtained
!                         for the current value of l from previous
!                         r2eta calculations
!               naccr   : accuracy of radial functions calculated for
!                         this value of l earlier using other methods.
!                         if this is the only method used the default
!                         value of minacc is used for naccr
!
!     output:   r2c     : characteristic of the prolate radial function
!                         of the second kind
!               ir2e    : exponent of the prolate radial function of the
!                         second kind
!               r2dc    : characteristic of the first derivative with
!                         respect to x of the prolate radial function
!                         of the second kind
!               ir2de   : exponent corresponding to r2dc
!               nacceta : estimated number of accurate decimal digits in
!                         r2 and r2d. computed from the wronskian if
!                         nacciop = 0; estimated from the series if
!                         nacciop = 1
!               nacciop : integer flag = 1 if the denominator in the
!                         expressions for r2 and r2d is computed using
!                         the theoretical wronskian and known values for
!                         r1 and r1d. nacciop = 0 otherwise
!               jeta    : maximum number of terms taken in the numerator
!                         sums for r2 and r2d
!               iopnee  : integer flag = 0 if none of the values used
!                         for eta for the present l has led to an
!                         estimated accuracy nacceta of at least 3
!                         decimal digits (else iopnee = 1) and the
!                         subtraction error in the numerator series
!                         for both r2 and r2d is not greater than
!                         ndec-minacc digits (else iopnee =2)
!               neemark : index for the eta value in the array storing
!                         eta values in the main program that
!                         corresponds to the last eta value used for
!                         which iopnee = 0
!               naccd   : estimated accuracy of denominator sum
!               naccn   : estimated accuracy of numerator sums
!               naccns  : larger of the subtraction errors in calculating
!                         r2 and r2d
!
!     input/
!     output:   naccnmax: maximum accuracy (in decimal digits) obtained
!                         for the numerator series for the current value
!                         of l from all previous r2eta calculations
!                         (input) and including the curent r2eta
!                         calculation (output)
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, dcon, dconi, dec, denom, dnew, dnewd1, dnewd2, dnewsum, &
         dnewdsum1, dnewdsum2, dold, doldd1, doldd2, eta, etas, &
         factor, pcoefn, pdcoefn, reld12, rm, rm2, r1c, r1dc, r2c, r2dc, &
         r2dcoef1, r2dcoef2, r2dtemp, r2dtemp1, r2dtemp2, r2est, &
         r2temp, r2test, sr2temp, sr2dtemp1, sr2dtemp2, sumcoef, &
         sumdnp1, sumdnp2, sumdp, sumnp, ten, test, testdm1, testdm2, &
         testd1, testd2, testm, teste, testeo, tx, txd1, txd2, wm, &
         wronc, wronca, wroncb, wroncm, wront, xet, xets, x1
    real(knd) enr(maxd), sneudr(maxlp), sneun(maxlp), pratb(maxp), &
         pratt(maxp), pdratt(maxp), sneuf(maxn), sneudf(maxn)
!
!  integer arrays
    dimension ineue(maxlp)
!
    ten = 10.0e0_knd
    nfac = nex / 2
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** (nfac)
    testeo = 1.0e0_knd / teste
    iscale = 0
    dec = ten ** (-ndec - 2)
    rm = m
    dcon = ten ** (-ndec - 2)
    dconi = ten ** (ndec + 2)
    etas = eta * eta
    xet = sqrt(x1 * (x1 + 2.0e0_knd) + etas)
    xets = xet * xet
    factor = 1.0e0_knd
     if(-ir1de - ineue(l + 1) - ipcoefn > nex - ndec - 30) then
     iscale = nex - ndec - 30
     factor = ten ** (-iscale)
     sumcoef = ten ** (-ir1de - ineue(l + 1) - ipcoefn - iscale)/ &
         (c * x1 * (x1 + 2.0e0_knd) * r1dc * sneun(l + 1) * pcoefn)
     else
     sumcoef = (ten ** (-ir1de - ineue(l + 1) - ipcoefn))/ &
        (c * x1 * (x1 + 2.0e0_knd) * r1dc * sneun(l + 1) * pcoefn)
     end if
    r2dcoef1 = -eta * wm / (xets * xet)
    r2dcoef2 = c * (x1 + 1.0e0_knd) / xet
    reld12 = (r2dcoef2 / r2dcoef1) * sneudr(l + 1) * (pcoefn/ &
        pdcoefn) * ten ** (ipcoefn - ipdcoefn)
    rm2 = m + m
    lm2 = (l - m) / 2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
    ix = l - m - 2 * lm2
    lim = limeta / 2 - ix
!
!  compute radial function of the second kind and its first derivative
!
!  backward series for denominator
    denom = 1.0e0_knd
    sumdp = 1.0e0_knd
    if (lm2 < 1) go to 20
    dold = 1.0e0_knd
     do 10 j = lm2, 1,-1
     jj = j + j + ix
     dnew = dold / (pratb(jj + 1) * enr(j))
     denom = denom + dnew
     if(dnew > 00.0e0_knd) sumdp = sumdp + dnew
     if(abs(dnew / denom) < dec) go to 20
     dold = dnew
10    continue
20   continue
!
!  forward series for denominator
    dold = 1.0e0_knd
     do 30 j = lm2 + 1, lim
     jj = j + j + ix
     dnew = dold * enr(j) * pratb(jj + 1)
     denom = denom + dnew
     if(dnew > 00.0e0_knd) sumdp = sumdp + dnew
     if(abs(dnew / denom) < dec) go to 40
     dold = dnew
30    continue
40   continue
    jden = j
    numsub = 0
    if(sumdp /= 0.0e0_knd) numsub = int(log10(abs(sumdp / denom)))
    if(numsub < 0) numsub = 0
    if(numsub > ndec) numsub = ndec
    naccd = ndec - max(2, int(log10(c))) - numsub
    if(naccd < 0) naccd = 0
    r2est = abs(sumcoef * denom)
    r2test = r2est * dconi
!
!  backward series for numerator
    dold = factor
    doldd1 = factor
    doldd2 = factor * reld12
    r2temp = dold
    sumnp = r2temp
    r2dtemp1 = 0.0e0_knd
    sumdnp1 = 0.0e0_knd
    r2dtemp2 = doldd2
    sumdnp2 = 0.0e0_knd
    if(doldd2 > 0.0e0_knd) sumdnp2 = doldd2
    if(l /= 0) r2dtemp1 = factor
    if(l /= 0) sumdnp1 = r2dtemp1
    if(lm2 == 0) go to 60
     do 50 j = lm2, 1,-1
     jj = j + j + ix
     dnew = -dold / (sneuf(jj + m) * pratt(jj + 1) * enr(j))
     dnewd1 = -doldd1 / (sneuf(jj + m) * pdratt(jj + 1) * enr(j))
     dnewd2 = -doldd2 / (sneudf(jj + m) * pratt(jj + 1) * enr(j))
     r2temp = r2temp + dnew
     r2dtemp1 = r2dtemp1 + dnewd1
     r2dtemp2 = r2dtemp2 + dnewd2
     if(dnew > 0.0e0_knd) sumnp = sumnp + dnew
     if(dnewd1 > 0.0e0_knd) sumdnp1 = sumdnp1 + dnewd1
     if(dnewd2 > 0.0e0_knd) sumdnp2 = sumdnp2 + dnewd2
     if(abs(dnew / r2temp) + abs(dnewd1 / r2dtemp1) + abs(dnewd2 / r2dtemp2) &
        < dcon) go to 60
     dold = dnew
     doldd1 = dnewd1
     doldd2 = dnewd2
50    continue
60   continue
     if(m == 0 .and. jj == 2) then
     r2dtemp1 = r2dtemp1 - dnewd1
     if(dnewd1 > 0.0e0_knd) sumdnp1 = sumdnp1 - dnewd1
     end if
!
!  forward series for numerator
    dold = factor
    doldd1 = factor
    doldd2 = factor * reld12
    test = 1.0e0_knd
    testd1 = 1.0e0_knd
    testd2 = 1.0e0_knd
    testm = 1.0e0_knd
    testdm1 = 1.0e0_knd
    testdm2 = 1.0e0_knd
    js = lm2
    jds1 = lm2
    jds2 = lm2
    sr2temp = r2temp
    sr2dtemp1 = r2dtemp1
    sr2dtemp2 = r2dtemp2
    tx = sumnp
    txd1 = sumdnp1
    txd2 = sumdnp2
    kount = 0
    kountd1 = 0
    kountd2 = 0
    dnewsum = 0.0e0_knd
    dnewdsum1 = 0.0e0_knd
    dnewdsum2 = 0.0e0_knd
     do 110 j = lm2 + 1, lim - 1
     jj = j + j + ix
     kount = kount + 1
     kountd1 = kountd1 + 1
     kountd2 = kountd2 + 1
     dnew = -dold * enr(j) * sneuf(jj + m) * pratt(jj + 1)
     dnewd1 = -doldd1 * enr(j) * sneuf(jj + m) * pdratt(jj + 1)
     dnewd2 = -doldd2 * enr(j) * sneudf(jj + m) * pratt(jj + 1)
     if((dnew / dold) <= 0.0e0_knd .or. kount == 100) go to 70
     dnewsum = dnewsum + dnew
     go to 80
70    r2temp = r2temp + dnewsum
     kount = 0
     if(abs(r2temp) > r2test) go to 120
     if(dnewsum > 0.0e0_knd) sumnp = sumnp + dnewsum
     if(dnewsum /= 0.0e0_knd) test = abs(dnewsum / r2temp)
     dnewsum = dnew
     if(test >= testm) go to 80
     testm = test
     sr2temp = r2temp
     js = j
     tx = sumnp
80    if((dnewd1 / doldd1 <= 0.0e0_knd) .or. (kountd1 == 100)) &
         go to 85
     dnewdsum1 = dnewdsum1 + dnewd1
     go to 90
85    r2dtemp1 = r2dtemp1 + dnewdsum1
     kountd1 = 0
     if(dnewdsum1 > 0.0e0_knd) sumdnp1 = sumdnp1 + dnewdsum1
     if(dnewdsum1 /= 0.0e0_knd) testd1 = abs(dnewdsum1 / r2dtemp1)
     dnewdsum1 = dnewd1
     if(testd1 >= testdm1) go to 90
     testdm1 = testd1
     sr2dtemp1 = r2dtemp1
     jds1 = j
     txd1 = sumdnp1
90    if(((dnewd2 / doldd2) <= 0.0e0_knd) .or. (kountd2 == 100)) &
         go to 95
     dnewdsum2 = dnewdsum2 + dnewd2
     go to 100
95    r2dtemp2 = r2dtemp2 + dnewdsum2
     kountd2 = 0
     if(dnewdsum2 > 0.0e0_knd) sumdnp2 = sumdnp2 + dnewdsum2
     if(dnewdsum2 /= 0.0e0_knd) testd2 = abs(dnewdsum2 / r2dtemp2)
     dnewdsum2 = dnewd2
     if(testd2 >= testdm2) go to 100
     testdm2 = testd2
     sr2dtemp2 = r2dtemp2
     jds2 = j
     txd2 = sumdnp2
100    if(test + testd1 + testd2 < dcon) go to 130
      if(abs(r2temp) > teste) then
      r2temp = r2temp * testeo
      r2dtemp1 = r2dtemp1 * testeo
      r2dtemp2 = r2dtemp2 * testeo
      sr2temp = sr2temp * testeo
      sr2dtemp1 = sr2dtemp1 * testeo
      sr2dtemp2 = sr2dtemp2 * testeo
      tx = tx * testeo
      txd1 = txd1 * testeo
      txd2 = txd2 * testeo
      sumnp = sumnp * testeo
      sumdnp1 = sumdnp1 * testeo
      sumdnp2 = sumdnp2 * testeo
      dnew = dnew * testeo
      dnewd1 = dnewd1 * testeo
      dnewd2 = dnewd2 * testeo
      dnewsum = dnewsum * testeo
      dnewdsum1 = dnewdsum1 * testeo
      dnewdsum2 = dnewdsum2 * testeo
      iscale = iscale + nfac
      r2test = r2test * testeo
      end if
     dold = dnew
     doldd1 = dnewd1
     doldd2 = dnewd2
110    continue
120   r2temp = sr2temp
    r2dtemp1 = sr2dtemp1
    r2dtemp2 = sr2dtemp2
    go to 135
130   tx = sumnp
    txd1 = sumdnp1
    txd2 = sumdnp2
    js = j
    jds1 = j
    jds2 = j
135   continue
    jmax = j
    jeta = max(js, jds1, jds2, jden)
    jds = max(jds1, jds2)
    jtestm = ndec
    if(testm /= 0.0e0_knd) jtestm = -int(log10(testm))
    if(jtestm < 0) jtestm = 0
    if(jtestm > ndec) jtestm = ndec
    jtestdm1 = ndec
     if(testdm1 /= 0.0e0_knd) then
     jtestdm1 = -int(log10(testdm1))
     end if
    if(jtestdm1 < 0) jtestdm1 = 0
    if(jtestdm1 > ndec) jtestdm1 = ndec
    jtestdm2 = ndec
     if(testdm2 /= 0.0e0_knd) then
     jtestdm2 = -int(log10(testdm2))
     end if
    if(jtestdm2 < 0) jtestdm2 = 0
    if(jtestdm2 > ndec) jtestdm2 = ndec
    naccns1 = 0
     if(abs(tx * r2temp) /= 0.0e0_knd) then
     naccns1 = int(log10(abs(tx / r2temp)))
     end if
    if(naccns1 < 0) naccns1 = 0
    if(naccns1 > ndec) naccns1 = ndec
    r2dtemp = r2dtemp1 + r2dtemp2
    naccns2 = 0
     if(abs((txd1 + txd2) * r2dtemp) /= 0.0e0_knd) then
     naccns2 = int(log10(abs((txd1 + txd2) / r2dtemp)))
     end if
    if(naccns2 < 0) naccns2 = 0
    if(naccns2 > ndec) naccns2 = ndec
     if(abs(r2dtemp1 * r2dtemp2) /= 0.0e0_knd) then
     icord = int(log10(abs(r2dtemp1 / r2dtemp2)))
     if(icord > 0) jtestdm2 = jtestdm2 + icord
     if(icord < 0) jtestdm1 = jtestdm1 - icord
     jtestdm = min(jtestdm1, jtestdm2)
     end if
    if(abs(r2dtemp1) == 0.0e0_knd) jtestdm = jtestdm2
    if(abs(r2dtemp2) == 0.0e0_knd) jtestdm = jtestdm1
    ncorr = max(0,-int(log10(x1) - 0.001e0_knd)) + 1
    if(x1 >= 0.05e0_knd) ncorr = 1
     if(x1 <= 0.02e0_knd) then
     naccn1 = min(jtestm - 1, ndec - ncorr, ndec - naccns1)
     naccn2 = min(jtestdm - 1, ndec - ncorr, ndec - naccns2)
     else
     naccn1 = min(jtestm - 1, ndec - ncorr, ndec - naccns1 - 1)
     naccn2 = min(jtestdm - 1, ndec - ncorr, ndec - naccns2 - 1)
     end if
    naccn = min(naccn1, naccn2)
    if(naccn > ndec - 2) naccn = ndec - 2
    if(naccn < 0) naccn = 0
    naccns = max(naccns1, naccns2)
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
    r2c = r2temp * sneun(l + 1) * pcoefn / denom
    iterm = 0
    if(r2c /= 0.0e0_knd) iterm = int(log10(abs(r2c)))
    ir2e = ineue(l + 1) + ipcoefn + iterm + iscale
    r2c = r2c * 10.0e0_knd ** (-iterm)
    r2dc = r2dcoef1 * r2dtemp * sneun(l + 1) * pdcoefn / denom
    iterm = 0
    if(r2dc /= 0.0e0_knd) iterm = int(log10(abs(r2dc)))
    ir2de = ineue(l + 1) + ipdcoefn + iterm + iscale
    r2dc = r2dc * 10.0e0_knd ** (-iterm)
if (debug) then
    write(40, 140) jmax, jden, lim, js, jtestm, naccns1, jds, jtestdm, &
           naccns2, naccn, naccd
140   format(8x,'r2eta: numerator, denominator converged in ', &
        i6,' ,',i6,' terms; ',i6,' terms available.',/, &
        15x,'best r2 at ',i6,' terms with convergence to',i3, &
        ' digits;',i3,' digits subtr. error.',/,15x, &
        'best r2d at ',i6,' terms with convergence to ',i3, &
        ' digits;',i3,' digits subtr. error.',/,15x, &
        'estimated numerator and denominator accuracy is ',i4, &
        ' and',i4,' digits.')
end if
    wronca = r1c * r2dc * (10.0e0_knd ** (ir1e+ir2de))
    wroncb = r2c * r1dc * (10.0e0_knd ** (ir2e+ir1de))
    wronc = wronca - wroncb
    wront = 1.0e0_knd / (c * x1 * (x1 + 2.0e0_knd))
    wroncm = max(abs(wronca), abs(wroncb))
    nsubw =-int(log10(abs(wronc / wroncm) + dec))
    nacceta = -int(log10(abs((wronc - wront) / wront) + dec))
    if(nacceta < 0) nacceta = 0
    if(nacceta > ndec - 1) nacceta = ndec - 1
    nacciop = 0
    if(nacceta < minacc .and. naccn - nsubw > naccd .and. naccn > &
      nacceta .and. (naccn >= min(naccr, 4) .or. (jtestm >= ndec - 1 .and. &
      jtestdm >= ndec - 1))) nacciop = 1
    if(jtestm < 5 .or. jtestdm < 5) nacciop = 0
    if(naccn - nsubw <= naccr) nacciop = 0
    if(nacciop == 0) go to 160
if (debug) then
    write(40, 150)
150   format(15x,'denominator calculated using wronskian.')
end if
    nacceta = naccn - nsubw
    r2c = r2c * wront / wronc
    iterm = 0
    if(r2c /= 00.0e0_knd) iterm = int(log10(abs(r2c)))
    ir2e = ir2e + iterm
    r2c = r2c * 10.0e0_knd ** (-iterm)
    r2dc = r2dc * wront / wronc
    iterm = 0
    if(r2dc /= 00.0e0_knd) iterm = int(log10(abs(r2dc)))
    ir2de = ir2de + iterm
    r2dc = r2dc * 10.0e0_knd ** (-iterm)
160   if(abs(r2c) >= 1.0e0_knd) go to 170
    r2c = r2c * 10.0e0_knd
    ir2e = ir2e - 1
170   continue
    if(abs(r2dc) >= 1.0e0_knd) go to 180
    r2dc = r2dc * 10.0e0_knd
    ir2de = ir2de - 1
180   continue
    naccnmaxp = naccnmax
    if(naccnmax > naccn .and. naccnmax >= 3) iopnee = 1
    if(naccn > naccnmax) naccnmax = naccn
    if(ndec - naccns < naccr + 2) iopnee = 2
    if(nacciop == 0 .and. naccn > nacceta) naccnmax = max(naccnmaxp, &
      nacceta)
     if(nacceta < 3 .and. naccmax == 0 .and. ndec - naccns >= &
       naccr) then
     iopnee = 0
     neemark = nee
     end if
    if(ndec - naccns >= naccr + 2) iopnee = 0
    if(jtestm == ndec .and. jtestdm == ndec .and. naccn <= max(naccr, 5)) &
      iopnee = 2
    if(naccns1 == ndec .or. naccns2 == ndec) iopnee = 2
    return
    end subroutine
!
!
    subroutine geteig (l, m, c, eig2, eig3, eig4, eig5, eigval)
!
!  purpose:     To calculate an estimate of the eigenvalue.
!
!  parameters:
!
!     input:    l        : l
!               m        : m
!               c        : c
!               eig2-eig5: previous eigenvalues for l-4, l-3, l-2, and
!                          l-1
!
!     output:   eigval   : estimate of the eigenvalue
!
!!    use param
!
!  real(knd) scalars
    real(knd) c, csq, eigval, eig2, eig3, eig4, eig5, lam1, lam2, &
         lam3, lam4, q, q2, q4, r, rl, rl2, rm, rm2, sm
!
    csq = c * c
    rl = l
    rl2 = l + l
    rm = m
    rm2 = m + m
    r = l - m
    q = 2 * (l - m) + 1
    q2 = q * q
    q4 = q2 * q2
    sm = m * m
!
!  special case for small m and large c
    if(l == (m + 4) .and. c > 8.0e0_knd .and. m < 3) go to 60
!
!  if previous values have been computed use these to determine
!  the next value
    if(l > (m + 3)) go to 10
!
!  use expansion in terms of c**2 for low c, and
!  expansion in terms of c for large c (c>8)
    if(c > 8.0e0_knd) go to 60
    if(c > 6.0e0_knd .and. m < 4) go to 60
    if(c > 3.0e0_knd .and. l == 0) go to 60
    if(c > 4.0e0_knd .and. l == 1) go to 60
    if(c > 5.0e0_knd .and. l == 2) go to 60
    if(c > 5.0e0_knd .and. l > (m + 1)) go to 10
!
!  compute coefficients for c**2 expansion
    lam1 = rl * (rl + 1.0e0_knd)
    lam2=.5e0_knd * (1.0e0_knd - (rm2 - 1.0e0_knd) * (rm2 + 1.0e0_knd)/ &
       ((rl2 - 1.0e0_knd) * (rl2 + 3.0e0_knd)))
    lam3=.5e0_knd * (rl - rm - 1.0e0_knd) * (rl - rm) * (rl + rm - 1.0e0_knd)* &
       (rl + rm) / ((rl2 - 3.0e0_knd) * (rl2 - 1.0e0_knd) * (rl2 - 1.0e0_knd) &
       *(rl2 - 1.0e0_knd) * (rl2 + 1.0e0_knd)) - 0.5e0_knd * (rl - rm+ &
       1.0e0_knd) * (rl - rm + 2.0e0_knd) * (rl + rm + 1.0e0_knd) * (rl+ &
       rm + 2.0e0_knd) / ((rl2 + 1.0e0_knd) * (rl2 + 3.0e0_knd) ** 3* &
       (rl2 + 5.0e0_knd))
    lam4 = (4.0e0_knd * rm * rm - 1.0e0_knd) * ((rl - rm + 1.0e0_knd)* &
       (rl - rm + 2.0e0_knd) * (rl + rm + 1.0e0_knd) * (rl + rm + 2.0e0_knd) &
       /((rl2 - 1.0e0_knd) * (rl2 + 1.0e0_knd) * (rl2 + 3.0e0_knd) ** 5* &
       (rl2 + 5.0e0_knd) * (rl2 + 7.0e0_knd)) - (rl - rm - 1.0e0_knd)* &
       (rl - rm) * (rl + rm - 1.0e0_knd) * (rl + rm) / ((rl2 - 5.0e0_knd)* &
       (rl2 - 3.0e0_knd) * (rl2 - 1.0e0_knd) ** 5 * (rl2 + 1.0e0_knd)* &
       (rl2 + 3.0e0_knd)))
    eigval = lam1 + csq * (lam2 + csq * (lam3 + lam4 * csq))
    return
10   continue
!
!  pick the correct expansion (first through third order)
    if(l > (3 + m)) go to 50
    if(l > (2 + m)) go to 30
!
!  first order
    eigval = 2.0e0_knd * eig5 - eig4
20   return
!
!  second order
30   eigval = 3.0e0_knd * eig5 - 3.0e0_knd * eig4 + eig3
40   return
!
!  third order
50   eigval = 4.0e0_knd * eig5 - 6.0e0_knd * eig4 + 4.0e0_knd * eig3 - eig2
    return
60   continue
    ic = c
!
!  if m>6 the eigenvalues are very regularly spaced
    if(m > 6 .and. l > m + 1) go to 10
!
!  n*(n+1) behavior reached?
    if(m < (10 + ic)) go to 70
    eigval = rl * (rl + 1.0e0_knd)
    return
!
!  compute eigenvalue estimate with asymptotic expansion
70   eigval = q * c + sm - 0.125e0_knd * (q2 + 5.0e0_knd) - q * (q2 - 32.0e0_knd * sm+ &
        11.0e0_knd) / (64.0e0_knd * c) - (5.0e0_knd * (q4 + 26.0e0_knd * q2+ &
        21.0e0_knd) - 384.0e0_knd * sm * (q2 + 1.0e0_knd)) / (1024.0e0_knd* &
        csq) - (q * ((33.0e0_knd * q4 + 1594.0e0_knd * q2 + 5621.0e0_knd)/ &
        (128.0e0_knd * 128.0e0_knd) - sm * (37.0e0_knd * q2 + 167.0e0_knd)/ &
        (128.0e0_knd) + sm * sm / 8.0e0_knd)) / (csq * c) - ((63.0e0_knd * q2* &
        q4 + 4940.0e0_knd * q4 + 43327.0e0_knd * q2 + 22470.0e0_knd)/ &
        (256.0e0_knd * 256.0e0_knd) - sm * (115.0e0_knd * q4 + 1310.0e0_knd &
        *q2 + 735.0e0_knd) / 512.0e0_knd + 3.0e0_knd * sm * sm * (q2+ &
        1.0e0_knd) / 8.0e0_knd) / (csq * csq)
    return
    end subroutine
!
!
    subroutine conver (l, m, c, limd, blist, glist, eig1, eig3, eig4, &
              ndec, maxd, eigval, eig5, enr, ienr)
!
!  purpose:     To determine a converged eigenvalue using the
!               boukwamp method.
!  parameters:
!
!     input:    l     : l
!               m     : m
!               c     : c
!               limd  : number of enr values computed
!               blist : array of coefficients used in recursion relation
!               glist : array of coefficients used in recursion relation
!               eig1  : previous eigenvalue for l-4
!               eig3  : previous eigenvalue for l-2
!               eig4  : previous eigenvalue for l-1
!               ndec  : number of decimal digits for real(knd)
!               maxd  : dimension of enr,blist,glist arrays
!               eigval: estimated value of the eigenvalue
!
!     output:   eigval: converged eigenvalue
!               eig5  : set equal to eigval
!               enr   : array of scaled ratios of successive d
!                       coefficients
!               ienr  : index n of last d coefficient ratio used
!                       in computing first term in the denominator
!                       of the eigenvalue correction
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, cll, clu, cora, corb, csq, de, dec, dl, eig1, eig3, eig4, &
         eig5, eigdec, eigval, enrc, fl
    real(knd) blist(maxd), enr(maxd), glist(maxd)
!
    csq = c * c
    dec = 10.0e0_knd ** (-ndec - 1)
    eigdec = 10.0e0_knd ** (-ndec + 1)
    if(l == m .or. l == m + 1) imax = 1
!
!  set the original eigenvalue spacing estimate. this is arbitrary
    if(eigval < eig4) eigval = eig4
    if(l > m) cll = eig4
    if(l == (m + 2) .or. l == (m + 3)) clu = eigval + 0.5e0_knd * (eigval - eig3)
    if(l > (m + 3)) clu = eigval + 0.5e0_knd * (eig3 - eig1)
    lm2 = (l - m) / 2
    limdb = 2 * ienr + 50
    if(l == m .or. l == m + 1) limdb = 2 * ienr
    if(limdb > limd) limdb = limd
!
!  begin Bouwkamp procedure
    fl = eigval
    jnde = 0
    ix = l - m - 2 * lm2
    ifc = 1
    lim2 = limdb / 2 - ix
    iglim = lim2 + 1
    irio = lm2 + 1
    iw1 = lm2 + 2
40   enr(1) = eigval - glist(1)
    if(lm2 < 1) go to 60
!
!  evaluate the continued fraction
     do 50 i = 1, lm2
     enr(i + 1) = -blist(i) / enr(i) - glist(i + 1) + eigval
50    continue
60   enr(lim2) = -blist(lim2) / (glist(iglim) - eigval)
    iw15 = lim2 - 1
    ip = iw1 + iw15
    if(iw15 < iw1) go to 80
!
!  evaluate the continued fraction
     do 70 i = iw1, iw15
     ipi = ip - i
     enr(ipi) = -blist(ipi) / (glist(ipi + 1) - eigval + enr(ipi + 1))
70    continue
80   if(ifc == 50) go to 130
    enrc = -blist(irio) / (glist(irio + 1) - eigval + enr(irio + 1))
    de = enrc * enrc / blist(irio)
    corb = de
    if(lim2 < iw1) go to 100
!
!  compute first sum in the denominator of the correction
     do 90 i = iw1, lim2
     de = enr(i) * enr(i) / blist(i) * de
     corb = corb + de
     if(abs(de / corb) < dec .and. (l == m .or. l == m + 1 .or. &
       i > ienr - 20)) go to 100
90    continue
100   if((l == m .or. l == m + 1) .and. i > imax) imax = i
    if(l /= m .and. l /= m + 1 .and. i > ienr) ienr = i
    de = 1.0e0_knd
    cora = de
    if(lm2 < 1) go to 120
!
!  compute second term in the denominator of the correction
     do 110 i = 1, lm2
     de = blist(irio - i) / (enr(irio - i) * enr(irio - i)) * de
     cora = cora + de
     if(abs(de / cora) < dec) go to 120
110    continue
!
!  compute the correction to the eigenvalue
120   dl = (enrc - enr(irio)) / (cora + corb)
    eigval = dl + eigval
!
!  eigenvalue accurate enough?
    if(abs(dl / eigval) < eigdec) go to 130
    ifc = ifc + 1
    if(ifc < 50) go to 40
130   continue
!
!  is the eigenvalue the correct one
!  if not then modify the original guess
!
    if(l == m .or. l == (m + 1)) go to 180
    if(eigval > cll) go to 140
!
!  converged to next lower eigenvalue of the same parity
    cll = fl
    go to 170
140   if(eigval < clu) go to 180
!
!  converged to the next higher eigenvalue of the same parity
    clu = fl
!
!  eigenvalue is now somewhere within the ranges established
!  above, repeat entire Bouwkamp procedure for a value from this
!  range (the mid-point)
170   eigval = 0.5e0_knd * (cll + clu)
    fl = eigval
    ifc = 1
    jnde = jnde + 1
!
!  too many modifications are being made
!  error in the routine
    if(jnde == 50) go to 210
    go to 40
180   eig5 = eigval
    if(l == m .or. l == m + 1) ienr = imax
!
!  calculate the d coefficient ratios (enr)
    enr(1) = eigval - glist(1)
     do 190 i = 1, lm2 - 1
     enr(i + 1) = -blist(i) / enr(i) - glist(i + 1) + eigval
190    continue
    lim2 = limd / 2 - ix
    enr(lim2) = -blist(lim2) / (glist(lim2 + 1) - eigval)
    ilow = lm2 + 1
    iupp = lim2 - 1
    ip = ilow + iupp
     do 200 i = ilow, iupp
     ipi = ip - i
     enr(ipi) = -blist(ipi) / (glist(ipi + 1) - eigval + enr(ipi + 1))
200    continue
    return
!
!  error printout
210 continue
if (output) then
    write(20, 220) l, c, m
end if
if (warn) then
    write(60, 220) l, c, m
end if
220   format(1x,'error in eigenvalue at l= ',i5, 2x,'c = ',e23.14, &
        2x,'m= ',i5)
    return
    end subroutine
!
!
    subroutine dnorm (l, m, c, ndec, nex, limd, maxd, enr, sgn, d01, id01, &
             dmfnorm, idmfe, dmlmf, idmlmfe, dmsnorm, idmse, &
             dmlms, idmlmse, jnorm, jsub)
!
!  purpose:     To compute d coefficient ratios from n values and to
!               calculate the normalization of the d coefficients.
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent available in real(knd)
!               limd    : approximately twice the maximum number
!                         of terms available to be taken in the sum
!               maxd    : dimension of enr array
!               enr     : array of ratios of scaled d coefficients
!
!     output:   enr     : array of ratios of d coefficients.
!                         enr(i) = ratio of the d coefficient with
!                         subscript 2*i+ix to the d coefficient with
!                         subscript 2*(i-1)+ix. Here ix =0 when l-m is
!                         even and ix=1 when l-m is odd.
!                         If the user needs the d coefficent ratios,
!                         they are available below right before
!                         statement 20.
!               sgn     : sign of the d coefficient for n=l-m
!               d01     : characteristic of ratio of first d
!                         coefficient (either d0 or d1, depending on
!                         whether l-m is even or odd) to the d
!                         coefficient of order equal to l-m
!               id01    : exponent corresponding to d01
!               dmfnorm : characteristic of the Morse-Feshbach
!                         normalization factor of the d coefficients.
!                         equal to the reciprocal of the value of the
!                         d coefficient d(n = l - m) using this
!                         normalization for the angular functions
!               idmfe   : exponent associated with dmfnorm
!               dmlmf   : characteristic of the d coefficient with
!                         subscript l - m for the Morse and Feshbach
!                         normalization
!               idmlmfe : exponent associated with dmlmf
!               dmsnorm : characteristic of Meixner-Schafke normalization
!                         sum of the d coefficients. Equal to the
!                         reciprocal of the square of the d coefficient
!                         d(n = l - m) using this normalization for the
!                         angular functions
!               idmse   : exponent associated with dmsnorm
!               dmlms   : characteristic of the d coefficient with index
!                         l-m in the Meixner-Schafke normalization
!               idmlmse : exponent associated with dmlms
!               jnorm   : maximum index of enr required for convergence
!                         of dmfnorm and dmsnorm
!               jsub    : number of decimal digits of subtraction error
!                         incurred in calculating dmfnorm
!
!    use param
!
!  real(knd) scalars and array
    real(knd) aj, arr, c, coef, csq, dec, dmfnorm, dmsnorm, dmlmf, &
         dmlms, d01, ea, rm2, sgn, sump, ten, term, teste, testeo
    real(knd) enr(maxd)
!
    ten = 10.0e0_knd
    rm2 = m + m
    rm2m1 = m + m - 1
    rm2p1 = m + m + 1
    rm2m3 = m + m - 3
    dec = ten ** (-ndec - 1)
    nfac = nex / 3
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** nfac
    testeo = 1.0e0_knd / teste
    csq = c * c
    lm2 = (l - m) / 2
    ix = l - m - 2 * lm2
    lim2 = limd / 2 - ix
    sgn = 1.0e0_knd
     do 20 i = 1, lim2
     arr = ix + i + i
     ea = arr + arr + rm2
     if(i > lm2) go to 10
     if(enr(i) < (0.0e0_knd)) sgn = sgn * (-1.0e0_knd)
10    enr(i) = (ea - 1.0e0_knd) * (ea + 1.0e0_knd) * enr(i) / ((arr + rm2)* &
          (arr + rm2 - 1.0e0_knd) * csq)
20    continue
!
!  compute the Morse-Feshbach normalizing factor
    dmfnorm = 1.0e0_knd
    sump = 1.0e0_knd
    term = 1.0e0_knd
    jlow = l - m + 2
    jterm = lm2
    iflag = 0
    idmfe = 0
     do 30 j = jlow, limd, 2
     aj = j
     jterm = jterm + 1
     term = term * (aj + rm2) * enr(jterm) * (aj + rm2 - 1.0e0_knd)/ &
        (aj * (aj - 1.0e0_knd))
     if(term > 0.0e0_knd) sump = sump + term
     dmfnorm = dmfnorm + term
     if(abs(term / dmfnorm) < dec) go to 40
     if(abs(dmfnorm) < teste) go to 30
     dmfnorm = dmfnorm * testeo
     term = term * testeo
     sump = sump * testeo
     idmfe = idmfe + nfac
     iflag = 1
30    continue
40   jlow = l - m
    jmf = jterm
    if(jlow < 2 .or. iflag == 1) go to 60
    term = 1.0e0_knd
    jterm = lm2
     do 50 j = jlow, 2,-2
     aj = j
     term = term * aj * (aj - 1.0e0_knd) / ((aj + rm2 &
        -1.0e0_knd) * (aj + rm2) * enr(jterm))
     if(term > 0.0e0_knd) sump = sump + term
     jterm = jterm - 1
     dmfnorm = dmfnorm + term
     if(abs(term / dmfnorm) < dec) go to 60
50    continue
60   continue
    if(dmfnorm == 0.0e0_knd) dmfnorm = dec
    jsub = int(log10(abs(sump / dmfnorm) + dec))
    if(jsub < 0) jsub = 0
    if(jsub > ndec) jsub = ndec
    iterm = 0
    if(dmfnorm /= 0.0e0_knd) iterm = int(log10(abs(dmfnorm)))
    idmfe = idmfe + iterm
    dmfnorm = dmfnorm * (10.0e0_knd ** (-iterm))
    dmlmf = 1.0e0_knd / dmfnorm
    idmlmfe = -idmfe
if (debug) then
    write(40, 70) lim2, jmf, jsub
70   format(8x, i6,' "d" coefficients; mf norm. converged in ', &
        i6,' terms with ',i3,' digit subtr. error.')
end if
!
!  compute the d0(c|ml) or d1(c|ml)
    id01 = 0
    d01 = 1.0e0_knd
    if(lm2 == 0) go to 90
     do 80 kjl = 1, lm2
     kkjl = lm2 - kjl + 1
     d01 = d01 / enr(kkjl)
      if(abs(d01) > teste) then
      d01 = d01 * testeo
      id01 = id01 + nfac
      end if
      if(abs(d01) < testeo) then
      d01 = d01 * teste
      id01 = id01 - nfac
      end if
80    continue
    iterm = int(log10(abs(d01)))
    d01 = d01 * (10.0e0_knd ** (-iterm))
    id01 = id01 + iterm
90   continue
!
!  compute the Meixner-Schafke normalizing factor
    jflag = 0
    idmse = 0
    dmsnorm = 1.0e0_knd
    coef = 1.0e0_knd
    jlow = l - m + 2
    jterm = lm2
     do 150 j = jlow, limd, 2
     aj = j
     aj2 = aj + aj
     jterm = jterm + 1
     coef = coef * (aj + rm2) * enr(jterm) * (aj + rm2m1) * enr(jterm) &
        *(aj2 + rm2m3) / (aj * (aj - 1.0e0_knd) * (aj2 + rm2p1))
     dmsnorm = dmsnorm + coef
     if(abs(coef / dmsnorm) < dec) go to 160
      if(abs(dmsnorm) > teste) then
      dmsnorm = dmsnorm * testeo
      coef = coef * testeo
      idmse = idmse + nfac
      jflag = 1
      end if
150    continue
160   jlow = l - m
    jn = jterm
    jnorm = max(jmf, jn)
    if(jlow < 2 .or. jflag == 1) go to 180
    coef = 1.0e0_knd
    jterm = lm2
    j = jlow
     do 170 jj = 2, jlow, 2
     aj = j
     aj2 = aj + aj
     coef = coef * aj * (aj - 1.0e0_knd) * (aj2 + rm2p1) / ((aj2 + rm2m3)* &
         enr(jterm) * enr(jterm) * (aj + rm2) * (aj + rm2m1))
     jterm = jterm - 1
     j = j - 2
     dmsnorm = dmsnorm + coef
     if(abs(coef / dmsnorm) < dec) go to 180
170    continue
180   iterm = int(log10(dmsnorm))
    dmsnorm = dmsnorm * ten ** (-iterm)
    idmse = idmse + iterm
     if(2 * (idmse / 2) /= idmse) then
     idmse = idmse - 1
     dmsnorm = ten * dmsnorm
     end if
    dmlms = sgn / sqrt(dmsnorm)
    idmlmse = -idmse / 2
if (debug) then
    write(50, 190) jn, lim2
190   format(5x,' Meixner-Schafke normalization converged in ', &
        i6,' terms; ',i6,' terms available.')
end if
200   continue
    return
    end subroutine
!
!
    subroutine dalt (l, m, c, limdr, maxdr, maxmp, ndec, nex, ioppsum, &
             eigval, enrneg, drhor, dneg, idneg, nsdneg, nsdrho)
!
!  purpose:     To calculate d ratios with negative subscripts
!               and d-rho ratios.
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               limdr   : number of ratios of successive d-rho
!                         coefficients calculated
!               maxdr   : dimension of drhor array
!               maxmp   : dimension of enrneg array
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent for real(knd)
!               ioppsum : integer index = 0 if no d rho coefficients
!                         are calculated (psum not needed for r2leg)
!               eigval  : eigenvalue
!
!     output:   enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               dneg    : characteristic of the ratio of the d
!                         coefficient with index -2m+ix to the
!                         d coefficient with index ix, where
!                         ix = 0 if l-m is even and ix = 1 if
!                         l-m is odd
!               idneg   : exponent (base 10) of dneg
!               nsdneg  : subtaction error in calculating dneg and
!                         maximum error in calculating the enrneg
!                         array
!               nsdrho  : maximum subtraction error in calculating
!                         the drhor array
!
!    use param
!
!  real(knd) scalars
    real(knd) c, dneg, eigval, r, rm, rn, t, ten, teste, testeo, uterm, &
         vterm, wterm
    real(knd) amnsdrho, ansdneg, ansdrho, asub, bsub
!
!  real(knd) arrays
    real(knd) enrneg(maxmp), drhor(maxdr)
!
    ten = 10.0e0_knd
    nfac = nex / 3
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** (nfac)
    testeo = 1.0e0_knd / teste
!  if l-m is even, ix=0; if l-m is odd, ix=1
    ix = (l - m) - 2 * ((l - m) / 2)
    t = m + m - ix - ix
!
!  calculate ratios of d coefficients with negative subscripts
!
!       enrneg(k) = { d(-2m+2k-2)/d(-2m+2k), l-m even    }
!                   { d(-2m+1+2k-2)/d(-2m+1+2k), l-m odd }
!
    rm = m
     if(m == 0) then
     dneg = 1.0e0_knd
     idneg = 0
     go to 30
     end if
     do 10 i = 1, m + 1
     enrneg(i) = 0.0e0_knd
10    continue
!
!  first calculate enrneg(1)
    n = 2 - 2 * m + ix
    rn = n
    r = n + m + m
    uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r - 1.0e0_knd))
    r = n + m - 2
    vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
       1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
       1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (c * c)
!
!       calculations continue up to and including
!       enrneg(k=m) = { d(-2)/d(0), l-m even }
!                     { d(-1)/d(1), l-m odd  }
!
    enrneg(1) = -uterm / vterm
    dneg = enrneg(1)
    idneg = 0
    ansdneg = 0.0e0_knd
    nsdneg = 0
!
!  backward recursion beginning with enrneg(1) and
!  ending with enrneg(m)
!
     do i = 2, 2 * m - 2, 2
     ii = i - 2 * m + ix
     n = ii + 2
     j = i / 2
     rn = n
     r = n + m + m
     uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r- &
           1.0e0_knd))
     r = n + m - 2
     vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
        1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
        1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (c * c)
     r = n - 4
     wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm+ &
           3.0e0_knd) * (r + r + rm + rm + 1.0e0_knd))
     enrneg(j + 1) = -uterm / (wterm * enrneg(j) + vterm)
     dneg = dneg * enrneg(j + 1)
      if(wterm * enrneg(j) * vterm /= 0.0e0_knd) &
        then
      asub = log10(abs(vterm / (vterm+ &
               wterm * enrneg(j))))
      if(asub > 0.0e0_knd) ansdneg = ansdneg + asub
      bsub = log10(abs(vterm / (wterm * enrneg(j))))
      if(bsub > 0.0e0_knd) ansdneg = max(0.0e0_knd, ansdneg - bsub)
      if(int(ansdneg) > nsdneg) nsdneg = ansdneg
      end if
      if(abs(dneg) > teste) then
      dneg = dneg * testeo
      idneg = idneg + nfac
      end if
      if(abs(dneg) < testeo) then
      dneg = dneg * teste
      idneg = idneg - nfac
      end if
     end do
     if(nsdneg > ndec) nsdneg = ndec
    iterm = int(log10(abs(dneg)))
    dneg = dneg * (ten ** (-iterm))
    idneg = idneg + iterm
!
!  calculate ratios of d rho coefficients
!
!       drhor(k-m) = { d(rho|2k)/d(rh0|2k-2), l-m even  }
!                    { d(rho|2k-1)/d(rho|2k-3), l-m odd }
!
30   if(ioppsum == 0) go to 60
     ansdrho = 0.0e0_knd
     amnsdrho = 0.0e0_knd
     nsdrho = 0
     mnsdrho = 0
     do 40 i = 1, limdr
     drhor(i) = (0.0e0_knd, 0.0e0_knd)
40    continue
     do 50 i = 2 * limdr, 6,-2
     n = 4 - i + ix - m - m
     ii = (i - 2) / 2
     rn = n
     r = n + m + m
     uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r - 1.0e0_knd))
     r = n + m - 2
     vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
         1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
         1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (c * c)
     r = n - 4
     wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm + 3.0e0_knd)* &
        (r + r + rm + rm + 1.0e0_knd))
     drhor(ii) = -uterm / (wterm * drhor(ii + 1) + vterm)
      if(wterm * drhor(ii + 1) * vterm /= 0.0e0_knd) &
        then
      asub = log10(abs(vterm / (vterm + wterm * drhor(ii + 1))))
      if(asub > 0.0e0_knd) ansdrho = ansdrho + asub
      bsub = log10(abs(vterm / (wterm * drhor(ii + 1))))
      if(bsub > 0.0e0_knd) ansdrho = max(0.0e0_knd, ansdrho - bsub)
      if(ansdrho > amnsdrho) amnsdrho = ansdrho
      end if
50    continue
    n = -2 * m + ix
    r = n + m - 2
    vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm - 1.0e0_knd)/ &
       ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r - 1.0e0_knd))+ &
       (r * (r + 1.0e0_knd) - eigval) / (c * c)
    r = n - 4
    wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm + 3.0e0_knd)* &
       (r + r + rm + rm + 1.0e0_knd))
!
!       the final value of ii is 1;
!       drhor(1) has a special value:
!       drhor(1) = { d(rho|2m+2)/d(-2m), l-m even  }
!                  { d(rho|2m+1)/d(-2m+1), l-m odd }
!
    drhor(1) = 1.0e0_knd / ((t - 1.0e0_knd) * (t + 1.0e0_knd)* &
         (wterm * drhor(2) + vterm))
     if(wterm * drhor(2) * vterm /= (0.0e0_knd, 0.0e0)) &
        then
     asub = log10(abs(vterm / (vterm + wterm * drhor(2))))
     if(asub > 0.0e0_knd) ansdrho = ansdrho + asub
     bsub = log10(abs(vterm / (wterm * drhor(2))))
     if(bsub > 0.0e0_knd) ansdrho = max(0.0e0_knd, ansdrho - bsub)
     if(ansdrho > amnsdrho) amnsdrho = ansdrho
     end if
    nsdrho = int(amnsdrho) + 1
    if(ix == 1) drhor(1) = -drhor(1)
60   continue
    return
    end subroutine
!
!
    subroutine gauss (ndec, n, x, w)
!
!  purpose:     To evaluate the coordinates and weighting factors
!               for an nth order Gaussian quadrature
!
!  parameters:
!
!     input:    ndec: number of decimal digits in real(knd) arithmetic;
!                     usually equal to 33
!               n   : order of quadrature
!
!     output:   x   : coordinate values for quadrature
!               w   : weighting factors
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) delta, der, pi, s, t, test, u, v, z
    real(knd) x(n), w(n)
!
    test = 10.0e0_knd ** (-ndec)
    imax = (n + 1) / 2
    pi = acos(-1.0e0_knd)
     do 40 i = 1, imax
     z = cos(pi * (i - 0.25e0_knd) / (n + 0.5e0_knd))
      do 20 j = 1, 30
      u = 0.0e0_knd
      v = 1.0e0_knd
       do 10 k = 1, n
       t = u
       u = v
       v = ((k + k - 1) * z * u - (k - 1) * t) / k
10      continue
      s = z * z - 1.0e0_knd
      der = n * (z * v - u) / s
      delta = -v / der - 0.5e0_knd * v * v * ((n * n * s - n * z * z - n) * v+ &
          2.0e0_knd * n * z * u) / (der * der * der * s * s)
      z = z + delta
      if(abs(delta / z) < test) go to 30
20     continue
30    continue
     x(i) = -z
     x(n + 1 - i) = z
     w(i) = 2.0e0_knd / ((1.0e0_knd - z * z) * der * der)
     w(n + 1 - i) = w(i)
40    continue
    return
    end subroutine
!
!
    subroutine pleg (m, lim, maxp, limcsav, iopd, ndec, nex, barg, narg, &
             maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, &
             beta, gamma, coefa, coefb, coefc, coefd, coefe)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the first kind for given arguments barg.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lim    : two less than the number of associated Legendre
!                        function ratios calculated for given arguments
!               maxp   : dimension of alpha, beta, gamma, coefa, coefb,
!                        coefc, coefd, and coefe arrays and second
!                        dimension of pr and pdr arrays
!               limcsav: integer equal to the number of coefficients in
!                        each of the arrays alpha, beta, gamma, coefa,
!                        coefb, coefc, coefd, and coefe arrays that
!                        have already been calculated in earlier calls
!                        to pleg for this value of m and will not be
!                        calculated again. [Note that the minimum
!                        array index for the coefficients is 3 and
!                        the maximum array index is limcsav+2]
!               iopd   : integer that is set = 0 if derivatives of
!                        Legendre functions (i.e., their ratios)
!                        are not required when iopang = 1 and the
!                        first derivatives of the angular functions
!                        are not requested.
!                        iopd is set = 1 when iopang = 2 and pleg is
!                        also being used to obtain ratios of first
!                        derivatives of Legendre functions for use in
!                        computing the first derivatives of the angular
!                        functions.
!                        iopd is set = 2 when pleg is being used to
!                        compute ratios of Legendre functions for use in
!                        the calculation of the denominator term used in
!                        calculating the radial functions of the second
!                        kind and their first derivatives in r2eta.
!                        iopd is set = 3 when pleg is being used to
!                        compute ratios of both the Legendre functions
!                        and their first derivatives for use in the
!                        calculation of the numerator terms used
!                        in r2eta and in r2leg to calculate the radial
!                        functions of the second kind and their first
!                        deriatives.
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent for real(knd)
!               barg   : array of narg values of eta for which Legendre
!                        functions are to be calculated
!               narg   : number of specified values of eta in barg array
!               maxt   : dimension of barg array
!
!     output:   pr     : array of ratios of successive first kind
!                        associated Legendre functions of the same
!                        parity
!               pdr    : array of ratios of successive derivatives of
!                        first kind associated Legendre functions of
!                        the same parity
!               pdnorm : array of characteristics of the first
!                        derivatives of associated Legendre functions
!                        of the first kind of order m and degree m
!               ipdnorm: array of exponents of the first derivatives
!                        of associated Legendre functions of the first
!                        kind of order m and degree m
!               pnorm  : array of characteristics of the associated
!                        Legendre functions of the first kind of order m
!                        and degree m
!               ipnorm : array of exponents of the associated Legendre
!                        functions of the first kind of order m and
!                        degree m
!
!     input/output:
!               alpha  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               beta   : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               gamma  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               coefa  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefb  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefc  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefd  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefe  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) adec, ajterm, am2p1, anden1, anden2, an2tnp1, bargs, den, rm, &
         rm2, temp1, temp2, temp3, ten, term, teste, testeo, ta, tb, tc, &
         t1, t2
    real(knd) alpha(maxp), barg(maxt), beta(maxp), coefa(maxp), &
         coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp), &
         gamma(maxp), pdnorm(maxt), pdr(maxt, maxp), pdr1(maxp), &
         pr(maxt, maxp), pnorm(maxt)
!
!  integer array
    dimension ipdnorm(maxt), ipnorm(maxt)
!
    ten = 10.0e0_knd
    adec = ten ** (-ndec + 2)
    nfac = nex / 2
    if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
    teste = ten ** nfac
    testeo = 1.0e0_knd / teste
    rm = m
    rm2 = m * m
    am2p1 = m + m + 1
    m2 = m + m
    m2m1 = m2 - 1
    mm1 = m - 1
    mm2 = m - 2
    msqp1 = 2 * m * m + 1
!
!  calculate the coefficients alpha(j), beta(j), and gamma(j) for
!  the three term recursion relating the Legendre function ratios
!
!              m                m
!   pr(k,j) = P    (barg(k)) / P  (barg(k))
!              m+j-1            m+j-3
!
!  and calculate the coefficients coefa(j), coefb(j), coefc(j),
!  coefd(j), and coefe(j) in the expression used to calculate
!  ratios of Legendre function derivatives
!
!               m                 m
!   pdr(k,j) = P'    (barg(k)) / P'  (barg(k))
!               m+j-1             m+j-3
!
!  Note that pr(k,1) and pr(k,2) are not ratios but actually equal to
!   m      m                                              m       m
!  P  and P   . Also, pdr(k,1) and pdr(k,2) are equal to P'  and P' .
!   m      m+1                                            m       m+1
!
     if(limcsav >= lim) go to 30
     do 10 j = limcsav + 3, lim + 2
     n = m + j - 3
     n2 = n + n
     n2p3 = n2 + 3
     n2p1 = n2 + 1
     n2m1 = n2 - 1
     nmmm2 = n - mm2
     nmmm1 = n - mm1
     npmm1 = n + mm1
     npm = n + m
     npmm1 = n + mm1
     npmp1 = n + m + 1
     npmp2 = npmp1 + 1
     an2tnp1 = real(2 * n, knd) * real((n + 1), knd)
     anden1 = real(nmmm2, knd) * real(nmmm1, knd)
     anden2 = real(n2m1, knd) * anden1
     alpha(j) = real(n2p3, knd) * real(n2p1, knd) / anden1
     beta(j) = real(n2p1, knd) * (real(msqp1, knd) - an2tnp1) / anden2
     gamma(j) = -real(n2p3, knd) * real(npm, knd) * real(npmm1, knd) / anden2
     coefa(j) = -real(npmp2, knd) / real(nmmm1, knd)
     coefb(j) = real(n2p3, knd) * real((n + 2), knd) / anden1
     coefc(j) = -real(npmp1, knd) * real(npmp2, knd) / anden1
     coefd(j) = real(npmp1, knd) / real(nmmm2, knd)
     coefe(j) = -real((n + 1), knd) * real(n2p3, knd) / anden1
10    continue
    gamma(3) = 0.0e0_knd
    gamma(4) = 0.0e0_knd
    term = 1.0e0_knd
    iterm = 0
    if(m < 2) go to 30
     do jm = 2, m
     term = (jm + jm - 1) * term
      if(term > teste) then
      term = term * testeo
      iterm = iterm + nfac
      end if
     end do
    jterm = int(log10(term))
    term = term * (ten ** (-jterm))
    iterm = iterm + jterm
30   continue
!
!   calculate the ratios of successive Legendre functions of the same
!   parity using the three term recursion relationship
!
!   pr(k,j) = alpha(j)*barg(k)*barg(k) + beta(j) + gamma(j)/pr(k,j-2)
!
     do 140 k = 1, narg
     pnorm(k) = term
     ipnorm(k) = iterm
     pdnorm(k) = term
     ipdnorm(k) = iterm
!
!   define the first two ratios equal to unity and (2m+1)*barg(k)
     pr(k, 1) = 1.0e0_knd
     pr(k, 2) = am2p1 * barg(k)
     jdelta = 1
     if(abs(barg(k)) < adec) jdelta = 2
     bargs = barg(k) * barg(k)
      do 40 j = 3, lim + 2, jdelta
      pr(k, j) = alpha(j) * bargs + beta(j) + (gamma(j) / pr(k, j - 2))
40     continue
!
!   calculate the corresponding ratios of first derviatives of
!   successive Legendre functions of the same parity using the
!   following relationship (except for (1) eta equal to zero or unity,
!   where special expressions are used and (2) when the magnitude of the
!   argument barg is <= 0.1 or abs((m+1)*barg*barg - 1) < 0.01, where
!   recursion on the ratios of successive first derivatives of the same
!   parity is used instead)
!
!              (coefa(j)+coefb(j)*barg(k)*barg(k))*pr(k,j)+coefc(j)
!   pdr(k,j) = ----------------------------------------------------
!                  pr(k,j)+coefd(j)+coefe(j)*barg(k)*barg(k)
!
     if(iopd == 0 .or. iopd == 2) go to 120
     pdr(k, 1) = 1.0e0_knd
     pdr(k, 2) = 1.0e0_knd
     if(abs(barg(k)) >= adec) go to 50
     pdr(k, 2) = am2p1
      do j = 4, lim + 2, 2
      pdr(k, j) = -real(m2m1 + j, knd) / real(j - 2, knd)
      end do
     go to 140
50    if(abs(abs(barg(k)) - 1.0e0_knd) >= adec) go to 70
     if(m == 0) go to 60
     if(m /= 2) go to 130
     pdr(k, 1) = -2.0e0_knd * barg(k)
     go to 80
60    temp1 = 1.0e0_knd
     temp2 = 3.0e0_knd
     pdr(k, 2) = 1.0e0_knd
     pdr(k, 3) = 3.0e0_knd * barg(k)
      do j = 4, lim + 2
      temp3 = temp2 + j - 1
      pdr(k, j) = temp3 / temp1
      temp1 = temp2
      temp2 = temp3
      end do
     go to 140
70    if(m /= 0) go to 80
     pdr(k, 1) = 1.0e0_knd
     pdr(k, 2) = 1.0e0_knd
     pdr(k, 3) = 3.0e0_knd * barg(k)
     jlow = 4
     go to 90
80    pdr(k, 2) = am2p1 * ((rm + 1.0e0_knd) * bargs - 1.0e0_knd) / (rm * barg(k))
     if(pdr(k, 2) == 0.0e0_knd) pdr(k, 2) = ten ** (-ndec)
     jlow = 3
     if(abs((rm + 1.0e0_knd) * bargs - 1.0e0_knd) < 0.01e0_knd) go to 110
90    continue
     if(abs(barg(k)) <= 0.1e0_knd) go to 110
      do 100 j = jlow, lim + 2
      den = (pr(k, j) + coefd(j) + coefe(j) * bargs)
      if(den == 0.0e0_knd) den = ten ** (-ndec)
      pdr(k, j) = ((coefa(j) + coefb(j) * bargs) * pr(k, j) + coefc(j)) / den
100     continue
     go to 120
110    continue
     if(m /= 0) pdr1(1) = pdr(k, 2)
     if(m == 0) pdr1(2) = pdr(k, 3)
      do j = jlow - 1, lim + 1
      n = j + m - 1
      t1 = bargs - 1.0e0_knd
      t2 = (n * n * t1 + rm2)
      ta = j * t2
      tb = (n + n + 1) * barg(k) * (t2 + n * t1)
      tc = -(n + m) * (t2 + (n + n + 1) * t1)
      pdr1(j) = (tb + tc / pdr1(j - 1)) / ta
      end do
      do j = jlow, lim + 2
      pdr(k, j) = pdr1(j - 2) * pdr1(j - 1)
      end do
120    if(m == 0 .or. iopd == 2 .or. iopd == 3) go to 140
     if(abs(abs(barg(k)) - 1.0e0_knd) < adec) go to 130
     ajterm = rm * log10(1.0e0_knd - bargs) / 2.0e0_knd
     jterm = int(ajterm)
     ipnorm(k) = ipnorm(k) + jterm
     pnorm(k) = pnorm(k) * (ten ** (ajterm - jterm))
     if(iopd == 0) go to 130
     ajterm = log10(rm * abs(barg(k))) + (rm - 2.0e0_knd)* &
         log10(1.0e0_knd - bargs) / 2.0e0_knd
     jterm = int(ajterm)
     ipdnorm(k) = ipdnorm(k) + jterm
     pdnorm(k) = -pdnorm(k) * (ten ** (ajterm - jterm))
     if(barg(k) < 0.0e0_knd) pdnorm(k) = -pdnorm(k)
     go to 140
130    pnorm(k) = 0.0e0_knd
     ipnorm(k) = 0
     if(m /= 2) pdnorm(k) = 0.0e0_knd
     if(m /= 2) ipdnorm(k) = 0
140    continue
    return
    end subroutine
!
!
    subroutine qleg (m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, &
             iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the second kind for given c,x, and m.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the second
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lnum   : number of l values desired (=lmax+1);
!                        also equal to the dimension of the arrays
!                        ql, iql, qdl, and iqdl
!               limq   : the number of associated Legendre function
!                        ratios calculated for given m,lnum,c,ndec,
!                        and x1
!               maxq   : dimension of qr and qdr arrays
!               x1     : x - 10.0e0_knd
!               ndec   : number of decimal digits for real(knd)
!
!     output:   qdr    : ratios of derivatives of successive
!                        associated Legendre functions of the second
!                        kind
!               qdml   : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m-1
!               iqdml  : exponent corresponding to qdml
!               qdl    : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m
!               iqdl   : exponent corresponding to qdl
!               qr     : ratios of successive associated Legendre
!                        functions of the second kind
!               qml    : characteristic of the associated Legendre
!                        function of the second kind with order m
!                        and degree m-1
!               iqml   : exponent corresponding to qml
!               ql     : characteristic of the associated Legendre
!                        function of the second kind with order m and
!                                                           -m/2
!                        degree m, scaled by (2m-1)!!(x*x-1)
!               iql    : exponent corresponding to ql
!               termpq : characteristic of the relative size of the
!                        maximum terms in the positive degree q series
!                        and the p series used to calculate r2 and r2d
!                        in subroutine r2leg
!               itermpq: exponent corresponding to termpq
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) ajm, dec, qdml, qlow, qml, qupp, q00, q11, rin, rm, &
         term, termpq, tjm, tm, tmr, x, x1, x1d, xsqr
    real(knd) qdl(lnum), qdr(maxq), ql(lnum), qr(maxq)
!
!  integer arrays
    dimension iqdl(lnum), iql(lnum)
!
    dec = 10.0e0_knd ** (-ndec)
    rm = m
    tm = rm + rm
    tmr = tm / (tm + 1.0e0_knd)
    x = x1 + 1.0e0_knd
    x1d = (x + 1.0e0_knd) * x1
    xsqr = sqrt(x1d)
    mxqrest = limq + ndec * int((1.0e0_knd - 1.0e0_knd / log10(x - xsqr)))
    if(m == 0) mlimq = 50000 * ndec + limq
    if(m == 1) mlimq = 12000 * ndec + limq
    if(m == 2) mlimq = 5000 * ndec + limq
    if(m == 3) mlimq = 600 * ndec + limq
    if(m >= 4) mlimq = 100 * ndec + limq
    if(m == 1 .and. x1 < 1.0e-9_knd) mlimq = 50000 * ndec + limq
    mxqr = min(mxqrest, mlimq)
    mxqrpm = mxqr + m
if (debug) then
    write(40, 5) mxqrpm
5    format(15x,'used backward recursion to calculate ratios of q', &
        ' functions starting at order',i8)
end if
!
!                              m    m
!  calculate ratios qr(k+m) = q  / q    ,k=m+limq to k=m+1 using
!                              k    k-1
!
!                                         m       m               1/2
!  backward recursion from qr(maxm+2m) = q     / q      =x-(x*x-1),
!                                         mxqr+m  mxqr+m-1
!
!  where the last quantity is the asymptotic limit of the ratio as
!  mxqr approaches infinity
!
    qupp = x - xsqr
     do jn = mxqr + m, limq + m + 1,-1
     rin = jn
     qlow = (rin + rm - 1.0e0_knd) / (x * (rin + rin - 1.0e0_knd) &
        -(rin - rm) * qupp)
     qupp = qlow
     end do
    qr(limq + m + m) = qupp
     do 10 jn = limq + m, m + 2,-1
     rin = jn
     qr(jn + m - 1) = (rin + rm - 1.0e0_knd) / (x * (rin + rin - 1.0e0_knd) &
           -(rin - rm) * qr(jn + m))
10    continue
!
!                              m     m
!  calculate ratios qr(k+m) = q   / q   ,k=m-1 to k=-m+1 using
!                              k-1   k
!
!                                       m      m
!  backward recursion from qr(m-1+m) = q    / q     = x
!                                       m-2    m-1
!
20   if(m == 0) go to 100
    qr(m + m - 1) = x
    if(m == 1) go to 40
     do 30 jn = m - 1, 2 - m,-1
     rin = jn
     qr(jn + m - 1) = (x * (rin + rin - 1.0e0_knd) &
           -((rin - rm) / qr(jn + m))) / (rin + rm - 1.0e0_knd)
30    continue
40   continue
!
!                  m
!  calculation of q    , m > 0 by forward division of qr ratios
!                  m-1
!
!                 m
!  starting with q  calculated from its closed form expression,
!                 0
!
!                           -m/2
!  scaled by (2m-1)!!(x*x-1).
!
    qml = rm * log10(x + 1.0e0_knd) - log10(2.0e0_knd)
    iqterm = int(rm * log10(x1 / (x + 1.0e0_knd)))
    if(iqterm < -ndec) go to 50
    qml = qml + log10(1.0e0_knd - ((x1 / (x + 1.0e0_knd)) ** m))
50   continue
    term = 1.0e0_knd
    iterm = 0
    if(m < 2) go to 70
     do jm = 2, m
     ajm = jm
     term = term * (ajm - 1.0e0_knd) / (ajm + ajm - 1.0e0_knd)
     if(term > dec) go to 60
     term = term / dec
     iterm = iterm - ndec
60    continue
     end do
70   term = log10(term)
    qml = qml + term + iterm
    iqml = int(qml)
    qml = 10.0e0_knd ** (qml - iqml)
    if(2 * (m / 2) /= m) qml = -qml
    if(m < 2) go to 90
     do jm = 1, m - 1
     qml = qml / qr(jm + m)
     if(abs(qml) > dec) go to 80
     qml = qml / dec
     iqml = iqml - ndec
80    continue
     end do
90   continue
    iterm = int(log10(abs(qml)))
    qml = qml * 10.0e0_knd ** (-iterm)
    iqml = iqml + iterm
!
!                  m
!  calculation of q   by forward recursion in m starting with values
!                  m
!       0       1
!  for q   and q  obtained from their closed form expressions, scaled
!       0       1
!                    -m/2
!  by (2m-1)!!(x*x-1).
!
100   q00 = 0.5e0_knd * log((x + 1.0e0_knd) / x1)
    if(m /= 0) go to 110
    ql(1) = q00
    go to 130
110   q11 = x1d * q00 - x
    if(m /= 1) go to 120
    ql(1) = q11
    go to 130
120   qlow = q00
    qupp = q11
     do jm = 1, m - 1
     tjm = real((jm + jm), knd) / (jm + jm + 1)
     ql(1) = (x1d - tjm) * qupp + tjm * x1d * qlow
     qlow = qupp
     qupp = ql(1)
     end do
130   iql(1) = int(log10(abs(ql(1))))
    ql(1) = ql(1) * (10.0e0_knd ** (-iql(1)))
!
!  calculation of ratios of the first derivatives of q with respect
!  to x, using the relationships:
!
!                  m    m      [kx]qr(k+m)-(k+m)
!     qdr(k+m) = q'  / q'   =  ----------------- , k=m+lim to k=m+1
!                  k    k-1    [(k-m)]qr(k+m)-kx
!
!                  m      m    [(k-m)]qr(k+m)-kx
!     qdr(k+m) = q'   / q'  =  ----------------- , k=m-1 to k=-m+1
!                  k-1    k    [kx]qr(k+m)-(k+m)
!
     do jm = m + 1, m + limq
     ajm = jm
     qdr(jm + m) = (ajm * x * qr(jm + m) - (ajm + rm)) / ((ajm - rm) * qr(jm + m) - ajm * x)
     end do
!
    if(m == 0) go to 140
     do jm = 1 - m, m - 1
     ajm = jm
     qdr(jm + m) = (ajm * x * qr(jm + m) - (ajm - rm)) / ((ajm + rm) * qr(jm + m) - ajm * x)
     end do
140   continue
!
!                   m         m                      m        m
!  calculation of q'    and q'  from the values for q    and q  .
!                   m-1       m                      m-1      m
!
    if(m > 0) go to 150
    qdl(1) = -1.0e0_knd / x1d
    iqdl(1) = 0
    go to 160
150   qdml = -rm * x * qml / x1d
    iterm = int(log10(abs(qdml)))
    qdml = qdml * (10.0e0_knd ** (-iterm))
    iqdml = iqml + iterm
    qdl(1) = rm * (x * ql(1) - 2.0e0_knd * qml * (10.0e0_knd ** (iqml - iql(1)))) &
        /x1d
    iqdl(1) = iql(1)
160   continue
    m2m1 = m + m - 1
     do jl = 2, lnum
     ql(jl) = ql(jl - 1) * qr(m2m1 + jl)
     iql(jl) = iql(jl - 1)
     if(abs(ql(jl)) > 1.0e-10_knd) go to 170
     ql(jl) = ql(jl) * 1.0e10_knd
     iql(jl) = iql(jl) - 10
170    qdl(jl) = qdl(jl - 1) * qdr(m2m1 + jl)
     iqdl(jl) = iqdl(jl - 1)
     if(abs(qdl(jl)) > 1.0e-10_knd) go to 180
     qdl(jl) = qdl(jl) * 1.0e10_knd
     iqdl(jl) = iqdl(jl) - 10
180    end do
    termpq = rm * log10(xsqr)
    itermpq = int(termpq)
    termpq = 10.0e0_knd ** (termpq - itermpq)
    return
    end subroutine
!
!
    subroutine pint (c, m, lnum, x1, limint, maxint, maxlp, maxmp, ndec, &
            wg, xg, ngau, ngqs, rpint1, rpint2, pint1, pint2, &
            pint3, pint4, norme, pnorm, ipnorm, coefme, coefmo)
!
!  purpose:     To calculate integrals of the product of associated
!               Legendre functions and kernels containing spherical
!               Neumann functions and a window function. Four
!               different kernel functions are involved leading to
!               integrals of four different types. The integrals are
!               calculated using gaussian quadrature.
!
!  parameters:
!
!     input:    c      : c
!               m      : m
!               lnum   : number of l values desired
!               x1     : x - 10.0e0_knd
!               limint : number of integrals of each of the four types
!                        required
!               maxint : dimension of the integral arrays
!               maxlp  : dimension of characteristic and exponent
!                        arrays of integrals
!               maxmp  : dimension of the spherical Neumann function
!                        array
!               ndec   : number of decimal digits for real(knd)
!               wg     : gaussian quadrature weighting factors
!               xg     : gaussian quadrature arguments
!               ngau   : order of gaussian quadrature
!               ngqs   : number of gaussian quadrature steps the
!                        integrand is divided into
!
!     output:   rpint1 : array of ratios of successive integrals of
!                        the same parity of the first type (l-m even)
!                        or of the third type (l-m odd)
!               rpint2 : array of ratios of successive integrals of
!                        the same parity of the second type (l-m even)
!                        or of the fourth type (l-m odd)
!               pint1  : array of scaled values for the integrals of
!                        the first type
!               pint2  : array of scaled values for the integrals of
!                        the second type
!               pint3  : array of scaled values for the integrals of
!                        the third type
!               pint4  : array of scaled values for the integrals of
!                        the fourth type
!               norme  : array of scaling exponents for the spherical
!                        Neumann functions
!               pnorm  : array of characteristics for the scaling
!                        factors used for the associated Legendre
!                        functions
!               ipnorm : array of exponents for the scaling factors
!                        used for the associated Legendre functions
!               coefme : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m even)
!               coefmo : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m odd)
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) ak, amo2, an, arg, argb, arn, bn, c, coef, coefme, &
         coefmo, coefo, coef1, coef2, coef4, darg, dec, etai, &
         etaim1, etais, etal, etau, etcoef1, etcoef2, rm, rn, sargb, &
         scal, sneu1, sneu2, sneu3, stemp0, step, step0, step1, step2, &
         twom, twomi, x1, x1sm1
    real(knd) alpha(maxint), beta(maxint), p(maxint), pint1(maxint), &
         pint2(maxint), pint3(maxint), pint4(maxint), &
         pnorm(maxlp), rpint1(maxint), rpint2(maxint), &
         sneun(maxmp), wg(ngau), xg(ngau), ynormn(maxmp)
!
!  integer arrays
    dimension norme(maxmp), ipnorm(maxlp)
!
    x1sm1 = x1 * (x1 + 2.0e0_knd)
    dec = 10.0e0_knd ** (-ndec - 1)
    amo2 = 0.5e0_knd * m
    lim = limint
    step0 = 1.0e0_knd / ngqs
    coefme = (m * (x1 + 1.0e0_knd)) / x1sm1
    coefmo = (m * (x1 + 1.0e0_knd) ** 2 + x1sm1) / ((x1 + 1.0e0_knd) * x1sm1)
    rm = m
    if(x1 >= 0.1e0_knd) go to 10
    step1 = 1.0e0_knd / (ngqs * (1.0e0_knd - 3.0e0_knd * log10(x1)* &
       (1.0e0_knd + 3.0e0_knd * log10(rm + 1.0e0_knd))))
    if(step1 > step0) step1 = step0
    go to 20
10   step1 = step0
20   step2 = (1.0e0_knd - step1) / (ngqs - 1)
if (debug) then
    write(40, 30) ngau, step1
30   format(11x,'order of gauss quadrature =',i4,'. first step' &
        ' size = 'f10.6,'.')
end if
!
!  calculation of scaling factors for the associated Legendre functions
    pnorm(1) = 1.0e0_knd
    pnorm(2) = 1.0e0_knd
    ipnorm(1) = 0
    ipnorm(2) = 0
    if(m == 0) go to 50
     do 40 n = 1, m
     an = n + n
     bn = n + n + 1
     pnorm(1) = pnorm(1) * an
     pnorm(2) = pnorm(2) * bn
     iterm1 = int(log10(pnorm(1)))
     iterm2 = int(log10(pnorm(2)))
     pnorm(1) = pnorm(1) * 10.0e0_knd ** (-iterm1)
     pnorm(2) = pnorm(2) * 10.0e0_knd ** (-iterm2)
     ipnorm(1) = ipnorm(1) + iterm1
     ipnorm(2) = ipnorm(2) + iterm2
40    continue
50   twom = m + m
    pnorm(3) = pnorm(1) * (twom + 2) / 2
    iterm3 = int(log10(pnorm(3)))
    pnorm(3) = pnorm(3) * 10.0e0_knd ** (-iterm3)
    ipnorm(3) = iterm3 + ipnorm(1)
    if(lnum < 4) go to 70
     do 60 il = 4, lnum, 2
     pnorm(il) = pnorm(il - 2) * (twom + il - 1) / (il - 1)
     pnorm(il + 1) = pnorm(il - 1) * (twom + il) / (il)
     iterm1 = int(log10(pnorm(il)))
     iterm2 = int(log10(pnorm(il + 1)))
     ipnorm(il) = ipnorm(il - 2) + iterm1
     ipnorm(il + 1) = ipnorm(il - 1) + iterm2
     pnorm(il) = pnorm(il) * 10.0e0_knd ** (-iterm1)
     pnorm(il + 1) = pnorm(il + 1) * 10.0e0_knd ** (-iterm2)
60    continue
70   continue
!
!  calculation of the coefficients in the recursion relation used
!  for the scaled associated Legendre functions
    alpha(1) = (twom + 1.0e0_knd) * pnorm(1) / pnorm(2)
    alpha(1) = alpha(1) * 10.0e0_knd ** (ipnorm(1) - ipnorm(2))
    beta(1) = 0.0e0_knd
    alpha(2) = (twom + 3.0e0_knd) * pnorm(2) / (pnorm(3) * 2.0e0_knd)
    alpha(2) = alpha(2) * 10.0e0_knd ** (ipnorm(2) - ipnorm(3))
    beta(2) = -(twom + 1.0e0_knd) / (twom + 2.0e0_knd)
     do 80 il = 3, lim + 2
     alpha(il) = alpha(il - 2) * (twom + il - 1) * (twom + il + il - 1)* &
     (il - 2) / ((il - 1) * (twom + il) * (twom + il + il - 5))
     beta(il) = -(twom + il - 1) / (twom + il)
80    continue
!
     do 90 il = 1, lim + 2, 2
     pint1(il) = 0.0e0_knd
     pint2(il) = 0.0e0_knd
     pint3(il + 1) = 0.0e0_knd
     pint4(il + 1) = 0.0e0_knd
90    continue
!
!  calculation of the scaling exponents for the spherical Bessel
!  functions required for the four types of integrals
    twomi = 1.0e0_knd
    if(m == 0) go to 110
     do n = 1, m
       twomi = twomi * (n + n - 1) / (n + n)
     end do
110   continue
    arg = c * sqrt(x1 * (x1 + 2.0e0_knd))
    darg = 1.0e0_knd / arg
    stemp0 = -cos(arg) * darg
    norme(1) = int(log10(abs(stemp0)))
    ynormn(1) = stemp0 * 10.0e0_knd ** (-norme(1))
    ynormn(2) = stemp0 * darg - sin(arg) * darg
    norme(2) = norme(1)
    ynormn(2) = ynormn(2) * 10.0e0_knd ** (-norme(1))
!
     do 120 n = 3, m + 3, 2
     rn = n + n - 3
     arn = n + n - 1
     ynormn(n) = -ynormn(n - 2) + darg * rn * ynormn(n - 1)
     ynormn(n + 1) = -ynormn(n - 1) + darg * arn * ynormn(n)
     norme(n + 1) = log10(abs(ynormn(n + 1)))
     scal = 10.0e0_knd ** (-norme(n + 1))
     ynormn(n + 1) = ynormn(n + 1) * scal
     ynormn(n) = ynormn(n) * scal
     norme(n + 1) = norme(n + 1) + norme(n - 1)
     norme(n) = norme(n + 1)
120    continue
!
!  gaussian quadrature integration loops. first dividing integrand
!  into ngqs steps
     do 180 k = 1, ngqs
     ak = k
     if(k /= 1) go to 130
     etal = 0.0e0_knd
     etau = step1
     step = step1
     go to 140
130    etal = step1 + (ak - 2.0e0_knd) * step2
     etau = etal + step2
     step = step2
140    etcoef1 = (etau + etal) / 2.0e0_knd
     etcoef2 = (etau - etal) / 2.0e0_knd
!
!  gaussian quadrature integration over each step
      do 170 i = 1, ngau
      etai = etcoef1 + xg(i) * etcoef2
      etais = etai * etai
      etaim1 = 1.0e0_knd - etais
      argb = x1sm1 + etais
      coef = step * ((x1sm1 * etaim1 / argb) ** amo2)
      if(coef < dec) go to 190
      sargb = sqrt(argb)
      coefo = 1.0e0_knd / sargb
      arg = c * sargb
      stemp0 = -cos(arg) / arg
      sneun(1) = stemp0
      sneun(1) = sneun(1) * 10.0e0_knd ** (-norme(1))
      sneun(2) = stemp0 / arg - sin(arg) / arg
      sneun(2) = sneun(2) * 10.0e0_knd ** (-norme(2))
      darg = 1.0e0_knd / arg
       do 150 n = 3, m + 3, 2
       rn = n + n - 3
       arn = n + n - 1
       sneun(n) = -sneun(n - 2) + darg * rn * sneun(n - 1)
       sneun(n + 1) = -sneun(n - 1) + darg * arn * sneun(n)
       sneun(n + 1) = sneun(n + 1) * 10.0e0_knd ** (-norme(n + 1) + norme(n - 1))
       sneun(n) = sneun(n) * 10.0e0_knd ** (-norme(n) + norme(n - 1))
150      continue
      sneu1 = sneun(m + 1)
      sneu2 = sneun(m + 2) * 10.0e0_knd ** (norme(m + 2) - norme(m + 1))
      sneu3 = sneun(m + 3) * 10.0e0_knd ** (norme(m + 3) - norme(m + 1))
      p(1) = twomi * (etaim1 ** amo2)
      p(2) = alpha(1) * etai * p(1)
      coef1 = coef * sneu1 * wg(i)
      coef2 = coef * coefo * sneu2 * wg(i)
      coef4 = coef * coefo * coefo * etai * sneu3 * wg(i)
      pint1(1) = pint1(1) + coef1 * p(1)
      pint2(1) = pint2(1) + coef2 * p(1)
      pint4(2) = pint4(2) + coef4 * p(2)
!
       do 160 il = 2, lim + 1, 2
       p(il + 1) = alpha(il) * etai * p(il) + beta(il) * p(il - 1)
       p(il + 2) = alpha(il + 1) * etai * p(il + 1) + beta(il + 1) * p(il)
       pint1(il + 1) = pint1(il + 1) + coef1 * p(il + 1)
       pint2(il + 1) = pint2(il + 1) + coef2 * p(il + 1)
       pint4(il + 2) = pint4(il + 2) + coef4 * p(il + 2)
160      continue
170     continue
180    continue
190   continue
     do 200 il = 1, lim, 2
     pint3(il + 1) = (pint2(il + 2) - beta(il + 1) * pint2(il)) &
           /alpha(il + 1)
200    continue
!
!  calculation of ratios of integrals for ease in compution of r2 and
!  r2d in subroutine r2int
    rpint1(1) = 0.0e0_knd
    rpint1(2) = 0.0e0_knd
    rpint2(1) = 0.0e0_knd
    rpint2(2) = 0.0e0_knd
     do 210 il = 3, lim, 2
     rpint1(il) = pint1(il) * (twom + il - 1) / (pint1(il - 2) * (il - 1))
     rpint2(il) = pint2(il) * (twom + il - 1) / (pint2(il - 2) * (il - 1))
     rpint1(il + 1) = pint3(il + 1) * (twom + il) / (pint3(il - 1) * (il))
     rpint2(il + 1) = pint4(il + 1) * (twom + il) / (pint4(il - 1) * (il))
210    continue
    return
    end subroutine
!
!
    subroutine sphbes (c, x, limj, maxj, maxlp, sbesf, sbesdf, sbesn, &
              ibese, sbesdr)
!
!  purpose:     To calculate ratios of successive first kind spherical
!               Bessel functions of the same parity for given c and x.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Bessel functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    c      : c
!               x      : x
!               limj   : the number of spherical Bessel function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxj   : dimension of sbesf vector
!               maxlp  : the number of scale factors
!                        that are calculated
!
!     output:   sbesf  : ratios of successive first kind spherical
!                        Bessel functions of the same parity
!               sbesdf : ratios of first derivatives of successive
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : characteristics for the spherical
!                        Bessel functions
!               ibese  : exponents for the spherical
!                        Bessel functions
!               sbesdr : ratios of first derivatives of spherical Bessel
!                        functions to the corresponding spherical
!                        spherical functions
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, cx, rn, stemp0, stemp1, x
    real(knd) sbesdf(maxj), sbesdr(maxlp), sbesf(maxj), sbesn(maxlp)
!
!  integer array
    dimension ibese(maxlp)
!
    cx = c * x
    lim1 = cx + cx + 20
!
!  compute first kind Bessel function ratios
!        sbesf(k)= j(n=k,c*x)/j(n=k-1,c*x)
!        sbesn(k)= (j(n=k-1),c*x))*10.0e0_knd**(-ibese(k))
!
    if (cx < limj) go to 20
!
!  for c*x >= lim, use forward recursion to
!  get fcn. ratios:
!       j(n+1,c*x)/j(n,c*x)=(2*n+1)/(c*x)-1/(j(n,c*x)/j(n-1,c*x))
!
    stemp0 = sin(cx) / cx
    sbesf(1) = (stemp0 / cx - cos(cx) / cx) / stemp0
     do 10 n = 1, limj - 1
     rn = n + n + 1
     sbesf(n + 1) = (rn / cx) - (1.0e0_knd / sbesf(n))
10    continue
    go to 60
20   continue
!
!  for c*x < lim, use backward recursion to
!  get fcn. ratios:
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
    stemp0 = 0.0e0_knd
    if(lim1 <= limj) go to 40
     do 30 n = lim1, limj,-1
     rn = n + n + 1
     stemp1 = 1.0e0_knd / (rn / cx - stemp0)
     stemp0 = stemp1
30    continue
40   sbesf(limj) = stemp0
     do 50 n = limj - 1, 1,-1
     rn = n + n + 1
     sbesf(n) = 1.0e0_knd / (rn / cx - sbesf(n + 1))
50    continue
60   continue
!
!  for all c*x, calculate the amplitude and sign scale
!  factors by forward operation on the Bessel function
!  ratios.
    stemp0 = sin(cx) / cx
    stemp1 = stemp0 / cx - cos(cx) / cx
    ibese(1) = int(log10(abs(stemp0)))
    sbesn(1) = stemp0 * 10.0e0_knd ** (-ibese(1))
    if(abs(sin(cx)) < 0.5e0_knd .and. cx > 1.0e0_knd) go to 70
    sbesn(2) = sbesn(1) * sbesf(1)
    ibese(2) = int(log10(abs(sbesn(2))))
    sbesn(2) = sbesn(2) * 10.0e0_knd ** (-ibese(2))
    ibese(2) = ibese(2) + ibese(1)
    go to 80
70   ibese(2) = int(log10(abs(stemp1)))
    sbesn(2) = stemp1 * 10.0e0_knd ** (-ibese(2))
    sbesf(1) = stemp1 / stemp0
80   continue
     do 90 n = 3, maxlp
     sbesn(n) = sbesn(n - 1) * sbesf(n - 1)
     ibese(n) = log10(abs(sbesn(n)))
     sbesn(n) = sbesn(n) * 10.0e0_knd ** (-ibese(n))
     ibese(n) = ibese(n) + ibese(n - 1)
90    continue
!
!  calculate the ratios of the first derivatives of successive
!  Bessel functions using corresponding function ratios
     do 100 n = 1, limj
     rn = n - 1
     sbesdf(n) = (cx - (rn + 2.0e0_knd) * sbesf(n)) / (rn - cx * sbesf(n))
100    continue
!
!  calculate the ratios of the first derivative to the corresponding
!  spherical Bessel function
     do 110 n = 1, maxlp
     rn = n - 1
     sbesdr(n) = (rn / cx) - sbesf(n)
110    continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
     do 120 n = limj, 2,-1
     sbesf(n) = sbesf(n - 1) * sbesf(n)
     sbesdf(n) = sbesdf(n - 1) * sbesdf(n)
120    continue
    return
    end subroutine
!
!
    subroutine sphneu (c, x, limn, maxn, maxlp, sneuf, sneun, ineue, sneudf, &
              sneudr)
!
!  purpose:     To calculate ratios of spherical Neumann functions
!               and ratios of their first derivatives for given c and x.
!               to calculate the Neumann function characteristics
!               and exponents. To calculate ratios of the first
!               derivatives to corresponding functions.
!
!  parameters:
!
!     input:    c      : c
!               x      : x
!               limn   : the number of spherical Neumann function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxn   : dimension of sneuf and sneudf arrays
!               maxlp  : the number of values of scale factors
!                        that are calculated
!
!     output:   sneuf  : ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : characteristic for the spherical
!                        Neumann functions
!               ineue  : exponent for the spherical
!                        Neumann functions
!               sneudf : ratios of first derivatives of successive
!                        spherical Neumann functions of the same parity
!               sneudr : ratios of first derivatives of spherical
!                        Neumann functions to the corresponding
!                        function
!
!    use param
!
!  real(knd) scalars and arrays
    real(knd) c, cx, rn, rnn, stemp0, stemp1, x
    real(knd) sneudf(maxn), sneudr(maxlp), sneuf(maxn), sneun(maxlp)
!
!  integer arrays
    dimension ineue(maxlp)
!
!  compute first kind ratios of Neumann functions and ratios
!  of their first derivatives
!
!        sneuf(k)=y(n=k,c*x)/y(n=k-2,c*x)
!        sneun(k)=(y(n=k-1),c*x)*10.0e0_knd**(-ineue(k))
!        sneudf(k)=y'(n=k,c*x)/y'(n=k-2,c*x)
!        sneudr(k)=(y'(n=k-1),c*x)/y(n=k-1),c*x))
!
!  use forward recursion to compute function ratios
!
!       y(n+1,c*x)/y(n,c*x)=(2*n+1)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
!
!  compute derivative ratios at same time using function ratios.
    cx = c * x
    stemp0 = -cos(cx) / cx
    stemp1 = stemp0 / cx - sin(cx) / cx
    sneuf(1) = stemp1 / stemp0
    sneudf(1) = -(cx - 2.0e0_knd * sneuf(1)) / (cx * sneuf(1))
     do 10 n = 1, limn - 1
      rn = n
      rnn = n + n + 1
      sneuf(n + 1) = rnn / cx - 1.0e0_knd / sneuf(n)
      sneudf(n + 1) = (cx - (rn + 2.0e0_knd) * sneuf(n + 1)) / (rn- &
            cx * sneuf(n + 1))
10     continue
     sneuf(limn + 1) = 0.0e0_knd
     sneuf(limn + 2) = 0.0e0_knd
     sneudf(limn + 1) = 0.0e0_knd
     sneudf(limn + 2) = 0.0e0_knd
!
!  calculate the characteristic and exponent
!  by forward operation on the Neumann function ratios:
    ineue(1) = int(log10(abs(stemp0)))
    sneun(1) = stemp0 * 10.0e0_knd ** (-ineue(1))
    ineue(2) = int(log10(abs(stemp1)))
    sneun(2) = stemp1 * 10.0e0_knd ** (-ineue(2))
    nlimit = maxlp
    if(maxlp > limn + 1) nlimit = limn + 1
     do 20 n = 3, nlimit
     sneun(n) = sneun(n - 1) * sneuf(n - 1)
     ineue(n) = int(log10(abs(sneun(n))))
     sneun(n) = sneun(n) * 10.0e0_knd ** (-ineue(n))
     ineue(n) = ineue(n) + ineue(n - 1)
20    continue
!
!  calculate the ratios of the first derivatives to the corresponding
!  spherical Neumann functions
     do 30 n = 1, nlimit
     rn = n - 1
     sneudr(n) = (rn / cx) - sneuf(n)
30    continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
     do 40 n = limn, 2,-1
     sneuf(n) = sneuf(n - 1) * sneuf(n)
     sneudf(n) = sneudf(n - 1) * sneudf(n)
40    continue
    return
    end subroutine
! end module prolate_swf

! Link to C
! Double precision
!subroutine profcn_cpp_interface(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
!           r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
!           s1c, is1e, s1dc, is1de, naccs) bind(C)
!
!    use prolate_swf
!!    use param
!    implicit none
!
!    ! Types must match Fortran and C
!    real(knd), intent(in) :: c, x1, arg(*)
!    integer, intent(in) :: m, lnum, ioprad, iopang, iopnorm, narg
!    real(knd), intent(out) :: r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
!                              s1c(lnum, narg), s1dc(lnum, narg)
!    integer, intent(out) :: ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
!                            is1e(lnum, narg), is1de(lnum, narg), naccr(lnum), naccs(lnum, narg)
!
!    call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
!                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
!                s1c, is1e, s1dc, is1de, naccs)
!
!end subroutine profcn_cpp_interface


! C interface - changes name based on USE_QUAD flag
#ifdef USE_QUAD
subroutine profcn_cpp_interface_quad(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
#else
subroutine profcn_cpp_interface(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
#endif
           r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
           s1c, is1e, s1dc, is1de, naccs) bind(C)

    use iso_c_binding
    implicit none

    ! Thin ISO C binding wrapper used by psms.cpp for single-order requests.
    real(knd), intent(in) :: c, x1, arg(*)
    integer, intent(in) :: m, lnum, ioprad, iopang, iopnorm, narg
    real(knd), intent(out) :: r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                              s1c(lnum, narg), s1dc(lnum, narg)
    integer, intent(out) :: ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                            is1e(lnum, narg), is1de(lnum, narg), naccr(lnum), naccs(lnum, narg)

    call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs)

#ifdef USE_QUAD
end subroutine profcn_cpp_interface_quad
#else
end subroutine profcn_cpp_interface
#endif

! Batched C interface
#ifdef USE_QUAD
subroutine profcn_cpp_interface_batch_quad(c, m_start, m_count, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
#else
subroutine profcn_cpp_interface_batch(c, m_start, m_count, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
#endif
           r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
           s1c, is1e, s1dc, is1de, naccs) bind(C)

    use iso_c_binding
    implicit none

    ! Thin ISO C binding wrapper used by psms.cpp for batched m-block requests.
    real(knd), intent(in) :: c, x1, arg(*)
    integer, intent(in) :: m_start, m_count, lnum, ioprad, iopang, iopnorm, narg
    real(knd), intent(out) :: r1c(lnum, m_count), r1dc(lnum, m_count), &
                              r2c(lnum, m_count), r2dc(lnum, m_count), &
                              s1c(lnum, narg, m_count), s1dc(lnum, narg, m_count)
    integer, intent(out) :: ir1e(lnum, m_count), ir1de(lnum, m_count), &
                            ir2e(lnum, m_count), ir2de(lnum, m_count), naccr(lnum, m_count), &
                            is1e(lnum, narg, m_count), is1de(lnum, narg, m_count), &
                            naccs(lnum, narg, m_count)

    call profcn_batch(c, m_start, m_count, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                      r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                      s1c, is1e, s1dc, is1de, naccs)

#ifdef USE_QUAD
end subroutine profcn_cpp_interface_batch_quad
#else
end subroutine profcn_cpp_interface_batch
#endif

#ifdef USE_QUAD
end module prolate_swf_quad
#else  
end module prolate_swf
#endif
