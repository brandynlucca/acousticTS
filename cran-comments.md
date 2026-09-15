## Resubmission

This is a resubmission (version 2.0.6) addressing the latest review feedback:

* Removed the `\dontrun{}` wrapper from the short `plot(krill)` example. Removed
  the two optional quad precision examples, which depend on native support
  in the package build, while retaining runnable double precision examples
  for Smn and Rmn. The precision argument documentation explains this build
  requirement. No `\dontrun{}` or `\donttest{}` wrappers remain in the help examples.
* Added Arnie Lee Van Buren `[ctb, cph]` and Jeffrey E. Boisvert `[ctb]` to
  `Authors@R`, with comments identifying them specifically as authors of the
  upstream `prolate_swf` code adapted for this package. They did not participate
  directly in development of `acousticTS`. This upstream code was heavily modified for the package. The upstream copyright and MIT license
  notices remain preserved in `inst/COPYRIGHTS` and the source attribution
  headers. The `inst/COPYRIGHTS` file also clarifies the scope of the attribution.
* Simplified License to GPL-3 and removed the redundant top-level `LICENSE`
  file. The package's GPL-3 licensing is unchanged, and the upstream MIT
  license notice is retained in `inst/COPYRIGHTS`.
* Added a compiler-capability check for LLVM Flang builds that diagnose R's global `-Wall` setting as an unused command-line argument. When that exact behavior is detected, the generated package Makevars applies the supported diagnostic-suppression flag to the package's Fortran compilation only.

The words reported as possibly misspelled by the CRAN incoming check are correctly spelled author surnames (Lucca, MacLennan, and Simmonds), Latin citation terms (et al.), and accepted fisheries-acoustics terminology (backscatter and scatterer).

The previous correction replacing two relative links in the installed
boundary-conditions vignette with absolute pkgdown URLs is retained.

## Validation for this resubmission

0 errors | 0 warnings | 0 notes

* Local Windows 11 x64 (R 4.5.2), from the source tarball with vignettes   enabled (`R CMD check --no-manual`).
* Local WSL (Ubuntu 26.04 LTS, R 4.5.2)
* Remote RStudio Server (Linux x86_64, R 4.5.2)
* Google Cloud Workstation (Linux container environment, R 4.5.2)
* GitHub Actions CI
    * ubuntu-latest (release)
    * ubuntu-clang (release)
    * ubuntu-latest (oldrel-1)
    * ubuntu-latest (no-suggests)
    * ubuntu-latest (devel)
    * macos-latest (release)
    * windows-latest (release)
* R-hub v2:
    * atlas
    * c23
    * clang16
    * clang17
    * clang18
    * clang19
    * clang20
    * clang21
    * clang22    
    * clang-asan
    * clang-ubsan
    * donttest
    * gcc-asan
    * gcc13
    * gcc14
    * gcc15
    * gcc16
    * intel
    * linux (R-devel)
    * lto
    * m1-san (R-devel)
    * macos (R-devel)
    * macos-arm64 (R-devel)
    * mkl
    * nold
    * noremap
    * ubuntu-clang
    * ubuntu-gcc12
    * ubuntu-next
    * ubuntu-release
    * valgrind   
    * vnu
    * windows (R-devel)   

## Previous submission: R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Previous submission: test environments

* Local Windows 11 x64 (R 4.5.2)
* Local WSL (Ubuntu 26.04 LTS, R 4.5.2)
* Remote RStudio Server (Linux x86_64, R 4.5.2)
* Google Cloud Workstation (Linux container environment, R 4.5.2)
* GitHub Actions CI
    * ubuntu-latest (release)
    * ubuntu-clang (release)
    * ubuntu-latest (oldrel-1)
    * ubuntu-latest (no-suggests)
    * ubuntu-latest (devel)
    * macos-latest (release)
    * windows-latest (release)
* R-hub v2:
    * atlas
    * c23
    * clang16
    * clang17
    * clang18
    * clang19
    * clang20
    * clang21
    * clang22    
    * clang-asan
    * clang-ubsan
    * donttest
    * gcc-asan
    * gcc13
    * gcc14
    * gcc15
    * gcc16
    * intel
    * linux (R-devel)
    * lto
    * m1-san (R-devel)
    * macos (R-devel)
    * macos-arm64 (R-devel)
    * mkl
    * nold
    * noremap
    * ubuntu-clang
    * ubuntu-gcc12
    * ubuntu-next
    * ubuntu-release
    * valgrind   
    * vnu
    * windows (R-devel)    

## Previous submission: note on R-hub `rchk`

The R-hub `rchk` job exits non-zero because its wrapper treats protection-balance diagnostics emitted for Rcpp's `Armor` and `Shield` helper headers as fatal. These diagnostics point to the installed Rcpp headers, not the package-authored C or C++ code. The acousticTS package contains no direct `PROTECT()` or `UNPROTECT()` calls. The additional "too many states" and `objdump` diagnostics are documented by `rchk` as ignorable. No package-level change is available or appropriate for these upstream diagnostics.

## Reverse dependencies

There are currently no downstream dependencies for this package.
