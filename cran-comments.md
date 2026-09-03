## Resubmission

This is a resubmission. In this version I have:

* Replaced the two relative links in the installed boundary-conditions
  vignette with absolute pkgdown URLs, resolving the invalid file URI note
  reported by CRAN.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Test environments

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

## Note on R-hub `rchk`

The R-hub `rchk` job exits non-zero because its wrapper treats protection-balance diagnostics emitted for Rcpp's `Armor` and `Shield` helper headers as fatal. These diagnostics point to the installed Rcpp headers, not the package-authored C or C++ code. The acousticTS package contains no direct `PROTECT()` or `UNPROTECT()` calls. The additional "too many states" and `objdump` diagnostics are documented by `rchk` as ignorable. No package-level change is available or appropriate for these upstream diagnostics.

## Reverse dependencies

There are currently no downstream dependencies for this package.
