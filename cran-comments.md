## Submission

This is a maintenance update to acousticTS, version 2.0.8

This release corrects numerical edge cases in Bessel functions, spheroidal angular function indexing, and SVD solvers. It also improves numerical stability and compiler portability, and expands regression test coverage.

## R CMD check results

0 errors | 0 warnings | 0 notes

Other validation checks:

* Local Windows 11 x64 (R 4.5.2), from the source tarball with vignettes enabled (`R CMD check --no-manual`).
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

## Reverse dependencies

There are currently no downstream dependencies for this package. Checked on 2026-10-08.
