## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new submission.

## Test environments

* Local Windows 11 x64 (R 4.5.2)
    * `R CMD build`
    * `R CMD check --as-cran`
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
    * clang-asan
    * clang-ubsan
    * clang16
    * clang17
    * clang18
    * clang19
    * clang20
    * donttest
    * gcc-asan
    * gcc13
    * gcc14
    * gcc15
    * intel
    * mkl
    * nold
    * noremap
    * ubuntu-clang
    * ubuntu-gcc12
    * ubuntu-next
    * ubuntu-release
    * valgrind
    * linux (R-devel)
    * macos (R-devel)
    * macos-arm64 (R-devel)
    * windows (R-devel)
    * m1-san (R-devel)

## Note on R-devel
* Some additional experimental R-devel environments encountered upstream dependency installation failures under the current R-devel toolchain, not failures in acousticTS itself.

## Reverse dependencies

There are currently no downstream dependencies for this package.