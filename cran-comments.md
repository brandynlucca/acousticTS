## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* Local Windows 11 x64, R 4.5.2
    * `R CMD build`
    * `R CMD check --as-cran`
* GitHub Actions CI
    * ubuntu-latest (release)
    * ubuntu-clang (release)
    * ubuntu-latest (oldrel-1)
    * ubuntu-latest (no-suggests)
    * ubuntu-latest (devel)
    * macos-latest (release)
    * windows-latest (release)
* R-hub v2:
    * ubuntu-release
    * gcc15

## Note on R-devel
* Some additional experimental R-devel environments encountered upstream dependency installation failures under the current R-devel toolchain, not failures in `acousticTS` itself.

## Reverse dependencies

There are currently no downstream dependencies for this package.