# Model library

## Introduction

The package now contains enough model families that the site works
better when the model pages are treated as a library rather than as a
flat list in the top navigation. This page is the main entry point for
that library.

### Model-status policy

Status tags are derived from the package’s internal validation registry
and used conservatively: - Benchmarked means the family has a documented
comparison against a canonical benchmark ladder or stored benchmark
values. - Validated means the family has a documented comparison against
at least one external implementation or software package. - Experimental
means the family is available to use, but its interface or validation
scope should still be treated as provisional. - Unvalidated means the
package site does not yet document either benchmark evidence or an
external implementation comparison for that family.

These statuses are not all mutually exclusive: - `Benchmarked` and
`Validated` are evidence tags and can appear together. - `Experimental`
is a lifecycle tag and can coexist with either evidence tag. -
`Unvalidated` is reserved for families lacking both evidence tags, so it
should not be combined with other evidence tags.

## Modal-series families

### [SPHMS](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md)

Benchmarked Validated

Spherical modal-series solution for canonical spherical targets.

### [FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md)

Benchmarked Validated

Finite-cylinder modal-series solution for straight cylindrical targets.

### [PSMS](https://brandynlucca.github.io/acousticTS/articles/psms/index.md)

Benchmarked Validated

Prolate-spheroidal modal-series solution for smooth elongated canonical
bodies.

### [SOEMS](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md)

Benchmarked Validated

Solid elastic spherical model used mainly for calibration spheres.

### [ESSMS](https://brandynlucca.github.io/acousticTS/articles/essms/index.md)

Unvalidated

Elastic-shelled spherical family for layered shell targets.

### [BCMS](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md)

Experimental Unvalidated

Bent-cylinder modal-series family for straight and uniformly bent
cylinders.

### [ECMS](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md)

Experimental Unvalidated

Elastic-cylinder modal-series family for fully elastic solid cylinders.

## Approximation and ray-based families

### [DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)

Benchmarked Validated

Weak-scattering elongated-body approximation for fluid-like targets.

### [SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)

Benchmarked Validated

Stochastic DWBA family for unresolved phase variability.

### [KRM](https://brandynlucca.github.io/acousticTS/articles/krm/index.md)

Benchmarked Validated

Kirchhoff-ray mode model for segmented fish-like body-plus-inclusion
targets.

### [HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)

Benchmarked Validated

High-pass approximation for compact asymptotic screening.

### [TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md)

Benchmarked

Two-ray cylindrical family for high-frequency locally cylindrical
targets.

### [PCDWBA](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)

Validated Experimental

Phase-compensated DWBA for bent weakly scattering targets.

## Composite and emerging families

### [BBFM](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)

Experimental Unvalidated

Composite flesh-plus-backbone family for swimbladder-less fish-like
targets.

### [VESM](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md)

Validated Experimental

Viscous-elastic layered-sphere family for gas-core, shell, and
viscous-layer targets.

### [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)

Benchmarked Validated Experimental

Single-target transition-matrix family for retained monostatic and
angle-dependent scattering products across supported canonical shapes.

## How to use this page

Each family page is intended to be the local navigation hub for that
model family. Use the family page to move between overview, theory,
implementation, and related package pages, rather than relying on a long
global theory or implementation menu.
