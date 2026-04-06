# Choosing a model

## Introduction

Model selection in this package mirrors the broader literature: match
geometry, material contrast, and anatomical complexity before worrying
about computational detail ([**jech_etal_2015?**](#ref-jech_etal_2015);
[**clay_horne_1994?**](#ref-clay_horne_1994);
[**stanton_acoustic_1996?**](#ref-stanton_acoustic_1996)).

There is no single deterministic rule that maps every scatterer to
exactly one model. In practice, model choice depends on at least five
questions:

1.  What is the gross geometry of the target?
2.  Is the target fluid-like, gas-filled, rigid, or elastic-shelled?
3.  Is a geometry-matched exact or modal-series solution available?
4.  Is the main goal physical fidelity, computational speed, or
    comparative screening?
5.  Is unresolved variability important enough that a fully coherent
    deterministic model is too sharp?

The flowchart below is therefore intended as a screening guide rather
than as a strict decision rule. It identifies which models are typically
the most natural starting points and where multiple models may be
reasonable for the same target. Its purpose is to narrow the candidate
set while keeping the physical tradeoffs visible.

The most important principle is that model selection should follow the
target and the scientific question, not just model familiarity. A model
can be mathematically sophisticated and still be the wrong starting
point if its geometry, boundary assumptions, or approximation regime do
not match the target. Conversely, a simpler model can be the better
first choice if it captures the dominant physics of the question being
asked.

![Revised general model-selection
flowchart](model-selection-flowchart.png)[](https://brandynlucca.github.io/acousticTS/articles/krm/index.md "KRM overview")[](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md "BBFM overview")[](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md "HPA overview")[](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md "SPHMS overview")[](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md "Calibration overview")[](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md "FCMS overview")[](https://brandynlucca.github.io/acousticTS/articles/essms/index.md "ESSMS overview")[](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md "TRCM overview")[](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md "BCMS overview")[](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md "ECMS overview")[](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md "VESM overview")[](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md "SDWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md "DWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/psms/index.md "PSMS overview")[](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md "PCDWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")

The chart separates composite internal structure, canonical geometry,
and asymptotic screening into distinct decisions rather than folding
them into a single forced branch. That leaves room for families such as
`BBFM`, `BCMS`, `ECMS`, `PCDWBA`, `VESM`, and the shape-specific `TMM`
branches, and it also makes clearer that `HPA`, `SDWBA`, and `TRCM` are
often neighboring screening or regime-specific alternatives rather than
absolute replacements for the more geometry-matched families.
Experimental families are shown in a separate color so they remain
visible without being confused with the longer-established model paths,
while gray boxes identify future-development families that have been
discussed but are not public package models.

Boundary condition is often as decisive as geometry in this decision
process. A spherical pressure-release cavity, a fluid-filled sphere, and
an elastic shell may look similar at gross scale while belonging to
different model families because the interface physics is different in
each case.

## Core principle

When a geometry-matched exact or modal-series model exists and its
assumptions are acceptable, that is usually the best first choice for
mechanistic interpretation. Approximate models become most useful when
the target geometry is not one of the canonical separable shapes, when
the target is segmented or anatomically composite, when the acoustic
regime favors asymptotic rather than modal physics, or when a faster
screening model is more useful than a full boundary-value solution.

This principle can be restated more concretely. Model choice usually
benefits from asking four questions in order:

1.  What geometry is actually defensible?
2.  What boundary or material interpretation is physically important?
3.  Is there a geometry-matched model that represents that boundary
    physics directly?
4.  If not, which approximation family best captures the dominant
    scattering mechanism?

That sequence helps prevent a common mistake, namely choosing a model
because its name sounds close to the target rather than because its
physical assumptions match the problem.

## Geometry-first screening

### Sphere family

If the target is well represented as a sphere, the sphere-family models
are the natural starting point.
[SPHMS](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-theory.md)
is the appropriate first choice when the object is spherical and the
relevant boundary condition is fluid-filled, rigid, or pressure-release.
[ESSMS](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md)
becomes appropriate when the object is a fluid-filled elastic shell and
shell elasticity is physically important.
[calibration/SOEMS](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-theory.md)
is the natural choice for a solid calibration sphere rather than a
biological or fluid sphere.
[TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md)
is relevant when the user wants to stay inside a transition-matrix
workflow for later angular or orientation-based post-processing, even
though `SPHMS` remains the more direct exact spherical model.
[HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
belongs to this family only as a compact asymptotic alternative when a
full modal solution is not needed.

The main point in the spherical branch is that “sphere” alone is not
enough to choose a model. A gas-like or fluid-like spherical target, an
elastic shell, and a solid calibration sphere may all share a gross
geometry while requiring very different physics. Readers should
therefore treat geometry as the first screen, not as the full answer.

### Prolate spheroid

If the target is elongated but still smooth and well approximated by a
prolate spheroid,
[PSMS](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md)
is usually the most geometry-matched exact choice. It should be
preferred when the prolate spheroidal geometry is a meaningful physical
approximation and the boundary condition lies within the model’s scope.
[TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md)
is the relevant alternative when a retained transition-matrix
representation is the point of the calculation rather than only a
monostatic spectrum.
[HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
can be useful for faster trend-level approximation when a full
spheroidal modal calculation is not required.

This branch is most useful when the reader is confident that the target
is smoothly elongated and that a spheroidal description is more faithful
than a cylindrical or arbitrary segmented one. If the gross elongation
is right but the ends, curvature, or internal organization matter
strongly, then a more approximate but more physically relevant model
family may be preferable to a mathematically elegant but poorly matched
canonical geometry.

### Oblate spheroid

If the target is flattened rather than elongated and is better
represented as an oblate spheroid, the package path is
[TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md).
That branch is still experimental, but it is the active shape-matched
route for single-target oblate scattering in the package. The chart also
shows `OSMS` in gray because a true oblate modal-series companion model
remains a future-development target rather than a public package model.

This oblate branch matters because it separates two different ideas that
are easy to conflate. The first is “can the target be represented as an
oblate spheroid?” The second is “does the package already have a
dedicated oblate modal-series family?” The answer is yes to the first
question and no to the second, which is why the flowchart routes the
public choice through `TMM` and reserves `OSMS` for future development.

### Finite cylinder

If the target is better approximated by a straight finite cylinder than
by a sphere or spheroid,
[FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md)
is the natural modal-series choice. It is most appropriate when the
cross-section is circular, the cylinder is reasonably straight, and the
separated finite-cylinder reduction is defensible.
[TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md)
is shown for cylinders as a guarded experimental branch rather than as a
co-equal default, because its documented supported scope is narrower
than the sphere, oblate, and prolate `TMM` branches.
[HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
is useful when only a compact approximate response is needed, while
[TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md)
is more appropriate when the dominant physics is high-frequency specular
and through-body ray interference rather than modal structure.

The cylinder branch often involves an especially important judgment
call. A reader has to decide whether they care most about modal
boundary-condition physics, broad asymptotic trends, or a high-frequency
ray picture. Those are not minor stylistic differences. They correspond
to different dominant scattering mechanisms.

## Approximation-family screening

### Weakly scattering elongated fluid-like bodies

If the target is elongated, fluid-like, and weakly contrasting relative
to the surrounding water, the
[DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
family is usually the general-purpose starting point.
[DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
is appropriate when the body can be represented as a weakly scattering
centerline plus local cross-sections and deterministic coherence is
acceptable.
[SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md)
should be preferred when the same weak-scattering geometry is
appropriate but unresolved posture, roughness, or internal variability
is expected to smear the interference pattern.

This is one of the most important non-unique branches in the package.
For many zooplankton-like or krill-like bodies, `DWBA` and `SDWBA` are
both defensible. The difference is not geometry alone. It is whether
phase variability is a real part of the scientific question. If the
reader wants the deterministic response of a specified body posture,
`DWBA` is often the more direct choice. If the reader wants the response
expected after unresolved variability has already blurred the fine
interference pattern, `SDWBA` is often closer to the quantity of
interest.

### High-frequency locally cylindrical targets

If the target is acoustically large enough that specular reflection
dominates and a locally cylindrical high-frequency description is
reasonable,
[TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md)
becomes relevant. It is appropriate when the dominant physics is well
described by two coherent paths, a near-interface reflection and a
through-body return.
[DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
should be preferred when weak-scattering volume integration is more
physically appropriate than a two-ray asymptotic picture, while
[FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md)
should be preferred when a geometry-matched cylindrical modal solution
is both available and desirable.

This branch is often where acoustic size matters most. A target may have
cylindrical geometry at all frequencies, but the most useful model still
changes depending on whether the reader wants weak-scattering behavior,
modal behavior, or a high-frequency specular description.

### Fish-like bodies with a swimbladder

If the object is a segmented fish body with a gas-filled swimbladder,
[KRM](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md)
is usually the most natural model family in the current library. It
should be used when the target has a fish-like body-plus-bladder
structure and a hybrid Kirchhoff-plus-modal treatment is physically
sensible.
[DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
or
[SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md)
may be more appropriate when the body is weakly scattering and the
bladder is not the dominant feature.
[TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md)
or
[HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
are better treated as simplified trend-level or high-frequency
approximations, not as replacements when swimbladder physics is central.

This is another case where internal structure matters as much as outer
shape. A fish body may be elongated enough to invite comparison with
other elongated-body models, but if the scientific question is driven
primarily by the interaction between weakly scattering tissue and a
strong internal inclusion, KRM is often the more physically transparent
first choice.

## When several models are reasonable

### Sphere or near-sphere

For a sphere or near-sphere, the most common overlap is between `SPHMS`
and `HPA`, where `SPHMS` should be preferred for mechanistic
boundary-condition fidelity and `HPA` for a faster asymptotic
approximation. Another overlap occurs between `SPHMS` and `ESSMS`, but
`ESSMS` is only appropriate when shell elasticity genuinely matters.

### Elongated fluid-like bodies

For elongated fluid-like bodies, the most common comparison is `DWBA`
versus `SDWBA`, where `SDWBA` should be preferred when unresolved phase
variability is part of the expected physics. `DWBA` versus `TRCM` is a
choice between weak-scattering volume physics and high-frequency two-ray
specular physics. `DWBA` versus `HPA` is a choice between
geometry-specific body structure and a compact trend-level
approximation.

### Canonical elongated shapes

For canonical elongated shapes, `PSMS` versus `HPA` is a choice between
a true prolate spheroidal model and a rougher asymptotic descriptor.
`FCMS` versus `TRCM` is a choice between modal boundary-condition
physics and high-frequency ray-interference behavior. `FCMS` versus
`HPA` is a choice between resolving modal structure and retaining only
the broad scaling and orientation trend.

Whenever several models are reasonable, the most informative approach is
often comparison rather than forced selection. Running two defensible
models on the same geometry, material properties, and frequency grid can
reveal whether the scientific conclusion is robust to model structure or
whether it depends strongly on the assumed scattering mechanism. That is
often more valuable than trying to declare a winner before the
comparison is made.

## Practical rule set

If a user wants a quick screening workflow, a simple rule set is usually
enough. Start with the closest canonical geometry you can defend
physically. If an exact or modal-series model exists for that geometry,
use it first unless speed or simplification is the main goal. If the
target is elongated and weakly scattering, start with `DWBA`. If the
same target also shows unresolved phase variability, compare against
`SDWBA`. If the target is high-frequency and specular with locally
cylindrical physics, test `TRCM`. If the target is fish-like and
swimbladder-dominated, test `KRM`. Use `HPA` when a broad asymptotic
approximation is more useful than a geometry-matched exact model.

That rule set is intentionally conservative. It aims to get the reader
to a defensible first model, not to guarantee that only one model should
ever be run. In many cases, the better habit is to begin with the
best-matched first model and then run one neighboring alternative to see
how sensitive the interpretation is to the choice of model family.

## What this flowchart does not decide

The flowchart does not resolve every scientific choice. It does not
decide whether the available morphology is good enough for a given
geometry simplification, whether a weak-scattering assumption is
quantitatively justified, whether a shell should be treated as rigid,
fluid, or elastic, or whether the question is about absolute prediction,
comparative trend, or sensitivity analysis.

It also does not decide what level of agreement is good enough once
multiple models have been run. That part of the workflow still requires
judgment about whether the comparison should be made in `TS`, in
`sigma_bs`, over a narrow band, over a broad band, or against particular
features such as peaks, nulls, or orientation trends.

Those are still modeling judgments. The purpose of the decision tree is
to narrow the candidate set to the models whose assumptions are closest
to the target and to make clear where a user should compare more than
one model rather than commit to only one.
