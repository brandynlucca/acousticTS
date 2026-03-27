# Plotting and Inspecting Results

## Introduction

Plot inspection is most useful when curves are read against known
benchmark shapes, resonances, and asymptotic regimes from the published
scattering literature ([Foote 1990](#ref-foote_spheres_1990); [Jech et
al. 2015](#ref-jech_etal_2015)).

The package stores geometry, parameters, and model outputs in nested
scatterer objects. That design is useful only if the contents can be
inspected and visualized efficiently. This article is the companion
guide to those object-level operations.

Inspection is more than a cosmetic last step. In practice, it is one of
the main safeguards against carrying an incorrect target, an unintended
model choice, or a misleading output representation into later analysis.
A clean line plot is not automatically a trustworthy result. Readers
often learn more from a short inspection pass than from immediately
producing a polished figure.

## Main tools

The main inspection tools are
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for shape and
model views, `show()` for quick object summaries, and
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
for pulling nested geometry, parameters, and model outputs. Together,
these functions provide a short inspection loop that can be repeated
whenever geometry, material properties, or model settings change.

These tools answer different questions.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is usually the
fastest way to judge whether something looks physically sensible at a
glance. `show()` gives a quick structural summary of what the object
contains.
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
is the most useful tool when the goal shifts from visual checking to
careful comparison, tabulation, or custom analysis. A robust workflow
usually uses all three rather than relying on only one.

## Shape plots versus model plots

The `plot.Scatterer` method can be used for both geometry and modeled
output, but those plots answer different questions. Shape plots answer
whether the object geometry and orientation are sensible. Model plots
answer whether the predicted frequency or acoustic-size response is
sensible. It is usually worth checking both before treating a model
result as trustworthy, because a well-behaved target-strength curve does
not rescue a physically incorrect shape, and a carefully constructed
shape does not guarantee that the stored model output is the one the
user intended to inspect.

Shape plots are especially useful for catching problems that are obvious
visually but hard to detect numerically. A rotated body may be pointing
the wrong way, a segmented object may be coarser than intended, or a
swimbladder may be positioned implausibly relative to the body. Those
issues often reveal themselves immediately in the geometry view.

Model plots answer a different set of questions. Is the frequency range
correct? Are the returned values in the expected domain? Does the curve
vary smoothly where it should, or show sharp structure where the chosen
model would make that plausible? Even here, the goal is not just to
admire the curve. It is to check whether the plotted behavior is
consistent with the target, the model family, and the expected acoustic
regime.

## Why `extract()` matters

The
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
accessor is one of the most important workflow utilities in the package
because the modeled outputs are nested within the object rather than
returned as flat tables by default. It is the most direct way to
retrieve stored quantities for custom comparison, plotting, or
downstream analysis. In practice,
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
is often the point where a user moves from package-level plotting into
project-specific summaries, tables, and comparative graphics.

It also matters because the package stores more than one kind of
information inside the same object. Geometry, material properties, model
parameters, and model outputs are not interchangeable, and
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
makes it possible to pull exactly the part that is relevant to the
current question. That reduces the risk of treating a quick plot as the
only inspection step when a direct look at the stored data frame or
nested component would reveal something important.

## A practical inspection sequence

The most useful inspection sequence is usually to begin with
`plot(..., type = "shape")`, because geometry and orientation mistakes
are easiest to catch visually. The next step is to use `show()` or
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
to confirm that the object contains the expected model entries and
parameter fields. Only after those checks is it worth building a final
model plot or exporting the result for external analysis.

This order matters because the nested object structure is designed to
preserve context. Geometry, material properties, and results stay
attached to the same target record. Inspection works best when it
follows that same structure rather than jumping directly to a final line
plot.

A practical sequence often looks like this:

1.  inspect the shape plot to confirm geometry, orientation, and
    component placement,
2.  inspect the object summary to confirm that the expected target class
    and stored components are present,
3.  extract the relevant model output to confirm that the requested
    model, frequency grid, and returned variables are actually there,
4.  plot the model response in the domain most relevant to the question,
    and
5.  only then move to comparison, export, or custom downstream graphics.

That sequence is deliberately simple, but it catches a large share of
avoidable workflow mistakes. Many apparent model problems turn out to be
geometry problems, parameterization problems, or output-selection
problems that become obvious during a careful inspection pass.

## What a good inspection actually checks

A good inspection step is not only about whether a figure renders
without error. It should help answer a few practical questions:

1.  Is this the object I think it is?
2.  Is this the model I meant to run?
3.  Are the returned variables the ones I need for the next step?
4.  Is the result being viewed in the right domain for the scientific
    question?

Those questions matter because a workflow can look smooth while still
answering the wrong question. A result may be correct in `TS` but
inconvenient for coherent addition that should be done in `sigma_bs`. A
plot may look plausible while hiding the fact that the object
orientation was not what the user intended. Inspection is the step that
keeps those mismatches visible.

## Inspection before comparison

Inspection is especially important before comparing models or targets.
If two curves disagree, the first task is to confirm that they were
generated from comparable objects, comparable frequency grids, and
comparable output domains. Without that check, a model comparison can
become a comparison of inconsistent bookkeeping instead of a comparison
of physics.

This is one reason
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
is so central. It allows readers to compare exactly what was stored,
rather than only what was drawn by a default plotting method. For
careful work, looking directly at the extracted results is often just as
important as looking at the plot.

## Related workflow pages

- [Running target-strength
  models](https://brandynlucca.github.io/acousticTS/articles/running-models/running-models.md)
- [Comparing models on the same
  target](https://brandynlucca.github.io/acousticTS/articles/comparing-models/comparing-models.md)
- [Working with real example
  data](https://brandynlucca.github.io/acousticTS/articles/example-data/example-data.md)

## References

Foote, K. G. 1990. “Spheres for Calibrating an Eleven-Frequency Acoustic
Measurement System.” *ICES Journal of Marine Science* 46 (3): 284–86.
<https://doi.org/10.1093/icesjms/46.3.284>.

Jech, J. Michael, John K. Horne, Dezhang Chu, David A. Demer, David T.
I. Francis, Natalia Gorska, Benjamin Jones, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.
