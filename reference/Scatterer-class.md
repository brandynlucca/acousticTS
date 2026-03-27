# Scatterer-class object for target strength estimation

The `acousticTS` package uses a variety of defined S4-class objects
comprising different types of scatterers, such as fish with gas-filled
swimbladders
([SBF](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md))
and fluid-like crustaceans
([FLS](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)).

## Slots

- `metadata`:

  List containing relevant metadata

- `model_parameters`:

  Model parameters necessary for predicting TS

## Data Organization

- `metadata`::

  A `list` comprising any identifying information associated with the
  scatterer. The default metadata entry includes `ID` that uses a
  default `character` value of `"UID"` (i.e. `ID = "UID"`). This can
  otherwise be formatted in any manner for book keeping purposes.

- `body/bladder`::

  A `list` that includes information relevant to the scatterer's
  position vector, material properties, tilt/orientation, etc. For some
  scatterers, this may only include `body`, but other targets may have
  an additional parameter such as `bladder`. Generally, each entry
  includes:

  - rpos: the relevant position vector (r₀) that includes axes such as
    `x`, `y`, `z`, etc., that depend on the type of scatterer being
    used.

  - `radius`: in some cases this includes the radius measurements for
    each cylinders depending on the type of scatterer object.

  - `theta`: the orientation of the scatterer relative to the
    transmitting transducer or sound source (θ_(animal)) that can be
    represented either by degrees or radians, although all functions
    require radians.

  - `g, h`: material properties that represent the density and sound
    speed contrasts (g and h, respectively) relative to the
    ambient/surrounding fluid. Some targets may instead have standard
    sound speed (c_(animal), m s⁻¹) and density (ρ_(animal), kg m³).

- `shape_parameters`::

  A `list` that includes metadata pertaining to the shape of the
  scatterer and any other features of interest (e.g. gas-filled
  swimbladder). Generally, each entry includes: overall body length, the
  number of discrete cylinders that make up the shape (if applicable),
  and the units related to both θ_(animal) (e.g. rad, °) and length
  (e.g. mm, m).

- `model_parameters`::

  A `list` that contains relevant model parameterization once an object
  has been initialized for modeling σ_(bs). This is typically broken up
  into three categories:

  - `parameters`: A `list` that includes information such as frequency
    (Hz), acoustic wavenumber (i.e. k), etc.

  - `medium`: A `data.frame` including information such as the material
    properties of the ambient medium.

  - `scatterer`: A `list` containing summarized information used to
    parameterize certain scattering models.

- `model`::

  A `list` that collects model results from one or more models in the
  linear domain (i.e. σ_(bs)).

## Supported Scatterers

- `Elastic-based scatterers`
  ([ELA](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md))

- `Composite multi-component scatterers`
  ([CSC](https://brandynlucca.github.io/acousticTS/reference/CSC-class.md))

- `Backboned fish`
  ([BBF](https://brandynlucca.github.io/acousticTS/reference/BBF-class.md))

- `Calibration spheres`
  ([CAL](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md))

- `Elastic-shelled scatterers`
  ([ESS](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md))

- `Fluid-like scatterers`
  ([FLS](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md))

- `Gas-filled scatterers`
  ([GAS](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md))

- `Swimbladdered fish`
  ([SBF](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md))
