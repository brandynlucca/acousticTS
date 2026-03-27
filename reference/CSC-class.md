# Composite scatterer (CSC) object/class.

A shared parent S4 class for composite scatterers that retain a primary
body component together with one or more additional anatomically
distinct sub-components. The `CSC`-class is intended to provide modular
scaffolding for multi-component biological targets, such as
swimbladder-bearing fish and swimbladder-less fish represented as flesh
plus backbone. It stores the primary body, a general component list, the
shared model output slot, and shape metadata without implying any one
specific resonant or elastic inclusion type.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md),
[SBF](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md),
[BBF](https://brandynlucca.github.io/acousticTS/reference/BBF-class.md)
