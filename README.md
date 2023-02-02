# acousticTS

Acoustic target strength (TS) represents the intensity of an echo returning from an individual scatterer such as bubbles, fish, or zooplankton. TS can be used to convert integrated or volumetric backscatter collected from fisheries acoustic surveys into units of number density (e.g. animals per m^3), abundance (e.g. number of animals), and biomass (e.g. kg). This parameter can also be used to aid in classifying backscatter, such as separating likely echoes of large predatory fish (e.g. adult cod) from smaller prey (e.g. shrimp). One way to estimate TS is to use physics-based models to calculate theoretical TS that comprise exact and approximate solutions as well as analytical approaches. The models provided can help provide TS estimates over broad statistical distirbutions of model parameters.

_General DOI_
https://doi.org/10.5281/zenodo.7600659

_Latest release DOI_
[![DOI](https://zenodo.org/badge/161965429.svg)](https://zenodo.org/badge/latestdoi/161965429)

## Installation

You can install the development version of acousticTS like so:

``` r
devtools::install_github("brandynlucca/acousticTS@test-branch")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(acousticTS)
## Let's create a calibration sphere 
cal_sphere <- cal_generate()
## The default inputs here are a 38.1 mm diameter and a tungsten carbide (WC) material.
## Let's define frequency
frequency <- c(38e3, 70e3, 120e3, 200e3)
# Calculate TS; update original CAL object
cal_sphere <- target_strength(object = cal_sphere,
                              frequency = frequency,
                              model = "calibration")
# Extract model results
model_results <- extract(cal_sphere, "model")
# Print the results
print(model_results)
```
