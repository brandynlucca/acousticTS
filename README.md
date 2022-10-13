---
output:
  html_document:
    keep_md: true
---
# acousticTS

Acoustic backscatter from a single target or organism is expressed as the 
intensity of an echo typically denoted as the *backscattering cross-section* 
(&sigma;~bs~, m^2^). Target strength (TS, dB re. 1 m^2^) is the logarithmic 
representation of &sigma;~bs~ where: $TS = 10~log_{10}(\sigma_{bs})$. TS can 
be used to convert integrated (e.g. nautical area scattering coefficient, S~A~,
dB re. 1(m^2^ nmi^-2^) or volumetric backscatter (e.g. S~v~, dB re. 1 m^-1^) 
collected from fisheries acoustic surveys into units of number density, such as 
the volumetric density of a fish school (i.e. animals m^-3^). This parameter can 
also aid in classifying backscatter based on the multifrequency response of 
targets, such as separating likely echoes of large predatory fish from smaller
prey. While there are several approaches for estimating TS, one common method is
to apply physics-based models to predict theoretical TS that comprise exact and
approximate solutions. The models provided in the `acousticTS` package can help
provide TS estimates parameterized using broad statitsical distributions of 
inputs. This package is in a constant state of development with updates to 
the available model library, computational efficiency, and quality-of-life 
improvements.

## Installation

You can install the current released version of acousticTS via: 

```r
devtools::install_github("brandynlucca/acousticTS")
```

Or you can install the development version of acousticTS like so:

```r
devtools::install_github("brandynlucca/acousticTS@test-branch")
```

## Examples

Below are two examples with predicting TS for a tungsten carbide calibration 
sphere and a sardine with a gas-filled swimbladder.

### Kirchoff Ray-Mode approximation for a Sardine with a gas-filled swimbladder


```r
### Call in the library
library(acousticTS)
### Call in the built-in sardine shape dataset
data(sardine)
### Inspect the object
print(sardine)
```

```
## SBF object 
##  Gas-filled swimbladdered scatterer 
##  ID: UID 
##  Body length: 0.21 (n = 379 cylinders) 
##  Bladder length: 0.085 (n = 154 cylinders) 
##  Body orientation (relative to transducer axis): 1.571 
##  Bladder orientation (relative to transducer axis): 1.571 
##  Material properties (body): density = 1070 (kg/m^3); sound speed = 1570 (m/s) 
##  Material properties (bladder): density = 1.24 (kg/m^3); sound speed = 345 (m/s)
```

```r
plot(sardine)
```

![](test_md_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
### We will now define a frequency range to predict TS over
frequency <- seq(1e3, 400e3, 1e3)
### And now we use the target_strength(...) function to model TS for this fish
sardine <- target_strength(sardine, 
                           frequency = frequency, 
                           model = "KRM")
### Extract model results
sardine_ts <- extract(sardine, "model")$KRM
### This prints out the summed linear backscatter from the fluid-like bodily
### tissues, soft tissue representing the swimbladder, the summation of these
### two tissue-types, and TS.
head(sardine_ts)
```

```
##                       f_fluid                    f_soft
## 1 -7.283091e-05-4.958492e-04i -0.004076481+0.002142353i
## 2 -1.756594e-04-6.856564e-04i -0.006913023+0.003296627i
## 3 -2.785704e-04-8.198628e-04i -0.009070826+0.003887462i
## 4 -3.741593e-04-9.262943e-04i -0.010796866+0.004113333i
## 5 -4.602035e-04-1.016880e-03i -0.012217241+0.004083154i
## 6 -5.362574e-04-1.097972e-03i -0.013402666+0.003863448i
##                        f_bs        TS
## 1 -0.004149312+0.001646504i -47.00541
## 2 -0.007088682+0.002610970i -42.43618
## 3 -0.009349396+0.003067599i -40.14029
## 4 -0.011171026+0.003187039i -38.69830
## 5 -0.012677445+0.003066274i -37.69246
## 6 -0.013938924+0.002765477i -36.94775
```

```r
### Plot results
plot(x = frequency * 1e-3,
     y = sardine_ts$TS,
     type = 'l',
     xlab = "Frequency (kHz)",
     ylab = expression(Target~strength~(dB~re.~1~m^2)))
```

![](test_md_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

### Calibration sphere


```r
### Let's create a calibration sphere 
### Default inputs here are a 38.1 mm diameter and a tungsten carbide 
### (WC) material properties.
cal_sphere <- cal_generate()
### We will use the same frequency range as the previous example
### Calculate TS
cal_sphere <- target_strength(object = cal_sphere,
                              frequency = frequency,
                              model = "calibration")
### Extract model results
calibration_ts <- extract(cal_sphere, "model")
### Plot the results
plot(x = frequency * 1e-3,
     y = calibration_ts$TS,
     type = 'l',
     xlab = "Frequency (kHz)",
     ylab = expression(Target~strength~(dB~re.~1~m^2)))
```

![](test_md_files/figure-html/unnamed-chunk-2-1.png)<!-- -->
