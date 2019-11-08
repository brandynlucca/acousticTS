# acousticTS
### Calibration sphere target strength model
Current version: v0.7 [08 November 2019]

Comprises multiple Bessel functions: Bessel functions of the first [ja(...)] and second [ya(...)] kind, spherical Bessel functions of the first [jl(...)] and second [yl(...)] kind, first derivative of the spherical Bessel functions of the first [jd(...)] and second [yd(...)] kind, and the second derivative of the spherical Bessel function of the first kind [jdd(...)]. 

There are three target strength model functions for calibration spheres: 1) sphere.ts(...) which calculates the TS at a single frequency, 2) sphere.cw(...) which calculates the TS of a FM sweep at a given center frequency, and 3) sphere.spec(...) which calculates the TS-frequency spectrum of a desired bandwidth. 

There are a few support and wrapper functions as well. The sphere.spec_plot(...) function takes the output of sphere.spec(...) and plots the TS-frequency spectrum. Another is sphere_param(...) which allows the user to toggle through preset material property measurements (i.e., transversal/longitudinal sound speeds, density) of different types of calibration spheres (i.e., tungsten carbide, copper, stainless steel, and aluminum). 

### Citations

Amos D.E., "AMOS, A portable package for Bessel functions of a complex argument and nonnegative order", http://netlib.org/amos. 

MacLennan D. N. (1981). The theory of solid spheres as sonar calibration targets. Scottish Fisheries Research No. 22, Department of Agriculture and Fisheries for Scotland. 



## **To Do**
- [X] Fix required package imports 
- [ ] Create full README.md file to provide overview and walkthrough of package
- [X] Flesh out documentation for all functions; separate out each of the Bessel functions and their respective derivatives
- [X] Replace Python-sourced Bessel functions to cut down on load-time 
- [ ] Add sound speed, density, and absorption coefficient calculator
