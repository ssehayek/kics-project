# Measuring diffusion and blinking simultaneously in complex systems using autocorrelations

## Motivation

In this project we successfully measured diffusion and photoblinking fluorescent probe dynamics using an [extended k-space image correlation spectroscopy (kICS) method](https://doi.org/10.1016/j.bpr.2021.100015) in live cells. Our model assumes one immobile and one freely diffusing population of particles, both with the same photophysical properties. 

Application of image correlation methods generally requires homogeneity within the region of interest (ROI) being analyzed. On the other hand, our technique is able to analyze ROIs with non-uniform distributions of immobile particles that are suubject to photophysical processes. 

To summarize the main advantages: 

  * Rapid way to determine apparent diffusion coefficient, photoblinking rates, and fraction of diffusing particles within an ROI
  * Fit is solely dependent on these parameters
  * Able to analyze ROIs with high particle densities
  * Able to analyze ROIs with an inhomogeneous distribution of particles
  
## Features

  * [Demo code and simulations](kics-project/demo)
  * [Simulator](kics-project/simulation-files)
  * [kICS autocorrelation code](kics-tools)
  * [Fitting functions](kics-project/kics-fitting)
  
## Usage

Please see the [demo](kics-project/demo) folder for an example of how to use the code in this repository. Please also extract "simulation-files.zip" to run the [demo code](kics-project/demo).

## References

Please refer to our [published manuscript](https://doi.org/10.1016/j.bpr.2021.100015) for further information on this method.
