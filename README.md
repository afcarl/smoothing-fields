# smoothing-fields
Implementations from "Smoothing Fields of Weighted Collections with Applications to Diffusion MRI Processing" by Sigurdsson and Prince

This folder containts scripts used to create some of the experiments in 'Smoothing Fields of Weighted Collections with Applications to Diffusion MRI Processing' which appeared in SPIE Medical Imaging 2014

coeffs.m is a function that calculates coefficients for the data
smooth1D.m uses said coefficients to do smoothing
collectionsmooth*.m uses these functions to smooth a slice/volume

idealdemo.m is a script that creates a toy phantom and does smoothing
fibertracking.m is a script that does simple fiber tracking before and after smoothing
comparisonwithgaussian.m is a script used to compare with smoothing in the measurement domain. This script is missing steps where the reconstruction algorithm CFARI was used (implemented in the Java Imaging Science Toolkit)


Gunnar Atli Sigurdsson, 2013
