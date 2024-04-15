# SANDIX
Matlab code to fit models of diffusion with exchange in brain gray matter
As described in:
Jonas L. Olesen, Leif Østergaard, Noam Shemesh, Sune N. Jespersen,
Diffusion time dependence, power-law scaling, and exchange in gray matter,
NeuroImage, Volume 251, 2022, 118976
https://doi.org/10.1016/j.neuroimage.2022.118976 
Usage is free, but please cite this publication.

The main function is fitwrapper.m, which allows to fit a handfuld of model variants, SMEX, SMEX-dot, SANDIX and eSANDIX. A good place to start is the demo in fittingexample.m, which showcases the usage with a small example data set. 

The function plotlandscape can be used to visualize the fitting landscape (cost function) for a given signal. See landscape_example to learn how to use it.
