# SANDIX
Matlab code to fit models of diffusion with exchange in brain gray matter
As described in:
Jonas L. Olesen, Leif Østergaard, Noam Shemesh, Sune N. Jespersen,
Diffusion time dependence, power-law scaling, and exchange in gray matter,
NeuroImage, Volume 251, 2022, 118976
https://doi.org/10.1016/j.neuroimage.2022.118976 
Usage is free, but please cite this publication.

The main function is fitwrapper.m, which allows to fit a handfuld of model variants, SMEX, SMEX-dot, SANDIX and eSANDIX. Notably, the scripts take into account a finite diffusion gradient pulse width. A good place to start is the demo in fittingexample.m, which showcases the usage with a small example data set. 

The function plotlandscape can be used to visualize the fitting landscape (cost function) for a given signal. See landscape_example to learn how to use it.

![1-s2 0-S1053811922001057-gr1](https://github.com/Neurophysics-CFIN/SANDIX/assets/20967007/a2b32d91-7d00-4d6e-80a5-3c4dd82c1caa)

(From https://doi.org/10.1016/j.neuroimage.2022.118976 under Creative Commons License, https://creativecommons.org/licenses/by/4.0/)

![SMEX_SDE_NN](https://github.com/user-attachments/assets/3e56cc7f-647b-4af0-8ddb-b8a397541419)

The code can use a Neural Net for speed - the accuracy is shown here compared to the ordinary differential equation solver for 100,000 randomly generated kernels ([1 3 3 1]'.*rand(4,1) + [0 0 0 0]') and random protocols in the range of b = 0...100 ms/um2, Delta = 10...50 ms, delta = 1..10.
