# MRS_Apodization
This repository contains 2 files. First 'Apodization_in_TD.m' and second 'lsqnonlin_fit.m'.

The first file takes input from .mat files (raw signal and processing parameter values) and performs apodization, plots them and also compares different apodization.

Zero and First Order Phase correction is also done and values for the parameters are fetched from the .mat files as mentioned earlier. In our case, .mat file was generated from the MRS processing software (KALPANA).

Baseline correction is also performed. In this case, we used Asymmetric Least Squares algorithm for the baseline correction. This could be replaced by any algorithm you deem fit.

A comparative study has NOT been done by us comparing the effects and accuracy in calculation of metabolites when different baseline correction algorithms are used.



The second file is a semi-automatic curve fitting technique. It takes a few inputs from the user like Total number of peaks to fit as "Number of metabolites"
Initial guess of the Amplitude and the ppm position of the peaks as "Estimated peaks" and initial guess of FWHM as "Estimated FWHMs"

A gaussian curve is fit to the signal with certain deviation (iterates between a lower and an upper bound to find the best fit) allowed from the estimated guess. The Amplitude has been given a larger leeway because the convolution affects the maximum amplitude of different metabolites. (Allows multipeak fitting)

Individual gaussians are fit to each signal and then a sum is calculated. From the summed gaussian and the processed Time Domain signal, a residual signal is calculated as well.

Using this, CRLB and Fiterror is calculated.


Should you notice any mistakes or want me to add something else, please feel free to write to me at zoheb95@gmail.com
