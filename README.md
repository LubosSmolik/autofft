## autofft - A time-frequency analyser for Matlab
The Signal Processing Toolbox™ provides several functions to power spectrum estimation, including `pspectrum`, `pwelch` and `stft`. Although these functions are of high quality and well-documented, they might be cumbersome in some engineering applications. These applications include estimating magnitudes of vibration, noise and other discrete-time signals in engineering units or comparing theoretical results with measurements.

In such applications, you can use `autofft` to estimate the discrete Fourier transform (DFT), which mimics the properties of the Brüel & Kjaer FFT analysers. Based on your input, autofft segments signal, applies window functions and performs spectral averaging. The resulting averaged spectrum, also called _modified periodogram_, can be returned in various engineering spectral units, including decibels, magnitude, root mean square (RMS), peak-to-peak and power spectral density (PSD). autofft can also estimate spectral derivation or spectral integral of DFT and perform the short-time Fourier transform (STFT).

### Capabilities of autofft
- Does not require the Signal Processing Toolbox™.
- Performs significantly better than `pwelch` and `pspectrum` in computationally intensive problems.
- Estimates magnitudes of components in your data in engineering units, e.g. dB, V or Pa.
- Provides control over the setup of the frequency analyser, which is impossible with library functions.
- Can apply high-pass and frequency weighting filters, spectral derivation or spectral integration.

### What's new in v1.5?
- __v1.5.2:__ _Changed functionality_: The package no longer requires the Signal Processing Toolbox™.
- __v1.5.2:__ _Changed functionality_: A first-order Butterworth digital filter is now used for high-pass filtering rather than a first-order elliptic filter.
- __v1.5.2:__ _New functions_: The package is now distributed with functions that can construct Blackman-Harris, flat-top, Hamming, Hann, Kaiser and uniform windows and can design an n-th order Butterworth digital filter. These functions can be found in `+utilities` directory.
- __v1.5.1:__ _Code optimisation_: The STFT is now computed more efficiently.
- __v1.5.1:__ _Bug fix_: In same cases, times for the STFT were evaluated more than once. This has been fixed.
- _New functionality_: The output spectra can be returned in decibel scale using the `'dbReference'` parameter.
- _New function_: A `freqWeight` function, which applies frequency weighting filters to the power spectrum, is now included in the package.
- _Documentation_: New example added.
- _Code optimisation_: Times at which the STFT is evaluated are computed more efficiently.

### Getting started

For more information read [user manual](https://github.com/CarlistRieekan/autofft/blob/master/user_manual.pdf) or visit  [![View Frequency and time-frequency analysis in Matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/69534-frequency-and-time-frequency-analysis-in-matlab)
