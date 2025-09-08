## autofft - A time-frequency analyser for Matlab

`autofft` is a package of MATLAB® functions for automated frequency analysis of discrete-time signals, inspired by the operation of Brüel & Kjær FFT analysers.
It automatically filters, segments, weights and averages input signals, returning spectra in a wide range of engineering units (dB, RMS, peak-to-peak, PSD, etc.).
It works without the Signal Processing Toolbox™ and is optimised for performance.

### Introduction
Signal Processing Toolbox™ provides several functions for power spectrum estimation, including `pspectrum`, `pwelch` and `stft`. While these functions are high-quality and well-documented, they can be cumbersome for certain engineering applications. These applications include estimating magnitudes of vibration, noise and other discrete-time signals in engineering units or comparing theoretical results with measurements.

The `autofft` package is designed to address these needs by mimicking the operation of Brüel & Kjær FFT analysers. Based on user input, `autofft` automatically filters the signal, divides it into segments, applies window functions, and performs spectral averaging. The resulting averaged spectrum, also called modified periodogram, can be returned in various engineering spectral units, including decibels, root mean square (RMS) and peak-to-peak magnitudes, and power spectral density (PSD). In addition, `autofft` can compute spectral derivatives and integrals, as well as perform short-time Fourier transform (STFT) analysis.

### Key features
- __No toolbox dependency:__ works without Signal Processing Toolbox™.
- __High performance:__ performs significantly better than `pwelch` and `pspectrum` in computationally intensive problems.
- __Works with engineering units:__ works not only with relative units (dB) but also with engineering units (EU), such as RMS, peak-to-peak or power spectral density.
- __Configurable analyser setup:__ provides greater control over the setup of the frequency analyser than built-in functions.
- __Built-in filtering:__ can apply high-pass filters in the time domain.
- __Advanced post-processing:__ can apply frequency weighting filters per IEC 61672-1:2013, spectral derivation or spectral integration.

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](https://github.com/LubosSmolik/autofft/blob/main/%2Bvalidation/autofft_cpu_time.png)
__Figure:__ CPU times required to run `autoFFT` and similar Matlab functions measured on a Code Ocean using `timeit` function, see this [reproducible capsule](https://codeocean.com/capsule/9899368/tree). Each run includes a predefined number of segments weighted using a 4096-point Hann window.

### What's new?
- __v1.5.5:__ _Bug fix_: Uncommented text preventing use of custom filters has been removed.
- __v1.5.5:__ _Bug fix:_ Erroneous coefficient that caused bandstop filters to be generated incorrectly using autoButter function has been corrected.
- __v1.5.5:__ _Code optimisation:_ Generation of windows with a very large number of samples has been optimised.
- __v1.5.4:__ _Bug fix_: Automatic plotter no longer uses white backgroud color in R2025a dark mode.
- __v1.5.4:__ _Bug fix_: Plotting error when the user selected tiled layout for time-frequency analysis results from only one channel has been fixed.
- __v1.5.4:__ _Code optimisation_: Code optimisation: Evaluation of magnitude and spectral unit optimised. Other minor optimisations implemented and code refactored. Implemented changes reduced CPU time by up to 35 %.
- __v1.5.4:__ _Code optimisation_: Error handling during filtering has been improved.
- __v1.5.4:__ _Documentation_: Introduction has been revised.

### Getting started

1. Download the latest release from the [Releases page](https://github.com/LubosSmolik/autofft/releases), extract it and add the folder to your MATLAB path:
   `addpath('path_to_autofft_folder');
   savepath;`

3. [![View Frequency and time-frequency analysis in Matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/69534-frequency-and-time-frequency-analysis-in-matlab) By clicking this link, you can install the package directly from MATLAB File Exchange.

For further information read [user manual](https://github.com/LubosSmolik/autofft/blob/master/user_manual.pdf).
