## autofft - A frequency analyser for Matlab
Matlab Signal Processing Toolbox™ provides several functions to power spectrum estimation, including `pspectrum`, `pwelch` and `stft`. Although these functions are of high quality and well-documented, they might be cumbersome in some engineering applications. These applications include estimating magnitudes of vibration, noise and other discrete-time signals in engineering units or comparing theoretical results with measurements.

In such applications, you can use `autofft` to estimate the discrete Fourier transform (DFT), which mimics the properties of the Brüel & Kjaer FFT analysers. Based on your input, `autofft` segments signal, applies window functions and performs spectral averaging. The resulting DFT can be returned in various engineering spectral units, including decibels, magnitude, root mean square (RMS), peak-to-peak and power spectral density (PSD). autofft can also estimate spectral derivation or spectral integral of DFT and perform the short-time Fourier transform (STFT).

### Capabilities of autofft
- It performs significantly better than `pwelch` in computationally intensive problems.
- Estimates magnitudes of components in your data in engineering units, e.g. dB, V or Pa.
- Provides control over the setup of the frequency analyser, which is impossible with library functions.
- Can apply high-pass and frequency weighting filters, spectral derivation or spectral integration.

For more information read [user manual](https://github.com/CarlistRieekan/autofft/blob/master/user_manual.pdf) or visit  [![View Frequency and time-frequency analysis in Matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/69534-frequency-and-time-frequency-analysis-in-matlab)
