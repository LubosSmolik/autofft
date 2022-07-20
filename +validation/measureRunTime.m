% AutoFFT - Performance measurement
% Authors: Luboš Smolík, Jan Rendl and Roman Pašek
% 2022, doi: 10.24433/CO.3004443.v1
% 
% This code measures CPU times of autofft and Matlab built-in functions. It is also
% published as reproducible capsule at https://codeocean.com/capsule/1932803/tree/v1
% 
% Requires the Signal Processing Toolbox
%
% Parameters of test signal
fs   = 4096;                        % Sampling frequency (Hz)
dt   = 1/fs;                        % Sampling period (s)
nseg = round(logspace(0, 5, 21));   % Number of segments in a signal (-)

% Arrays to store results
tautofft   = zeros(1, length(nseg));
tpspectrum = zeros(1, length(nseg));
tpwelch    = zeros(1, length(nseg));

stautofft     = zeros(1, length(nseg));
stpspectrum   = zeros(1, length(nseg));
ststft        = zeros(1, length(nseg));
stspectrogram = zeros(1, length(nseg));

% DFT setup
nfft = fs;          % Number of DFT points
df   = 1.46832;     % Target frequency resolution for pspectrum
win  = hann(nfft);  % Window
leak = 0.85;        % Leakage for pspectrum - approximates the Hann window
olap = 75;          % Overlap percentage
olal = olap * fs / 100; % Overlap length

% Note: The frequency resolution of pspectrum is calculated using the sampling
%       frequency, number of samples in a segment and leakage of the used window.
%       Here, we use the frequency resolution which yields the number of segments
%       equal to other methods.

% Specify autofft setup
setupPSp = struct("FFTLength",     nfft, ...
                  "Window",        win, ...
                  "OverlapLength", olal );

setupSTFT = struct("FFTLength",     nfft, ...
                   "Window",        win, ...
                   "OverlapLength", olal, ...
                   "Averaging",     "none");


for i = 1:length(nseg)
    % Generate test signal
    t = dt:dt:(nseg(i) * (nfft - olal) + olal) / fs;
    s = sin(2*pi*50*t) + 0.5 * sin(2*pi*150*t);
    x = s + 2*randn(size(t));
    
    % Averaged DFT
    f = @() autofft(x, fs, setupPSp);
    tautofft(i) = timeit(f, 2);

    f = @() pspectrum(x, fs, FrequencyResolution = df, Leakage = leak);
    tpspectrum(i) = timeit(f, 2);

    f = @() pwelch(x, win, olal, fs);
    tpwelch(i) = timeit(f, 2);

    % STFT
    f = @() autofft(x, fs, setupSTFT);
    stautofft(i) = timeit(f, 4);

    f = @() pspectrum(x, fs, "spectrogram", FrequencyResolution = df, ...
                      Leakage = leak, OverlapPercent = olap);
    stpspectrum(i) = timeit(f, 3);

    f = @() validation.timeit_stft(x, fs, win, olal);
    ststft(i) = timeit(f, 3);

    f = @() spectrogram(x, win, olal, fs);
    stspectrogram(i) = timeit(f, 4);
end

% Clear t, s, x to limit the size of the results
clear t s x

% Save results
save("../results/results.mat")

% Plot results
figure;
axes(Box = "on", NextPlot = "add", XScale = "log", YScale = "log");
plot(nseg, tautofft, DisplayName = "autofft", LineWidth = 1);
plot(nseg, tpspectrum, DisplayName = "psepctrum", LineWidth = 1);
plot(nseg, tpwelch, DisplayName = "pwelch", LineWidth = 1);
legend;
saveas(gcf, "../results/powerspectrum.fig");

figure;
axes(Box = "on", NextPlot = "add", XScale = "log", YScale = "log");
plot(nseg, stautofft, DisplayName = "autofft", LineWidth = 1);
plot(nseg, stpspectrum, DisplayName = "psepctrum", LineWidth = 1);
plot(nseg, ststft, DisplayName = "stft", LineWidth = 1);
plot(nseg, stspectrogram, DisplayName = "spectrogram", LineWidth = 1);
legend;
saveas(gcf, "../results/stft.fig");
