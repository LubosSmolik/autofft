% Parameters of test signal
fs   = 4096;                        % Sampling frequency (Hz)
dt   = 1/fs;                        % Sampling period (s)
nseg = round(logspace(0, 4, 16));   % Number of segments in a signal (-)

stautofft     = zeros(1, length(nseg));
stautofft2    = zeros(1, length(nseg));
stautofft3    = zeros(1, length(nseg));

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

setupSTFT = struct("FFTLength",     nfft, ...
                   "Window",        win, ...
                   "OverlapLength", olal, ...
                   "Averaging",     "none");


for i = 1:length(nseg)
    i

    % Generate test signal
    t = dt:dt:(nseg(i) * (nfft - olal) + olal) / fs;
    s = sin(2*pi*50*t) + 0.5 * sin(2*pi*150*t);
    x = s + 2*randn(size(t));
    
    % STFT
    f = @() autofft(x, fs, setupSTFT);
    stautofft(i) = timeit(f, 4);

    f = @() autofft2(x, fs, setupSTFT);
    stautofft2(i) = timeit(f, 4);

    f = @() autofft3(x, fs, setupSTFT);
    stautofft3(i) = timeit(f, 4);
end

% Clear t, s, x to limit the size of the results
clear t s x

% Plot results
figure;
axes(Box = "on", NextPlot = "add", XScale = "log", YScale = "log");
plot(nseg, stautofft, DisplayName = "autofft", LineWidth = 1);
plot(nseg, stautofft2, DisplayName = "autofft2", LineWidth = 1);
plot(nseg, stautofft3, DisplayName = "autofft2", LineWidth = 1);
legend;