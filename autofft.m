function [spectrum, freq, varargout] = autofft(xs, ts, userSetup)
% AUTOFFT Evaluates a frequency spectrum of a signal using wFFT algorithm
%
%  Copyright (c) 2017-2022         Lubos Smolik, University of West Bohemia
% v1.5.2 (build 7. 7. 2022)        e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
% s = autofft(xs, fs)
% s = autofft(xs, ts)
% s = autofft(___, setup)
% [s, f] = autofft(___)
% [s, f, setup] = autofft(___)
% [s, f, t, setup] = autofft(___)
%
% s = autofft(xs, fs) returns the DFT or STFT s of xs using sampling
%   frequency fs (Hz). xs can be either a vector or an array consisting of
%   column vectors.
%
% s = autofft(xs, ts) returns the DFT or STFT s of xs using a vector of
%   time stamps ts (s). Also returns the frequencies f at which the DFT or
%   STFT s is evaluated.
%
% [___] = autofft(___, setup) performs the DFT or STFT name-value pair
%   arguments specified in structure setup.
%
% [s, f] = autofft(___) returns the frequencies f at which the DFT or STFT
%    s is evaluated.
%
% [s, f, setup] = autofft(___) returns the setup of the FFT analyser.
%
% [s, f, t, setup] = autofft(___) returns the times t at which the STFT is
%     evaluated.
%
%       
% Construction of setup:
%	setup = struct('param', 'value', ...);
%
%	List of parameters (all parameters and values are case insensitive):
%   - 'FFTLength' - [ positive integer {length(xs) or size(xs, 2)} ]
%     The length of each segment in samples, i.e. number of DFT points.
%
%   - 'TimeResolution' - [ real positive scalar ]
%     The duration of each segment in seconds.
%
%   - 'FrequencyResolution' - [ real positive scalar ]
%     The desired frequency resolution of the analyser in Hz.
%
%     Note: Default values of the three parameters above are set so that
%       the input data are treated as one segment, i.e. no spectral averag-
%       ing is performed. If more than one of these three parameters is
%       specified, 'FFTLength', 'TimeResolution' and 'FrequencyResolution'
%       have the highest, middle and lowest priority, respectively.
%
%   - 'Mode' - [ {'onesided'} | 'twosided' ]
%     Selects an algorithm to construct of the freqeuncy spectrum:
%     - 'onesided' - computes one-sided spectrum estimates over [0, maxF]
%     - 'twosided' - computes centered, two-sided spectrum estimates over
%       [-maxF, maxF], where maxF is given by 'LowPassFrequency'
%
%   - 'HighPassFrequency' - [ {NaN} | real scalar ]
%     Specifies the passband frequency of an elliptic highpass filter.
%     The filter has a slope -20 dB/dec from the passband frequency.
%     The exact filtering process depends on the parameter value:
%     - NaN             - no highpass filtering
%     - 0               - DC filtering with the use of detrend function
%     - positive scalar - DC filtering and subsequent highpass filtering
%     - negative scalar - highpass filtering without DC filtering
%
%   - 'LowPassFrequency' - [ positive scalar in (0, fs/2) {fs/2} ]
%     Specifies the maximum frequency in Hz over which autofft computes
%     estimates. Recommended value is fs/2.56.
%
%   - 'Window' - [ character {'u'} | string | vector ]
%     Specifies a window function which is than multiplied with each
%     segment. The window function can be specified either as a vector
%     or using a character or a string from the list below:
%     - 'b' - Blackmann-Harris window
%     - 'f' - flat-top window
%     - 'h' - Hann window
%     - 'm' - Hamming window
%     - 'k' - Kaiser window with shape factor beta = 1.6
%     - 'kA.A' - Kaiser window with shape factor beta = A.A
%     - 'u' - uniform (rectangular) window                        {default}
%     Alternatively, a string containing a name of the desired window can
%     be used instead of the character.
%
%   - 'OverlapLength' - [ 50 % of 'FFTLength' | integer ] 
%      Number of overlapped samples of two successive segments. If both
%      'OverlapLength' and 'OverlapPercentage' are specified, only 
%      'OverlapLength' is considered.
%
%   - 'OverlapPercentage' - [ real scalar in [0, 100] {50} ] 
%      Overlap percentage of two successive segments. Use value 100 for the
%      maximum possible overlap, i.e. 'FFTLength' - 1.
%
%   - 'Averaging' - [ energy | {linear} | max | median | min | var | none]
%     Specifies the spectral averaging mode from: 
%     - 'energy', 'rms'     - energy (rms) averaging
%     - 'linear', 'lin'     - linear averaging                    {default}
%     - 'max', 'pk', 'peak' - maximum-hold aveargaing
%     - 'median'            - returns median specified spectral unit
%     - 'min'               - minimum-hold aveargaing
%     - 'var'               - returns variance of specified spectral unit
%     - 'none'              - returns the STFT with no spectral averaging
%
%   - 'jwWeigthing' - [ 1/jw2 | 1/jw | {1} | jw | jw2 ]
%     Use this parameter to apply a frequency-domain post-weighting to the
%     output spectra e.g., to estimate displacement from acceleration. 
%     The parameter can be specified as follows:
%     - '1/jw2' - double integration
%     - '1/jw'  - single integration
%     - '1'     - as input                                        {default}
%     - 'jw'    - single differentiation
%     - 'jw2'   - double differentiation
%
%   - 'SpectralUnit' - [ {pow} | rms | pk | pp | psd | rsd ]
%     Specifies absolute unit used to compute estimates from:
%     - 'pow', 'power' - autospectrum (square of rms magnitudes)  {default}
%     - 'rms'          - linear spectrum with rms magnitude
%     - 'pk', '0-pk'   - linear spectrum with 0-peak magnitude
%     - 'pp', 'pk-pk'  - linear spectrum with peak-peak magnitude
%     - 'psd'          - power spectral density
%     - 'rsd','rmssd'  - root mean square of power spectral density 
%
%   - 'dbReference' - [ {NaN} | 0 | real positive scalar ]
%     Specifies the reference value to calculate the decibel scale.
%     - NaN - The output spectrum is not expressed in dB.
%     - 0   - The output spectrum is expressed in dB. The reference value
%             is selected automatically so that 0 dB is the maximum.
%     - positive scalar - The output spectrum is expressed in dB with the
%             reference value specified by the user. 
%
% What's new in v1.5?
% v1.5.1: Code optimisation: The STFT is now computed more efficiently.
% v1.5.1: Bug fix: In same cases, times for the STFT were evaluated more
%   than once. This has been fixed.
% New functionality: The output spectra can be returned in decibel scale
%   using the 'dbReference' parameter.
% New function: A freqWeight function, which applies frequency weighting
%   filters to the power spectrum, is now included in the release.
% Documentation: New example added.
% Code optimisation: Times at which the STFT is evaluated are computed more
%   efficiently.

%% Validate number of input arguments
narginchk(2, 3);

%% Convert row vectors to column vectors if needed
if size(xs, 1) == 1     % Samples
    xs = xs(:);
end

if length(ts) == 1      % User-specified sampling frequency 
    fs = ts;
    ts = 0.5/fs;
else                    % Compute sampling frequency from time stamps
    fs = (length(ts) - 1) / (ts(end) - ts(1)); 
end

%% Generate a structure array containing a default analyser setup
setup = struct("SamplingFrequency",    fs, ...
               "DataDuration",         size(xs, 1) / fs, ...
               "DataLength",           size(xs, 1), ...
               "FFTLength",            NaN, ...
               "TimeResolution",       NaN, ...
               "FrequencyResolution",  NaN, ...
               "Mode",                 "onesided", ...
               "HighPassFrequency",    NaN, ...
               "LowPassFrequency",     fs / 2, ...
               "Window",               "uniform", ...
               "WindowNoiseBandwidth", 1, ...
               "OverlapLength",        NaN, ...
               "OverlapPercentage",    50, ...
               "Averaging",            "linear", ...
               "NumberOfAverages",     NaN, ...
               "jwWeigthing",          "none", ...
               "SpectralUnit",         "power", ...
               "dbReference",          NaN);
setupFields = fieldnames(setup);

% Set user-specified parameters
if nargin == 3
    % Convert fftset to setup (due to compatibility with v1.1 and v1.2)
    userFields = fieldnames(userSetup);
    oldFields = ["nwin" "twin" "df" "highpass" "lowpass" "overlap" "jw" "unit"];
    newFields = [4 5 6 8 9 13 16 17];

    for i = 1:size(userFields, 1)
        ind = strcmpi(userFields{i}, oldFields);
        if any(ind)
            userSetup.(setupFields{newFields(ind)}) = userSetup.(userFields{i});
        end
    end

    % Refresh user-specified field names (these might changed above)
    userFields = fieldnames(userSetup);

    % Merge the user-specified setup with the default analyser setup
    for i = 4:size(setupFields, 1)
        ind = strcmpi(setupFields{i}, userFields);
        if any(ind)
            setup.(setupFields{i}) = userSetup.(userFields{ind});
        end
    end
end

%% Generate remaining parameters in the analyser setup
% Check if there is user-defined FFT length
if ~isnan(setup.FFTLength)
    if setup.FFTLength > setup.DataLength
        setup.FFTLength = setup.DataLength;
        warning("Specified FFT length is too high. " + ...
                "FFT length has been changed to " + ...
                num2str(setup.FFTLength, '%d') + " samples.");
    end
% Check if there is user-defined time resolution
elseif ~isnan(setup.TimeResolution)
    if round(setup.TimeResolution * fs) > setup.DataLength
        setup.FFTLength = setup.DataLength;
        warning("Specified time resolution is too long. " + ...
                "Time resolution has been changed to " + ...
                num2str(setup.FFTLength / fs, '%.2g') + " s.");
    else
        setup.FFTLength = round(setup.TimeResolution * fs);
    end
% Check if there is user-defined frequency resolution
elseif ~isnan(setup.FrequencyResolution)
    if round(fs / setup.FrequencyResolution) > setup.DataLength
        setup.FFTLength = setup.DataLength;
        warning("Specified frequency resolution cannot be reached. " + ...
                "Frequency resolution has been changed to " + ...
                num2str(fs / setup.FFTLength, '%.2g') + " Hz.");
    else
        setup.FFTLength = round(fs / setup.FrequencyResolution);
    end
% Do if none of FFT length, time res. and frequency res. is defined 
else
    setup.FFTLength = setup.DataLength;
end

% Compute exact time resolution and frequency resolution
setup.FrequencyResolution = fs / setup.FFTLength;
setup.TimeResolution      = setup.FFTLength / fs;

% Check if the low pass frequency is not higher than the Nyquist frequency
if setup.LowPassFrequency > fs
    setup.LowPassFrequency = fs / 2;
    warning("Maximum frequency has been changed to " + fs + " Hz.");
end

% Generate frequency vector
switch lower(setup.Mode)
    case "twosided" % Two-sided spectrum
        freq = [-flip(0:setup.FrequencyResolution:setup.LowPassFrequency) ...
                setup.FrequencyResolution:setup.FrequencyResolution:setup.LowPassFrequency].';
        setup.Mode = "twosided";

        % Vector of indices for manipulation with the DFT
        sc  = 1;                            % Scaling constant
        st  = (size(freq, 1) - 1) / 2 + 1;  % Index of the static component
        ind = [setup.FFTLength - st + 2:setup.FFTLength, 1:st].';
    otherwise       % One-sided spectrum 
        freq = (0:setup.FrequencyResolution:setup.LowPassFrequency).';
        setup.Mode = "onesided";

        % Vector of indices for manipulation with the DFT
        sc  = 2;                            % Scaling constant
        st  = 1;                            % Index of the static component
        ind = (1:size(freq, 1)).';
end

% Calculate number of overlaping samples
if isnan(setup.OverlapLength)
	setup.OverlapLength = round(setup.FFTLength * setup.OverlapPercentage / 100);
end

% Check whether the overlap length is permissible
if setup.OverlapLength >= setup.FFTLength
    setup.OverlapLength = setup.FFTLength - 1;
end
    
% Recalculate overlap percentage due to rounding
setup.OverlapPercentage = 100 * setup.OverlapLength / setup.FFTLength;

% Calculate the number of segments and segment indices
setup.NumberOfAverages = floor((setup.DataLength - setup.OverlapLength) / ...
                               (setup.FFTLength - setup.OverlapLength));
seg = zeros(setup.NumberOfAverages, 2);         % matrix of segment indices
ni  = 1;                                        % pointer
for i = 1:setup.NumberOfAverages                % cycle through segments
    seg(i, 1) = ni; 
    ni = ni + setup.FFTLength - 1;
    seg(i, 2) = ni;
    ni = ni - setup.OverlapLength + 1;
end

% Set a window function
if isnumeric(setup.Window)
    % Convert a window function to a column vector
    setup.Window = setup.Window(:);
    
    % Interpolate missing values in the window function if necessary
    if  size(setup.Window, 1) ~= setup.FFTLength
        % Use the modified Akima interpolation (Matlab v2019b or newer)
        try   
            setup.Window = makima(linspace(0, 1, size(setup.Window, 1)), ...
                             setup.Window, linspace(0,1,setup.FFTLength).');
        % Use the spline interpolation (Matlab v2019a or older)
        catch 
            setup.Window = spline(linspace(0, 1, size(setup.Window, 1)), ...
                             setup.Window, linspace(0,1,setup.FFTLength).');
        end
        warning("The window function has been interpolated to " + ...
                num2str(setup.FFTLength, "%d") + " samples.");
    end
else
    % Generate the window function internally
    switch lower(extractBefore(setup.Window, 2))
        case "b"    % Blackmann-Harris
            setup.Window = blackmanharris(setup.FFTLength);
            windowName   = "Blackmann-Harris";
        case "f"    % flat-top
            setup.Window = flattopwin(setup.FFTLength);
            windowName   = "flat-top";
        case "h"    % Hann or Hamming
            if strncmpi(setup.Window, "ham", 3) % Hamming
                setup.Window = hamming(setup.FFTLength);
                windowName   = "Hamming";
            else    % Hann
                setup.Window = hann(setup.FFTLength);
                windowName   = "Hann";
            end
        case "k"    % Kaiser
            % Try to calculate parameter beta
            beta = str2double(regexprep(setup.Window,"[a-zA-Z,=\s]",""));

            if isnan(beta) % Use the default Kaiser window (beta = 1.6)
                setup.Window = kaiser(setup.FFTLength, 1.6);
                windowName   = "Kaiser, b = 0.5";
            else
                setup.Window = kaiser(setup.FFTLength, beta);
                windowName   = "Kaiser, b = " + string(beta);
            end
        case "m"    % Hamming
            setup.Window = hamming(setup.FFTLength);
            windowName   = "Hamming";
        otherwise   % uniform
            setup.Window = rectwin(setup.FFTLength);
            windowName   = "uniform";
    end
end

% Normalise the window function for correct magnitude estimation
setup.Window = setup.Window ./ mean(setup.Window);

% Calculate the noise power bandwidth of the window function in Hz
% exploiting the fact that mean(setup.Window) == 1
setup.WindowNoiseBandwidth = fs * mean(setup.Window.^2) / setup.FFTLength;

%% Filtering
if ~isnan(setup.HighPassFrequency)
    % Remove the mean value from xs
    if setup.HighPassFrequency >= 0
        xs = detrend(xs, 0);
    end
    
    % Construct and apply an elliptic high-pass filter
    if setup.HighPassFrequency ~= 0
        n     = 1;    % Filter order
        aPass = 0.1;  % Passband ripple (dB)
        aStop = 20;   % Stopband attenuation (dB)

        % Construct an FDESIGN object (5x faster than designfilt)
        h  = fdesign.highpass("N,Fp,Ast,Ap", n, abs(setup.HighPassFrequency), ...
                                             aStop, aPass, fs);
        hp = design(h, "ellip");
        
        % Filter xs
        xs = filter(hp, xs);
    end
end

%% FFT
% Preallocate an array for the DFT of individual segments
tSegments = zeros(setup.FFTLength, setup.NumberOfAverages, size(xs, 2));

% Fast Fourier transform of the time-weighted segments
for i = 1:size(xs, 2)
    for j = 1:setup.NumberOfAverages
        % FFT of the individual segment
        tSegments(:, j, i) = xs(seg(j,1):seg(j,2), i);
    end
end

% Apply time weighting
tSegments = setup.Window .* tSegments;

% Perform FFT
tSpectrum = fft(tSegments, [], 1);

% Application of the jw weigthing and computation of unscaled magnitudes
switch lower(setup.jwWeigthing)
    case {"1/jw2", "double integration"}
        tSpectrum(ind,:,:) = abs((-1 ./ (4 * pi^2 * freq.^2)) .* ...
                                 tSpectrum(ind,:,:));
        setup.jwWeigthing  = "double integration";
    case {"1/jw", "single integration"}
        tSpectrum(ind,:,:) = abs((1 ./ (2i * pi * freq)) .* ...
                                 tSpectrum(ind,:,:));
        setup.jwWeigthing  = "single integration";
    case {"jw", "single differentiation"}
        tSpectrum(ind,:,:) = abs(2i * pi * freq .* tSpectrum(ind,:,:));
        setup.jwWeigthing  = "single differentiation";
    case {"jw2", "double differentiation"}
        tSpectrum(ind,:,:) = abs(-4 * pi^2 * freq.^2 .* tSpectrum(ind,:,:));
        setup.jwWeigthing  = "double differentiation"; 
    otherwise
        tSpectrum(ind,:,:) = abs(tSpectrum(ind,:,:));
end

% Spectral averaging and the limitation of the maximum frequency
switch lower(setup.Averaging)
    case {"energy", "rms"}                  % Energy averaging
        spectrum = rms(tSpectrum(ind,:,:), 2);
        setup.Averaging = "energy";       
    case {"maximum", "max", "pk", "peak"}   % Maximum-hold averaging
        spectrum = max(tSpectrum(ind,:,:), [], 2);
        setup.Averaging = "maximum";
    case {"median", "med"}                  % Median-hold averaging
        spectrum = median(tSpectrum(ind,:,:), 2);
        setup.Averaging = "median";
    case {"minimum", "min"}                 % Minimum-hold averaging
        spectrum = min(tSpectrum(ind,:,:), [], 2);
        setup.Averaging = "minimum";
    case "none"                             % No averaging
        spectrum = tSpectrum(ind,:,:);
    case {"variance", "var"}                % Variance-hold averaging
        spectrum = var(tSpectrum(ind,:,:), 0, 2);
        setup.Averaging = "variance";
    otherwise                               % Linear averaging
        spectrum = mean(tSpectrum(ind,:,:), 2);
        setup.Averaging = "linear";
end

% Evaluation of spectral unit
switch lower(setup.SpectralUnit)
    case "rms"                          % RMS magnitude
        spectrum(st,:)        = spectrum(st,:) / setup.FFTLength;
        spectrum(1:end~=st,:) = (sc/sqrt(2)) * spectrum(1:end~=st,:) / setup.FFTLength;
        setup.SpectralUnit    = "RMS";
    case {"pk", "0-pk", "peak"}         % 0-peak magnitude
        spectrum(st,:)      = spectrum(st,:) / setup.FFTLength;
        spectrum(1:end~=st,:)  = sc * spectrum(1:end~=st,:) / setup.FFTLength;
        setup.SpectralUnit = "0-pk";
    case {"pp", "pk-pk", "peak2peak"}   % Peak-peak magnitude
        spectrum(st,:)      = spectrum(st,:) / setup.FFTLength;
        spectrum(1:end~=st,:) = 2 * sc * spectrum(1:end~=st,:) / setup.FFTLength;
        setup.SpectralUnit = "pk-pk";
    case {"asd", "psd"}                 % Power spectral density
        spectrum(st,:)      = (spectrum(st,:) / setup.FFTLength).^2 / ...
                             setup.WindowNoiseBandwidth;
        spectrum(1:end~=st,:)  = sc * (spectrum(1:end~=st,:) / setup.FFTLength).^2 / ...
                             setup.WindowNoiseBandwidth;
        setup.SpectralUnit = "PSD";
    case {"rsd", "rmssd"}               % Root mean square spectral density
        spectrum(st,:)        = spectrum(st,:) / ...
                            (setup.WindowNoiseBandwidth * setup.FFTLength);
        spectrum(1:end~=st,:) = (sc/sqrt(2)) * spectrum(1:end~=st,:) / ...
                            (setup.WindowNoiseBandwidth * setup.FFTLength);
        setup.SpectralUnit = "RMSSD";
    otherwise                           % Autospectrum / power spectrum
        spectrum(st,:)        = (spectrum(st,:) / setup.FFTLength).^2;
        spectrum(1:end~=st,:) = sc * (spectrum(1:end~=st,:) / setup.FFTLength).^2;
        setup.SpectralUnit = "power";
end

% Squeeze the resulting spectrum
spectrum = squeeze(spectrum);

% Transform resulting spectra to decibels
if ~isnan(setup.dbReference)
    % Set the reference level automatically in not specified by the user
    if setup.dbReference == 0
        if setup.SpectralUnit == "PSD" || setup.SpectralUnit == "power"
            setup.dbReference = sqrt(max(spectrum(:)));
        else
            setup.dbReference = max(spectrum(:));
        end
    end

    % Calculate decibels
    if setup.SpectralUnit == "PSD" || setup.SpectralUnit == "power"
        % Calculate decibels for powerspectra, i.e. 20 log (s_i / s_ref)
        spectrum = 10 * log10(spectrum / (setup.dbReference^2));
    else
        % Calculate decibels for linear spectra, i.e. 10 log (s_i / s_ref)
        spectrum = 10 * log10(spectrum / setup.dbReference);
    end  
end

%% Return the analyser setup
% Return the analyser setup only
if nargout == 3
    % Convert numeric representation of the window to string
    if exist("windowName", "var")
        setup.Window = windowName;
    end

    varargout{1} = setup;

% Return segment times and the analyser setup
elseif nargout > 3
    % Compute the relative segment times
    tshift = (setup.FFTLength - setup.OverlapLength) / fs;
    tseg   = (0:tshift:(setup.NumberOfAverages - 1) * tshift).';
    
    % Compute the absolute segment times        
    varargout{1} = mean(setup.Window .* (0:setup.FFTLength-1).') / fs + ts(1) + tseg;

    % Convert numeric representation of the window to string
    if exist("windowName", "var")
        setup.Window = windowName;
    end

    varargout{2} = setup;
end

% End of main fucntion
end