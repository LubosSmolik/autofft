function [spectrum, freq, varargout] = autofft(xs, ts, userSetup)
% AUTOFFT Evaluates a frequency spectrum of a signal using wFFT algorithm
%
% Copyright (c) 2017-2025              Luboš Smolík, Jan Rendl, Roman Pašek
% v1.5.4 (build 7. 8. 2025)            e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
%
% autofft(xs, fs)
% autofft(xs, ts)
% autofft(___, setup)
% s = autofft(___)
% [s, f] = autofft(___)
% [s, f, setup] = autofft(___)
% [s, f, t] = autofft(___)
% [s, f, t, setup] = autofft(___)
%
%
% autofft(xs, fs) computes and diplays the autospectrum or spectrogram of
%   input xs using sampling frequency fs in Hz. xs can be either a vector
%   or an array consisting of column vectors. In the latter case, each
%   spectrum or spectrogram is shown in a separate figure.
%
% autofft(xs, ts) computes and diplays the autospectrum or spectrogram of
%   input xs using a vector of time stamps ts.
% 
% autofft(___, setup) computes and diplays the autospectrum or spectrogram
%   of input xs using name-value pair arguments specified in a structure
%   array setup.
% 
% s = autofft(___) computes and returns autospectrum or the short-time
%   Fourier transform (STFT) of input xs. Does not visualise the results
%   unless specified in setup.
%
% [s, f] = autofft(___) also returns the vector of frequencies f at which
%   the autospectrum or spectrogram s are evaluated.
%
% [s, f, setup] = autofft(___) also returns the setup of the FFT analyser
%   if s is the autospectrum.
%
% [s, f, t] = autofft(___) also returns the vector of times t at which the
%   STFT slices are evaluated if the STFT is performed.
%
% [s, f, t, setup] = autofft(___) returns all information available.
%
%       
% Construction of setup:
%	setup = struct('param', 'value', ...);
%
% List of valid name-value pair arguments (all parameter names and values
% are case insensitive):
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
%   - 'HighPassFrequency' - [ {NaN} | real scalar | cell of vectors [b] and
%     [a] | filter object ]
%     Specifies the passband frequency of an elliptic highpass filter.
%     The filter has a slope -20 dB/dec from the passband frequency.
%     The exact filtering process depends on the parameter value:
%     - NaN             - no highpass filtering                   {default}
%     - 0               - DC filtering with the use of detrend function
%     - positive scalar - DC filtering and subsequent highpass filtering
%     - negative scalar - highpass filtering without DC filtering
%     - {[b],[a]}       - [b] and [a] are vectors specifiying the numerator
%                         and denominator coefficients of a rational 
%                         transfer function used to filter xs. Internally,
%                         filter(b,a,xs) is used to filter the input data.                  
%     - filter object   - Object, used to filter the input data. Since the
%                         filter validation depends on the DSP toolbox, the
%                         the object is not validated, which may cause an
%                         unexpected behaviour.
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
%   - 'jwWeigthing' - [ '1/jw2' | '1/jw' | {'1'} | 'jw' | 'jw2' ]
%     Use this parameter to apply a frequency-domain post-weighting to the
%     output spectra e.g., to estimate displacement from acceleration. 
%     The parameter can be specified as follows:
%     - '1/jw2' - double integration
%     - '1/jw'  - single integration
%     - '1'     - as input                                        {default}
%     - 'jw'    - single differentiation
%     - 'jw2'   - double differentiation
%
%   - 'SpectralUnit' - [ {'pow'} | 'rms' | 'pk' | 'pp' | 'psd' | 'rsd' ]
%     Specifies absolute unit used to compute estimates from:
%     - 'pow', 'power' - autospectrum (square of rms magnitudes)  {default}
%     - 'rms'          - linear spectrum with rms magnitude
%     - 'pk', '0-pk'   - linear spectrum with 0-peak magnitude
%     - 'pp', 'pk-pk'  - linear spectrum with peak-peak magnitude
%     - 'psd'          - power spectral density
%     - 'rsd','rmssd'  - root mean square of power spectral density 
%
%   - 'EngineeringUnit' - [ character {'EU'} | string ]
%     Specifies a unit of measure of xs. This parameter is used to generate
%     labels for data visualisation.
%
%   - 'dbReference' - [ {NaN} | 0 | real positive scalar ]
%     Specifies the reference value to calculate the decibel scale.
%     - NaN - The output spectrum is not expressed in dB.
%     - 0   - The output spectrum is expressed in dB. The reference value
%             is selected automatically so that 0 dB is the maximum.
%     - positive scalar - The output spectrum is expressed in dB with the
%             reference value specified by the user. 
%
%   - 'PlotLayout' - [ 'none' | 'separated' | 'stacked' | 'tiled' ]
%      Specifies, which layout is used to visualise results.
%      - 'none'      - Does not visualise any results. It is a {default}
%                      option if auttofft is called with no output.
%      - 'separated' - Visualise each spectrum or spectrogram in a separate
%                      figure. It is a {default} option if auttofft is
%                      called with at least one output variable.
%      - 'stacked'   - Plots spectra into single axes. This value cannot be
%                      used to visualise spectrograms.
%      - 'tiled'     - Creates a panel for each spectrum or spectrogram in
%                      a single figure.
%
% What's new in v1.5.4?
% v1.5.4: Bug fix: Plotting error when the user selected tiled layout for
%    time-frequency analysis results from only one channel has been fixed.
% v1.5.4: Code optimisation: Evaluation of spectral unit optimised. Minor 
%   code optimisations and refactoring reduced CPU time by 1-2 %.
% v1.5.4: Code optimisation: Error handling during filtering has been
%   improved.
%

%% Validate number of input and output arguments
narginchk(2, 3);
nargoutchk(0, 4);

%% Convert row vectors to column vectors if needed
if size(xs, 1) == 1
    xs = xs(:);
end

if isscalar(ts)
    % User-specified sampling frequency 
    fs = ts;
    ts = 0.5 / fs;
else
    % Compute sampling frequency from time stamps
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
               "EngineeringUnit",      "EU", ...
               "dbReference",          NaN, ...
               "PlotLayout",           "none");
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
        conSc = 1;                           % Scaling constant
        indSt = (size(freq, 1) - 1) / 2 + 1; % Index of the static component
        ind   = [setup.FFTLength - indSt + 2:setup.FFTLength, 1:indSt].';
    otherwise       % One-sided spectrum 
        freq = (0:setup.FrequencyResolution:setup.LowPassFrequency).';
        setup.Mode = "onesided";

        % Vector of indices for manipulation with the DFT
        conSc = 2;                           % Scaling constant
        indSt = 1;                           % Index of the static component
        ind   = (1:size(freq, 1)).';
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
            setup.Window = makima(0:1/(size(setup.Window, 1)-1):1, ...
                                  setup.Window, (0:1/(setup.FFTLength-1):1).');
        % Use the spline interpolation (Matlab v2019a or older)
        catch 
            setup.Window = spline(0:1/(size(setup.Window, 1)-1):1, ...
                                  setup.Window, (0:1/(setup.FFTLength-1):1).');
        end

        % Throw warning
        warning("The window function has been interpolated to " + ...
                num2str(setup.FFTLength, "%d") + " samples.");
    end
else
    % Generate the window function internally
    switch lower(extractBefore(setup.Window, 2))
        case "b"    % Blackmann-Harris
            setup.Window = utilities.autoBlackmanHarris(setup.FFTLength);
            windowName   = "Blackmann-Harris";
        case "f"    % flat-top
            setup.Window = utilities.autoFlatTop(setup.FFTLength);
            windowName   = "flat-top";
        case "h"    % Hann or Hamming
            if strncmpi(setup.Window, "ham", 3) % Hamming
                setup.Window = utilities.autoHamming(setup.FFTLength);
                windowName   = "Hamming";
            else    % Hann
                setup.Window = utilities.autoHann(setup.FFTLength);
                windowName   = "Hann";
            end
        case "k"    % Kaiser
            % Try to calculate parameter beta
            beta = str2double(regexprep(setup.Window,"[a-zA-Z,=\s]",""));

            if isnan(beta) % Use the default Kaiser window (beta = 1.6)
                setup.Window = utilities.autoKaiser(setup.FFTLength, 1.6);
                windowName   = "Kaiser, b = 0.5";
            else
                setup.Window = utilities.autoKaiser(setup.FFTLength, beta);
                windowName   = "Kaiser, b = " + string(beta);
            end
        case "m"    % Hamming
            setup.Window = utilities.autoHamming(setup.FFTLength);
            windowName   = "Hamming";
        otherwise   % uniform
            setup.Window = utilities.autoUniform(setup.FFTLength);
            windowName   = "uniform";
    end
end

% Normalise the window function for correct magnitude estimation
winMean = sum(setup.Window) / setup.FFTLength;
setup.Window = setup.Window ./ winMean;

% Calculate the noise power bandwidth of the window function in Hz
% exploiting the fact that mean(setup.Window) == 1
setup.WindowNoiseBandwidth = fs * sum(setup.Window.^2) / setup.FFTLength.^2;

%% Filtering
if ~isnan(setup.HighPassFrequency)
    % setup.HighPassFrequency is a single number
    if isnumeric(setup.HighPassFrequency)
        % Remove the mean value from xs
        if setup.HighPassFrequency >= 0
            xs = detrend(xs, 0);
        end
        
        % Construct and apply a Butterworth high-pass filter
        if setup.HighPassFrequency ~= 0
            n  = 1;                                 % Filter order
            fc = abs(setup.HighPassFrequency) / 10; % Cutoff frequency (-3 dB)
    
            % Design filter using an internal function
            [b, a] = utilities.autoButter(n, fc, fs, 'high');
            
            % Filter xs
            xs = filter(b, a, xs);
        end

    % setup.HighPassFrequency is not a single number
    else
        try
            numerator and denominator coefficients b and a
            % Try to filter xs
            if iscell(setup.HighPassFrequency) && all(cellfun(@isnumeric,c))
                % setup.HighPassFrequency is a cell array consisting of
                % vectors b and a (i.e. {[b], [a]}). Vectors b and a
                % contain numerator and denominator coefficients of a
                % rational transfer function
                temp = filter(setup.HighPassFrequency{1}(:), ...
                              setup.HighPassFrequency{2}(:), xs);
            else
                % Try to apply object stored in setup.HighPassFrequency
                % Note: This is a crude solution, but Matlab can recognise
                %   filters only if the DSP Toolbox is installed using e.g.
                %   isfir and issos functions
                temp = filter(setup.HighPassFrequency, xs);
            end

            % Validate size of the filtered output
                if all(size(xs) == size(temp))
                    xs = temp;
                else
                    warning("Filtered data has unexpected size." + newline + ...
                            "Frequency analysis will be performed on unfiltered data.");
                end
        catch ME
            % Throw warning if the filtering fails
            warning(ME.identifier, "Filtering failed with error: %s" + newline + ...
                    "Frequency analysis will be performed on unfiltered data.", ...
                    ME.message);
        end
    end
end

%% FFT
% Preallocate an array for the FFT of individual segments
tSegments = zeros(setup.FFTLength, setup.NumberOfAverages, size(xs, 2));

% Apply time weighting using partial vectorisation
for i = 1:setup.NumberOfAverages
    tSegments(:, i, :) = setup.Window .* xs(seg(i,1):seg(i,2), :);
end

% Fast Fourier transform of the time-weighted segments
tSpectrum = fft(tSegments, [], 1);

% Application of the jw weigthing and computation of unscaled magnitudes
switch lower(setup.jwWeigthing)
    case {"1/jw2", "double integration"}
        % Apply double integration
        conjw = -4 * pi^2;
        tSpectrum(ind,:,:) = abs((1 ./ (conjw .* freq.^2)) .* ...
                                 tSpectrum(ind,:,:));
        setup.jwWeigthing  = "double integration";      % Update setup

    case {"1/jw", "single integration"}
        % Apply single integration
        conjw = 2i * pi;
        tSpectrum(ind,:,:) = abs((1 ./ (conjw .* freq)) .* ...
                                 tSpectrum(ind,:,:));
        setup.jwWeigthing  = "single integration";      % Update setup

    case {"jw", "single differentiation"}
        % Apply single differentiation
        conjw = 2i * pi;
        tSpectrum(ind,:,:) = abs(conjw .* freq .* tSpectrum(ind,:,:));
        setup.jwWeigthing  = "single differentiation"; % Update setup

    case {"jw2", "double differentiation"}
        % Apply double differentiation
        conjw = -4 * pi^2;
        tSpectrum(ind,:,:) = abs(conjw .* freq.^2 .* tSpectrum(ind,:,:));
        setup.jwWeigthing  = "double differentiation"; % Update setup

    otherwise
        % Compute the absolute value of spectrum
        tSpectrum(ind,:,:) = abs(tSpectrum(ind,:,:));
        setup.jwWeigthing  = "none";                   % Update setup
end

% Spectral averaging and the limitation of the maximum frequency
switch lower(setup.Averaging)
    case {"energy", "rms"}
        % Perform energy averaging
        spectrum = rms(tSpectrum(ind,:,:), 2);
        setup.Averaging = "energy";     % Update setup

    case {"maximum", "max", "pk", "peak"}
        % Permorm maximum-hold averaging
        spectrum = max(tSpectrum(ind,:,:), [], 2);
        setup.Averaging = "maximum";    % Update setup

    case {"median", "med"}
        % Perform median-hold averaging
        spectrum = median(tSpectrum(ind,:,:), 2);
        setup.Averaging = "median";     % Update setup

    case {"minimum", "min"}
        % Perform minimum-hold averaging
        spectrum = min(tSpectrum(ind,:,:), [], 2);
        setup.Averaging = "minimum";    % Update setup

    case "none" 
        % No averaging
        spectrum = tSpectrum(ind,:,:);

    case {"variance", "var"}
        % Compute variance at each spectral line
        spectrum = var(tSpectrum(ind,:,:), 0, 2);
        setup.Averaging = "variance";   % Update setup

    otherwise
        % Perfotm linear averaging
        spectrum = mean(tSpectrum(ind,:,:), 2);
        setup.Averaging = "linear";     % Update setup
end

% Evaluation of spectral unit
% Prepare indices 
indDyn = [1:indSt-1, indSt+1:size(spectrum, 1)].';

switch lower(setup.SpectralUnit)
    case "rms"
        % RMS magnitude
        conSc = conSc / (sqrt(2) * setup.FFTLength);
        spectrum(indSt,:)  = spectrum(indSt,:) ./ setup.FFTLength;   % 0 Hz
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:);
        setup.SpectralUnit = "RMS";  % Update setup

    case {"pk", "0-pk", "peak"} 
        % 0-peak magnitude
        conSc = conSc / setup.FFTLength;
        spectrum(indSt,:)  = spectrum(indSt,:) ./ setup.FFTLength;   % 0 Hz
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:);
        setup.SpectralUnit = "0-pk";  % Update setup

    case {"pp", "pk-pk", "peak2peak"}
        % Peak-peak magnitude
        conSc = 2 * conSc / setup.FFTLength;
        spectrum(indSt,:)  = spectrum(indSt,:) ./ setup.FFTLength;   % 0 Hz
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:);
        setup.SpectralUnit = "pk-pk"; % Update setup

    case {"asd", "psd"}
        % Power spectral density
        conSc = conSc / (setup.WindowNoiseBandwidth * setup.FFTLength.^2);
        spectrum(indSt,:)  = spectrum(indSt,:).^2 ./ ...             % 0 Hz
                         (setup.WindowNoiseBandwidth * setup.FFTLength.^2);
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:).^2;
        setup.SpectralUnit = "PSD";   % Update setup

    case {"rsd", "rmssd"}
        % Root mean square spectral density
        conSc = conSc / (sqrt(2) * setup.WindowNoiseBandwidth * setup.FFTLength);
        spectrum(indSt,:)  = spectrum(indSt,:) / ...                 % 0 Hz
                            (setup.WindowNoiseBandwidth * setup.FFTLength);
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:);
        setup.SpectralUnit = "RMSSD"; % Update setup

    otherwise
        % Autospectrum / power spectrum
        conSc = conSc / setup.FFTLength.^2;
        spectrum(indSt,:)  = spectrum(indSt,:).^2 ./ setup.FFTLength.^2; % 0 Hz
        spectrum(indDyn,:) = conSc .* spectrum(indDyn,:).^2;
        setup.SpectralUnit = "power"; % Update setup
end

% Squeeze the resulting spectrum
spectrum = squeeze(spectrum);

% Transform resulting spectra to decibels
if ~isnan(setup.dbReference)
    % Set the reference level automatically in not specified by the user
    if setup.dbReference == 0
        if setup.SpectralUnit == "PSD" || setup.SpectralUnit == "power"
            setup.dbReference = sqrt(max(spectrum,[],"all"));
        else
            setup.dbReference = max(spectrum,[],"all");
        end
    end

    % Calculate decibels
    if setup.SpectralUnit == "PSD" || setup.SpectralUnit == "power"
        % Calculate decibels for powerspectra, i.e. 20 log (s_i / s_ref)
        % Note that the spectrum array contains squared values, i.e. s_i^2
        spectrum = 10 * log10(spectrum / (setup.dbReference^2));
    else
        % Calculate decibels for linear spectra, i.e. 10 log (s_i / s_ref)
        spectrum = 10 * log10(spectrum / setup.dbReference);
    end  
end

%% Display results
% Compute segment times
if nargout > 3 || (setup.NumberOfAverages > 1 && setup.Averaging == "none")
    % Compute the relative segment times
    tshift = (setup.FFTLength - setup.OverlapLength) / fs;
    tseg   = (0:tshift:(setup.NumberOfAverages - 1) * tshift).';

    % Compute the absolute segment times        
    tseg = sum(setup.Window .* (0:setup.FFTLength-1).') / (setup.FFTLength * fs) + ts(1) + tseg;
end

% Display results
if nargout == 0 || setup.PlotLayout ~= "none"
    % Change the PlotLayout parameter if the user expects graphical output
    if setup.PlotLayout == "none"
        setup.PlotLayout = "separated";
    end

    % Display spectra
    if setup.NumberOfAverages == 1 || setup.Averaging ~= "none"
        utilities.autoPlot(setup, spectrum, freq);

    % Display spectrograms
    else
        utilities.autoPlot(setup, spectrum, freq, tseg);
    end
end
 
% Arrange outputs
% Return the analyser setup only
if nargout == 3 && (setup.NumberOfAverages == 1 || setup.Averaging ~= "none")
    % Convert numeric representation of the window to string
    if exist("windowName", "var")
        setup.Window = windowName;
    end
   
    varargout{1} = setup;
 
% Return the segment times only
elseif nargout == 3 && setup.NumberOfAverages > 1 && setup.Averaging == "none"
    varargout{1} = tseg;

% Return segment times and the analyser setup
elseif nargout > 3
    % Convert numeric representation of the window to string
    if exist("windowName", "var")
        setup.Window = windowName;
    end

    varargout{1} = tseg;
    varargout{2} = setup;
end

% End of main fucntion
end