function [spectrum, freq, setup] = autofft(xs, ts, fftset)
% AUTOFFT Evaluates a frequency spectrum of a signal using wFFT algorithm
%
%  Copyright (c) 2017-2020         Lubos Smolik, University of West Bohemia
% v1.2.5a (build 30. 1. 2020)            e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-2-Clause License.
%
% [s, f] = autofft(xs, fs)
% [s, f] = autofft(xs, ts)
% [___]  = autofft(___, options)
% [___, setup] = autofft(___)
%
% [s, f] = autofft(xs, fs) returns the DFT or STFT s of xs using sampling
%   frequency fs (Hz). Also returns the frequencies f at which the DFT or
%   STFT s is evaluated. xs can be either a vector or an array consisting
%   of column vectors.
% [s, f] = autofft(xs, ts) returns the DFT or STFT s of xs using a vector
%   of time stamps ts (s). Also returns the frequencies f at which the DFT
%   or STFT s is evaluated. 
% [___]  = autofft(___, options) performs the DFT or STFT name-value pair
%   arguments specified in structure options.
% [s, f, setup] = autofft(___) also returns the setup of the FFT analyser.
%       
% Construction of options:
%	options = struct('param', 'value', ...);
%
%	List of parameters:
%     - 'nwin' - [ positive integer {length(xs) or size(xs, 2)} ]
%       The length of each segment in samples, i.e. number of DFT points.
%
%     - 'twin' - [ real positive scalar ]
%       The duration of each segment in seconds.
%
%     - 'df' - [ real positive scalar ]
%       The desired frequency resolution of the analyser in Hz.
%
%       Note: Default values of the three parameters above are set so that
%             the input data are treated as one segment, i.e. no spectral 
%             averaging is performed.
%             If more than one of these three parameters is specified,
%             'nwin' is preferred over 'twin' and 'twin' is preferred over
%             'df'.
%
%     - 'highpass' - [ {NaN} | real scalar ]
%       Specifies the passband frequency of an elliptic highpass filter.
%       The filter has a slope -20 dB/dec from the passband frequency.
%       The exact filtering process is implemented as follows:
%       - NaN             - No highpass filtering.
%       - 0               - DC filtering with the use of detrend function.
%       - positive scalar - DC filtering and subsequent highpass filtering.
%       - negative scalar - Highpass filtering without DC filtering.
%
%     - 'lowpass' - [ positive scalar in (0, fs/2) {fs/2} ]
%       Specifies the maximum frequency in Hz over which autofft computes
%       estimates. Recommended value is fs/2.56.
%
%     - 'window' - [ character {'u'} | vector ]
%       Specifies a window function which is than multiplied with each
%       segment. The window function can be specified either as a vector
%       or using a character from the list below:
%       - 'b' - Blackmann-Harris window
%       - 'f' - flat-top window
%       - 'h' - Hann window
%       - 'm' - Hamming window
%       - 'k' - Kaiser window with shape factor beta = 0.5
%       - 'kA.A' - Kaiser window with shape factor beta = A.A
%       - 'u' - uniform (rectangular) window                      {default}
%
%     - 'overlap' - [ real scalar in [0, 100) {50} ] 
%        Overlap percentage of two successive segments.
%
%     - 'averaging' - [ energy | {linear} | max | median | min | var | none]
%       Specifies the spectral averaging mode from: 
%       - 'energy', 'rms'     - energy (rms) averaging
%       - 'linear', 'lin'     - linear averaging                  {default}
%       - 'max', 'pk', 'peak' - maximum-hold aveargaing
%       - 'median'            - returns median specified spectral unit
%       - 'min'               - minimum-hold aveargaing
%       - 'var'               - returns variance of specified spectral unit
%       - 'none'              - returns the STFT with no spectral averaging
%
%     - 'jw' - [ 1/jw2 | 1/jw | {1} | jw | jw2 ]
%       Use this parameter to apply a frequency-domain post-weighting to
%       the output spectra, for example to estimate displacement from
%       acceleration. The parameter can be specified as follows:
%          - '1/jw2' - double integration
%          - '1/jw'  - single integration
%          - '1'     - as input                                   {default}
%          - 'jw'    - single differentiation
%          - 'jw2'   - double differentiation
%
%     - 'unit' - [ {pow} | rms | pk | pp | psd | rsd ]
%       Specifies absolute unit used to compute estimates from:
%        - 'pow', 'power' - autospectrum (square of rms magnitudes)  {def.}
%        - 'rms'          - linear spectrum with rms magnitude
%        - 'pk', '0-pk'   - linear spectrum with 0-peak magnitude
%        - 'pp', 'pk-pk'  - linear spectrum with peak-peak magnitude
%        - 'psd'          - power spectral density
%        - 'rsd','rmssd'  - root mean square of power spectral density 
%
% Changelist
% v1.2.5a- Outputs which appeared twice in setup now appear only once.
% v1.2.5 - An optional highpass filtering has been implemented.
% v1.2.4b- Overlap length (in samples) returned in setup is now correct.
% v1.2.4a- Error occuring during estimation of autospectrum has been fixed.
% v1.2.4 - Documentation has been improved.
%        - Results of the STFT of multiple signals are now returned as 3D
%           array rather than cell array of 2D arrays.
%        - New types of averaging: median filter and variance of spectral unit. 
%        - Performance has been improved significantly.
%        - Accuracy of PSD and RMSSD estimates has been slightly improved.
%        - Dealing with a content at the Nyquist frequency has been improved.
% v1.2.3 - User can define frequency resolution of the analyser.
%       - The window function can be directly specified as a vector.
%       - The analyser setup can be returned as an output variable.
%       - Warning messages are now displayed.
% v1.2.2a- Error occuring during peak hold averaging has been fixed.
% v1.2.2 - The Kaiser-Bessel window parameter (beta) can now be specified.
%       - Autospectrum is now properly square of RMS rather than 0-Pk.
%       - Relations for evaluation of PSD and RMSPSD now consider the noise
%         power bandwidth of the used window function. 
% v1.2.1 - New function - low-pass filtering
%        - New types of averaging - no averaging     ('none')
%                                 - energy averaging ('energy' or 'rms')
%                                 - minimum value    ('min')
%        - Options for 'unit' and 'peak' has been merged (into 'unit').
%        - Relations for evaluation of PSD amd RMSPSD have been fixed.
%        - Dealing with a content at the Nyquist frequency has been fixed.
% v1.2.0 - Input parameters are now specified in a structured variable.
%        - v1.2 is not compatible with v1.12 and older versions!
% v1.1.2 - Input can now be an array.
% v1.1.1 - Performance optimization
%        - Handling of input vectors with the even number of samples has
%          been fixed.
% v1.1.0 - Handling of non-uniform window functions has been fixed.
%        - Parameters nwin and overlap can now be skipped by user.
%
%% nargin check
if nargin < 2
    error('Not enough input arguments.');
elseif nargin > 3
    error('Too many input arguments.');
end
%
%% Convert row vectors to column vectors if needed
if size(xs, 1) == 1         % samples
    xs = xs(:);                     
end
if size(ts(:), 1) == 1      % sampling frequency
	fs = ts;                    
else
    fs = 1 / (ts(2) - ts(1)); 
end
%
%% Specify the default setup
defset = struct('nwin', size(xs, 1), ...
                'overlap', 50, ...
                'highpass', NaN, ...
                'lowpass', fs/2, ...
                'window', 'u', ...
                'averaging', 'lin', ...
                'jw', '1', ...
                'unit', 'pow');
deffields = fieldnames(defset);
%
%% Set analyser parameters
if nargin == 2  % use default fftset
    fftset = defset;
else            % use user-defined fftset  
    % Check if there is user-defined 'nwin' parameter
    if isfield(fftset, 'nwin') && fftset.nwin > size(xs, 1)
        fftset.nwin = size(xs, 1);
        warning("Window function has more samples than input data. " + ...
                "Length of the window function has been changed to " + ...
                num2str(fftset.nwin, '%d') + " samples.");
    % Check if there is user-defined 'twin' parameter        
    elseif isfield(fftset, 'twin')
        if round(fftset.twin * fs) > size(xs, 1)
            fftset.nwin = size(xs, 1);
            warning("window function is longer than input data. " + ...
                    "Length of the window function has been changed " + ...
                    "to " + num2str(fftset.nwin / fs, '%.2g') + " s.");           
        else
            fftset.nwin = round(fftset.twin * fs);
        end
        
    % Check if there is user-defined 'df' parameter
    elseif isfield(fftset, 'df')
        if round(fs / fftset.df) > size(xs, 1)
            fftset.nwin = size(xs, 1);
            warning("Specified frequency resolution cannot be reached. " + ...
                    "The frequency resolution has been changed to " + ...
                    num2str(fs / fftset.nwin, '%.2g') + " Hz.");           
        else
            fftset.nwin = round(fs / fftset.df);
        end
    end
    
    % Set unspecified parameters to default
    for i = 1:numel(deffields)        
        if isfield(fftset, deffields{i}) == 0
            fftset.(deffields{i}) = defset.(deffields{i});
        end
    end
end
% Generate frequency vector
freq = (fs * (0:(fftset.nwin/2)) / fftset.nwin)';
% Set allowed frequencies for the limitation of the maximum freqeuncy
freq = freq(freq <= fftset.lowpass);
maxf = size(freq, 1);
% Calculate number of overlaping samples
fftset.overlap = round(fftset.nwin * fftset.overlap / 100);
% Set indices for the signal segmentation
imax = floor((size(xs, 1)-fftset.overlap) / (fftset.nwin-fftset.overlap));
                                                % number of windows
ind = zeros(imax,2);                            % matrix of indices
ni = 1;                                         % pointer
for i = 1:imax                                  % cycle through windows
    ind(i,1) = ni; 
    ni = ni + fftset.nwin - 1;
    ind(i,2) = ni;
    ni = ni - fftset.overlap + 1;
end
% Generate a structure array containing the analyser setup
setup = struct("SamplingFrequency",    fs, ...
               "DataDuration",         size(xs, 1) / fs, ...
               "DataLength",           size(xs, 1), ...
               "HighPassFrequency",    fftset.highpass, ...
               "LowPassFrequency",     fftset.lowpass, ...
               "FFTLength",            fftset.nwin, ...
               "FrequencyResolution",  fs / fftset.nwin, ...
               "TimeResolution",       fftset.nwin / fs, ...
               "Window",               "uniform", ...
               "WindowNoiseBandwidth", 1, ...
               "OverlapLength",        fftset.overlap, ...
               "OverlapPercentage",    100 * fftset.overlap / fftset.nwin, ...
               "Averaging",            "none", ...
               "NumberOfAverages",     imax, ...
               "SpectralUnit",         "power", ...
               "jwWeigthing",          "none");
% Set a window function
if isnumeric(fftset.window)
    % Use the custom user-specified window function
    if length(fftset.window) == fftset.nwin
        % Add info to the structure array containing the analyser setup
        setup.Window = "custom";
    else
        % Interpolate missing values in the window function if necessary
        try   % Use the modified Akima interpolation (v2019b or newer)
            fftset.window = makima(linspace(0,1,length(fftset.window)), ...
                                 fftset.window, linspace(0,1,fftset.nwin));
        catch % Use the spline interpolation (v2019a or older)
            fftset.window = spline(linspace(0,1,length(fftset.window)), ...
                                 fftset.window, linspace(0,1,fftset.nwin));
        end
        
        % Print warning and add info to the information resource
        warning("The window function has been reinterpolated to " + ...
                num2str(fftset.nwin, "%d") + " samples.");
        % Add info to the structure array containing the analyser setup
        setup.Window = "custom reinterpolated to " + ...
                       num2str(fftset.nwin, "%d") + " samples ";
    end
else
    % Generate the window function internally
    [fftset.window, setup.Window] = windowfunc(fftset.window, fftset.nwin);
end
% Normalise the window function for correct magnitude estimation
fftset.window = fftset.window / mean(fftset.window);
% Calculate the noise power bandwidth of the window function in Hz
fftset.noiseband = enbw(fftset.window, fs);
% Add info to the structure array containing the analyser setup
setup.WindowNoiseBandwidth = fftset.noiseband;

%% Filtering
if ~isnan(fftset.highpass)
    % Remove the mean value from xs
    if fftset.highpass >= 0
        xs = detrend(xs, 0);
    end
    
    % Apply highpass filtering
    if fftset.highpass ~= 0
        n     = 1;    % Filter order
        aPass = 0.1;  % Passband ripple (dB)
        aStop = 20;   % Stopband attenuation (dB)

        % Construct an FDESIGN object and call its ELLIP method.
        h  = fdesign.highpass('N,Fp,Ast,Ap', n, abs(fftset.highpass), ...
                                             aStop, aPass, fs);
        Hd = design(h, 'ellip');
        
        % Filter xs
        xs = filter(Hd, xs);
    end
end

%% FFT
% Preallocate an array for the DFT of individual segments
tSpectrum = zeros(fftset.nwin, size(xs, 2), imax);
% Fast Fourier transformation of the time-weighted segments
for i = 1:size(xs, 2)
    for j = 1:imax
        tSpectrum(:, i, j) = fft(fftset.window .* xs(ind(j,1):ind(j,2), i), ...
                                 fftset.nwin);
    end
end
% Scaling and application of the jw weigthing
switch lower(fftset.jw)
    case '1/jw2'
        tSpectrum(1:maxf, :, :) = (- 1 ./ (4 * pi^2 * freq.^2)) .* ...
                                  (tSpectrum(1:maxf, :, :) / fftset.nwin);
        setup.jwWeigthing = "double integration";
    case '1/jw'
        tSpectrum(1:maxf, :, :) = (1 ./ (2i * pi * freq)) .* ...
                                  (tSpectrum(1:maxf, :, :) / fftset.nwin);
        setup.jwWeigthing = "single integration";
    case 'jw'
        tSpectrum(1:maxf, :, :) = 2i * pi * freq .* ...
                                  tSpectrum(1:maxf, :, :) / fftset.nwin;
        setup.jwWeigthing = "single differentiation"; 
    case 'jw2'
        tSpectrum(1:maxf, :, :) = - 4 * pi^2 * freq.^2 .* ...
                                  tSpectrum(1:maxf, :, :) / fftset.nwin;
        setup.jwWeigthing = "double differentiation"; 
    otherwise
        tSpectrum(1:maxf, :, :) = tSpectrum(1:maxf, :, :) / fftset.nwin;
end
% Evaluation of spectral unit
switch lower(fftset.unit)
    case 'rms'           % Linear spectrum with rms magnitude
        tSpectrum(1,:,:)      = abs(tSpectrum(1,:,:));
        tSpectrum(2:maxf,:,:) = (2/sqrt(2)) * abs(tSpectrum(2:maxf,:,:));
        setup.SpectralUnit = "RMS";
    case 'pk'            % Linear spectrum with 0-peak magnitude
        tSpectrum(1,:,:)      = abs(tSpectrum(1,:,:));
        tSpectrum(2:maxf,:,:) = 2 * abs(tSpectrum(2:maxf,:,:));
        setup.SpectralUnit    = "0-pk";
    case 'pp'            % Linear spectrum with peak-peak magnitude
        tSpectrum(1,:,:)      = abs(tSpectrum(1,:,:));
        tSpectrum(2:maxf,:,:) = 4 * abs(tSpectrum(2:maxf,:,:));        
        setup.SpectralUnit = "pk-pk";
    case {'asd','psd'}   % Power spectral density
        tSpectrum(1,:,:)      = tSpectrum(1,:,:) .* conj(tSpectrum(1,:,:)) ...
                                 / fftset.noiseband;
        tSpectrum(2:maxf,:,:) = 2 * tSpectrum(2:maxf,:,:) .* conj( ...
                                tSpectrum(2:maxf,:,:)) / fftset.noiseband;
        setup.SpectralUnit = "PSD";
    case {'rsd','rmssd'} % Root mean square spectral density (RMS of PSD)
        tSpectrum(1,:,:)      = abs(tSpectrum(1,:,:)) / fftset.noiseband;
        tSpectrum(2:maxf,:,:) = (2/sqrt(2)) * abs(tSpectrum(2:maxf,:,:)) ...
                                 / fftset.noiseband;
        setup.SpectralUnit = "RMSSD";
    otherwise            % Autospectrum
        tSpectrum(1,:,:)      = tSpectrum(1,:,:) .* conj(tSpectrum(1,:,:));
        tSpectrum(2:maxf,:,:) = 2 * tSpectrum(2:maxf,:,:) .* conj(tSpectrum(2:maxf,:,:));
end
% Spectral averaging and the limitation of the maximum frequency
switch lower(fftset.averaging)
    case {'energy', 'rms'}     % Energy averaging
        spectrum = rms(tSpectrum(1:maxf,:,:), 3);
        setup.Averaging = "energy";       
        
    case {'max', 'pk', 'peak'}  % Maximum peak hold averaging
        spectrum = max(tSpectrum(1:maxf,:,:), [], 3);
        setup.Averaging = "maximum";
        
    case 'median'               % Maximum peak hold averaging
        spectrum = median(tSpectrum(1:maxf,:,:), 3);
        setup.Averaging = "median";
    case 'min'                  % Minimum peak hold averaging
        spectrum = min(tSpectrum(1:maxf,:,:), [], 3);
        setup.Averaging = "minimum";
    case 'none'                 % No averaging
        spectrum = squeeze(tSpectrum(1:maxf,:,:));
        
    case 'var'                  % Variance of averaging
        spectrum = var(tSpectrum(1:maxf,:,:), 0, 3);
        setup.Averaging = "variance";
    otherwise                   % Linear averaging
        spectrum = mean(tSpectrum(1:maxf,:,:), 3);
        setup.Averaging = "linear";
end
% End of main fucntion
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction windowfunc generates the window function
%
% Input
%   - sym - symbol for the window function
%   - n   - length of the window function (samples)
%
% Output
%   - window - a discrete representation of the window function
%  
function [window, windowName] = windowfunc(sym, n)
    % Generate the specified window function and its name
    switch sym(1)
        case 'b'    % Blackmann-Harris
            window     = blackmanharris(n);
            windowName = "Blackmann-Harris";
        case 'f'    % flat-top
            window = flattopwin(n);
            windowName = "flat-top";
        case 'h'    % Hann
            window = hann(n);
            windowName = "Hann";
        case 'k'    % Kaiser-Bessel
            if length(sym) == 1
                % window with default beta = 0.5
                window = kaiser(n, 0.5);
                windowName = "Kaiser, b = 0.5";
            else
                % window with user specified beta
                window = kaiser(n, str2double(sym(2:end)));
                windowName = "Kaiser, b = " + string(sym(2:end));
            end
        case 'm'    % Hamming
            window = hamming(n);
            windowName = "Hamming";
        otherwise   % uniform
            window = rectwin(n);
            windowName = "uniform";
    end
end
