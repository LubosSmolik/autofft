function [ws, setup] = freqWeight(s, f, wType, setup)
% FREQWEIGHT Applies frequency weighting filter to the input power spectrum
%
% Copyright (c) 2022-2025           Lubos Smolik, University of West Bohemia
% v1.0.1 (build 2. 8. 2025)         e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
% ws = autofft(s, f, wType)
% [ws, setup] = autofft(s, f, wType, setup)
%
% ws = autofft(s, f, wType) applies the frequency weighting wType to power
%   spectrum s computed at frequencies f and returns weighted power
%   spectrum ws. Both s and ws are given in dB.
% 
% [ws, setup] = autofft(s, f, wType, setup) also insert the information
%   about the applied frequency weighting to the setup of the FFT analyser.
% 
% Input variable wType can be selected from the following list:
% "A" - Applies A-type frequency weighting specified in IEC 61672-1:2013.
%       The A-type frequency weighting accounts for the relative loudness
%       perceived by the human ear which is less sensitive to sound at low
%       and high frequencies. This weighting is intented for noise
%       measurements.
% "B" - Applies B-type frequency weighting specified in IEC 60651:1979
%       (withdrawn). It was developed because the A-type frequency
%       weighting does not properly describe the relative loudness at low
%       to intermediate levels.
% "C" - Applies B-type frequency weighting specified in IEC 61672-1:2003.
%       Similarly to the A-type frequency weighting, the C-type frequency
%       weighting attenuates components at low frequencies. This weighting
%       is designed for peak sound pressure measurements.
% "D" - Applies D-type frequency weighting specified in IEC 60651:1979
%       (withdrawn). This weighting was designed for measuring high-level
%       aircraft noise in accordance with IEC 537.
% "M" - Applies M-type frequency weighting as specified in ISO 21727:2016.
%       The M-type frequency weighting is intended for measurement of short
%       duration motion-picture sound and allows assessment of the
%       subjective loudness and annoyance in rooms calibrated in accordance
%       to ISO 2969. This weighting emulates ITU-R 468 noise weighting but
%       places 0 dB offset to 2 kHz rather than 1 kHz. 
% "Z" - Applies Z-type frequency weighting specified in IEC 61672-1:2003.
%       This option returns the input power spectrum unchanged.

% CHANGELOG
% v1.0.1 - Terminology complies with ISO 21727:2016 and IEC 61672-1:2013.
%        - Weighting functions now complies with IEC 61672-1:2013.

%% Validate input arguments
narginchk(3, 4);

% Convert wType to lower case
wType = lower(wType);

%% Apply the frequency weighting filter
% Specify the reference frequency for the offset value
if wType == "m"
    fref = 2000;
else
    fref = 1000;
end

% Compute offset value
off = 20 * log10(weighting(fref, wType));   

% Compute array of weights (in dB)
w   = 20 * log10(weighting(f, wType)) - off;
w   = repmat(w(:), size(s, 2), size(s, 3));

% Apply the weights
ws  = s + w;

%% Store information about the frequency weighting to setup input argument
if nargin == 4
    switch wType
        case "a"
            setup.FrequencyWeighting = "A-weighting";
        case "b"
            setup.FrequencyWeighting = "B-weighting";
        case "c"
            setup.FrequencyWeighting = "C-weighting";
        case "d"
            setup.FrequencyWeighting = "D-weighting";
        case "m"
            setup.FrequencyWeighting = "M-weighting";
        otherwise
            setup.FrequencyWeighting = "Z-weighting";
    end
end

%% Return vector of weights at frequencies f
%
function wSpec = weighting(f, weightType)
    switch weightType
        case "a"
            % A-type frequency weighting per IEC 61672-1:2013
            % Denominators
            den   = [20.6, 107.7, 737.9, 12194].^2;

            % Compute weighting function
            wSpec = (den(4) * f.^4) ./ ((f.^2 + den(1)) .* sqrt((f.^2 + den(2)) .* (f.^2 + den(3))) .* (f.^2 + den(4)));

        case "b"
            % B-type frequency weighting per IEC 60651:1979 (withdrawn)
            % Denominators
            den   = [20.6, 158.5, 12194].^2;
            
            % Compute weighting function
            wSpec = (den(3) * f.^3) ./ ((f.^2 + den(1)) .* sqrt((f.^2 + den(2))) .* (f.^2 + den(3)));

        case "c"
            % C-type frequency weighting per IEC 61672-1:2013
            % Denominators
            den   = [20.6, 12194].^2;

            % Compute weighting function
            wSpec = (den(2) * f.^2) ./ ((f.^2 + den(1)) .* (f.^2 + den(2)));

        case "d"
            % D-type frequency weighting per IEC 60651:1979 (withdrawn)
            % Coefficients and denominators
            cfs = [1037918.48, 1080768.16, 9837328, 11723776];
            den = [6.8966888496476e-5, 79919.29, 1345600];
            
            % Auxiliary response function
            hf  = ((cfs(1) - f.^2).^2 + cfs(2) * f.^2) ./ ((cfs(3) - f.^2).^2 + cfs(4) * f.^2);

            % Compute weighting function
            wSpec = (f / den(1)) .* sqrt(hf ./ ((f.^2 + den(2)) .* (f.^2 + den(3))));

        case "m"
            % M-type frequency weighting per ISO 21727:2016
            % Coefficients
            cfs1 = [-4.737338981378384e-24,  2.043828333606125e-15, -1.363894795463638e-7];
            cfs2 = [ 1.306612257412824e-19, -2.118150887518656e-11,  5.559488023498642e-4];

            % Auxiliary response functions
            h1f  = cfs1(1) * f.^6 + cfs1(2) * f.^4 + cfs1(3) * f.^2 + 1;
            h2f  = cfs2(1) * f.^5 + cfs2(2) * f.^3 + cfs2(3) * f;
            
            % Compute weighting function
            wSpec = (1.246332637532143e-4 * f) ./ sqrt(h1f.^2  + h2f.^2);

        otherwise
            % Z-type frequency weighting per IEC 61672-1:2013
            wSpec = zeros(size(f));
    end
end

% End of main function
end