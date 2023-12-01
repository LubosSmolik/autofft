function [ws, setup] = freqWeight(s, f, wType, setup)
% FREQWEIGHT Applies frequency weighting filter to the input power spectrum
%
%  Copyright (c) 2022              Lubos Smolik, University of West Bohemia
% v1.0.0 (build 27. 4. 2022)       e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
% ws = autofft(s, f, wType)
% [ws, setup] = autofft(s, f, wType, setup)
%
% ws = autofft(s, f, wType) applies the frequency weighting filter wType to
%   power spectrum s computed at frequencies f and retursn weighted power
%   spectrum ws. Both s and ws are given in dB.
% 
% [ws, setup] = autofft(s, f, wType, setup) also insert the information
%   about the applied filter to the setup of the FFT analyser.
% 
% Input variable wType can be selected from the following list:
% "A" - Applies A-weighting as specified in IEC 61672:2003. This weighting
%       accounts for the relative loudness perceived by the human ear since
%       the ear is less sensitive to sound at low and high frequencies.
%       The A-weighting is intented for noise measurements.
% "B" - Applies B-weighting as specified in withdrawn IEC 60651. This
%       weighting was developed because the A-weighting does not properly
%       describe the relative loudness at low to intermediate levels.
% "C" - Applies C-weighting as specified in IEC 61672:2003. Similarly to
%       the A-weighting, the C-weighting attenuates components at low
%       frequencies. The C-weighting is designed for peak sound pressure
%       measurements.
% "D" - Applies D-weighting as specified in withdrawn IEC 60651. The
%       B-weighting was designed for measuring high-level aircraft noise in
%       accordance with IEC 537.
% "M" - Applies M-weighting as specified in ISO 21727. The M-weighting
%       emulates ITU-R 468 noise weighting but places 0 dB offset to 2 kHz
%       rather than 1 kHz. The M-weighting is intended to gauge loudness of
%       cinema soundtracks.
% "Z" - Applies Z-weighting as specified in IEC 61672:2003. This option
%       returns the input power spectrum unchanged.

%% nargin check
narginchk(3, 4);

%% Apply the frequency weighting filter
% Specify the frequency for the offset value
if lower(wType) == "m"
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

%% Store information about the frequency weighting
if nargin == 4
    switch lower(wType)
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
    switch lower(weightType)
        case "a"
            % Denominators
            den   = [20.6, 107.7, 737.9, 12194].^2;

            % Compute weighting function (A-weighting per IEC 61672:2003)
            wSpec = (den(4) * f.^4) ./ ((f.^2 + den(1)) .* sqrt((f.^2 + den(2)) .* (f.^2 + den(3))) .* (f.^2 + den(4)));

        case "b"
            % Denominators
            den   = [20.6, 158.5, 12194].^2;
            
            % Compute weighting function (B-weighting per IEC 60651 - now withdrawn)
            wSpec = (den(3) * f.^3) ./ ((f.^2 + den(1)) .* sqrt((f.^2 + den(2))) .* (f.^2 + den(3)));

        case "c"
            % Denominators
            den   = [20.6, 12194].^2;

            % Compute weighting function (C-weighting per IEC 61672:2003)
            wSpec = (den(2) * f.^2) ./ ((f.^2 + den(1)) .* (f.^2 + den(2)));

        case "d"
            % Coefficients and denominators
            cfs = [1037918.48, 1080768.16, 9837328, 11723776];
            den = [6.8966888496476e-5, 79919.29, 1345600];
            
            % Auxiliary response function
            hf  = ((cfs(1) - f.^2).^2 + cfs(2) * f.^2) ./ ((cfs(3) - f.^2).^2 + cfs(4) * f.^2);

            % Compute weighting function (D-weighting per IEC 60651 - now withdrawn)
            wSpec = (f / den(1)) .* sqrt(hf ./ ((f.^2 + den(2)) .* (f.^2 + den(3))));

        case "m"
            % Coefficients
            cfs1 = [-4.737338981378384e-24,  2.043828333606125e-15, -1.363894795463638e-7];
            cfs2 = [ 1.306612257412824e-19, -2.118150887518656e-11,  5.559488023498642e-4];

            % Auxiliary response functions
            h1f  = cfs1(1) * f.^6 + cfs1(2) * f.^4 + cfs1(3) * f.^2 + 1;
            h2f  = cfs2(1) * f.^5 + cfs2(2) * f.^3 + cfs2(3) * f;
            
            % Compute weighting function (M-weighting - emulation of ITU-R 468 noise weighting)
            wSpec = (1.246332637532143e-4 * f) ./ sqrt(h1f.^2  + h2f.^2);

        otherwise
            % Compute weighting function (Z-weighting per IEC 61672:2003)
            wSpec = zeros(size(f));
    end
end
% End of main function
end