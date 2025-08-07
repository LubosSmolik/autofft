% AutoFFT - Performance measurement (supplement function)
% Authors: Luboš Smolík, Jan Rendl and Roman Pašek
% 2022, doi: 10.24433/CO.3004443.v1
% 
% This code measures CPU times of autofft and Matlab built-in functions. It is also
% published as reproducible capsule at https://codeocean.com/capsule/1932803/tree/v1
% 
% Requires the Signal Processing Toolbox

function [s, f, t] = timeit_stft(x, fs, window, overlap)

% Perform STFT
[s, f, t] = stft(x, fs, Window = window, OverlapLength = overlap, ...
                 FrequencyRange = "onesided");

% s is a complex array and has to be converted to a power spectrum
s = (2 / length(window).^2) .* s .* conj(s);

% Note: For a simplicity sake, the static component is two times higher
%       than it should be.

end
