function [Z, P, G] = autoButter(n, wn, varargin)
% AUTOBUTTER Butterworth digital filter design
%
% Copyright (c) 2014                Jan Simon - original code
% Copyright (c) 2022                Lubos Smolik - validation, revisions
% v1.1.0 (build 14. 7. 2022)        e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
% [b, a] = autoButter(n, wn)
% [b, a] = autoButter(n, wn, ftype)
% [b, a] = autoButter(n, fc, fs, ftype)
% [z, p, g] = autoButter(...)
%
% [b,a] = autoButter(n,wn) designs an nth order lowpass digital Butterworth
%   filter with 3 dB attenuation at cutoff frequency wn. Cutoff frequency
%   wn must be 0 < wn < 1, with 1 being the Nyquist frequency, i.e. half of
%   the sampling frequency. Vectors b and a contain filter numerator and
%   denominator coefficients, respectively. The coefficients are listed in
%   descending powers of z. 
%   If wn is a two-element vector, i.e. wn = [w1 w2], autoButter returns a
%   bandpass filter with passband w1 < w < w2.
%
% [b,a] = autoButter(n,wn,ftype) allows for the definition of the filter
%   type using ftype parameter:
%   [b,a] = autoButter(n,wn,'high') designs a highpass filter.
%   [b,a] = autoButter(n,wn,'low') designs a lowpass filter.
%   [b,a] = autoButter(n,wn,'stop') is a bandstop filter if wn = [w1 w2].
%   [b,a] = autoButter(n,wn,'bandpass') designs a bandpass filter if wn =
%                                                                  [w1 w2].
%   
% [b,a] = autoButter(n,fc,fs,ftype) designs a digital Butterworth filter
%   with cutoff frequency or frequencies fc and sampling frequency fs.
%   Parameter ftype specifies the filter type.
%
% [z,p,g] = autoButter(...) returns vectors containing zeros z and poles p,
%   and scalar g corresponding to the filter gain.

% Validate number of input arguments and input parameters
narginchk(2, 4);
nargoutchk(2, 3);

validateattributes(n,{'numeric'},{'scalar','integer','positive'},'autoButter','n');
validateattributes(wn,{'numeric'},{'vector','real','positive','finite'},'autoButter','wn');
switch nargin
    case 2  % Use a default value for the filter type if not specified
        if length(Ws) == 1
            varargin{1} = 'low';
        else
            varargin{1} = 'bandpass';
        end
    case 4 % Validate sample frequency
        validateattributes(varargin{1},{'numeric'},{'scalar','real','positive','finite'},'autoButter','fs');
        % Compute normalized cutoff frequency wn from fc and fs
        wn = 2 * wn / varargin{1};
end

% Validate filter type - accept also partial strings
ftype = validatestring(varargin{end},{'low','high','bandpass','stop'},'autoButter','ftype');

% Show warning if the filter order is too high
if n > 15
    warning(['Use of autoButter for filter orders higher than 15 is ' ...
             'not recommended due to limited accuracy'])
end

% Construct the n-th order analog lowpass prototype
V = tan(0.5 * pi * wn);

% Temporary poles located  on the unit circle in the left-half plane
Q = exp((0.5i * pi / n) * ((2 + n - 1):2:(3 * n - 1)));
nQ = size(Q, 2);

% Transform the analog prototype to state-space
switch ftype
   case 'stop'
      Sg = 1 / prod(-Q);
      c  = -V(1) * V(2);
      b  = (V(2) - V(1)) * 0.5 ./ Q;
      d  = sqrt(b .* b + c);
      Sp = [b + d, b - d];
      Sz = sqrt(c) * (-1) .^ (0:2 * nQ - 1);
   case 'bandpass'
      Sg = (V(2) - V(1)) ^ nQ;
      b  = (V(2) - V(1)) * 0.5 * Q;
      d  = sqrt(b .* b - V(1) * V(2));
      Sp = [b + d, b - d];
      Sz = zeros(1, nQ);
   case 'high'
      Sg = 1 ./ prod(-Q);
      Sp = V ./ Q;
      Sz = zeros(1, nQ);
   case 'low'
      Sg = V ^ nQ;
      Sp = V * Q;
      Sz = [];
end

% Bilinear transform
P = (1 + Sp) ./ (1 - Sp);
Z = repmat(-1, size(P));
if isempty(Sz)
   G = real(Sg / prod(1 - Sp));
else
   G = real(Sg * prod(1 - Sz) / prod(1 - Sp));
   Z(1:length(Sz)) = (1 + Sz) ./ (1 - Sz);
end

% From zeros, ples and gain to b (numerator) and a (denominator):
if nargout == 2
   Z = G * real(poly(Z'));  % b
   P = real(poly(P));       % a
end
