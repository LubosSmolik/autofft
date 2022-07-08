function [Z, P, G] = autoButter(n, Wn, varargin)
% AUTOBUTTER Butterworth-like digital filter design
%
%  Copyright (c) 2014              Jan Simon - original code
%  Copyright (c) 2022              Lubos Smolik - validation, revisions
% v1.1.0 (build 7. 7. 2022)        e-mail: carlist{at}ntis.zcu.cz
%
% This code is published under BSD-3-Clause License.
%
% Docs - TO DO

% Validate number of input arguments and input parameters
narginchk(2, 4)

validateattributes(n,{'numeric'},{'scalar','integer','positive'},'autoButter','N');
validateattributes(Wn,{'numeric'},{'vector','real','finite'},'autoButter','Wn');
switch nargin
    case 2  % Use a default value for the filter type if not specified
        varargin{1} = 'low';
    case 3  % Validate filter type
        validateattributes(varargin{1},{'char','string'},{'scalartext'},'autoButter','ftype');
    case 4 % Validate sample frequency and filter type
        validateattributes(varargin{1},{'numeric'},{'scalar','real','positive','finite'},'autoButter','fs');
        validateattributes(varargin{2},{'char','string'},{'scalartext'},'autoButter','ftype');
        % Compute normalized cutoff frequency Wn from fc and fs
        Wn = 2 * Wn / varargin{1};
end

% Show warning if the filter order is too high
if n > 15
    warning(['Use of autoButter for filter orders higher than 15 is ' ...
             'not recommended due to limited accuracy'])
end

% Construct the n-th order analog lowpass prototype
V = tan(Wn * 1.5707963267948966);

% Temporary poles located  on the unit circle in the left-half plane
Q = exp((1.5707963267948966i / n) * ((2 + n - 1):2:(3 * n - 1)));
nQ = size(Q, 2);

% Transform the analog prototype to state-space
switch lower(varargin{end})
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
   otherwise
      error('user:myButter:badFilter', 'Unknown filter type: %s', pass)
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

% From Zeros, Poles and Gain to B (numerator) and A (denominator):
if nargout == 2
   Z = G * real(poly(Z'));  % B
   P = real(poly(P));       % A
end