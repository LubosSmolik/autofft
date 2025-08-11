function win = autoKaiser(n, beta, flag)
%AUTOKAISER Generates a Kaiser window
%
% Copyright (c) 2022-2025, Lubos Smolik, Jan Rendl
% v1.1.0 (build 11. 8. 2025)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoKaiser(n)
% win = autoKaiser(n, beta)
% win = autoKaiser(n, beta, flag)
%
% win = autoKaiser(n) returns an N-point symmetric Kaiser window with shape
%   parameter beta = 0.5. The window is returned as a column vector.
%
% win = autoKaiser(n, beta) returns an N-point symmetric Kaiser window with
%   shape parameter beta. The window is returned as a column vector.
%
% win = autoKaiser(n, beta, flag) returns an N-point Kaiser window with
%   shape parameter beta. The window can be either symmetric if the flag
%   equals to 'symmetric' or periodic if it equals to 'periodic'.


% CHANGELOG
% v1.1.0 - Added support for periodic Kaiser windows
%        - Computation of the Kaiser window refactored
% v1.0.1 - Input validation and performance has been improved
%        - Function description has been corrected

% Validate number of inputs and outputs
narginchk(1,3);
nargoutchk(0,1);

% Validate n and check if the window is trivial
[n, win, istrivial] = utilities.validateN(n);
if istrivial
    return;
end

% Set beta if not specified by user
if nargin < 2
    beta = 0.5;
else
    validateattributes(beta, 'numeric', {'scalar','real','nonnegative', ...
                       'finite'}, '', 'beta', 2);
    % Covert beta to double to enforce precision
    beta = double(beta);
end

% Set flag to symmetric if not specified and validate it
if nargin < 3
    flag = 'symmetric';
else
    flag = validatestring(flag, {'symmetric','periodic'}, '', 'flag', 3);
end

% Generate the window
switch flag
    case 'periodic'
        win = generateWindow(n + 1, 1);        
    case 'symmetric'
        win = generateWindow(n, 0);
end

% Calculate the window
function win = generateWindow(n, ind)
    % Generate a sequence of N samples
    x0 = (n-1) / 2;
    xi = sqrt(1 - ((transpose(0:n-1) - x0) ./ x0).^2);
    win = besseli(0, beta .* xi) ./ besseli(0, beta);
    win = win(1:end-ind);
end

% End of main function
end