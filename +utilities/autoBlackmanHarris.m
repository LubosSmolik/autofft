function win = autoBlackmanHarris(n, flag)
%AUTOBLACKMANHARRIS Generates a four-term Blackmann-Harris window
%
% Copyright (c) 2022-2025, Lubos Smolik, Jan Rendl
% v1.0.2 (build 12. 8. 2025)
%
% This code is published under BSD-3-Clause License.
%
% win = autoBlackmanHarris(n)
% win = autoBlackmanHarris(n, flag)
%
% win = autoBlackmanHarris(n) returns an N-point symmetric Blackmann-Harris
%   window in a column vector.
%
% win = autoBlackmanHarris(n, flag) returns an N-point Blackmann-Harris
%   window. The window can be either symmetric if the flag equals to
%   'symmetric' or periodic if it equals to 'periodic'.

% CHANGELOG
% v1.0.2 - Multiplications of vectors by scalars now use element-wise
%          operators, i.e. '.*'
%        - Fixed erroneous description.
% v1.0.1 - Input validation has been improved
%        - Function description has been corrected

% Validate number of inputs and outputs
narginchk(1,2);
nargoutchk(0,1);

% Validate n and check if the window is trivial
[n, win, istrivial] = utilities.validateN(n);
if istrivial
    return;
end

% Set flag to symmetric if not specified and validate it
if nargin == 1
    flag = 'symmetric';
else
    flag = validatestring(flag, {'symmetric','periodic'}, '', 'flag', 2);
end

% Generate the window
switch flag
    case 'periodic'
        win = generateWindow(n, 1);        
    case 'symmetric'
        win = generateWindow(n - 1, 0);
end

%GENERATEWINDOW Calculates the window function
function win = generateWindow(n, ind)
    % Define window coefficients
    a0 = 0.35875;
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;
    
    x = (pi / n) .* transpose(0:n);
    
    % Generate a sequence of N + 1 samples
    win = a0 - a1 .* cos(2 .* x) + a2 .* cos(4 .* x) - a3 .* cos(6 .* x);
    win = win(1:end-ind);
end

% End of main function
end