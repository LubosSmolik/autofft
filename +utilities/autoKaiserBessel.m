function win = autoKaiserBessel(n, flag)
%AUTOKAISERBESSEL Generates a four-term Kaiser-Bessel window
%
% Copyright (c) 2025, Lubos Smolik
% v1.0.0 (build 12. 8. 2025)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoKaiserBessel(n)
% win = autoKaiserBessel(n, flag)
%
% win = autoKaiserBessel(n) returns an N-point symmetric Kaiser-Bessel
%   window in a column vector.
%
% win = autoKaiserBessel(n, flag) returns an N-point Kaiser-Bessel window.
%   The window can be either symmetric if the flag equals to 'symmetric'
%   or periodic if it equals to 'periodic'.

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
    a0 = 1;
    a1 = 1.24;
    a2 = 0.244;
    a3 = 0.00305;
    
    x = (pi / n) .* transpose(0:n);
    
    % Generate a sequence of N + 1 samples
    win = a0 - a1 * cos(2 * x) + a2 * cos(4 * x) - a3 * cos(6 * x);
    win = win(1:end-ind);
end

% End of main function
end