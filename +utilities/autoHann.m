function win = autoHann(n, flag)
%AUTOHANN Generates a Hann window
%
% Copyright (c) 2022-2025, Lubos Smolik, Jan Rendl
% v1.0.2 (build 12. 8. 2025)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoHann(n)
% win = autoHann(n, flag)
%
% win = autoHann(n) returns an N-point symmetric Hann window in a column
%   vector.
%
% win = autoHann(n,flag) returns an N-point Hann window. The window can be
%   either symmetric if the flag equals to 'symmetric' or periodic if it
%   equals to 'periodic'.

% CHANGELOG
% v1.0.2 - Multiplications of vectors by scalars now use element-wise
%          operators, i.e. '.*'
% v1.0.1 - Input validation has been improved

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
    % Generate a sequence of N + 1 samples
    win = 0.5 - 0.5 .* cos((2 * pi / n) .* transpose(0:n));
    win = win(1:end-ind);
end

% End of main function
end