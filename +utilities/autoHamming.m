function win = autoHamming(n, flag)
%AUTOHAMMING Generates a Hamming window
%
% Copyright (c) 2022-2024, Lubos Smolik, Jan Rendl
% v1.0.2 (build 12. 9. 2024)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoHamming(n)
% win = autoHamming(n, flag)
%
% win = autoHamming(n) returns an N-point symmetric Hamming window in a
%   column vector.
%
% win = autoHamming(n,flag) returns an N-point Hamming window. The window
%   can be either symmetric using 'symmetric' flag or periodic using
%   'periodic' flag.

% CHANGELOG
% v1.0.2 - Coefficients a_0 and a_1 are now optimal in the equiripple sense
%        - Input validation has been improved
%        - Function description has been corrected
% v1.0.1 - Coefficients a_0 and a_1 now contain four significant digits
%          instead of two.


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

% Calculate the window 
function win = generateWindow(n, ind)
    % Generate a sequence of N + 1 samples
    win = 0.53836 - 0.46164 * cos((2 * pi / n) * transpose(0:n));
    win = win(1:end-ind);
end

% End of main function
end