function win = autoHamming(n, flag)
%AUTOHANN Generates a Hamming window
%
% Copyright (c) 2022, Lubos Smolik, Jan Rendl
% v1.0.0 (build 19. 7. 2022)  
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
    flag = validatestring(flag,{'symmetric','periodic'},'autoHamming','flag');
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
    win = 0.54 - 0.46 * cos((2 * pi / n) * transpose(0:n));
    win = win(1:end-ind);
end

% End of main function
end