function win = autoFlatTop(n, flag)
%AUTOHANN Generates a flat top window
%
% Copyright (c) 2022, Lubos Smolik, Jan Rendl
% v1.0.0 (build 19. 7. 2022)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoFlatTop(n)
% win = autoFlatTop(n, flag)
%
% win = autoFlatTop(n) returns an N-point symmetric flat top window in a
%   column vector.
%
% win = autoFlatTop(n,flag) returns an N-point flat top window. The window
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
    flag = validatestring(flag,{'symmetric','periodic'},'autoFlatTop','flag');
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
    a0 = 0.21557895;
    a1 = 0.41663158;
    a2 = 0.277263158;
    a3 = 0.083578947;
    a4 = 0.006947368;
    
    x = (pi / n) * transpose(0:n);
    
    % Generate a sequence of N + 1 samples
    win =   a0 - a1 * cos(2 * x) + ...
            a2 * cos(4 * x) - ...
            a3 * cos(6 * x) + ...
            a4 * cos(8 * x);
    win = win(1:end-ind);
end

% End of main function
end