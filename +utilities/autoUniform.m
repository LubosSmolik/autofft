function win = autoUniform(n)
%AUTOHANN Generates a uniform (rectangular) window
%
% Copyright (c) 2022, Lubos Smolik, Jan Rendl
% v1.0.0 (build 19. 7. 2022)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoUniform(n)
%
% win = autoUniform(n) returns an N-point uniform window in a column vector.

% Validate number of inputs and outputs
narginchk(1,1);
nargoutchk(0,1);

% Validate n and check if the window is trivial
[n, win, istrivial] = utilities.validateN(n);
if istrivial
    return;
end

% Generate the window
win = ones(n, 1);

% End of main function
end