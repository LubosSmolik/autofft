function win = autoUniform(n)
%AUTOUNIFORM Generates a uniform (rectangular) window
%
% Copyright (c) 2022-2024, Lubos Smolik, Jan Rendl
% v1.0.1 (build 12. 9. 2024)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoUniform(n)
%
% win = autoUniform(n) returns an N-point uniform window in a column vector.

% CHANGELOG
% v1.0.1 - Function description has been corrected

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