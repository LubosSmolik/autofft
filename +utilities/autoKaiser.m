function win = autoKaiser(n, beta)
%AUTOHANN Generates a Kaiser window
%
% Copyright (c) 2022, Lubos Smolik, Jan Rendl
% v1.0.0 (build 19. 7. 2022)  
%
% This code is published under BSD-3-Clause License.
%
% win = autoKaiser(n)
% win = autoKaiser(n, beta)
%
% win = autoKaiser(n) returns an N-point symmetric Kaiser window with shape
%   parameter beta = 0.5. The window is returned as a column vector.
%
% win = autoKaiser(n, beta) returns an N-point Kaiser window with shape
%   parameter beta. The window is returned as a column vector.

% Validate number of inputs and outputs
narginchk(1,2);
nargoutchk(0,1);

% Validate n and check if the window is trivial
[n, win, istrivial] = utilities.validateN(n);
if istrivial
    return;
end

% Set beta if not specified by user
if nargin == 1
    beta = 0.5;
else
    validateattributes(beta,{'numeric'},{'scalar','real','nonnegative','finite'},'autoKaiser','beta');
    % Covert beta to double to enforce precision
    beta = double(beta);
end

% Generate a sequence of N + 1 samples
win = besseli(0, beta*sqrt(1-((transpose(0:n-1)-(n-1)/2)/((n-1)/2)).^2))...
                 / besseli(0,beta);

% End of main function
end