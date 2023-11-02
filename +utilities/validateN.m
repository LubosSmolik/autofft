function [nout, win, istrivial] = validateN(nin)
%VALIDATEN Checks input is integer and generates a trivial window if needed
%
% Copyright (c) 2022, Lubos Smolik, Jan Rendl
% v1.0.0 (build 19. 7. 2022)  
%
% This code is published under BSD-3-Clause License.
%
% [nout, win, istrivial] = validateN(nin) rounds nin to the nearest integer
%   if necessary. Note that nout is reutrned with the double precision. In
%   special cases (nin is [], 0, 1) validateN returns a trivial window win
%   and set istrivial to true.

% Check if nin is empty
if isempty(nin)
    nout = 0;
    win  = zeros(nout, 1);
    istrivial = true;
    return;    
else
    % Validate nin
    validateattributes(nin,{'numeric'},{'scalar','finite','real','nonnegative'},'validateN','n');
end

% Check if nin is integer
if nin == floor(nin)
    % Convert to double to enforce precision
    nout = double(nin);
else
    % Round and convert to double to enforce precision
    nout = double(round(nin));
    warning('N has been rounded to the nearest integer');
end

% Check if nin is 0 or 1
if nin == 0
    win = zeros(nout, 1); % Empty matrix: 0-by-1
    istrivial = true;
elseif nin == 1
    win = 1;
    istrivial = true;
else
    win = 0;
    istrivial = false;
end
% End of main function
end