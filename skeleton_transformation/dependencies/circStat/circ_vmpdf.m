function [p alpha] = circ_vmpdf(alpha, thetahat, kappa)

% [p alpha] = circ_vmpdf(alpha, w, p)
%   Computes the circular von Mises pdf with preferred direction thetahat 
%   and concentration kappa at each of the angles in alpha
%
%   The vmpdf is given by f(phi) =
%   (1/(2pi*I0(kappa))*exp(kappa*cos(phi-thetahat)
%
%   Input:
%     alpha     angles to evaluate pdf at, if empty alphas are chosen to
%               100 uniformly spaced points around the circle
%     [thetahat preferred direction, default is 0]
%     [kappa    width, default is 1]
%
%   Output:
%     p         von Mises pdf evaluated at alpha
%     alpha     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, Fisher
%
% Circular Statistics Toolbox for Matlab

% By Marc J. Velasco, 2009
% velasco@ccs.fau.edu
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/

if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

%intervals on the circle
dd = abs(circ_dist(alpha, [alpha(2:end); alpha(1)]));

C = 1/(2*pi*bessi0(kappa));
p = C*exp(kappa*cos(alpha-thetahat));
%integrate
p = p.*dd;
