function [r mu] = circ_axialmean(alphas, m);
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [m		axial correction (2,3,4,...)]
%
%   Output:
%     r		mean resultant length
%     mu		mean direction
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/


if nargin < 2
    m = 1;
end

zbarm = mean(exp(i*alphas*m));

r = abs(zbarm);
mu = angle(zbarm)/m;

