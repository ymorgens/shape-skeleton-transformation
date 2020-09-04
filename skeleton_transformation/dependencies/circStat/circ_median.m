function md = circ_median(alpha)
%
% mu = circ_median(alpha, w)
%   Computes the median direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%
%   Output:
%     mu		mean direction
%
% PHB 3/19/2009
%
% References:
%   Biostatistical Analysis, J. H. Zar (26.6)
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/

if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end
alpha = mod(alpha,2*pi);
n = length(alpha);

done = false;
sz = ang2rad(1);
while ~done
  sz = sz/10;
  dg = 0:sz:pi;
  m = zeros(size(dg));
  for i=1:length(dg)
    m(i) = sum((alpha > dg(i) & alpha < pi + dg(i)));    
  end
  [mn idx] = min(abs(2*m-n));
  if (mn == 0 && mod(n,2)==0)  || ...
      (mn == 1 && mod(n,2)==1)
    done = true;
  else 
    continue
  end
  
  rawmd = dg(idx);
  
  if mod(n,2)==1
    [foo idx2] = min(abs(circ_dist(alpha,rawmd)));
    md = alpha(idx2);
  else
    [foo idx2] = sort(abs(circ_dist(alpha,rawmd)),1,'ascend');
    md = circ_mean(alpha(idx2(1:2)));
  end
  
  if circ_dist(circ_mean(alpha),md)>circ_dist(circ_mean(alpha),md+pi)
    md = mod(md+pi,2*pi);
  end
end
  