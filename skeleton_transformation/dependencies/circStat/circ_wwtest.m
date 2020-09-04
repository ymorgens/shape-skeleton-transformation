function [pval mu F] = circ_wwtest(varargin)
%
% [pval, F] = circ_wwtest(alpha, idx, [w])
% [pval, F] = circ_wwtest(alpha1, alpha2, [w1, w2])
%   Parametric Watson-Williams multi-sample test for equal means. Can be
%   used as a one-way ANOVA test for circular data.
%
%   H0: the s populations have equal means
%   HA: the s populations have unequal means
%
%   Note: 
%   Use with binned data is only advisable if binning is finer than 10 deg.
%   In this case, alpha is assumed to correspond
%   to bin centers.
%
%   The Watson-Williams two-sample test assumes underlying von-Mises
%   distributrions. All groups are assumed to have a common concentration
%   parameter k.
%
%   Input:
%     alpha   angles in radians
%     idx     indicates which population the respective angle in alpha
%             comes from, 1:s
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     pval    p-value of the Watson-Williams multi-sample test. Discard H0 if
%             pval is small.
%     mu      best estimate of shared population mean if H0 is not
%             discarded at the 0.05 level and NaN otherwise.
%     F       test statistic of the Watson-Williams test.
%
%
% PHB 3/19/2009
%
% References:
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/

[alpha, idx, w] = processInput(varargin{:});

% number of groups
u = unique(idx);
s = length(u);

% number of samples
n = length(idx);

% compute relevant quantitites
pn = zeros(s,1); pr = pn;
for t=1:s
  pidx = idx == u(t);
  pn(t) = sum(pidx);
  pr(t) = circ_r(alpha(pidx),w(pidx));
end

r = circ_r(alpha,w);
rw = sum(pn.*pr)/n;

% make sure assumptions are satisfied
checkAssumption(rw,mean(pn))

% test statistic
kk = circ_r2kappa(rw);
beta = 1+3/(8*kk);    % correction factor
A = sum(pr.*pn) - r*n;
B = n - sum(pr.*pn);

F = beta * (n-s) * A / (s-1) / B;
pval = 1-fcdf(F,s-1,n-s);

% compute estimate of population mean
if pval > 0.05
  mu = circ_mean(alpha,w);
else 
  mu = NaN;
end



function checkAssumption(rw,n)

if n > 10 && rw<.45
  warning('Test not applicable. Average resultant vector length < 0.45.') %#ok<WNTAG>
elseif n > 6 && rw<.5
  warning('Test not applicable. Average number of samples per population < 11 and average resultant vector length < 0.5.') %#ok<WNTAG>
elseif n >=5 && rw<.55
  warning('Test not applicable. Average number of samples per population < 7 and average resultant vector length < 0.55.') %#ok<WNTAG>
elseif n < 5
  warning('Test not applicable. Average number of samples per population < 5.') %#ok<WNTAG>
end


function [alpha, idx, w] = processInput(varargin)

if nargin == 4
  alpha1 = varargin{1}(:);
  alpha2 = varargin{2}(:);
  w1 = varargin{3}(:);
  w2 = varargin{4}(:);
  alpha = [alpha1; alpha2];
  idx = [ones(size(alpha1)); ones(size(alpha2))];
  w = [w1; w2];
elseif nargin==2 && sum(abs(round(varargin{2})-varargin{2}))>1e-5
  alpha1 = varargin{1}(:);
  alpha2 = varargin{2}(:);
  alpha = [alpha1; alpha2];
  idx = [ones(size(alpha1)); 2*ones(size(alpha2))];
  w = ones(size(alpha));
elseif nargin==2
  alpha = varargin{1}(:);
  idx = varargin{2}(:);
  if ~(size(idx,1)==size(alpha,1))
    error('Input dimensions do not match.')
  end
  w = ones(size(alpha));  
elseif nargin==3
  alpha = varargin{1}(:);
  idx = varargin{2}(:);
  w = varargin{3}(:);
  if ~(size(idx,1)==size(alpha,1))
    error('Input dimensions do not match.')
  end
  if ~(size(w,1)==size(alpha,1))
    error('Input dimensions do not match.')
  end  
else
  error('Invalid use of circ_wwtest. Type help circ_wwtest.')
end

