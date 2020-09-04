function [h U UC] = circ_raotest(alpha, sig)

% [h U UC] = circ_raotest(alpha, [sig])
%   Calculates Rao's spacing test by comparing distances between points on
%   a circle to those expected from a uniform distribution.
%
%   H0: Data is distributed uniformly around the circle. 
%   H1: Data is not uniformly distributed around the circle.
%
%   Alternative to the Rayleigh test and the Omnibus test. Less powerful
%   than the Rayleigh test when the distribution is unimodal on a global
%   scale but uniform locally.
%
%   Input:
%     alpha     sample of angles
%     [sig      significance level, can be any of 0.01, 0.05 and 0.1, 
%               default is 0.05]
%
%   Output:
%     h         result of hypothesis test at sig-level
%     U         computed value of the test-statistic u
%     UC        critical value of the test statistic at sig-level
%
%
%   References:
%     Batschelet, 1981, Sec 4.6
%
% Circular Statistics Toolbox for Matlab

% By Marc J. Velasco, 2009
% velasco@ccs.fau.edu
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/


if nargin < 2
   sig = 0.05; 
end

if ~ismember(sig, [.01 .05 .10])
   error('Significance level chosen is not supported: sig is restricted to (.01, .05, .10).');
end

% for the purpose of the test, convert to angles
alpha = rad2ang(alpha);
n = length(alpha);
alpha = sort(alpha);

% compute test statistic
U = 0;
lambda = 360/n;
for j = 1:n-1
    ti = alpha(j+1) - alpha(j);
    U = U + abs(ti - lambda);
end

tn = (360 - alpha(n) + alpha(1));
U = U + abs(tn-lambda);

U = (1/2)*U;

% get critical value from table
table = gettable();

nn = table(:,1);
[easy row] = ismember(n, nn);
if ~easy
   % find closest value if no entry is present)
   row = length(nn) - sum(n<nn); 
   if row == 0
       error('N too small.');
   else
      warning('N=%d not found in table, using closest N=%d present.',n,nn(row)) %#ok<WNTAG>
   end
end

% find column corresponding to significance level
if sig == .01, col = 2; 
elseif sig == .10, col = 4; 
else col = 3;
end
    
UC = table(row, col);

% compare test statistic to critical value
h = U > UC;

function table = gettable()

% table L from Batschelet, 1981, p. 339
% alpha = .01, .05, .10

table = [...
    4   221.0   186.5   171.7  ; ...
    5   212.0   183.6   168.8  ; ...
    6   206.0   180.7   166.3  ; ...
    7   202.7   177.8   164.9  ; ...
    8   198.4   175.7   163.4  ; ...
    9   195.1   173.5   162.4  ; ...
        ...
    10  192.2   172.1   161.3  ; ...
    11  189.7   170.3   160.2  ; ...
    12  187.6   169.2   159.2  ; ...
    13  185.8   167.8   158.4  ; ...
    14  184.0   166.7   157.7  ; ...
        ...
    15  182.2   165.6   157.0  ; ...
    16  180.7   164.9   156.6  ; ...
    17  179.6   164.2   155.9  ; ...
    18  178.2   163.1   155.2  ; ...
    19  177.1   162.4   154.8  ; ...
        ...
    20  176.0   161.9   154.4  ; ...
    25  171.9   158.9   152.7  ; ...
    30  168.8   156.7   151.4  ; ...
    35  166.4   155.0   150.3  ; ...
    40  164.4   153.6   149.5  ; ...
    45  162.7   152.4   148.7  ; ...
    50  161.2   151.4   148.1  ; ...
        ...
    100 152.8   146.8   143.7  ; ...
    200 146.8   142.6   140.4  ; ...
    
    ];

