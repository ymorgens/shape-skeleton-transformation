function [d j len] = curvedist(c,crit);

% calculate length of curve by sequentially finding the distance between 
% neighbouring points (x,y) of a curve and then summing distances.
% d is the sum, 
% c is the contour, 

if nargin<2
    crit = NaN;
end

d = sum(sqrt(sum((c(1:end-1,:) - c(2:end,:)).^2,2)));

% sometimes will need to find point on contour after some length
if any(crit)
    if nargin<2
        crit = NaN;
    end
    
    distfun = @(x,y) sum(sqrt(sum((x-y).^2,2)));
    
    for j = 1:(size(c,1))-1;
        pt1 = c(j,:);
        % if j<size(c,1)
        pt2 = c(j+1,:);
        %     else
        %         pt2 = c(1,:);
        %     end
        len(j) = distfun(pt1, pt2);
        d = sum(len);
        
        
        if any(crit)
            if d>=crit| abs(d - crit)<0.0001
                j = j+1; % include pt2 in the total
                break
            end
            
        end
    end
end