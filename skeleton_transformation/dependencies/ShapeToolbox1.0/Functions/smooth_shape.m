function smoothed_shape = smooth_shape(shape,mask,mask_type)
% smoothed_shape = smooth_shape(shape,mask,mask_type)
% 
% Smooth shape using mask of size 'mask'. 
% mask_type = 'box','gaussian' determines shape of mask.
%
% With 'box', mask is size of box.
% With 'gaussian', mask is standard deviation, and the neighborhood used is three times as big.
%   This means that points beyond 3sd are ignored.
%
% Default is 'box'.
%................
if nargin == 2
	mask_type = 'box';
end;
switch mask_type
	case 'box' 
		neighborhood_size = mask;
	case 'gaussian' 
		if mask == 1 % Special case: if requested mask is 1, by convention this means "no smoothing", regardless.
			neighborhood_size = 1;
		else % normal case:
			neighborhood_size = 3 * mask;
		end;
end;
%................
smoothed_shape = [];
n = size(shape,1);
for i=1:n
    neighborhood = [];
    for c = [1:neighborhood_size] - round(neighborhood_size/2)  
        if i + c < 1
            ii = i + c + n;
        elseif i + c > n
            ii = i + c - n;
        else
            ii = i + c;
        end;
        neighborhood = [neighborhood; shape(ii,:)];  %collects points in neighborhood
    end;
    %neighborhood is now set of points to average over
	% ... Make mask: ............
	switch mask_type
		case 'box'
			mask_filter  = (1/neighborhood_size) * ones(1,neighborhood_size);  % equal values at each point, summing to 1
		case 'gaussian'
			if mask == 1
				mask_filter = [1];
			else
				mask_filter  = normpdf(-round(neighborhood_size/2-1):round(neighborhood_size/2)-1,0,mask); % a guassian mask, centered at middle, sd = mask, summing to 1
				mask_filter = (1/sum(mask_filter)) * mask_filter; % normalize mask filter to sum to 1
			end;
			% If mask size is even, we need to add a zero at the beginning (left) to make the mask_filter the right size
			if isequal(neighborhood_size/2,round(neighborhood_size/2))
				mask_filter = [0 mask_filter];
			end;
		otherwise
			error('unknown mask type [smooth_shape]');
	end;
	new_point_x = dot(mask_filter,neighborhood(:,1));
	new_point_y = dot(mask_filter,neighborhood(:,2));
	smoothed_shape = [smoothed_shape; new_point_x new_point_y];
end;