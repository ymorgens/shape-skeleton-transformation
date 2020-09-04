function [fixed_shape,reversed] = insideright_shape(shape)
% [fixed_shape,reversed] = insideright_shape(shape)
%
% Correct the figural parity of shape so that its normals point OUT.
%
% Accomplishes this by summing angle clockwise around the
% shape. If the sum is 2pi, it's clockwise; -2pi, counterclockwise. 
% In practice numerical errors are substantial, so we simply test if the 
% sum is positive or negative. 
%
% Note that if the shape has a loop, this won't work, because topologically
% it is now equivalent to a loop of 2pi and another of -2pi -- no matter
% how small one of the loops is. In practice this happens sometimes when
% the shape has been "fractalized". We minimize this by smoothing first,
% but it can always happen if the random loops is larger than the smoothing
% mask. 

original_shape = shape;
shape = smooth_shape(shape,round(size(shape,1)/4));  % Use very smoothed version of shape.
   % If the shape is very noisy, it can contain loops. Loops completely
   % reverse the direction of the angular sum, so for example 1 loop leads
   % to a toal angular sum of approximately zero (whose sign is thus noise
   % based on rounding errors). So we smooth to hopefully eliminate loops.
   % Based on the analysis o the smoothed version, we either reverse the
   % (original) shape or don't. 
sum = 0;
if isequal(shape(1,:),shape(end,:))  % if first and last points are the same
    shape = shape(1:end-1,:);        % delete the last point
end;
% Augment the shape to include the first two points at the end again:
shape = [shape; shape(1,:); shape(2,:)];
n = size(shape,1);
for i=3:n
    p1 = shape(i-2,:);
    p2 = shape(i-1,:);
    p3 = shape(i,:);
    v1 = p2-p1;
    v2 = p3-p2;
    turning_angle =  angle_from_to(v1,v2);
    sum = sum + turning_angle;
end;
if sum > 0   % backwards
	fixed_shape = original_shape(end:-1:1,:);
    reversed = true;  % flag indicating we had to reverse the shape
else
	fixed_shape = original_shape; % Ok in the first place
    reversed = false;
end;

%------------------------------------------------------
function turning_angle = angle_from_to(v1,v2)
% Compute the angle from vector 1 to vector 2, including correct sign.
% Positive if clockwise. 

magnitude = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
if ~isreal(magnitude) % happens if argument is > 1 -- numerical precison problem
    magnitude = 1;
end;
cp = cross([v1,0],[v2,0]);
turning_angle = sign(cp(3))*magnitude; % alpha