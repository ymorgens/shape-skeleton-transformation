function [shape,skeleton] = random_axial_shape(n_axes,style)
% [shape,skeleton] = random_axial_shape(n_axes,style)
%
% Generates a shape from a randomly generated skeleton with n_axes.
%
% Procecdure is quite klugy. 

if nargin == 1 % no style specified
    style = 'nonparametric';
end;
if n_axes < 1
    error('Zero or fewer axes requested from random_axial_shape');
end;

%---------------------------------------------
global n_points_per_axis;
center = [50 50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch style % parametric or nonparametric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'nonparametric' % generate axes point by point
%---------------------------------------------
% First, make an initial axis
initial_orientation = 2 * pi *rand;
[Xs Ys] = random_curve(initial_orientation,center,n_points_per_axis);
root_axis = [Xs Ys];
skeleton(1).contour = root_axis;
skeleton(1).parent = -1;
skeleton(1).index = 1;
%draw_skeleton(skeleton);
%---------------------------------------------
% Next, make additional axes off the root
for a= 2:n_axes
    perm = randperm(a-1);
    theparent = perm(1); % random number between 1 and a-1 
    skeleton(a).parent = theparent;
    skeleton(a).index = a;
    parentaxis = skeleton(theparent).contour;
    parentnormals = jfnormals(parentaxis(:,1),parentaxis(:,2));
    perm = randperm(n_points_per_axis);
    start_index = perm(1);
    start_location = parentaxis(start_index,:);
    d = rand; sign = (d-0.5)/abs(d-0.5); % set random sign (kluge)
    start_vector = sign* parentnormals(start_index,:);
    %....
    start_direction = atan2(start_vector(2),start_vector(1));
    [Xs Ys] = random_curve(start_direction,start_location,n_points_per_axis);
    skeleton(a).contour = [Xs Ys];
end;
%skeleton(2).contour
%---------------------------------------------
widths = 5 + 10 * rand(n_axes,1);
shape = smooth_shape(skeleton2shape(skeleton,'constant',widths),3); % Smooth in order to minimize self-intersections at curves. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'parametric' % parameters for spline construction of axes
        % not implemented  yet
end;




          
          
