function [skeleton,description_length] = parameters2skeleton(parameters,axis_hierarchy,use_spacing_to_match_shape)
% [skeleton,description_length] = parameters2skeleton(parameters,axis_hierarchy,use_spacing_to_match_shape)
%
% Reconstructs a skeleton from parameters, given axis_hierarchy (in
% practice often stored globally).
% 
% This version has arbitrary numbers of parameters (and knotpoints) per axis, so a lot of key information is stored in the axis_hierarchy.

% Each row of axis_hierarchy relates to one axis, containing:
% [index  parent_index parent_knotpoint first_paramter last_parameter]

% NOTE The parmaeters for each axis refer to the locations of knotpoints 2:end, wrt the first. 
% The whole thing is then translated to the location of rootpoint, assumped equal to the first knot. 
% (Note thought there are situations where these are not equal -- NYI). 

global shape_spacing;
if nargin == 2 %if use_spacing_to_match_shape is NOT set : default
	use_spacing_to_match_shape = true; % Using this causes the number of points per contour to change over optimization, so it causes problems
end;
%----------------------------------------------------------
global current_axis_hierarchy;
if nargin == 1
    axis_hierarchy = current_axis_hierarchy;
	use_spacing_to_match_shape = true;
end;
%----------------------------------------------------------
n_axes = size(axis_hierarchy,1);
firstrootpoint = [parameters(1) * cos(parameters(2)) parameters(1) * sin(parameters(2))];
%----------------------------------------------------------
% The structure of the following procedure is all based around the fact that parents have to be made before their children,
% since each axis starts from a rootpoint determined wrt its parent.
skeleton_in_progress= [];
for depth = 1:length(unique(axis_hierarchy(:,1)))  % number of distinct parents-- ie max depth
    axes_this_depth = find_axes_this_depth(axis_hierarchy,depth); % find axes that have depth = depth; 
    for i=1:length(axes_this_depth)
        this_axis = axes_this_depth(i);
        this_axis_parent =                  axis_hierarchy(this_axis,2);
        this_axis_parent_knotpoint =        axis_hierarchy(this_axis,3);
        parameter_start =                   axis_hierarchy(this_axis,4);
        parameter_end =                     axis_hierarchy(this_axis,5);
        parameters_this_axis = parameters(parameter_start:parameter_end);
        %... Compute rootpoint:
        % Note that parent axis is assumed made already. Its knotpoints are stored in skeleton(parent).knotpoints.
        if this_axis_parent == -1 % root
            rootpoint = firstrootpoint; 
			previous_parent_knotpoint = [0 0];
        else % normal case
            parent_knotpoints = skeleton_in_progress(this_axis_parent).knotpoints;
            rootpoint = parent_knotpoints(this_axis_parent_knotpoint,:);
			if this_axis_parent_knotpoint > 1
				previous_parent_knotpoint = parent_knotpoints(this_axis_parent_knotpoint - 1,:);
			else % special case where parent knotpoint is FIRST on that axis; take one step "back" to make imaginary previous knot
				first_parent_knot  = parent_knotpoints(1,:);
				second_parent_knot = parent_knotpoints(2,:);
				previous_parent_knotpoint = first_parent_knot + first_parent_knot - second_parent_knot; % this matches what is done in the opposite direction
			end;
        end; % now rootpoint and previous_parent_knotpoint are set
        % Next, determine starting angle of branch
        if this_axis_parent == -1 % root
            starting_angle = parameters(2);
        else % not root
			starting_angle = atan2(rootpoint(2) - previous_parent_knotpoint(2), rootpoint(1) - previous_parent_knotpoint(1));
        end;
        [axis_contour,knotpoints] = parameters2axis(rootpoint,previous_parent_knotpoint,parameters_this_axis,starting_angle);
        %-----------------------------------------------------------------
        % "Nudge" axis_contour so that its rootpoint matches nearest point
        % on parent. 
        
        adjustment_method = 1; % nudge whole contour
        
        if this_axis_parent > 0
            ds = [];    
            n = size(axis_contour,1);
            parent_contour = skeleton_in_progress(this_axis_parent).contour;  % Not going to change this. 
            % Determine which point on the PARENT is closest to the root. 
            for jj=1:n
                ds = [ds norm(axis_contour(1,:) - parent_contour(jj,:))];
            end;
            [min_d,min_index] = min(ds);  % min_index is the index of the point on parent we want to nudge TOWARDS
            
            if adjustment_method == 1
                % New version. Instead of nudging rootpoint, we nudge the
                % WHOLE CONTOUR over so rootpoint matches. 
                adjustment = parent_contour(min_index,:) - axis_contour(1,:); 
                % Add adjustment to all points in contour
                for i=1:size(axis_contour,1)
                    axis_contour_(i,:) = axis_contour(i,:) + adjustment;
                end;
            else
            % Nudge rootpoint
                axis_contour(1,:) = parent_contour(min_index,:); 
            end;
        else
            axis_contour_ = axis_contour;
        end;
            %-----------------------------------------------------------------
        % Establish skeleton elements:
        skeleton_in_progress(this_axis).contour = axis_contour_;
        skeleton_in_progress(this_axis).index  = this_axis;
        skeleton_in_progress(this_axis).parent = this_axis_parent;
        % Finally, add supplental fields that are used by children axes in being constructed:
        skeleton_in_progress(this_axis).knotpoints = knotpoints; 
          % keeps a record of which vertices are knotpoints in the construction of this axis... for the construction of descendents.    
         

    end;
end;

use_spacing_to_match_shape = false;  
%--------------------------------------------------------------
% Next, finespace skeleton if use_spacing_to_match_shape is turned on.
skeleton = skeleton_in_progress; % Note this has some extra stuff in it, but that shouldn't cause a problem
if use_spacing_to_match_shape  % Flag that forces axis sampling to match shape sampling
	for a=1:length(skeleton) 	% finespace to match shape
		this_contour = skeleton(a).contour;
        first_point = this_contour(1,:);
        [x,y] = finespacedcontour(this_contour(:,1),this_contour(:,2),shape_spacing,0);
		evenlyspaced_version  = [x,y]; % this is guaranteed to include the original first point, which here is the one shared with the parent axis
		% Now we have to make sure the parent also contains the first branch point.
		if skeleton(a).parent > 0
			this_parent_contour = skeleton(skeleton(a).parent).contour;
			% Now find closest point in this_parent_contour and "nudge" it to match first point in branch
            %rootpoint might be first point, might be last point!
            first_point = this_contour(1,:);
            last_point  = this_contour(end,:);
            ds_to_first = [];
            ds_to_last  = [];
            for i=1:size(this_parent_contour,1)
				ds_to_first = [ds_to_first norm(this_parent_contour(i,:) - first_point)];
                ds_to_last = [ds_to_last norm(this_parent_contour(i,:) - last_point)];
			end;
            [min_d_first winning_index_first] = min(ds_to_first); % winning index is the point closest to first branch point
            [min_d_last winning_index_last] = min(ds_to_last); % winning index is the point closest to first branch point
            if min_d_last <  min_d_first % last point is root point
                this_contour = evenlyspaced_version(end:-1:1,:);  % invert contour, replacing it with evenly spaced version
                this_parent_contour(winning_index_last,:) = this_contour(1,:);  % Adjust closest point in parenet contour to match rootpoint
            else
                this_contour = evenlyspaced_version;              % replacie contour with evenly spaced version
                this_parent_contour(winning_index_first,:) = this_contour(1,:); % Adjust "
            end;
            skeleton(skeleton(a).parent).contour = this_parent_contour; % Replace this parent contour, with "nudged" point
                 
                
		end;
	end;
end;

%----------------------------------------------------------
function [axis_contour,knots] = parameters2axis(rootpoint,previous_parent_knotpoint,parameters_this_axis,starting_cum_angle)
% Constructs an axis_contour from a bunch of parameters; also returns the actual knotpoints from which the spline is derived.
% Generally, parameters refer to the bend & logstretch of the 2:end
% rootpoint is the point at which this axis starts (ie its first knot).
% previous_parent_knotpoint is the previous knot on the parent (or origin in the case of the root), with respect to which orientations are computed.
% Special case: when the knotpoint on the parent is the FIRST on the parent, the "previous knotpoint" is an imaginary knotpoint one step backwards
% from the first (distance and direction equal to the vetor from the first to the second).
global n_points_per_axis;   % set in shape_tool
global shape_spacing;
axis_contour = [];
n_knotpoints = (length(parameters_this_axis)/2) + 1;  % there are 2 parameters specifying each knotpoint from 2:end, plus the first
knots = zeros(n_knotpoints,2);
knots(1,:) = rootpoint;
cum_angle = starting_cum_angle;
previous_knot_x = knots(1,1);
previous_knot_y = knots(1,2);
previous_leg = norm(rootpoint - previous_parent_knotpoint);
for i=2:n_knotpoints
    this_bend =         parameters_this_axis(2*(i-1)-1);
    this_logstretch =   parameters_this_axis(2*(i-1));
    cum_angle = cum_angle + this_bend;
    this_leg = previous_leg * exp(this_logstretch);
    this_knot_x  = previous_knot_x + this_leg * cos(cum_angle);
    this_knot_y  = previous_knot_y + this_leg * sin(cum_angle);
    previous_knot_x = this_knot_x;
    previous_knot_y = this_knot_y;
    previous_leg = this_leg;
    knots(i,:) = [this_knot_x this_knot_y];
    %scatter(this_knot_x,this_knot_y);
end;
% Next, spline the knots. (This should be a standard function -- right now it's just the following lines, which are repeated in various places.)
t = 1:n_knotpoints;
tt = [1:(n_knotpoints-1)/(n_points_per_axis-1):n_knotpoints]; 
xx = spline(t,knots(:,1),tt);
yy = spline(t,knots(:,2),tt);
axis_contour = [xx' yy'];


%-------------------
function axes = find_axes_this_depth(axis_hierarchy,depth)
% find those axes at the given depth
axes = [];
for i=1:size(axis_hierarchy,1) % n axes
    if depth_of_axis(axis_hierarchy,i)==depth 
        axes = [axes i];
    end;
end;
%----------------------------------------------------------------
function depth = depth_of_axis(axis_hierarchy,axis)  % recursive
    if parent_of_axis(axis_hierarchy,axis) == -1
        depth = 1;
    else
        depth = depth_of_axis(axis_hierarchy,parent_of_axis(axis_hierarchy,axis)) + 1;
    end;
%----------------------------------------------------------------
function parentindex = parent_of_axis(axis_hierarchy,axis) 
    parentindex = axis_hierarchy(find(axis_hierarchy(:,1)==axis),2);
        
%----------------------------------------------------------------
function angle = find_direction_of_axis_at_point(parent_contour,point)
n_points = size(parent_contour,1);
distances = [];
for i=1:n_points
    this_point = parent_contour(i,:);
    this_distance = norm(this_point - point);
	distances = [distances this_distance];
	[winning_distance,winning_index] = min(distances);
end; % Now we now which index on parent_contour point is at
this_point = parent_contour(winning_index,:);
if winning_index == 1 % first point is a special case
    second_point = parent_contour(2,:);
    direction_on_parent = atan2(second_point(2) - this_point(2), second_point(1) - this_point(1));
elseif winning_index == size(parent_contour,1) % last point, also a special case
	second_to_last_point = parent_contour(end - 1,:);
	direction_on_parent = atan2(this_point(2) - second_to_last_point(2), this_point(1) - second_to_last_point(1));
else
    % Normal case-- we use previous and next point
    previous_point = parent_contour(winning_index-1,:);
	next_point = parent_contour(winning_index+1,:);
   direction_on_parent = atan2(next_point(2) - previous_point(2), next_point(1) - previous_point(1));
end;
angle = direction_on_parent - (pi/2);  % so it points to the right relative to parent

