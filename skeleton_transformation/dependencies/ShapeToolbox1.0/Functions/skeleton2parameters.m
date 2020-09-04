function [parameters,axis_hierarchy] = skeleton2parameters(skeleton,tolerance)
% [parameters,axis_hierarchy] = skeleton2parameters(skeleton,tolerance)
% 
% Given a skeleton, compute parameters that encode spline approximations of its axes. Also returns the axis_hierarchy, which contains all
% the "structural" information essential to reconstructing the skeleton, such as the parentage of each axis. Normally this is set globally
% in current_axis_hierarchy.
%
% Each row of axis_hierarchy relates to one axis, containing:
% [index  parent_index parent_knotpoint first_paramter last_parameter]
%
% parent_knotpoint indicates the number of the knotpoint on the parent (parent_index) at which this axis attaches. 
% The parameters for this axis are numbers first_parameter:last_parameter of parameters. Normally there 
% are 2*(k-1) parameters for each axis encoded by k knotpoints, and several (usually 2; (global) n_initial_parameters says how many)
% initial 'extra' parameters that relate to he entire skeleton, eg the position of the rootpoint of the 1st axis. 

% The two parameters for each knotpoint are:
%   bend logstretch

% Note: 
% Knot-point locations are at branch points, in order to ensure that when an initial skeleton is used,
% (eg a medial-axis estimate) the branches come out at the correct locations.
%------------------------------------------------
global current_axis_hierarchy;   % We're going to set this in this function, because the parameters we make can only be interepreted with this in place.
global current_axis_index;  % compute_parameters_of_axis uses these -- needs to know if we are on 1st axis, where on previous axis we were, etc.
global current_parentindex; % 
global current_parentvertex;
%------------------------------------------------
% Default tolerance
if nargin == 1
    tolerance = .1;
end;
%------------------------------------------------
% Make initial parameters:
% Find root:
for a = 1:length(skeleton)
	if skeleton(a).parent == -1
		rootindex = a;
	end;
end;
parameters = [norm(skeleton(rootindex).contour(1,:))  atan2(skeleton(rootindex).contour(1,2),skeleton(rootindex).contour(1,1))];  % These are the "special" first parameters, not connected with any one axis
%....... 
% Set up basic records:
knotpoint_locations = zeros(1,length(skeleton));  % this is going to be continually updated as we go through the axes.
%.... 
% Make axis parameters
for axis = 1:length(skeleton)
    parameter_counter = length(parameters); % ie, so far
    [skeleton,parameters_this_axis] = parameterize_this_axis(axis,skeleton,tolerance);
    parameter_start = parameter_counter + 1;
    parameter_end = parameter_counter + length(parameters_this_axis);
    skeleton(axis).parameter_locations = [parameter_start parameter_end];
    parameters = [parameters parameters_this_axis];
end; % Now, parameters and knotpoint_locations are complete.
% Finally, create axis hierarchy. This has all the "structual" information required to recreate the skeleton from the parameters.
axis_hierarchy = create_axis_hierarchy(skeleton);  % Note that skeleton now includes fields: knotpoints and parameter_locations
%=============================================================================================================================================================
% Local functions ============================================================================================================================================
function [skeleton,parameters] = parameterize_this_axis(axis,skeleton,tolerance)
% Return parameters for the axis-th axis of skeleton. 
% We place knotpoints at each child, and enough extra to fit the curve at tolerance less than some threshold ('tolerance'). 
% "Extra" knotpoints are placed at the max deviation from the previous attempt. 
% Also returns the skeleton, which is augmented by addition of the field "knotpoints" (ie vertices at which knot points have been chosen).
%tolerance = 10;
minimum_knot_separation = 2; % don't put knotpoints closer than this, otherwise artifacts ensue.
this_axis = skeleton(axis).contour;
% .. First, find children:
children_this_axis = [];
for potential_child = 1:length(skeleton)  % Collect children of this axis
    if skeleton(potential_child).parent == axis
        children_this_axis = [children_this_axis potential_child];
    end;
end;
% .. Next, determine the vertices on this axis where children attach (this may be first or last point on each child axis).
child_attachment_vertices = [];
for i=1:length(children_this_axis)
    this_child = children_this_axis(i);
    this_child_axis = skeleton(this_child).contour;
    ds = [];
    for v=1:size(this_axis,1)
        this_point = this_axis(v,:);
        d_top =    norm(this_point - this_child_axis(1,:));
        d_bottom = norm(this_point - this_child_axis(end,:));
        ds = [ds min([d_top d_bottom])];
    end; %
    [d,winning_index] = min(ds);
    child_attachment_vertices = [child_attachment_vertices winning_index];
end;
initial_knot_vertices = unique([1 length(this_axis) child_attachment_vertices]);  % 1, end, children
% now initial_knot_vertices is a list of vertices at which at least one child attaches, plus first and last, in order.
% .. Check whether child knotpoints are sufficent:
[worst_deviation,worst_index] = evaluate_these_knotpoints(this_axis,initial_knot_vertices,minimum_knot_separation);
% .. Next, if necessary, progressively add knots until deviation is below tolerance.
previous_knotpoints = initial_knot_vertices;
max_cycles = 20;
counter = 0;
done_adding_knots = false;
if worst_deviation > tolerance
    while ~done_adding_knots
     counter = counter + 1;
     % worst_index was worst on previous round
     % Check that potential new knotpoint is actually far enough ( > minimum_knot_distance) from its knot neighbors
     lower_neighbor = max(previous_knotpoints(find(previous_knotpoints < worst_index)));
     upper_neighbor = min(previous_knotpoints(find(previous_knotpoints > worst_index)));
     if ~isempty(lower_neighbor) ...
			 && norm(this_axis(worst_index,:) - this_axis(lower_neighbor,:)) < minimum_knot_separation   
         new_knotpoint_vertex = round(mean(lower_neighbor,worst_index));
         previous_knotpoints = setdiff(previous_knotpoints,lower_neighbor);  % remove neighbor because the average of it and the worst is about to be added
	 elseif ~isempty(upper_neighbor) && norm(this_axis(worst_index,:) - this_axis(upper_neighbor,:)) < minimum_knot_separation
         new_knotpoint_vertex = round(mean(upper_neighbor,worst_index));
         previous_knotpoints = setdiff(previous_knotpoints,upper_neighbor);
	 else 
		 new_knotpoint_vertex = worst_index;
	 end;
     new_knotpoints = unique([previous_knotpoints new_knotpoint_vertex]); % puts them in order
     [worst_deviation,worst_index] = evaluate_these_knotpoints(this_axis,new_knotpoints,minimum_knot_separation);
     previous_knotpoints = new_knotpoints;
     done_adding_knots = (worst_deviation < tolerance) | counter > max_cycles;
    end;
end; % Now, previous_knotpoints should be set to a set of v's that gives good enough fit (or has max_knotpoints).
final_knotpoints = previous_knotpoints;
skeleton(axis).knotpoints = final_knotpoints; % this keeps a record of WHICH VERTEX on this_axis are knotpoints-- for recovery later.
parent_previous_knot = find_parent_previous_knot(axis,skeleton); % this finds the tangent on the parent axis at the point the i-th axis sprouts.
parameters = knotpoints2parameters(this_axis(final_knotpoints,:),parent_previous_knot);
% This gives parameters-- but still need to finalize things like which knotpoints the children are at.

%.........................................................................................................
function [worst_deviation,worst_index] = evaluate_these_knotpoints(this_axis,vertices,minimum_knot_separation);
% Fits a spline at vertices, evaluating the mismatch at each vertex; returns worse mismatch and index at which it occurs.
%...Next, make spline in canonic coords-- should pass through above points
% Note: we need to make sure that we only evaluate points that are "far enough" from one of the existing knotpoints,
% because we never want to add a new knotpoint too close to an existing one.
knots = this_axis(vertices,:);
n_knotpoints_this_axis = length(knots);
n_points_desired = size(this_axis,1);
t = 1:n_knotpoints_this_axis;
%tt = [1:1/(n_points_desired-1):n_knotpoints_this_axis];
tt = [1:(n_knotpoints_this_axis-1)/(n_points_desired-1):n_knotpoints_this_axis]; 
xx = spline(t,knots(:,1),tt)';
yy = spline(t,knots(:,2),tt)';
% Detrmine points that are far enough from current knots:
legitimate_candidates = [];
for i = 1:n_points_desired
	distance_from_current_points =  min(abs(i - vertices));
	if distance_from_current_points > minimum_knot_separation
		legitimate_candidates = [legitimate_candidates i];
	end;
end;
%... now evaluate deviations: 
deviations = [];
for i=1:length(legitimate_candidates)
	this_candidate = legitimate_candidates(i);
	deviations = [deviations norm([xx(this_candidate) yy(this_candidate)] - this_axis(this_candidate,:))];
end;
[worst_deviation,worst_legitimate_index]  = max(deviations);
worst_index = legitimate_candidates(worst_legitimate_index);


%.........................................................................................................
function parameters = knotpoints2parameters(knots,parent_previous_knot)
% Compute parameters encoding the locations of these knotpoints. 
n_knots = size(knots,1);
parameters = [];
for i=2:n_knots
    if i==2  % so two points before is "virtual zeroth point"
        previous_previous = parent_previous_knot;
    else
        previous_previous = knots(i-2,:);
    end;
    previous_knot = knots(i-1,:);
    this_knot = knots(i,:);
    bend =       turning_angle(previous_previous,previous_knot,this_knot);  % 
    logstretch = log(norm(previous_knot - this_knot)/norm(previous_previous - previous_knot));
    parameters = [parameters bend logstretch];
end;

%........
function angle = turning_angle(a,b,c)
% Computes the turning angle from point a to b to c.
v1 = b - a;
v2 = c - b;
if norm(v1/norm(v1) -v2/norm(v2)) < 0.01  % check for collinearity -- if so atan gives imaginary output
	angle = 0;
else
	cos_angle = dot(v1,v2)/(norm(v1) * norm(v2));
	angle_magnitude =  acos(cos_angle);
	cp = cross([v1 0],[v2 0]);
	angle_sign = sign(cp(3));
	angle = angle_magnitude * angle_sign;
end;
%.........................................................................................................
function axis_hierarchy = create_axis_hierarchy(skeleton)
% Create axis hieararchy using skeleton (which should include additional field 'knotpoints').
% Each row of axis_hierarchy relates to one axis, containing:
% [index  parent_index parent_knotpoint first_paramter last_parameter]
axis_hierarchy = [];
for a=1:length(skeleton)
    this_index = skeleton(a).index;
    this_parent = skeleton(a).parent;
    this_first_parameter = skeleton(a).parameter_locations(1);
    this_last_parameter =  skeleton(a).parameter_locations(2);
    this_parent = skeleton(a).parent;
    % Finally, figure out which knotpoint on parent this guy is closest to.
    if this_parent == -1 % root
        this_parent_knotpoint = -1; % root is not attached to its parent anywhere (-1), since it doesn't have a parent (-1)
    else % normal case
        top =    skeleton(a).contour(1,:);
        bottom = skeleton(a).contour(end,:);
        ds = [];
        this_parent_knotvertices = skeleton(this_parent).knotpoints;
        this_parent_knotpoints = skeleton(this_parent).contour(this_parent_knotvertices,:);
        for v=1:size(this_parent_knotpoints,1)
            this_knotpoint = this_parent_knotpoints(v,:);
            ds = [ds min([norm(this_knotpoint - top) norm(this_knotpoint - bottom)])];
        end;
        [winning_distance winning_index] = min(ds);
        this_parent_knotpoint = winning_index;  % ie index WITHIN knotpoints (1... n_knotpoints) of knotpoint to which it belongs.
    % Finally, put it all together in axis_hierarchy
    end;
    this_hierarchy_row = ...
        [this_index this_parent this_parent_knotpoint this_first_parameter this_last_parameter];
    axis_hierarchy = [axis_hierarchy; this_hierarchy_row]; % add one row to hiearchy
end;

    
%..........................................................................
function parent_tangent = find_parent_direction(index,skeleton);
% finds the normal direction at the root of the i-th axis
this_index = skeleton(index).index;
this_parent = skeleton(index).parent;
if this_parent == -1 % root
    parent_tangent = [1 0];
else 
	distances = [];
    this_point = skeleton(index).contour(1,:);
    % Now find which point on the parent matches this point
    parent_contour = skeleton(this_parent).contour;
    n_parent_points = size(parent_contour,1);
    for i=1:n_parent_points
        this_distance = norm(parent_contour(i,:) - this_point);
		distances = [distances this_distance];
		[winning_distance,winning_index] = min(distances);
    end;
    if winning_index == 1 % first point, need to take NEXT link
        parent_direction = parent_contour(2,:) - parent_contour(1,:);
	elseif winning_index == n_parent_points
		parent_direction = parent_contour(winning_index - 1,:) - parent_contour(winning_index,:); % careful-- this means the head of the vector
		      % points INWARD along the contour regardless of which way the index numbers nominally point. This makes it direction invariant. 
	else 
	% NORMAL case: use previous and next
        parent_direction = parent_contour(winning_index+1,:) - parent_contour(winning_index - 1,:);
		  % Here the direction does depend on the numbering but value itself should be invariant to direction (because bidirectional). 
    end;
    parent_tangent = parent_direction/norm(parent_direction);
end;
    
%..........................................................................
function parent_previous_knot = find_parent_previous_knot(index,skeleton);
% finds the knot on the parent right "before" this branch point. Used for calculating parameters later. 
this_index = skeleton(index).index;
this_parent = skeleton(index).parent;
if this_parent == -1 % root
    parent_previous_knot = [0 0];
else 
	distances = [];
    this_point = skeleton(index).contour(1,:);
    % Now find which point on the parent matches this point
    parent_contour = skeleton(this_parent).contour;
	parent_contour_knotpoint_indices = skeleton(this_parent).knotpoints; % which indexs in parent_contour are the knots
	n_parent_knotpoints = length(parent_contour_knotpoint_indices);
    for i=1:n_parent_knotpoints
		this_parentknotpoint = parent_contour(parent_contour_knotpoint_indices(i),:);
        this_distance = norm(this_parentknotpoint - this_point);
		distances = [distances this_distance];
		[winning_distance,winning_index] = min(distances);
    end;
    if winning_index == 1 % first knotpoint, special case; invent a virtual previous point right before the first, using 1-2 vector in the opposite direction
		first_knot  = parent_contour(parent_contour_knotpoint_indices(1),:);
		second_knot = parent_contour(parent_contour_knotpoint_indices(2),:);
		parent_previous_knot = first_knot + first_knot - second_knot;
	else 
	% NORMAL case: simply use previous knot
	    parent_previous_knot = parent_contour(parent_contour_knotpoint_indices(winning_index - 1),:);
    end;
end;

