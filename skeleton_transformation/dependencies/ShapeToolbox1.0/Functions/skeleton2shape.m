function shape = skeleton2shape(skeleton,varargin)
% shape = skeleton2shape(skeleton,varargin)
%
% Generate a random shape given a skeleton, model, and parameters.

% Generates mean rib lengths from a negative exponential, with parameters determined based on generation of
% part, so that root is thickest, children are thinner, etc. 

%------------------------------------------------------------------------
global distance_decay;
n_points_around_end = 10;
%------------------------------------------------------------------------
% defaults:
root_width = 10;  % expected width of root; children have progressively smaller mean widths.
%------------------------------------------------------------------------
% parameters:
dummy = option_value(varargin,'root_width'); 
    if dummy ~= 'null' 
        root_width = dummy;
    end;
%------------------------------------------------------------------------
n_axes = length(skeleton);
% First, make "blobs"--- polygons corresponding to each axis
for a = 1:n_axes
    this_blob = [];
    this_axis = skeleton(a).contour;
    this_axis_normals = jfnormals(this_axis(:,1),this_axis(:,2));
    %this_mean_rib = parameters(a);
    this_mean_rib = compute_mean_rib_length(root_width,a,skeleton);
    for v=1:size(this_axis,1)  % first we go up one side
        this_point = this_axis(v,:) + this_mean_rib * this_axis_normals(v,:);
        this_blob = [this_blob; this_point];
    end;
    unit_step_around_end = pi/n_points_around_end;
    for i=0:n_points_around_end-1
        angle = -i * unit_step_around_end;
        rotation_matrix = [cos(angle) sin(angle); -sin(angle) cos(angle)];
        final_normal = this_axis_normals(end,:);
        this_normal = final_normal * rotation_matrix;
        this_point = this_axis(end,:) + this_mean_rib * this_normal;
        this_blob = [this_blob; this_point];
    end;
    for v=size(this_axis,1):-1:1 % then down the other
        this_point = this_axis(v,:) - this_mean_rib * this_axis_normals(v,:); % the other way
        this_blob = [this_blob; this_point];
    end;
    if a==1 % only for root:
        for i=n_points_around_end-1:-1:0
        angle = i * unit_step_around_end;
        rotation_matrix = [cos(angle) sin(angle); -sin(angle) cos(angle)];
        first_normal = this_axis_normals(1,:);
        this_normal = first_normal * rotation_matrix;
        this_point = this_axis(1,:) + this_mean_rib * this_normal;
        this_blob = [this_blob; this_point];
        end;
    end;
    skeleton(a).blob = this_blob;  % add these polygons to the skeleton structure
end;
%... Now join blobs together to make one big shape.
shape = skeleton(1).blob;
for a=2:length(skeleton)
%    shape = shape_union(shape,skeleton(a).blob); %  joing a-th blob with shape so far
     shape = combine_shapes(shape,skeleton(a).blob,'union'); %  joing a-th blob with shape so far
end; % now shape is entire n-way union

%------------------------------------------------------------------------
function this_mean_rib = compute_mean_rib_length(root_width,a,skeleton)
% determines mean rib based on generation of a within skeleton
reduction_fraction = 0.5; % reduction in width for each generation
generation = determine_generation(a,skeleton);
this_mean_rib  = root_width * reduction_fraction ^ generation;

%.........
function generation = determine_generation(a,skeleton);
if a == 1
    generation = 0;
else
    parents_generation = determine_generation(skeleton(a).parent,skeleton);
    generation = parents_generation + 1;
end;

%------------------------------------------------------------------------
function value = option_value(options,name)
% Given an arg cell array options, returns the value corresponding to option called 'name'

value = 'null'; % returns if 'name' is not set in options
for i=1:length(options)
    if isequal(options{i},name)
        value = options{i+1};
    end;
end;