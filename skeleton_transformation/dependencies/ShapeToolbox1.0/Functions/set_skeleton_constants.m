function set_skeleton_constants
% set_skeleton_constants
%
% Sets globals for the shape_toolbox.
% ----------------------------------------------------------------
% Declarations:
global shape_axis_handle;
global window_width;
global default_skeleton_drawing_style; 
global current_shape;
global current_skeleton;
global currently_collecting; % 1 or 0; determines whether mousemotions have any effect
global which_collecting;  % 'shape','skeleton'
                %if we are collecting, determines whether we're adding to
global last_point; % last point collected on whatever we're collecting
global min_distance_new_point;  % how far mouse must travel before new point is added
global current_skeleton_index;  % which axis are we currently adding
global current_axis_contour;  % buffer for the contour of the current axis we're working on
global current_axis;            % current axis we're working on
global current_parentindex;          % parent of axis contoru we're working on
global current_parentvertex;    %   ... and vertx on it.
global axis_handle_dls;         % handle for plot of dl's. 
global current_preskeleton_axis;  % does this need to be global?
global n_knotpoints_per_axis;       % n KNOT POINTS per axis -- default 3
global current_skeleton_parameters;  % contains parameters of current skeleton, n_parameters per axis
global n_points_per_axis;
global current_axis_index;
global current_axis_hierarchy; % vector of indices; the i-th one is the parent of the i-th axis
global shape_spacing;
global show_candidates;
global button_left;
global button_width;
global button_height;
global button_bottom;
global current_movie;
global rerib_frequency;  % how often (num opt steps) we recompute the ribs
global optimization_counter; % how many steps we've taken
global current_shape_handle;
global distance_decay;  % ... These two constants control both "responsibility" computation and likelihood computation. 
global vonMises_b;      % ..
global branching_penalty;
global line_width;
global n_attneave_axes;  % for random attneave shapes
global attneave_smoothing;
global label_column;
global current_image;
global require_twosided_explanation;
global optimization_tolerance;
global max_optimization_iterations;
global skeleton_approximation_tolerance;
global filming;
global show_progressive_computation; % Flag for showing (1) progressive pruning of branches and (2) gradient descent procedure
global show_ribs; 
global light_gray;
global usegradientdescent;
global optimization_scale_factor;

%---------------------------------------------------------------------------------------
% Initializations:
n_knotpoints_per_axis = 5;
require_twosided_explanation = true;
distance_decay = 10;
vonMises_b = 1;
branching_penalty = 50;
line_width = 2;  % use 2 for plots to export, 1 normal
filming = false;  % 1 or 0, determines if we are taking snapshots (time consuming)
min_distance_new_point = 3;
rerib_frequency = 2;  % how often (num opt steps) we recompute the ribs
optimization_counter = 0; % how many steps we've taken
default_skeleton_drawing_style = 'none'; 
n_points_per_axis = 30; % 20 is a typical value, 15 to simplify things
n_attneave_axes = 7;  % for random attneave shapes
attneave_smoothing = 5;
shape_spacing = 5; % distance from point to point along shape (note that in optimization we are taking one likelihood per point, so the
                    % larger this parameter, the fewer points, the fewer likelihoods, and the quicker the opimization can proceed
skeleton_approximation_tolerance = .1;                    
button_left = 40;
button_width = 85;
button_height = 20;
button_bottom = 70;
window_width = 300;
label_column = 120;
show_candidates = 1;  % whether to show each candidate skel as it's considired-- not real time (yet)
current_skeleton_index = 1;
current_preskeleton_axis = [];
current_skeleton_parameters = [];
current_skeleton = [];
current_axis_contour = [];
current_axis = [];
current_axis_index = 0;
currently_collecting = 0;   % will be turned on by first mouse click
which_collecting = 'none'; % default
optimization_tolerance = 0.1; % doesn't seem to work
max_optimization_iterations = 15;
light_gray = [.7 .7 .7];
show_ribs = 1;
usegradientdescent = 0; % Default is NOT to use gradient descent optimization step. 
show_progressive_computation = 1;
optimization_scale_factor = 4;

