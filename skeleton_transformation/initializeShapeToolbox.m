% start shapeToolbox

% Declarations:
global shape_axis_handle;
global animals;
global window_width;
global current_shape;
global default_skeleton_drawing_style; 
global part_model;
global current_skeleton;
global currently_collecting; % 1 or 0; determines whether mousemotions have any effect
global which_collecting;  % 'shape','skeleton'
                %if we are collecting, determines whether we're adding to
global last_point; % last point collected on whatever we're collecting
global min_distance_new_point;  % how far mouse must travel before new point is added
global current_skeleton_index;  % which axis are we currently adding
global current_axis_contour;  % buffer for the contour of the current axis we're working on
global current_axis;            % current axis nwe're working on
global current_parentindex;          % parent of axis contoru we're working on
global current_parentvertex;    %   ... and vertx on it.
global axis_handle_dls;         % handle for plot of dl's. 
global label_column;
global current_preskeleton_axis;  % does this need to be global?
%global n_initial_parameters;    % num of parameters preceding the first per-axis parameter set (now 1)
%global n_parameters; % number of parameters. Currently 4 but can be changed in future
global n_knotpoints_per_axis;       % n KNOT POINTS per axis -- default 3current_skeleton_parameters
global current_skeleton_parameters;  % contains parameters of current skeleton, n_parameters per axis
global n_points_per_axis;  % number of points (vertices, NOT knotpoints) on the splines that are fit to axis segments.
global current_axis_index;
global line_width;
global current_axis_hierarchy; % vector of indices; the i-th one is the parent of the i-th axis
global shape_spacing;
global distance_decay;  % ... These two constants control both "responsibility" computation and likelihood computation. 
global vonMises_b; 
global branching_penalty;  % cost of one new branch in skeleton
global show_candidates;
global rerib_frequency;  % how often (num opt steps) we recompute the ribs
global optimization_counter; % how many steps we've taken
global comet_trail;  % when optimizing, saves intermediate candidate skeleton parameters
global h_comet_trail;
global button_left;
global button_width;
global button_height;
global button_bottom;
global current_coribs;
global filming;
global n_attneave_axes;  % for random attneave shapes
global attneave_smoothing;
global current_movie;
global current_shape_normals;  % store these globally so compute_coribs can use them
global current_parent;     % used in parametric interface; the index of the axis that contains the
global current_image;
global optimization_tolerance; 
global max_optimization_iterations;
global h_figure;
global light_gray;
global parameter_label_altitude;
%global h_branching_penalty;  
global h_graphical_output;
global show_progressive_computation;
global show_ribs; % Determines whether display of ribs is included when skeletons are drawn.
global h_show_ribs; 
global h_usegradientdescent;
global usegradientdescent;
global optimization_scale_factor; % Controls step size in optimization proess
global current_augmented_skeleton; % Separate current skeleton, so we can keep the regular one "clean"

%addpath('Functions');
set_skeleton_constants;  % sets globals 
%require_twosided_explanation = false;
rand('state',sum(100*clock));   % initialize the random number generator