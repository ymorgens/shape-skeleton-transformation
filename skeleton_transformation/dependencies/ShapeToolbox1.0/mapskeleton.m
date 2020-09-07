function [skeleton,optimal_skeleton_pars,axis_hierarchy,final_description_length] = mapskeleton(shape,varargin)
% [skeleton,optimal_skeleton_pars,axis_hierarchy,final_description_length] = mapskeleton(shape,varargin)
%
% Estimate the MAP skeleton (see Feldman & Singh, PNAS 2006). 
%
% Input shape is an nx2 array of numbers, corresponding to the n points
% defining the boundary of a shape. 
% Output skeleton is a hierarchical set of axes, defined as a structure with 
% fields .index, .parent, and .contour. Each axis's parent is the axis to
% which it is attached; index 1 is the root axis. Each axis itself consists
% of an open curve defined by skeleton(i).contour. Use draw_skeleton to
% draw the skeleton. 
%
% Additional outputs include the final DL, which is the -log probability of
% the shape conditioned on the MAP skeleton, as well as several others used 
% internally by other functions in the toolbox and not intended to be
% interepreted directly. 
%
% Options:
%
% gradientdescent:  include optimization [EM] step (slower)
%
% graphics: Draw the shape and MAP skeleton. 
%
% show_progressive_computation: Graphically display the pruning process,
%    and (if optimization step is include) display the progressive
%    optimization process.  This is kind of cool and is the default in the GUI. 
%
% using_gui: Flag used only by GUI, which changes where the some parameters
% come from. 
%
% Default settings (no extra arguments) give fastest results, e.g.:
% skel = mapskeleton(shape);
% draw_skeleton(skel);
%
% For slower results more closely corresponding to the Feldman & Singh
% paper, include the gradient descent step, e.g.:
% skel = mapskeleton(shape,'gradientdescent');
% draw_skeleton(skel);
%
% For additional documentation, see README.txt.
%
%------------------------------------------------------------------------------------------:
%addpath('Functions');      % All required functions are located in this folder, except those available at the command line. 
% Set defaults:
graphics = false;           % If true, draws shape and skeleton while done.
show_progressive_computation = false; % If true, shows the pruning and gradient descent process graphically. 
gradientdescent = false;   % If true, include the "EM" stage. Skipping this saves a lot of time, so the default is FALSE.
if ismember('gradientdescent',varargin)
	gradientdescent = true;
end;
if ismember('graphics',varargin)
	graphics = true;
end;
if ismember('show_progressive_computation',varargin)
	show_progressive_computation = true;
end
if ismember('using_gui',varargin)
	using_gui = true; % This means we are inside the GUI, which changes a few minor things. 
else
    using_gui = false;
end;
optimal_skeleton_pars =[];
% Above need to set as global, get constants if NOT using GUI
%------------------------------------------------------------------------------------------
% Initialization
global current_shape;
global current_skeleton;
global current_coribs;
global current_shape_normals;
global current_skeleton_parameters;
global shape_axis_handle;
global current_axis_hierarchy;
global optimization_scale_factor;
global max_optimization_iterations;
global comet_trail;
comet_trail = [];    % This is the sequnce of parameters via the gradient descent: can be used to track or display the optimization.
%------------------------------------------------------------------------------------------
if ~using_gui % We only set the constants here if mapskeleton is being used stand-alone. Otherwise the GUI controls the constants. 
    set_skeleton_constants; % Sets lots of globals that are used down the line. 
end;
if isequal(shape(1,:),shape(end,:))  % if first and last points are the same
    shape = shape(1:end-1,:);       % delete the last point
end;
current_shape = shape;
current_shape = insideright_shape(current_shape);    % fix shape polarity 
current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));  % Install shape normals
% --------------------------------------------------------------------------------------	
% Prune skeleton
mat = make_medial_axis(smooth_shape(current_shape,10));  % We start with medial axis transform (of a slightly smoothed shape).
pruned_skeleton  = metaprune_skeleton(mat,show_progressive_computation);  % prune skeleton, showing progressive pruning if flag is TRUE
%------------------------------------------------------------------------------------------
% Gradient descent sequence
if gradientdescent
    final_skeleton = optimize_skeleton(pruned_skeleton,current_shape);
else
    final_skeleton = pruned_skeleton;
end;
% set final globals: .....................
current_skeleton = final_skeleton;  % because draw_coribs needs it
% Drawing:
if show_progressive_computation
    if using_gui
    	axes(shape_axis_handle);        
    end;
	cla;                           
	current_shape_handle = draw_shape(current_shape);  % redraw just the shape
	draw_skeleton(final_skeleton,'rainbow');
end;
%........................................
% Final return variables:
final_description_length = description_length(final_skeleton);
axis_hierarchy = current_axis_hierarchy;
skeleton = final_skeleton;
%........................................
if graphics
    cla;
    draw_shape(shape);
    draw_skeleton(skeleton);
end;
    