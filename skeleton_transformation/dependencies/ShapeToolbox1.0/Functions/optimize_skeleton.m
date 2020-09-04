function final_skeleton = optimize_skleton(starting_skeleton,shape)
% final_skeleton = optimize_skleton(starting_skeleton)
%
% Optimize skeleton with respect to shape, using fminsearch (Nelder-Mead).
% We use a (global) optimziation_scale_factor, which is a trick to get
% fminsearch to use different size (usually larger) steps. 

% Global declarations:
% Assumes set_skeleton_constants has been run.
global current_skeleton_parameters;
global current_axis_hierarchy;
global current_shape;
global current_skeleton;
global max_optimization_iterations;
global current_coribs;
global current_shape_handle;
global optimization_scale_factor;
global shape_axis_handle;
global comet_trail;
comet_trail = [];
if isequal(shape(1,:),shape(end,:))
    shape = shape(1:end-1,:);
end;
current_shape = shape; % Assign globally
[current_skeleton_parameters,current_axis_hierarchy] = skeleton2parameters(starting_skeleton);
fprintf('Optimizing posterior; please wait .... \n');
%fprintf('If you lose patience, interrupt (control-c) and then type "winner")\n');
initial_skeleton_pars = current_skeleton_parameters;
options = optimset('MaxIter',max_optimization_iterations);
current_coribs = compute_coribs(parameters2skeleton(initial_skeleton_pars,current_axis_hierarchy)); % set initial coribs
optimal_skeleton_pars = fminsearch(@eval_skeleton,optimization_scale_factor * initial_skeleton_pars',options); % optimization step
optimal_skeleton_pars = optimal_skeleton_pars/optimization_scale_factor;
%format_parameters(optimal_skeleton_pars);
winning_skeleton = parameters2skeleton(optimal_skeleton_pars,current_axis_hierarchy);
current_skeleton = winning_skeleton;  % because draw_coribs needs it
axes(shape_axis_handle);        % makes shape axis current
cla;                           % clear axes
current_shape_handle = draw_shape(current_shape);  % redraw just the shape
draw_skeleton(winning_skeleton,'rainbow');
draw_coribs(compute_coribs(winning_skeleton),'rainbow');  % or default_skeleton_drawing_style
