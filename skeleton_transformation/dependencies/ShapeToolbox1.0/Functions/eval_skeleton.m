function error = eval_skeleton(candidate_skeleton_parameters)
% error = eval_skeleton(candidate_skeleton_parameters)
% 
% Computes the DL of the skeleton corresponding to the given parameters. 
% Called by fminsearch as part of the gradient descent procedure in
% mapskeleton. 

% This version uses (global) optimziation_scale_factor.
% optimziation_scale_factor is a cheat to force fminsearch to use different step sizes. 
% We multiply all the parameters by this, then divide again before using
% them. fminsearch takes steps of size .05 x the parameter value, which can't
% easily be changed. So by multilying by a scale factor we increases or
% decrease the step size. 

% Note: See parameters2axis3 etc to understand parameterization.
% Basically, there is 1 initial paramter, and then n_parameters_per_axis per axis, 
% n_parameters_per_axis =  3 + (2 * (n_knotpoints_per_axis - 2)),
% eg for 3 knots points there are 5 params per axis, which are
%
% position            = parameters(1);
% initial_angle       = parameters(2);
% bend                = parameters(3);
% logstretch          = parameters(4);
% scale               = parameters(5);
% 
global show_candidates;
global comet_trail;
global part_model;
global current_coribs;  % For use if we don't recompute them every opt step
global rerib_frequency;  % how often (num opt steps) we recompute the ribs
global optimization_counter; % how many steps we've taken
global current_shape;
global current_movie;
global filming;
global current_skeleton;
global current_axis_hierarchy;
global optimization_scale_factor;

%------------------------------------------------------------------------------------------------------
candidate_skeleton_parameters = (1/optimization_scale_factor) * candidate_skeleton_parameters; % Undo adjustment
%------------------------------------------------------------------------------------------------------
skeleton = parameters2skeleton(candidate_skeleton_parameters,current_axis_hierarchy,0);
current_skeleton = skeleton;
optimization_counter = optimization_counter + 1;
% It would be nice to re-prune occasionally, but this actually conflicts with the use of fminsearch, which requires a constant parameterization!
if mod(optimization_counter,rerib_frequency) == 1
	coribs = compute_coribs(skeleton);   % re-rib
	current_coribs = coribs;
    cla;
    draw_shape(current_shape);
    comet_trail = [comet_trail; candidate_skeleton_parameters];  % This version, only adding comet trail on mod==1 steps
    draw_skeleton(skeleton,'rainbow');     % draw current estimate
    drawnow;
else
    coribs = current_coribs;              % if not re-ribbing-- inherit coribs from most recent re-ribbification
end;
comet_trail = [comet_trail; candidate_skeleton_parameters];
[error,newskeleton] = description_length(skeleton);
current_skeleton = newskeleton;
if mod(optimization_counter,rerib_frequency) == 1 
       fprintf('Description length: %f\n',error);
end;