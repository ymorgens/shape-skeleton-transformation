function complexity = skeleton_complexity(skeleton)
% complexity = skeleton_complexity(skeleton)
% 
% Compute cumulative suprisal summed over all axes, including branching
% complexity.

global branching_penalty;

surprisal_weighting = 10;

complexity_vector = [];
for axis = 1:length(skeleton)
    Xs = skeleton(axis).contour(:,1);
    Ys = skeleton(axis).contour(:,2);
    complexity_vector = [complexity_vector; cumulative_surprisal(Xs, Ys, 0, 1)];
end;
%---------------------------------------
branching_complexity = branching_penalty * length(skeleton); 
complexity = surprisal_weighting*sum(complexity_vector) + branching_complexity;