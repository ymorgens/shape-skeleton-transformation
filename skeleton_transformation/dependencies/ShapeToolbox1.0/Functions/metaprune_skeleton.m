function final_skeleton = metaprune_skeleton(initial_skeleton,show_progressive_pruning)
% final_skeleton = metaprune_skeleton(initial_skeleton,show_progressive_pruning)
%
% Repeatedly prune skeleton until n axes stabilizes. Called by mapskeleton.

global current_coribs;
global current_temp_coribs;
reparameterizing = 0;
if reparameterizing
    [p,h] = skeleton2parameters(initial_skeleton);
    skeleton = parameters2skeleton(p,h);  
else
    skeleton = initial_skeleton;
end;
done = false;
while ~done
	current_coribs = compute_coribs(skeleton); 
	new_coribs = current_coribs;
	n_initial = length(skeleton);
	[new_skeleton,n_axes_removed] = prune_skeleton(skeleton,new_coribs,show_progressive_pruning);
	n_next = length(new_skeleton);
	done = isequal(n_axes_removed,0);    
	skeleton = new_skeleton;
end;
current_coribs = current_temp_coribs; % inherited from prune_skeleton
final_skeleton = new_skeleton;
