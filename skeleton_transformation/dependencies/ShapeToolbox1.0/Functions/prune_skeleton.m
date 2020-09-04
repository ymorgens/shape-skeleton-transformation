function [pruned_skeleton,n_axes_removed] = prune_skeleton(original_skeleton,original_coribs,show_progressive_computation)
% [pruned_skeleton,n_axes_removed] = prune_skeleton5(original_skeleton,original_coribs,show_progressive_computation)
%
% Removes "nonsignificant" axes, using posterior ratio test as described in
% Feldman & Singh (2006). 

% Proceeds in "depth order", trimming leaves, then their parents, then
% their parents, etc. 
%
% We use "partial reribbing" in order to go faster, meaning that when we 
% consider removing an axis, we rerib ONLY those points explained by the axis 
% in questions, and leave everything else alone. This leads to
% a number of complexities about what index system we are using to refer to
% the axes.
%
% Some details on the reribbing process:
% All pruning, skeletons, etc are constructed using indexes from original
% "reference" skeleton. This keeps everything clean -- for example, you have
% a fixed depth chart, so you know when you are done. (IE indices have never changeded.)
% But: the problem is how to "rerib" efficiently, keeping a running tab of
% ribs, even though the progressively smaller skeletons have different
% numbering systems. 
%
% Solution: We keepk "parallel books" in two indexing systems. 
% Main loop is controlled by original depth chart using original indexing
% system (reference_skeleton). We also have a running current_temp_coribs and corresponnding 
% skeleton. Each time we consider an index in depth_chart, we figure out
% what index it corresponds to in current_temp-skeleton, and delete and
% rerib there. THis is necessary in order to have an efficient re-ribbing
% process, where we can just rerib the shape points corresponding to that
% index. If we keep the axis, skeleton indexes don't change. If we prune
% the axis, we recompute the skeleton with a compressed (changed) index
% system. THen repeat until lowest depth on original depth chart. 

%...................................................................................................
global current_temp_coribs;  % These two will hold current hypothesized skeleton and corresponding ribs. 
global current_temp_skeleton; 
global current_temp_tt;  % tt in force to refer to current_temp_skeleton (wrt reference)
global current_coribs;       % This is used by other functions (e.g. description_length) -- temp ones are not.
global current_shape;

%...................................................................................................
n_original = length(original_skeleton);
%...................................................................................................
% First, trim from roots:
if isempty(original_coribs)
	current_coribs = compute_coribs(original_skeleton);   % re-rib
else
	current_coribs = original_coribs;
end;
dl_original = description_length(original_skeleton);
possibly_derooting_skeletons = false;
if possibly_derooting_skeletons   % flag for allowing the root to be pruned-- requires an elaborate special procedure
	skeleton = deroot_skeleton(original_skeleton,dl_original);  % skeleton is now starting point for further trimming:
else
	skeleton = original_skeleton;
end;

%...................................................................................................
% Next, trim in reverse depth order, leaves to root. 

reference_skeleton = skeleton; % to which all indices in topologies refer. 
starting_axis = 1;
original_topology = skeleton2topology(reference_skeleton,starting_axis);
depth_chart = make_depth_chart(original_topology);  % Creates a list of depths for each axis
current_temp_skeleton = reference_skeleton;
current_temp_coribs   = current_coribs; % Axis indices refer to current_temp_skeleton (and current_shape). NOT reference_skeleton
current_temp_tt = 1:length(current_temp_skeleton);  % initialize tt. Skeleeton hasn't been pruned yet, so tt is just 1:n. 
current_topology = original_topology; % This will hold evolving topology as it gets trimmed.
  % Indexes always refer to original skeleton. 
max_depth = max(depth_chart); % highest depth in the depth chart, ie depth of most remote leaves
for depth=0:max_depth-1  %from leaves to root
    axes_at_this_depth = find(depth_chart==depth);  % All axes that have this depth   
    for i=1:length(axes_at_this_depth)
        axis_under_consideration = axes_at_this_depth(i);  % REF coordinates
        [skeleton_with,skeleton_without,ribs_with,ribs_without,pruned_topology,tt_without,tt_with] = find_relative_skeletons_coribs(axis_under_consideration,current_topology,reference_skeleton);
        current_coribs = ribs_with;  % set global for description_length
        bad_explainer = bad_axis(ribs_with,find(tt_with==axis_under_consideration));  
        keep = 0;
        if ~bad_explainer
            % If not bad, compute DLs
            
            DL_with = description_length(skeleton_with);
            current_coribs = ribs_without;  % set global for description_length
            DL_without = description_length(skeleton_without);
            keep = (DL_with < DL_without);
        end;
        % Prune topology: ---
        if ~keep 
            current_topology = pruned_topology;  % Remove this axis from evolving topology. 
            current_temp_skeleton = skeleton_without;
            current_temp_coribs = ribs_without;
            current_temp_tt = tt_without;
            if show_progressive_computation   % Flag for whether or now to show display graphical progress
                cla; 
                draw_shape(current_shape);
                draw_skeleton(current_temp_skeleton,'plain');
                drawnow;
            end;
        end;
    end;
end;
pruned_skeleton = current_temp_skeleton; 
n_axes_removed= n_original - length(pruned_skeleton);

%...................................................................................................
function [skel_with,skel_without,ribs_with,ribs_without,pruned_topology,tt_without,tt_with] = find_relative_skeletons_coribs(axis_under_consideration,current_topology,reference_skeleton)
% Find skeletons and ribs WITH this axis and then re-rib WITHOUT this axis. 
% The ribs are "relative" in the sense that the WITHOUT version of the ribs is computed
% by removing one axis (and its descendants) from the WITH version. However
% the two skeletons then DON'T have common indexing systems, so we need to
% give ribs for each of them in their respective indexing systems. 
%
% We use translation_tables to figure out which axis in the WITH skeleton
% corresponds to which in the WITHOUT skeleton. Coribs must alwasy be
% paired with skeletons. 

% axis_under_consideration refers to reference_skeleton and topology


global current_temp_coribs;  
global current_temp_skeleton;
global current_temp_tt; % alsways refers to current_temp_skeleton
pruned_topology = remove_from_topology(axis_under_consideration,current_topology);  % Ref coordinates
skel_with = current_temp_skeleton;   % WITH coordinates
ribs_with  = current_temp_coribs;    % WITH coordinates
tt_with = current_temp_tt;  % for translating from reference to ribs_with
%[skel_with,tt_with] =               topology2skeleton(current_topology,reference_skeleton); % These have NOVEL indicies
[skel_without,tt_without] =         topology2skeleton(pruned_topology,reference_skeleton); % These have NOVEL indicies
descendants = find_descendants(axis_under_consideration,current_topology);  % Need to trim these indices from skel_with_topology
  % NB: descendants are in REF coordinats. 
% Convert descendants into WITH coordinates:
for i=1:length(descendants)
    descendants(i) = find(descendants(i) == tt_with);
end;
points_explained = ribs_with(ismember(ribs_with(:,2),descendants),1); % these are the points explained by the axes to be trimmed
new_ribs = compute_coribs(skel_without,points_explained); % "rapid" recomputation of ribs just replacing the trimmed ones. 
% Convert ribs_with to WITHOUT coordinate system, 
% then choose the ribs in this system that explain the OTHER points
ribs_to_keep = ribs_with(~ismember(ribs_with(:,1),points_explained),:);  % Fibs that explain OTHER points, in WITH coordinate system
for i=1:size(ribs_to_keep,1)
    this_rib_with = ribs_to_keep(i,2);  % WITH coordinate system
    this_rib_ref = current_temp_tt(this_rib_with); % REF coordiante system   
    this_rib_without = find(tt_without == this_rib_ref); % % WITHOUT coordinate system
    ribs_to_keep(i,2) = this_rib_without;  % replace rib in ribs_to_keep with translated coord
end;
ribs_without = [ribs_to_keep; new_ribs]; % combine ribs on NONtrimmed points with new ribs for those points. 

%...................................................................................................
function depth_chart = make_depth_chart(topology)
% Create list of all axes in topology along with their depths
% 
% This version gives 0 to leaves, max depth to root
% Ie each branch's depth is its number of generations of descendants. 
n_axes = size(topology,1); %
depth_chart = zeros(n_axes,1); 
for i=1:n_axes
    depth_chart(i) = n_generations_hanging_from(i,topology);
end;
% ----
function n_generations = n_generations_hanging_from(index,topology)
immediate_children = topology(topology(:,1)==index,2);
if isempty(immediate_children)
    n_generations = 0;
else
    kids_depths = [];
    for i=1:length(immediate_children)
        kids_depths = [kids_depths n_generations_hanging_from(immediate_children(i),topology)];
    end;
    n_generations = 1+max(kids_depths);
end;


%...................................................................................................
function keep = evaluate_axis(index,current_topology,reference_skeleton)
global current_coribs;
global current_shape;

skeleton_with = topology2skeleton(current_topology,reference_skeleton);
current_coribs = compute_coribs(skeleton_with);  % In new reference frame
skeleton_without = topology2skeleton(remove_from_topology(index,current_topology),reference_skeleton);
current_coribs = compute_coribs(skeleton_with);  %
DL_with = description_length(skeleton_with);
%...
skeleton_without = topology2skeleton(remove_from_topology(index,current_topology),reference_skeleton);
current_coribs = compute_coribs(skeleton_without);
DL_without = description_length(skeleton_without);
keep = DL_with < DL_without; 
if show_progressive_computation
    if keep
        cla;
        draw_skeleton(skeleton_with,'plain');
        drawnow;
    else
        cla;
        draw_skeleton(skeleton_without,'plain');
        drawnow;
    end;
end;

%...................................................................................................
function trimmed_topology = remove_from_topology(index,topology)
% Removes  index and all its descendents from topology
descendants = find_descendants(index,topology);
trimmed_topology = topology(~ismember(topology(:,1),descendants),:);  
trimmed_topology = topology(~ismember(topology(:,2),descendants),:); 

%...................................................................................................
function descendants = find_descendants(index,topology)
% Returns indices of children, grandchildren, etc. Including itself.
descendants = [index];
immediate_children = topology(topology(:,1)==index,2);
for i=1:length(immediate_children)
    descendants = [descendants find_descendants(immediate_children(i),topology)];
end;

%...................................................................................................
function [skeleton,translation_table] = topology2skeleton(topology,reference_skeleton)
% Removes axes not mentioned in topology, then "compact" (renumber) axes 
% so there are no gaps. 
%
% First, find the axes that are included.
% Then, renumber everything so only they exist. 
%
% translation_table has one entry per axis of the outgoing skeleton
% ("skeleton"), giving the corresponding index in reference_skeleton. 
if isempty(topology)
    error('Can"t make a skeleton from an empty topology.');
end;
included_axes = setdiff(unique(topology(:)),-1);  % list of axes to include
n_axes = length(included_axes);
% Now we use included_axes as a translation table, ie. included_axes are renumbered as 1:n_axes
for i=1:n_axes
    skeleton(i).index = i;
    reference_index = included_axes(i);  % index that this axis has in the original skeleton 
    translation_table(i) = reference_index; % Enter index of this axis in reference_skeleton
    the_parent = find(included_axes == reference_skeleton(reference_index).parent); % Find "new index" of this skeleotn's parent
    if isempty(the_parent)
        skeleton(i).parent = -1;
    else
        skeleton(i).parent = the_parent;
    end;
    skeleton(i).contour = reference_skeleton(reference_index).contour;
end;

% %...................................................................................................
function bad = bad_axis(ribs,axis) 
% Checks if given axis explains from ONLY one point, which we don't like. 
set_of_vertices_explaining = unique(ribs(ribs(:,2)==axis,3));
bad = (length(set_of_vertices_explaining) < 2);