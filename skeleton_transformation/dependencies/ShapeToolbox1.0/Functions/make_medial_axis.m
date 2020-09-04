function skeleton = make_medial_axis(shape)
% skeleton = make_medial_axis(shape)
% 
% Computes Blum medial axis of shape, returning the skeleton in the format
% from Feldman skeleton code (e.g. as drawn by draw_skeleton). 
% Uses a voronoi procedure as suggested by Ogniewicz and Kubler, 1995.
%
% input shape is an nx2 vector of points. output is in skeleton format,
% i.e. has fields .index, .parent, and .contour.

% This procedure is much, much more complicated than strictly necessary to
% compute the MAT, because it is designed to be part of the computation of 
% the map skeleton, and has many bells and whistles associated with the larger
% project (e.g. creating a hierarchy of relatively smooth axes arranged in
% a hierarchy).

% By JF 10/05, revised (cleaned up) 11/15. 
%----------------------------------------------------------------------------------------------------------------
[all_points,polygons] = voronoin(unique(shape,'rows'));  
global counter; 
global points_used_so_far;  % will contain 1s where points are used
counter = 0;
points_used_so_far = [];


inside_indices = find(inpolygon(all_points(:,1),all_points(:,2),shape(:,1),shape(:,2)));
% NOTE: from here on, point indices will be into all_points, for consistency -- this is how polygons refers to them 
terminator_list = []; % --> points that end (ie, just one neighbor); indices into all_points
% The following finds the neighbors of each point of "all_points"--> left and
% right neighbors. If the point is outside of the obj, the neighbors of
% that point is empty. "p" is the index of a point in "all_points" array
for p=1:size(all_points,1)  % make a matrix of 
    if ~ismember(p,inside_indices)
        neighbors(p) = {[]}; % assign exterior points no neighbors, just as a place holder-- we don't really care about them
    else        
        neighbors(p) = {neighbors_of(p,polygons,inside_indices)};  % now neighbors{i} is a vector of p's neighbors -- inside ones only
    end;
    if length(neighbors{p}) == 1  % terminator
        terminator_list = [terminator_list p]; % add this guy to the list of terminators
    end;
end;
% Later, we will use last_terminator as the endpoint of the root.

indexskeleton_so_far = [];
raw_index_skeleton = add_to_skeleton(1,1,terminator_list(1),[],neighbors,indexskeleton_so_far,all_points); % recursive call, make progressive skeleton
adjoinments = make_adjoinments(raw_index_skeleton,all_points);  % makes the list of things that arec connected to each other
smoother_index_axis_set = merge_adjoinments(raw_index_skeleton,adjoinments,all_points); % merges the adjoinments as much as possible
root_index = choose_root(smoother_index_axis_set); % chooses the index with highest lenght (? and graph centrality)
index_skeleton = construct_skeleton(...
	root_index,... % starting from top, recursively make skeleton, BREAKING axes where necessary.
	smoother_index_axis_set,...
	raw_index_skeleton,...
	all_points);
for a =1:length(index_skeleton)
	contour = index_skeleton(a).contour;
	contour_wo = remove_duplicates(contour);
	index_skeleton(a).contour = contour_wo;
end;
skeleton = indexes2skeleton(all_points,index_skeleton);  % make the actual skeleton from index skeleton
skeleton = cleanup_skeleton(skeleton); % clean up a few details of how the skel is organized

%.....................................................................................................................
function cleaned_skeleton = cleanup_skeleton(skeleton)
% Given skeleton, make sure that each axis is affixed in the "right direction", ie that its 1st vertex is actually the one that attaches.
% This comes out wrong some times from earlier procedures, and though it is not visible it screws up the parameterization we use later.
global shape_spacing;
if isempty(shape_spacing) 
    shape_spacing = 5; % default value
end;
for a=1:length(skeleton)
	% check if contour is "reversed", and if so, invert it.
	if skeleton(a).parent > 0 % if not root
		this_contour = skeleton(a).contour;
		parent_contour = skeleton(skeleton(a).parent).contour;
		intersection_point = intersect(this_contour,parent_contour,'rows');
		if intersection_point == this_contour(end,:) % if intersection point is the LAST point on the contour,
			skeleton(a).contour = this_contour(end:-1:1,:);
		end;
	end;
end;
for a=1:length(skeleton) 	% finespace to match shape
	this_contour = skeleton(a).contour;
	[x,y] = finespacedcontour(this_contour(:,1),this_contour(:,2),shape_spacing,0);
	evenlyspaced_version  = [x,y]; % this is guaranteed to include the original first point, which here is the one shared with the parent axis
	% Now we have to make sure the parent also contains the first branch point.
	if skeleton(a).parent > 0
		this_parent_contour = skeleton(skeleton(a).parent).contour;
		% Now find closest point in this_parent_contour and "nudge" it to match first point in branch
		root_point = this_contour(1,:);
		ds = [];
		for i=1:size(this_parent_contour,1)
			ds = [ds norm(this_parent_contour(i,:) - root_point)];
		end;
		[min_d winning_index] = min(ds); % winning index is the point closest to first branch point
		this_parent_contour(winning_index,:) = root_point;
		skeleton(skeleton(a).parent).contour = this_parent_contour;
	end;
end;
cleaned_skeleton = skeleton;


%.....................................................................................................................
function index_skeleton = construct_skeleton(root_index,smoother_index_axis_set,raw_index_skeleton,all_points);
% Make an (index)skeleton starting from root.
% smoother_index_axis_set is a set of long, perhaps too long, axial (index) pieces.
% No neighbor relations have been made yet.
% Starting from root, make skeleton, making neighbor relations as we go, and BREAKING pieces as necessary.
% this is critical to avoid loops. When all are dealt with, we stop regardless of who is neighboring on who.
% First, add the root:
piece_to_add = smoother_index_axis_set(root_index).contour;
index_skeleton_so_far = [];
% Mark 1 as used:
used_so_far = zeros(1,length(smoother_index_axis_set));
used_so_far(root_index) = true;
% Next, start recursing.
index_skeleton = add_to_branching_skeleton_from(piece_to_add,-1,...  % add this "piece", ie collection of points indexes, to skeleton
	index_skeleton_so_far,smoother_index_axis_set,used_so_far);

%........................................................................................................
% Recursive call for construct_skeleton
% Calling this on piece_to_add means "add this piece as the next axis", with parent "parent". Then find its neighbors and recurse.
function [index_skeleton,used_so_far] = add_to_branching_skeleton_from(piece_to_add,its_parent,...
	index_skeleton_so_far,axis_set,used_so_far)
% First, add the piece:.....................
index = length(index_skeleton_so_far) + 1;
index_skeleton_so_far(index).contour = piece_to_add;
index_skeleton_so_far(index).parent = its_parent;
index_skeleton_so_far(index).index = index;
% Next, find its neighbors: We are choosing from among "axis_set".
neighbors_this_piece = find_neighbors2(piece_to_add,axis_set,used_so_far); % find neighbors, excluding any in used_so_far
% Next, go through its neighbors (indexes into axis_set) and either (a) recurse on them or (b) break them up and then recurse on them.
for i=1:length(neighbors_this_piece);
	this_neighbor = neighbors_this_piece(i); % an index that we know connects somehow
	% Mark this index as taken care of
	if used_so_far(this_neighbor) % If  already used, 
		% do nothing
	else
		% Otherwise, mark it as used, and proceed
		used_so_far(this_neighbor) = true; % mark it as dealt with 
		childish_axis_index = this_neighbor; % "child-ish" because we don't yet know if it connects properly, as a branch
		childish_contour = axis_set(childish_axis_index).contour;  
		intersection_point = intersect(piece_to_add,childish_contour);
		% intersection_point is the index of the point they share
		n_points_childish = length(childish_contour);
		intersection_point_location = find(childish_contour == intersection_point); % position along childish contour at which they intersect
		if (intersection_point_location == 1) | (intersection_point_location == n_points_childish)
			% childish intersects at one of ITS endpoints, ie correctly, as a branch.
			[index_skeleton_so_far,used_so_far] = ...
				add_to_branching_skeleton_from(childish_contour,index,index_skeleton_so_far,axis_set,used_so_far);
		else
		   % childish contour connects wrong, and needs to be broken up
		  new_piece1 = childish_contour(1:intersection_point_location);
		  new_piece2 = childish_contour(intersection_point_location:n_points_childish); 
		  [index_skeleton_so_far,used_so_far] = ....
			  add_to_branching_skeleton_from(new_piece1,index,index_skeleton_so_far,axis_set,used_so_far);
		  [index_skeleton_so_far,used_so_far] = ...
			  add_to_branching_skeleton_from(new_piece2,index,index_skeleton_so_far,axis_set,used_so_far);
		end;
	end;
end;
index_skeleton = index_skeleton_so_far;		

% %......................
% function neighbors = find_neighbors(index,neighbor_relations,parent_list)
% % collects everything that is a neighbor of index, and that isn't anywhere on the parent list
% % Includes only possible children, ie no stemming-into
% if ~isempty(neighbor_relations)
%     neighbors = ...
%         setdiff(...
%                 neighbor_relations(neighbor_relations(:,1)==index,2),parent_list(:));
% else 
%     neighbors = [];
% end;

%...............................................
function neighbors = find_neighbors2(point_set,axis_set,exclude)
% Given set of points (indexes), choose those axes that are neighbors, ie SHARE A POINT with it.
% Exclude those that are in "exclude". Typically these are ones we've used before and don't want to use again.
neighbors = [];
for a=1:length(axis_set)
	if length(intersect(point_set,axis_set(a).contour)) == 1 & not(exclude(a))  % == 1
		neighbors = [neighbors a];
	end;
end;
neighbors = unique(neighbors);

%.....................................................................................................................
function root = choose_root(axis_set)
% Choose an axis to be the "root". 
% Simply chooses longest. (Should choose based on graph centrality, but NYI.)

% axis set has fields axis_set(i).contour
contour_lengths = [];
for i= 1:length(axis_set)
	contour_lengths = [contour_lengths length(axis_set(i).contour)];
end;
[longest_length index] = max(contour_lengths);
root = index; 
	

%.....................................................................................................................
function skeleton = indexes2skeleton(points,indexskeleton)
% convert index skeleton to regular skeleton -- indices point into points
for a=1:length(indexskeleton)
    skeleton(a).index =  indexskeleton(a).index;
    skeleton(a).parent =  indexskeleton(a).parent;
    for v=1:length(indexskeleton(a).contour)
        skeleton(a).contour(v,1:2) = points(indexskeleton(a).contour(v),:);
    end;
end;
       
%.....................................................................................................................
function indexskeleton = add_to_skeleton(a_start,v_start,this_point,last_point,neighbors,indexskeleton_so_far,all_points)
% Recursively make a skeleton by adding THIS_POINT to the skeleton_so_far at axis a_start, vertex v_start
% last_point is the last point added before, which the function uses to figure out which way we're going
% returns skeleton and a,v being first "empty slot", ie last used axis with 
% Strategy:
% if this_point has 1 neighbor different from last_point, continue the axis
% if it has 2 or more, branch
% NOTE here we are dealing with indexskeletons, not regular skeletons

global counter;
global points_used_so_far;  % keeps track so we don't go in circles

points_used_so_far = [points_used_so_far this_point]; 
if isempty(last_point)  % no last point, ie this is start; we assume this_point is a terminator
    indexskeleton = []; % initialize
    indexskeleton(1).contour(1) = this_point;
    indexskeleton(1).index      = 1;
    indexskeleton(1).parent = -1;
    this_neighbor_set = neighbors{this_point}; 
    %... and continue from there:
    if length(this_neighbor_set)==1
        next_point = this_neighbor_set; 
    else
        error('Bad starting point for skeleton');
    end;
    indexskeleton=add_to_skeleton(1,2,next_point,this_point,neighbors,indexskeleton,all_points);
else % not the starting point:  
	counter = counter + 1;
	if counter > length(all_points)
		fprintf('Exceeded maximum possible recursions!');
		return;
	end;
	this_neighbor_set = setdiff(neighbors{this_point},last_point);  % this point's neighbors, excluding last_point 
	this_neighbor_set = setdiff(this_neighbor_set,points_used_so_far); % only points that haven't been
    switch length(this_neighbor_set)
        case 0 % terminator -- add this point to the end
                indexskeleton = indexskeleton_so_far;  % establish most of skel
                indexskeleton(a_start).contour(v_start) = this_point; % then add one point    
        case 1 % regular point -- continue axis
                indexskeleton = indexskeleton_so_far;
                indexskeleton(a_start).contour(v_start) = this_point;
                % and recursively continue skeleton:
                last_point = this_point;
                next_point = this_neighbor_set; % only 1 point
                % and recursively continue skeleton:
                indexskeleton = add_to_skeleton(a_start,v_start+1,next_point,this_point,neighbors,indexskeleton,all_points);
        otherwise  % branching point; we continue on 1st neighbor and start new branches on the rest.
            indexskeleton = indexskeleton_so_far;  % establish most of skel
            indexskeleton(a_start).contour(v_start) = this_point;  % add one point
            %indexskeleton = add_to_skeleton(a_start,v_start+1,this_neighbor_set(1),this_point,neighbors,indexskeleton,all_points);
            % Next, branch on other neighbors
            for j=1:length(this_neighbor_set)
                axis_next = length(indexskeleton)+1;
                indexskeleton(axis_next).contour(1) = this_point;
                indexskeleton(axis_next).parent = a_start;
                indexskeleton(axis_next).index  = axis_next;
                indexskeleton = add_to_skeleton(axis_next,2,this_neighbor_set(j),this_point,neighbors,indexskeleton,all_points); 
            end;
    end; %switch
end; % if/else
%.....................................................................................................................
function new_skeleton = clean_up_indexskeleton(old_skeleton)
new_skeleton = old_skeleton(1); % first axis preserved
for a=2:length(old_skeleton)
    this_parent = old_skeleton(a).parent;
    if ismember(old_skeleton(a).contour(1),old_skeleton(this_parent).contour)
        first_step = 2;
    else 
        first_step = 1;
    end;
    last_step = length(old_skeleton(a).contour);
    new_skeleton(a).contour = old_skeleton(a).contour(first_step:last_step);
    new_skeleton(a).parent = this_parent;
    new_skeleton(a).index =  old_skeleton(a).index;
end;
    
%.....................................................................................................................
function axis_index = find_axis_containing(indexskeleton_so_far,this_point)
axis_index = 0;
for i=1:length(indexskeleton_so_far)
     if ismember(this_point,indexskeleton_so_far(i).contour) axis_index = i; end;
end; % now axis_index contains index of axis that contains this_point
%.....................................................................................................................
function  is_neighbor = is_a_neighbor_of_p(q,p,points,polygons)
% yields true if q is a neighbor of p (p and q are indices into points, and neighbor relations are derived from polygons)
% Give a vector of size matching 
if isempty(q)   % q is the vector of points (indices) we're testing
    is_neighbor =  [];
else
    is_neighbor = [neighbor(p,q(1),point,polygons) is_a_neighbor_of_p(q(2:end),p,points,polygons)];  % apply test to 1st element, concat to results of applying it to rest
end;
%.....................................................................................................................
function neighbors = neighbors_of(p,polygons,inside_indices)  
% find neighbors of p, given polygons; use inside_indices to filter out exterior points
neighbors = [];
for i = 1:length(polygons)
    thispoly = polygons{i};
    index_of_p = find(thispoly == p); % index of p within this polygon -- or empty
    if ~isempty(index_of_p) % p is in poly-- we're assuming only once
        % set left_nieghbor
        if index_of_p > 1
            left_neighbor = thispoly(index_of_p-1);
        else
            if index_of_p == 1
                left_neighbor = thispoly(end);
            else
                left_neighbor = [];
            end

        end;
        % set right_neighbor
        if index_of_p < length(thispoly)
            right_neighbor = thispoly(index_of_p+1);
        else
            if index_of_p == length(thispoly)
                right_neighbor = thispoly(1);
            else
                right_neighbor = [];
            end
        end;
        
      
        if ~ismember(left_neighbor,inside_indices)  
            left_neighbor = []; 
        end; % remove left neighbor if it's outside
        
        if ~ismember(right_neighbor,inside_indices) 
            right_neighbor = []; 
        end; % remove right neighbor if it's outside
        
        neighbors = unique([neighbors left_neighbor right_neighbor]); % add neighbors of p, if they exist (and are inside), removing duplicates
    end; % if p has neighbors in this poly
end; %looping through polygons


%================================================================================================================
function adjoinments = make_adjoinments(raw_index_skeleton,all_points);
global already_used_indices;
already_used_indices = [];
%display(raw_index_skeleton);
% Procedure:
% 1. Go through all intersections, for each deciding which pair (if any) to join.
% 2. Join together axes according to above. For each joined axis. Port parentage over accordinly, joining parentage where appropriate.
% 3. Decide which axis should be root, based on a length-weighted measure of graph centrality, subject to the constraint that root must have no parents. 
%    (note parentage is not totally forced up to now-- only certain joints have forced parentage)
% 4. Starting from root, recursively make skeleton
% ....
% Raw_index_skeleton is the "list" of axial segments. This has parentage but we are going to ignore it, except to find intersections.
% Compute normals:
for a=1:length(raw_index_skeleton)
    this_axis = [];
    this_axis_indices = raw_index_skeleton(a).contour;
    for v=1:length(this_axis_indices)
        point_index = raw_index_skeleton(a).contour(v);
        this_axis = [this_axis; all_points(point_index,:)];
    end;
    this_axis_normals = jfnormals(this_axis(:,1),this_axis(:,2));
    raw_index_skeleton(a).normals = this_axis_normals;  % load up raw_index_skeleton with normals
end;
% 1. Find adjoinments. Put pairs of raw indices to be joined in curly-brace-list adjoinments{}.
adjoinments{1} = []; % init
for axis=1:length(raw_index_skeleton)
        children = [];
        for possible_child = 1:length(raw_index_skeleton)
            if raw_index_skeleton(possible_child).parent == axis
                children = [children possible_child];
            end;
        end; % now children are children of axis
    to_be_joined = determine_joint([axis children],raw_index_skeleton,all_points); % put the axis along with its children, and figure out which ones to paste together.
    adjoinments = add_to_joint_masterlist(adjoinments,to_be_joined); % adds it to master, figuring out whehter there is an existing chain to add it to
end; % now adjoinments contains all joints, but not the singletons (that are not mentioned in any joints), so add them:
used_so_far = [];
for adj = 1:length(adjoinments);
    used_so_far = [used_so_far adjoinments{adj}];
end; % Now used so far contains all indices used so far
for a=1:length(raw_index_skeleton)
    if ~ismember(a,used_so_far)
        top_adjoinments = length(adjoinments);
        adjoinments{top_adjoinments + 1} = [a];
    end;
end;  
% Now, clean out duplicates (successively repeated points) from the adjoinments.
for a = 1:length(adjoinments)
	adjoinment = adjoinments{a};
	cleaned_adjoinment = remove_duplicates(adjoinment);
	adjoinments{a} = cleaned_adjoinment;
end;


%....................................................................................................
function adjoinments = add_to_joint_masterlist(adjoinments,to_be_joined)
% add indices p/q to list of adjoinments. The trick is one of them may already be there, in which case add it to the existing set at one end.
% Note the list adjoinments{i} will be listed in connecting order, but each individual index may still need to be "flipped"; this will
% happen in merge_adjoinments
global already_used_indices;
switch length(to_be_joined)
	case 0 %............ 
		% do nothing
	case 1 %............ 
%if length(to_be_joined) == 1
    if ~ismember(to_be_joined,already_used_indices)  % if already used, we do nothing;
        %next_index = length(adjoinments) + 1;
        if isempty(adjoinments{1})
            next_index = 1;   %   
        else
            next_index = length(adjoinments)+1;
        end;
        adjoinments{next_index} = to_be_joined;
    end;  %............ 
	case 2
    p = to_be_joined(1);
    q = to_be_joined(2);
    
    dispensed_with = 0;
    for i = 1:length(adjoinments)
        this_set = adjoinments{i};
        if ~isempty(this_set)
            if this_set(1) == p % pq fits at the beginning
                adjoinments{i} = [q this_set];
                dispensed_with = 1;
            elseif this_set(1) == q
                adjoinments{i} = [p this_set];
                dispensed_with = 1;
            elseif this_set(end) == p % pq fits at the end
                adjoinments{i} = [this_set q];
                dispensed_with = 1;
            elseif this_set(end) == q 
                adjoinments{i} = [this_set p];
                dispensed_with = 1;
            end; % list of ways of joining pq to an existing set  
        end;
    end;
    if isempty(adjoinments{i})
        next_index = i;
    else
        next_index = i+1;
    end;
    if ~dispensed_with % pq wasn't joined to any existing axis group
        adjoinments{next_index} = to_be_joined; % put pq to a new slot
    end;
end; 
already_used_indices = [already_used_indices to_be_joined]; % now that they have been dealt with, add them to persistent list of ones used


%........
function smoother_index_axis_set = merge_adjoinments(raw_index_skeleton,adjoinments,all_points); 
% given adjointments{}, each field of which is a set of raw axis indices, make a the actual index axes.
% smoother_index_axis_set is the set of axes (in indices)
% neighbor_relations is an nx2 matrix of unordered pairwise neighbor relations among them, which later will turn into parentage
translation_table = []; % first column gives raw_index, second column gives that index's new name in the smoother_index_axis_set



for i=1:length(adjoinments)
    this_merged_axis = [];
    this_axis_index_set = adjoinments{i};
    % procedure for adding the vertices from the first axis component is special: we have to look ahead to the second one to get it flipped right.
    if length(this_axis_index_set) > 1 % for multi-axis collections, we need to correctly flip the first axis;
        first_point_of_this_axis = all_points(raw_index_skeleton(this_axis_index_set(1)).contour(1),:);
        last_point_of_this_axis = all_points(raw_index_skeleton(this_axis_index_set(1)).contour(end),:);
        first_point_of_second_axis = all_points(raw_index_skeleton(this_axis_index_set(2)).contour(1),:);
        last_point_of_second_axis = all_points(raw_index_skeleton(this_axis_index_set(2)).contour(end),:);
        distance_from_first_point = ...
           min([norm(first_point_of_this_axis - first_point_of_second_axis) ...
                norm(first_point_of_this_axis - last_point_of_second_axis)]);
        distance_from_last_point = ...
           min([norm(last_point_of_this_axis - first_point_of_second_axis) ...
                norm(last_point_of_this_axis - last_point_of_second_axis)]);
        if  distance_from_last_point < distance_from_first_point % add points IN ORDER
            start_vertex = 1;
            end_vertex = length(raw_index_skeleton(this_axis_index_set(1)).contour);
            %step = 1;
        else % add points BACKWARDS
            start_vertex = length(raw_index_skeleton(this_axis_index_set(1)).contour);
            end_vertex = 1;
            %step = -1;
        end;
        step = sign(end_vertex - start_vertex);
        % .. now add the first axis
        for v  = start_vertex:step:end_vertex % step through the contour, in the correct order (up or down the index numbers)
            this_merged_axis = [this_merged_axis raw_index_skeleton(this_axis_index_set(1)).contour(v)];
        end;   
        previous_point = all_points(this_merged_axis(end),:);
    else % this is the case of 1-index collections; just add the index in normal (arbitrary) order
        for v  = 1:length(raw_index_skeleton(this_axis_index_set).contour); % step through the contour
            this_merged_axis = [this_merged_axis raw_index_skeleton(this_axis_index_set).contour(v)];
        end;  
    end;
    for axis_component = 2:length(this_axis_index_set); % now handle the rest of the indices, if there are more than 1
        % add the next index in this set, figuring out whether to flip it or not wrt last one.
        index_axis_to_add = this_axis_index_set(axis_component);
        first_point = all_points(raw_index_skeleton(index_axis_to_add).contour(1),:);
        last_point =  all_points(raw_index_skeleton(index_axis_to_add).contour(end),:);
        % Now determine the order in which to add vertices 
        if norm(first_point-previous_point) < norm(last_point-previous_point) % add them in up order
            start_index = 1;  % Note though that sometimes the first point will be redundant with last point of previous axis.
			                  % We'll remove such repeated points later (at the end of make_adjoinments).
            end_index   = length(raw_index_skeleton(index_axis_to_add).contour);
            %step = 1;
        else
            top_step = length(raw_index_skeleton(index_axis_to_add).contour);
            if isequal(last_point,previous_point)
                start_index = top_step - 1;
            else
                start_index = top_step;
            end;        
            end_index = 2;
            %step = -1;
        end;
        step = sign(end_index - start_index);
        if step == 0 
            step = 1; 
        end;
        % .. Now actually add the contour 
        for v  = start_index:step:end_index % step through the contour, in the correct order (up or down the index numbers)
            this_merged_axis = [this_merged_axis raw_index_skeleton(index_axis_to_add).contour(v)];
            p = all_points(raw_index_skeleton(index_axis_to_add).contour(v),:);
            %scatter(p(1),p(2)); drawnow;
        end;
        previous_point = all_points(raw_index_skeleton(index_axis_to_add).contour(v),:); % last point added
    end;
    smoother_index_axis_set(i).contour = this_merged_axis;
end;
%.. Next, make neighbor_relations -- will contain i j if:
%    i and j  contain pieces that were neighbors in the original adjacency matrix, and
%    i doesn't stem into j.
neighbor_relations = []; % determine relationships among the indices i in smoother_index_axis_set
for i = 1:length(smoother_index_axis_set)  % smoother_index_axis_set(i) is the i-th joined axis, corresponding to the raw axes in adjoinments{i}
    for j = 1:length(smoother_index_axis_set)
        if i ~= j && ij_should_relate(raw_index_skeleton,adjoinments,i,j)  
			% i_should_relate means i contains a raw axis piece that is a neighbor of
			% one of j's raw axis pieces.
			if ~i_under_j(i,j,smoother_index_axis_set,all_points) % if i ISN'T stemming into J given the actual axes
				neighbor_relations = [neighbor_relations; i j];
            end;
        end;
    end;
end;
%............................................................
function inferiors = find_inferiors(axis,axis_set,neighbor_relations)
% returns indices of the axes that are immediate inferiors of axis; their neighbors; and their neighbors' neighbors; etc.
immediate_inferiors = [];
for i = 1:size(axis_set)
    if i_under_j(i,axis) % 
        immediate_inferiors = [immediate_inferiors i]; % add i to immediate inferiors
    end;
end;


%............................................................
function should_be_neighbors = ij_should_relate(raw_index_skeleton,adjoinments,i,j)
% Returns true if the i-th axis has a component that is parent or child of the j-th component in raw_index_skeleton
i_indices = adjoinments{i};
j_indices = adjoinments{j};
should_be_neighbors = false;
for a=1:length(raw_index_skeleton)
    index = raw_index_skeleton(a).index;
    parent = raw_index_skeleton(a).parent;
    if (ismember(index,i_indices) & ismember(parent,j_indices)) | (ismember(index,j_indices) & ismember(parent,i_indices))  
		%if i ~= j
			should_be_neighbors = true;
		%end;
    end;
end;

%..................................................................................................
function [i_stems_into_j,j_touch_index] = i_under_j(i,j,smoother_index_axis_skeleton,all_points)
% Return i_stems_into_j=true if one of i's endpoints is nearest to the INTERIOR of j (rather than one if its endpoints)
% Also optionally reports the actual index j_touch_index, which is sometimes needed
% (This rules out i being the parent of j.)
contour_i = smoother_index_axis_skeleton(i).contour;
contour_j = smoother_index_axis_skeleton(j).contour;
% We have previously determined that i and j are "neighbors", meaning that they are composed of pieces some of whome
% were neighbors in the original voronoi set. This means that SOME part of i touches SOME part of j. We first find those parts.
% First, we figure out how far tip and tail of i are from j, and where on j they touch
size_i = length(contour_i);   % this is the length of j-th contour-- we want to know if i meets j at its tip, tail (this index), or otherwise.
size_j = length(contour_j);   % this is the length of j-th contour-- we want to know if i meets j at its tip, tail (this index), or otherwise.

% First, figure out which END of i is closest to somewhere on j
i_tip = all_points(contour_i(1),:);  % 
i_tail = all_points(contour_i(end),:);
i_tip_distances = [];  % will hold distances from i's tip to all locations on j
i_tail_distances = []; % will hold distances from i's tail to all locations on j
for jj=1:size_j
		i_tip_distances =  [i_tip_distances norm(i_tip -   all_points(contour_j(jj),:))]; % distances from i-tip to EACH point in j
		i_tail_distances = [i_tail_distances norm(i_tail - all_points(contour_j(jj),:))]; % distances from i-tail to EACH point in j
end;
[i_tip_distance,i_tip_index]   = min(i_tip_distances); % shortest path from i-tip to j[tip_distance,tip_index]
[i_tail_distance,i_tail_index] = min(i_tail_distances);
% Next, figure out which END of j is closest to somewhere on i
j_tip = all_points(contour_j(1),:);  % 
j_tail = all_points(contour_j(end),:);
j_tip_distances = [];  % will hold distances from i's tip to all locations on j
j_tail_distances = []; % will hold distances from i's tail to all locations on j
for ii=1:size_i
		j_tip_distances =  [j_tip_distances norm(j_tip -   all_points(contour_i(ii),:))]; % distances from j-tip to EACH point in i
		j_tail_distances = [j_tail_distances norm(j_tail - all_points(contour_i(ii),:))]; % distances from j-tail to EACH point in i
end;
[j_tip_distance,j_tip_index]   = min(j_tip_distances); % shortest path from i-tip to j
[j_tail_distance,j_tail_index] = min(j_tail_distances);
% Now, figure out whether i's tip/tail is closer to j; or j's tip/tail is closer to i
%index_from_which_i_is_closer_to_j = 
if min([i_tip_distance,i_tail_distance]) > min([j_tip_distance,j_tail_distance]) %j's tip/tail is closer to i then viceversa, so i can't stem into j.
	i_stems_into_j = false;
else % else 
	% then i's tip or tail is closer to j, so potentially i stems into j. Whether it does depends on whether it is the tip or tail...
	if i_tip_distance < i_tail_distance
		j_touch_index = i_tip_index;
	else
		j_touch_index = i_tail_index;
	end;
	i_stems_into_j =  ... % Finally, i stems into j if the place where i touches j is not at one of j's endpoints.
		j_touch_index > 1 & j_touch_index < size_j;
end;
%j_touch_index is position on contour j at which contour i touches it.

%..................................................................................................
function root_index = determine_root2(axis_set,neighbor_relations,raw_index_skeleton);
% Choose a root index, by
%  1. choosing appropriate candidates-- can't be stems of H's; 
%  2. Among candidates, choose the one that maximizes the product of graph centrality and length

% For each axis, find all inferiors. This includes all axes that are stems or (recursive) neighbors of stems.
% Any axis that is on the list of inferiors for any axis is ruled out as a candidate root.
% After all axes have been considered, candidate roots are those that remain (if any). 
%
% Not currently used. 
if isempty(neighbor_relations)
    root_index = 1;
else
n_axes = length(axis_set);
candidate_status = ones(1,n_axes); % start by assuming they are all ok, then start ruling out
stems = [];
for a=1:n_axes
    for i=1:n_axes
        inferiors_of_a = neighbor_relations(neighbor_relations(:,1)==a,2);
        superiors_of_a = neighbor_relations(neighbor_relations(:,2)==a,1);
        missing_partners = setdiff(inferiors_of_a,superiors_of_a)'; % set of inferiors that are not also superiors
    end;
    stems = [stems missing_partners];
end; % Now, stems is the full set of axes that are someone's immediate inferior
ruleouts = stems; % initialize ruleout list
done = false;
ruleouts_old = unique(stems);
% Now expand the list of "stems". From here on, a stems is something that is a stem of a stem.
while ~done
    %old_ruleouts = ruleouts;
	for a = 1:length(ruleouts_old)% for each current stem, add its immediate neighbors (stems and edge-ons)
		this_stem = ruleouts_old(a);
		%left_neighbors_of_a = neighbor_relations(neighbor_relations(:,1)==a,2);
		right_neighbors_of_a = neighbor_relations(neighbor_relations(:,1)==a,2);
		ruleouts = unique([ruleouts right_neighbors_of_a']);
	end;
    done = isequal(ruleouts,ruleouts_old);  % if new ruleouts is same as old ruleouts, we're done.
	ruleouts_old = ruleouts;
end;
%... Next, choose longest of non-stems
candidate_list = setdiff(1:n_axes,ruleouts);
axis_lengths = zeros(size(candidate_list));
    %  else choose the longest
    for a = 1:length(candidate_list)
        axis_lengths(a) = length(axis_set(a).contour);
    end;
    [value,root_a] = max(axis_lengths);
    root_index = candidate_list(root_a);
end; 

%..................................................................................................
function root_index = determine_root(smoother_index_axis_set,neighbor_relations); 
candidate_list = [];
for a=1:length(smoother_index_axis_set)
    this_contour = smoother_index_axis_set(a).contour;
    % determine whether a stems into anyone
    inferiors_of_a = neighbor_relations(neighbor_relations(:,1)==a,2);
    superiors_of_a = neighbor_relations(neighbor_relations(:,2)==a,1);
    % a stems into someone if it has a superior that is not also an inferior:
    missing_partners = setdiff(superiors_of_a,inferiors_of_a);
    %a_stems_into_someone = ~isempty(missing_partners);
    if isempty(missing_partners);  % if no missing partners--- ie no superiors
        candidate_list = [candidate_list a];
    end;
end;
% .. Next, among candidates, choose the shortest.
if isempty(candidate_list)
    fprintf('WARNING: Bad skeleton -- no root axis found');
    root_index = 1;
else
    axis_lengths = zeros(size(candidate_list));
    %  else choose the longest
    for a = 1:length(candidate_list)
        axis_lengths(a) = length(smoother_index_axis_set(a).contour);
    end;
    [value,root_a] = max(axis_lengths);
    root_index = candidate_list(root_a);
end;

%..................................................................................................
function new_skeleton = reconstruct_skeleton(root_index,neighbor_relations,smoother_index_skeleton_set,all_points)
% Put the merged axes into a skeleton, with parentages determined wrt root
global current_axis_hierarchy;
% First, recursively make parent fields
parentage = collect_parents(root_index,neighbor_relations,[]);
% and add them to new_skeleton
for a=1:size(parentage,1)
    index =  parentage(a,1);
    theparent = parentage(a,2);
	new_skeleton(index).parent = theparent;
end;
new_skeleton(root_index).parent = -1;  
% Next, make contours, being sure to flip each axis in the right direction
for a=1:length(smoother_index_skeleton_set)
    parent_axis = [];
    if new_skeleton(a).parent >= 1 % parent is not -1
        parent_axis = smoother_index_skeleton_set(new_skeleton(a).parent).contour;
    end;
    d_fronts = [];
    d_backs = [];
    for v_parent = 1:length(parent_axis)
        d_fronts = [d_fronts norm(all_points(smoother_index_skeleton_set(a).contour(1),:) - all_points(smoother_index_skeleton_set(new_skeleton(a).parent).contour(v_parent),:))];
        d_backs =  [d_backs norm(all_points(smoother_index_skeleton_set(a).contour(end),:) - all_points(smoother_index_skeleton_set(new_skeleton(a).parent).contour(v_parent),:))];
    end;
    if min(d_fronts) < min(d_backs) | isempty(parent_axis) % add forwards
        new_skeleton(a).contour = smoother_index_skeleton_set(a).contour;
    else  % add backwards
        new_skeleton(a).contour = smoother_index_skeleton_set(a).contour(end:-1:1);
    end;
    new_skeleton(a).index = a;
end;
% Finally, make (global) current_axis_hierarchy
axis_hierarchy = [];
for a=1:length(new_skeleton)
	if ~isempty(new_skeleton(a).parent)
		% Only add a line for this index if it HAS a parent. This means that orphan axes are being dropped from the skeleton. No choice.
		axis_hierarchy = [axis_hierarchy; new_skeleton(a).index new_skeleton(a).parent];
	end;
end;
current_axis_hierarchy = axis_hierarchy;

%......................
function parent_list = collect_parents(index,neighbor_relations,parent_list)
% parentlist of the form [index parent; etc]
neighbors = find_neighbors2(index,neighbor_relations,parent_list);
if ~isempty(neighbors)
    for n = 1:length(neighbors)
        this_child_index = neighbors(n);
        parent_list = [parent_list; this_child_index index]; % indicates this_child has index as parent
        % Next, recurse on this child
        parent_list = collect_parents(this_child_index,neighbor_relations,parent_list);
    end; 
end;



%............................................................
function list_wo = remove_duplicates(list);
% Removes successive duplicate members without otherwise changing the order of the list
previous = list(1);
list_wo = [list(1)];
for i=2:length(list)
	if list(i) ~= previous
		list_wo  = [list_wo list(i)];
		previous = list(i);
	end;
end;


function to_be_joined = determine_joint(indices,raw_index_skeleton,all_points)
% Given axes given by (raw) indices, decide which ones are collinear enough to join.
% Returns a set of 2 or 0 (if none are sufficiently collinear)
% First, determine which end of each index axis is actually the one at the intersection.
% We assume they all intersect, so we find the intersection of the first two, then for each of the others
% choose the end closest to the intersection. Note the "intersection" is just a point closest to one end of each; 
% We don't assume they necessarily actually share a common point.

collinearity_threshold_deg = 60; % degrees
collinearity_threshold = collinearity_threshold_deg * pi/180; % radians

if length(indices) < 2
    to_be_joined = indices;
else % normal case:
% First, find the intersection of the first two.
%     a_first = all_points(raw_index_skeleton(indices(1)).contour(1),:);
%     a_last  = all_points(raw_index_skeleton(indices(1)).contour(end),:);
%     b_first = all_points(raw_index_skeleton(indices(2)).contour(1),:);
%     b_last  = all_points(raw_index_skeleton(indices(2)).contour(end),:);
%     distances = [norm(a_first - b_first) norm(a_first - b_last) norm(a_last - b_first) norm(a_last - b_last)];
%     [min_distance,winning_index] = min(distances);
%     start_indices = [0 0];
    parent_index = indices(1); % First index was the parent, as this is how the list was assembled; rest are chidren.
    
    first_child_index = indices(2);
    intersection_point_index = intersect(raw_index_skeleton(parent_index).contour,raw_index_skeleton(first_child_index).contour);
    % ----
    % Next, create a set of vectors, giving the direction of each axis at the
    % intersection point. By convention, each points AWAY from the (common)
    % intersection point. 
    % Note we don't treat the parent differently from the children, because
    % in principle we might decide to join two of the children without the
    % parent. In that case, later, the two children will become the parent
    % (otherwise it's an "H"), but that will be decided later. 

    vectors = []; %  One for each indix in indices;
    for i=1:length(indices)
        this_contour = raw_index_skeleton(indices(i)).contour;
        if isequal(intersection_point_index,this_contour(1)); % beginning
            v = all_points(this_contour(2),:) - all_points(this_contour(1),:);  
            
        
        else                                                  % end
            v = all_points(this_contour(end-1),:) - all_points(this_contour(end),:);
         

        end;
        v = v/norm(v);
        vectors = [vectors; v];
       
        
    end;  % Now vectors is a "start" with vectors pointing away from the intersection point. 
    
     % Next, determine collinearities and choose most collinear. 
     % Note all vectors are unit vectors.
     % Criteron: we look for the MAX of dot(v1,-v2). This will be 1
     % whendraw_in
     % they two vectors are collinear, but pointing away from each other. 
     
     collinearities = [];
     for i=1:length(indices)
        for j = 1:length(indices)
            if i < j
                this_dot = dot(vectors(i,:),-vectors(j,:));
                collinearities = [collinearities ; i j this_dot];
            end;
        end;
     end;
         
     [winning_dot,winning_index] = max(collinearities(:,3)); % collinearites, ie dot products, are in third column . Max is the one that is "straightest". 
    if winning_dot > cos(collinearity_threshold) % dot of least straight angle that is straight "enough"
        to_be_joined = indices(collinearities(winning_index,1:2)); % i and j columns of the winning row.
    else
		to_be_joined = [];
	end;  
end;


