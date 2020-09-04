function [points,neighbor_pairs] = medial_axis(shape,varargin)
% [points,neighbor_pairs] = medial_axis(shape,graphical_output)
%  
% Compute the Medial Axis Transform (Blum, 1973) of input 'shape', based
% on the voronoi diagram (see Ogniewicz and Kubler, 1995). 
% 
% Input shape should be an n x 2 array of points, corresponding to a 2D
% shape. Output is a set of connected points. 
% 
% If optional second argument graphical_output is true (or 1), also draws the MAT.


% Note: Some of this code is derived from make_medial_axis, which gives output in a
% hierarchically organized skeleton (suitable for use in map_skeleton), using a much
% more complicated procedure designed to produce an intuitively correct
% hierarchy. 

% ----------------------------------------------------------------------------------------------------------------
global points_used_so_far;  % will contain 1s where points are used
points_used_so_far = [];
% ----------------------------------------------------------------------------------------------------------------
% If second argument is 1 (true), draw the MAT. Otherwise, don't. 
if isempty(varargin) | isequal(varargin{1},0)
    graphical_output = false;
else
    graphical_output = true;
end;
    
% ----------------------------------------------------------------------------------------------------------------
if isequal(shape(1,:),shape(end,:)) % If first and last point duplicated -- often happens with closed shapes. 
    shape = shape(1:end-1,:);
end;
% ----------------------------------------------------------------------------------------------------------------
[all_points,polygons] = voronoin(shape);  % Voronoi diagram. 
inside_indices = find(inpolygon(all_points(:,1),all_points(:,2),shape(:,1),shape(:,2)));
% NOTE: from here on, point indices will be into all_points.
terminator_list = []; % points with just one neighbor
% -----------
% Find neighbors of each point, and make note of terminators. 
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
indexskeleton_so_far = [];

% ------------------------------------------------------
% Collect the interior points for output.
points = all_points(inside_indices,:);

% ------------------------------------------------------
% Collect the neighbor relations in a more accessible format:
neighbor_pairs = [];
for i=1:size(all_points,1)
    these_neighbors = neighbors{i};
    for j=1:length(these_neighbors);
        neighbor_pairs = [neighbor_pairs; i these_neighbors(j)];
    end;
end;

% ------------------------------------------------------
%  If graphical_output is requested, draw the skeleton.
%  We draw a line segment for each neighbor edge. 
if graphical_output
    hold on
    for i=1:size(neighbor_pairs,1)
        this_pair = neighbor_pairs(i,:);
        a = all_points(this_pair(1),:);
        b = all_points(this_pair(2),:);
        plot([a(1) b(1)],[a(2) b(2)],'k','LineWidth',1); % Draw edge, black, thin line
    end;
end;       


% .....................................................................................................................
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
end; % looping through polygons
