function coribs = compute_coribs(skeleton,index_range)
% coribs = compute_coribs(skeleton,index_range)
%
% Compute the ribs (sometimes called coribs) between given skeleton and
% (global) current_shape. 
%
% Attempts to find an explanation for every contour point, possibly mixed from 
% multiple axis points, and some explanation coming from both directions from every axis
% point -- basically maximizing along both rows and columns of large affinity matrix.
%
% Optional argument index_range:
% If this is absent, does the whole shape.
% If it is present, finds coribs for ONLY these index points of current_shape.
% This is useful for when we are pruning, when we only want to add
% explanations for certain points while leaving the rest untouched. 
%
% Output syntax:
% Coribs output is a matrix with each row of the form
% c a v s
%  Where c is a contour point, a v is an axis point, s is 1 or -1 (the direction of the normal to which c attaches)
% (So there are 2n axis locations at which to attach, where n is the number of actual axis points).
% Multiple lines can have the same c, and multiple lines can have the same (a v s).
% Later we will compute description lengh by adding up all the explanation of a given point c. 
%
% "Explanation" can be one-sided or two-sided. 
% One-sided means that each shape point must be explained by one axis-vertex-sign-- period.
% Two-sided means that, in addition, each axis-vertex that explains in ONE direction (sign) must also
% explain in the OTHER direction (sign).

% The major constraints are that:
% 1. Every c is explained by at least one a v s
% 2. Every a v s explains at least one c.
% Note this enforces 2-sideness.

% translation_table is of the form: translation_table(av,:) = 
% [a v sign]  
% sign is 1 or -1

%-------------------------------------------------------------------------------------------------------------
global current_shape;
global current_shape_normals;
global distance_decay;
global require_twosided_explanation;
global vonMises_b;
global current_coribs;  % Might be used by another program. 

logconstant = log(2*pi*besseli(0,vonMises_b));
if nargin == 2
	partial_explanation = true;  % flag that says we are only explaining part of the shape-- needed later.
 else
        partial_explanation = false;
end;

if nargin == 1
	shape_to_explain = current_shape(1:end,:); % normal case 
	index_range = 1:size(shape_to_explain,1);  % This is the set of indices we are trying to explain-- here simply 1 through n shape points (-1)
else
	shape_to_explain = current_shape(index_range,:);
	% Here, index_range is the one given to the function-- not all points
end;
	
%-------------------------------------------------------------------------------------------------------------
coribs = []; 
% make translation table: --- this puts all the axis points into one long matrix, with 2n rows (2 for each axis point) 
translation_table = []; % will contain TWO rows per axis point, one left and one right
axis_normal_table = []; % will contain TWO rows per axis point, ie one per axis location, containing info about that location
for a=1:length(skeleton)
    this_axis = skeleton(a).contour;
        %... First, establish normals for this axis.........................
        if size(this_axis,1) > 2 
         normals_this_axis = jfnormals(this_axis(:,1),this_axis(:,2));
        else
            normals_this_axis =  zeros(size(this_axis));
        end;
        %.......Normals are now in normals_this_axis. ...............................
    for v = 1:size(skeleton(a).contour,1)
        translation_table = [translation_table; a v 1; a v -1]; % make two rows of translation table
        axis_normal_table = [axis_normal_table; this_axis(v,:) normals_this_axis(v,:); this_axis(v,:) -normals_this_axis(v,:)]; % make two rows of normals table
    end;
end;
n_points_skeleton = size(translation_table,1); % total "locations" in skeleton (2n)
n_points_shape =    size(shape_to_explain,1);  % either entire current_shape or just part of it.
likelihood_matrix = zeros(n_points_shape,n_points_skeleton); % initialize likelihood matrix
for i=1:n_points_shape
	c = index_range(i); % the real contour index of the point are we currently explaining;
    for av = 1:n_points_skeleton
        axis_point = axis_normal_table(av,1:2); % first two fields of axis_normal_table
        axis_normal = axis_normal_table(av,3:4); % 3rd and 4th fields
        shape_point = shape_to_explain(i,:); % this is the c-th point on the ENTIRE shape, but the i-th point in shape_to_explain
        shape_normal = current_shape_normals(c,:);
        dl_due_to_distance =  distance_decay * norm(shape_point - axis_point);
        this_affinity = affinity(axis_point,axis_normal,shape_point,shape_normal,vonMises_b,logconstant);  % Note here we use a different (higher) value
        %this_affinity = 0;
        likelihood_matrix(i,av) = this_affinity + dl_due_to_distance; 
    end;
end;
% Now likelihood_matrix contains the whole grid of point-to-axis affinities-- ie how well each point (row) is explaine by each axis point (column)
%------------------------------------------------------------
[likelihoods,winning_axis_points] = min(likelihood_matrix,[],2);
[likelihoods,winning_contour_points] = min(likelihood_matrix,[],1);
coribs = ones(n_points_shape,4); % initialize. Columns are: c a v s
for i=1:n_points_shape % for each contour point, add winning axis point
    %winning_axis_points(c) is winner for this point, expressed as an av location
	actual_shape_vertex = index_range(i);
    coribs(i,:) = [actual_shape_vertex translation_table(winning_axis_points(i),:)];
%    coribs = [coribs; translation_table(winning_axis_points_right(i),:); translation_table(winning_axis_points_left(i),:)];
end;

% require_twosided_explanation is a global flag, that, if on,
% forces every axis point that explains ANY shape points to explain one from the other side as well. 
if require_twosided_explanation & ~partial_explanation  % if partial explanation is on, we IGNORE the requirement of two-sided explanation.
	singleton_mates = make_singleton_set(winning_axis_points); % returns the OPPOSITE av's of the singletons
	for i=1:length(singleton_mates) % all axis points that were a winner already for SOME contour point
		av = singleton_mates(i);
		coribs = [coribs ; winning_contour_points(av) translation_table(av,:)];   % add one rib
	end;
end;
current_coribs = coribs;		

%------------------------------------------------------------
function singletons_mates = make_singleton_set(winning_axis_points) % returns the OPPOSITE av's of the singletons
% returns OPPOSITE locations from points that won, but have empty opposites
singletons_mates = [];
for i=1:length(winning_axis_points)
    if winning_axis_points(i)/2==round(winning_axis_points(i)/2); % if even
        opposite = winning_axis_points(i)-1; % if even, opposite guy is 1 before
    else
        opposite = winning_axis_points(i)+1;  % if odd, opposite is 1 after
    end;
    if ~ismember(opposite,winning_axis_points)
        singletons_mates = [singletons_mates opposite];
    end;
end;
