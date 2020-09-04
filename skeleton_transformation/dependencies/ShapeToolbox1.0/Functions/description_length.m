function [dl_total,augmented_skeleton] = description_length(skeleton)
% [dl_total,augmented_skeleton] = description_length(skeleton)

% Compute the description length (-log probability) of the given skeleton
% with respect to (global) current_shape. 
%
% Optional output augmented_skeleton adds certain fields that are used
% internally in shape_tool GUI. 

global current_coribs;  % For use if we don't recompute them every opt step
global current_shape; 
global distance_decay;
global mean_rib_lengths;  % This is global so it can be  used by prototype_shape
huge_penalty = 10000; % imposed if axis point is outside shape
endpoint_penalty = exp(15);
burbeck_fraction = 0.1;
%---------------------------------------
% Install shape normals
shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));

%---------------------------------------------------------------------------------------------
% First, compute rib length matrix for each axis
for a = 1:length(skeleton)
    %... While we're here, install normals for this axis
    this_contour = skeleton(a).contour;
    axis_normals = jfnormals(this_contour(:,1),this_contour(:,2));
    skeleton(a).axisnormals = axis_normals; % tack the normals on "skeleton" so it's available the rest of this function
	this_axis_riblengths = [];
    %.... Now, collect the ribs connecting to for each point on the shape, axis by axis; store mean for each axis point a v
	for v=1:size(this_contour,1) % for each axis point in this_contour
		coribs_this_axis_point = current_coribs(current_coribs(:,2) == a & current_coribs(:,3) == v,:);
		rib_lengths_this_axispoint = [];
		if ~isempty(coribs_this_axis_point)
			for j = 1:size(coribs_this_axis_point,1)
				this_corib = coribs_this_axis_point(j,:);  % [c a v s] [contour_point axis_index vertex_on_axis sign]
				c = this_corib(1);
				%a = this_corib(2);
				%v = this_corib(3);
				%s = this_corib(4);
				axis_point = skeleton(a).contour(v,:);
				shape_point = current_shape(c,:);
				this_distance = norm(axis_point - shape_point);
				rib_lengths_this_axispoint = [rib_lengths_this_axispoint this_distance]; % add this distance to the vector of riblenghts for this axis point
			end;
		else
			rib_lengths_this_axispoint = [0]; % no ribs from this axis point-- assign 0 as placeholder
		end;
		this_mean_riblength = mean(rib_lengths_this_axispoint);
		this_axis_riblengths = [this_axis_riblengths; this_mean_riblength];
	end;
	skeleton(a).axisriblengths = this_axis_riblengths; % add new field to skeleton called "axisriblengths"
	% Now, skeleton is loaded with vectors of axisriblengths (mean at each axis point);
end;
%---------------------------------------------------------------------------------------------
% Next, install SMOOTHED riblengths
for a=1:length(skeleton)
	this_axis_contour = skeleton(a).contour;
	this_axis_means  = [];
	for v=1:size(this_axis_contour,1)
		this_mean = local_average_riblength(skeleton,a,v);
		this_axis_means = [this_axis_means; this_mean];
	end;
	skeleton(a).smoothedaxisriblengths = this_axis_means; % installs a new field on skeleton
end;  
		    
%---------------------------------------------------------------------------------------------
% Next, compute DL(Shape|Skel)
likelihood_vector = []; % will contain one likelihood per shapepoint
for c = 1:size(current_shape,1)-1  % for each shape point
    coribs_this_shapepoint = current_coribs(current_coribs(:,1) == c,:); % pulls out just the rows that attach to this shape point
    n_coribs_this_shapepoint = size(coribs_this_shapepoint,1); % number of "explanations" for this shape point; we will weight them all equally
    affinities_this_shapepoint = zeros(1,n_coribs_this_shapepoint); % initialize
    riblengthlikelihoods_this_shapepoint = zeros(1,n_coribs_this_shapepoint);
    for i= 1:n_coribs_this_shapepoint % for all coribs attached to this shape point
        this_corib = coribs_this_shapepoint(i,:); % [c a v s]
        %... Collect the necessary parameters:
        c = this_corib(1);
        a = this_corib(2);
        v = this_corib(3);
        s = this_corib(4);
        axis_point = skeleton(a).contour(v,:);
        axis_normal = s * skeleton(a).axisnormals(v,:);
        shape_point = current_shape(c,:);
        shape_normal = shape_normals(c,:);
        %.. compute affinity;
        affinities_this_shapepoint(i) = affinity(axis_point,axis_normal,shape_point,shape_normal);
            %if isempty(mean_rib_lengths) printf('NO RIBS'); end;
		local_mean = local_average_riblength(skeleton,a,v); % gives the mean in the K-neighborhood of axis a, vertex v, using axisriblengths. 
        sigma = 1;
        %sigma = burbeck_fraction * local_mean;
		
		normal_constant = -log(1/(sigma * sqrt(2 * pi)));
		if inpolygon(axis_point(1),axis_point(2),current_shape(:,1),current_shape(:,2)) % if the axis point is INSIDE the shape, as we want
			riblengthlikelihoods_this_shapepoint(i) = normal_constant + ((norm(axis_point - shape_point) - local_mean)/sigma)^2;  
		else % impose a huge penalty
			riblengthlikelihoods_this_shapepoint(i) = huge_penalty;
		end;

    end;
	% Now combine affinities and riblikelihoods by probability weighting:
	summed_likelihoods_this_shapepoint = affinities_this_shapepoint + riblengthlikelihoods_this_shapepoint; % Summing loglikehoods, ie multiplying probs, ie assuming independence
	ps_this_shapepoint = exp(-summed_likelihoods_this_shapepoint);
	if sum(ps_this_shapepoint) == 0% very unlikely point
		normalized_ps_this_shapepoint = huge_penalty * ones(1,n_coribs_this_shapepoint);
	else
		normalized_ps_this_shapepoint = ps_this_shapepoint / sum(ps_this_shapepoint);
	end;% normalize so they sum to one like probabilities
	% Finally, "probability-weight" the likelihoods. Big ones, ie low probability, get low weights.	
	likelihood_this_shapepoint = ...
		   dot(normalized_ps_this_shapepoint,affinities_this_shapepoint) ...
		 + dot(normalized_ps_this_shapepoint,riblengthlikelihoods_this_shapepoint);
    likelihood_vector = [likelihood_vector likelihood_this_shapepoint];
end;
dl_likelihood = sum(likelihood_vector);
%---------------------------------------------------------------------------------------------
% % Next, compute DL(Skel)
dl_skeleton = skeleton_complexity(skeleton);
%---------------------------------------------------------------------------------------------
dl_total = dl_skeleton + dl_likelihood + distance_decay * sum(mean_rib_lengths);
augmented_skeleton = skeleton; % this is (optionally) returned as well; it is the orignal plus new field smoothedaxisriblengths

%---------------------------------------------------------------------------------------------
% Function that computes mean in a neighborhood near an axis point
function mean_rib = local_average_riblength(skeleton,a,v)
K = 15; % continuity mask width
riblengths = skeleton(a).axisriblengths;
n = length(riblengths);
start_vertex = max([1 v - round(K/2)]); % either v - halfmask, or 0, whichever is higher)
end_vertex   = min([n v + round(K/2)]); 
mean_rib = mean(riblengths(start_vertex:end_vertex));

 



