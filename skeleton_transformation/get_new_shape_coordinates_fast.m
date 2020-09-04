function [new_shape] = get_new_shape_coordinates_fast(new_skeleton, parms)    

% sets coribs from original skeleton to new skelton with proper scale and
% relative angular position as original
% new shape then becomes the points that connect the new corib ends.
% "fast" is faster by replacing loop between object points to loop between object
% contours by taking advantage of matrix multiplication.

winning_skeleton = parms.winning_skeleton;
winning_coribs  =  parms.winning_coribs;
original_shape =parms.original_shape;
riblength = parms.rib_scale; %rand(length(new_skeleton)); % varies the size of the rib length along each skeletal index
missing_skel = parms.missing_skel_part;

[C,pt_ids,ic] = unique(winning_coribs(:,1),'rows','sorted');
distfun = @(x,y) sqrt(sum((x-y).^2,2));
% corib on winning_skeleton
cs = winning_coribs(pt_ids,1);
as = winning_coribs(pt_ids,2);
vs = winning_coribs(pt_ids,3);
%ss = winning_coribs(pt_ids,4);
as_un = unique(as);
new_shape = zeros(size(original_shape));
%new_coribs = zeros(size(winning_coribs));
for con_num = 1:length(as_un);
    a = as_un(con_num);
    id = find(ismember(as, a));
    v = vs(id);
    c = cs(id);
    
    
    skeleton_points = winning_skeleton(a).contour(v,:);
    shape_points = original_shape(c,:);
    
    % **need to transform shape point so that it's distance from skelepoint is scaled by mult
    
    % for first and last point on contour, use translation
    init_dist = distfun(shape_points,skeleton_points);
    new_dist = riblength(a)*parms.scale_val(a)*init_dist;
    
    new_skeleton_points = new_skeleton(a).contour(v,:);
    
    v_id =  [find(v==1); find(v==length(winning_skeleton(a).contour))];
    
    trans_shape_points =  tr_pnts_w_simTriROP(skeleton_points(v_id,:), shape_points(v_id,:), init_dist(v_id), new_dist(v_id)); % scale shape point by mult
    
    tr = skeleton_points(v_id,:) - new_skeleton_points(v_id,:); % find where to translate new shape point
    new_shape_points = bsxfun(@minus, trans_shape_points,tr);% trans_shape_points - tr; % translate new shape point
    new_shape(id(v_id),:) = new_shape_points;
    
    % for inner contour points, make sure relative angle between rib and contour are the same for original and new contour
    v_id = find(ismember(sum([double(v==1)  double(v==length(winning_skeleton(a).contour))],2),0));
    
    % [for original shape] compute vectors for local points around current skeletal point
    % for the original shape
    
    p1 = winning_skeleton(a).contour(v(v_id)-1,:);
    p3 = winning_skeleton(a).contour(v(v_id)+1,:);
    
    vskel = shape_points(v_id,:) - skeleton_points(v_id,:);
    uvskel = bsxfun(@rdivide, vskel, sqrt(sum(vskel.^2,2)));%vskel./sqrt(sum(vskel.^2,2));
    v1 = p1 - skeleton_points(v_id,:); uv1= bsxfun(@rdivide, v1,sqrt(sum(v1.^2,2)));
    v2 = p3 - skeleton_points(v_id,:);  uv2= bsxfun(@rdivide, v2,sqrt(sum(v2.^2,2)));
    
    % [for original shape] compute angles between vector going towards rib and vectors going to local neighbouring points
    theta1 = acosd(dot(uv1',uvskel'));
    theta2 = acosd(dot(uv2', uvskel'));
    theta_sum = theta1+theta2;
    angrat = theta1./(theta_sum);
    
    % [for new shape] compute vectors in neighbouring directions around transformed skeletal point
    pn1 = new_skeleton(a).contour(v(v_id)-1,:);
    pn3 = new_skeleton(a).contour(v(v_id)+1,:);
    vn1 = pn1 - new_skeleton_points(v_id,:); %uvn1 = vn1./sqrt(sum(vn1.^2));
    uvn1= bsxfun(@rdivide, vn1,sqrt(sum(vn1.^2,2)));
    
    vn2 = pn3 - new_skeleton_points(v_id,:);% uvn2 = vn2./sqrt(sum(vn2.^2));
    uvn2= bsxfun(@rdivide, vn2,sqrt(sum(vn2.^2,2)));
    
    % find angle between these vectors
    new_theta_sum = acosd(dot(uvn1',uvn2'));
    % where the angles are the same, set new_theta_sum to zero
    [C, ia, ib] = intersect(uvn2, uvn1,'rows');
    
    new_theta_sum(ia) = 0;
    
    % make sure that the angle chosen between new_theta1 and
    % new_theta2 is the same as theta1 and theta2
    if any(new_theta_sum>180 & theta_sum<180)
        
        new_theta_sum(new_theta_sum>180 & theta_sum<180) = 360 - new_theta_sum(new_theta_sum>180 & theta_sum<180);
    elseif  any(new_theta_sum<180 & theta_sum>180)
        new_theta_sum(new_theta_sum<180 & theta_sum>180) = 360 - new_theta_sum(new_theta_sum<180 & theta_sum>180);
    end
    
    % estimate scaled angle between vector pointing to neighbouring point
    % and vector pointing towards shape point
    
    new_theta1 = new_theta_sum.*angrat;
    
    % compute vector in the direction of angle new_theta1
    CO = cosd(new_theta1); SI = sind(new_theta1);
    % there are two possible vector solutions
    %R1 = [CO -SI;SI CO]; R2 = [CO SI;-SI CO];
    
    
    vnskel1 = [CO.*vn1(:,1)' - SI.*vn1(:,2)';  SI.*vn1(:,1)' + CO.*vn1(:,2)']';
    vnskel2 = [CO.*vn1(:,1)' + SI.*vn1(:,2)';  -SI.*vn1(:,1)' + CO.*vn1(:,2)']';
    
    % vnskel2 = R2*vn1';
    uvnskel1 = bsxfun(@rdivide, vnskel1, sqrt(sum(vnskel1.^2,2)));
    uvnskel2 = bsxfun(@rdivide, vnskel2, sqrt(sum(vnskel2.^2,2)));
    % uvnskel1 = [vnskel1./sqrt(sum(vnskel1.^2,2))]';
    % uvnskel2 = [vnskel2./sqrt(sum(vnskel2.^2,2))]';
    
    
    % select solution closest to uvskel
    uvnskel = zeros(size(uvnskel1));
    uvnskel([distfun(uvnskel1, uvskel)<distfun(uvnskel2,uvskel)],:) = uvnskel1([distfun(uvnskel1, uvskel)<distfun(uvnskel2,uvskel)],:);
    uvnskel([distfun(uvnskel1, uvskel)>=distfun(uvnskel2,uvskel)],:) = uvnskel2([distfun(uvnskel1, uvskel)>=distfun(uvnskel2,uvskel)],:);
    
    % calculate shape point from new_skelton point given above data
    new_vector_points = new_skeleton_points(v_id,:)+100*uvnskel;
    new_vector_points_dist = distfun(new_vector_points,new_skeleton_points(v_id,:));
    new_shape_points =  tr_pnts_w_simTriROP(new_skeleton_points(v_id,:), new_vector_points, new_vector_points_dist, new_dist(v_id));
    
    new_shape(id(v_id),:) = new_shape_points;
   % new_coribs(id(v_id),:) = [id(v_id) as(v_id) vs(v_id) ss(v_id)];
    
    
end