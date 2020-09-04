function [new_shape, new_coribs] = get_new_shape_coordinates(new_skeleton, parms)    

% sets coribs from original skeleton to new skelton with proper scale and
% relative angular position as original
% new shape then becomes the points that connect the new corib ends.

winning_skeleton = parms.winning_skeleton;
winning_coribs  =  parms.winning_coribs;
original_shape =parms.original_shape;
riblength = parms.rib_scale; %rand(length(new_skeleton)); % varies the size of the rib length along each skeletal index
missing_skel = parms.missing_skel_part;
count = 1;
for j=1:length(original_shape);

%     
    if j ~= length(original_shape)
        pt_id = find(winning_coribs(:,1)==j,1,'first');
        % pt_id = pt_id(1);
    else
        pt_id = find(winning_coribs(:,1)==1,1,'first');
       % pt_id = pt_id(1);
    end
    %pt_id = j;
    
    % get original shape points
    
    % corib on winning_skeleton
    c = winning_coribs(pt_id,1);
    a = winning_coribs(pt_id,2);
    v = winning_coribs(pt_id,3);
    s = winning_coribs(pt_id,4);
    if ~ismember(a,missing_skel) % remove skeletal branch 'a'' from new skeleton
        
        skeleton_point = winning_skeleton(a).contour(v,:);
        shape_point = original_shape(c,:);
        
        
        % **need to transform shape point so that it's distance from skelepoint is scaled by mult
        init_dist = pdist2(shape_point,skeleton_point);
        
        new_dist = riblength(a)*parms.scale_val(a)*init_dist;
        
        new_skeleton_point = new_skeleton(a).contour(v,:);
      
        % compute corib on new skeleton
        if v==1 | v==length(winning_skeleton(a).contour)% for first and last point on contour, use translation
            
            trans_shape_point =  tr_pnts_w_simTriROP(skeleton_point, shape_point, init_dist, new_dist); % scale shape point by mult
            
            tr = skeleton_point - new_skeleton_point; % find where to translate new shape point
            % trans_shape_point = shape_point;%parms.mult(a).*shape_point;% + parms.delta(a);
            new_shape_point = trans_shape_point - tr; % translate new shape point
        else % for inner contour points, make sure relative angle between rib and contour are the same for original and new contour
            
            % [for original shape] compute vectors for local points around current skeletal point
            % for the original shape
            p1 = winning_skeleton(a).contour(v-1,:);
            p3 = winning_skeleton(a).contour(v+1,:);
            vskel = shape_point - skeleton_point; uvskel =  vskel./sqrt(sum(vskel.^2));
            v1 = p1 - skeleton_point; uv1 = v1./sqrt(sum(v1.^2));
            v2 = p3 - skeleton_point; uv2 = v2./sqrt(sum(v2.^2));
            
            % [for original shape] compute angles between vector going towards rib and vectors going to local neighbouring points
            theta1 = acosd(dot(uv1,uvskel));
            theta2 = acosd(dot(uv2, uvskel));
            theta_sum = theta1+theta2;
            angrat = theta1./(theta_sum);
            
            % [for new shape] compute vectors in neighbouring directions around transformed skeletal point
            pn1 = new_skeleton(a).contour(v-1,:);
            pn3 = new_skeleton(a).contour(v+1,:);
            vn1 = pn1 - new_skeleton_point; uvn1 = vn1./sqrt(sum(vn1.^2));
            vn2 = pn3 - new_skeleton_point; uvn2 = vn2./sqrt(sum(vn2.^2));
            
            % find angle between these vectors
            if issame(uvn1,uvn2)
                new_theta_sum = 0;
            else
            new_theta_sum = acosd(dot(uvn1,uvn2));
            end
            % make sure that the angle chosen between new_theta1 and
            % new_theta2 is the same as theta1 and theta2
            if (new_theta_sum>180 & theta_sum<180) | (new_theta_sum<180 & theta_sum>180)
                new_theta_sum = 360 - new_theta_sum;
            end
            
            % estimate scaled angle between vector pointing to neighbouring point
            % and vector pointing towards shape point
            new_theta1 = new_theta_sum*angrat;
            
            
            % compute vector in the direction of angle new_theta1
            CO = cosd(new_theta1); SI = sind(new_theta1);
            % there are two possible vector solutions
            R1 = [CO -SI;SI CO]; R2 = [CO SI;-SI CO];
            vnskel1 = R1*vn1'; vnskel2 = R2*vn1';
            uvnskel1 = [vnskel1./sqrt(sum(vnskel1.^2))]';
            uvnskel2 = [vnskel2./sqrt(sum(vnskel2.^2))]';
            
            % select solution closest to uvskel
            if pdist2(uvnskel1, uvskel)<pdist2(uvnskel2,uvskel)
                uvnskel = uvnskel1;
            else
                uvnskel = uvnskel2;
            end
            
            % calculate shape point from new_skelton point given above data
            new_vector_point = new_skeleton_point+100*uvnskel;
            new_vector_point_dist = pdist2(new_vector_point,new_skeleton_point);
            new_shape_point =  tr_pnts_w_simTriROP(new_skeleton_point, new_vector_point, new_vector_point_dist, new_dist);
            
        end
        new_shape(count,:) = new_shape_point;
        
          new_coribs(count,:) = [count a v s];
        count = count+1;
        
        
%         if a == 4
%             subplot(1,2,1);
%             plot(winning_skeleton(a).contour(:,1), winning_skeleton(a).contour(:,2),'k'); hold on;
%             plot(skeleton_point(1), skeleton_point(2),'k>'); hold on;
%             plot(original_shape(1:count-1,1), original_shape(1:count-1,2),'bx'); hold on;
%             plot(shape_point(1), shape_point(2),'ro');
%             pdist2(shape_point, skeleton_point)
%             subplot(1,2,2);
%             plot(new_skeleton(a).contour(:,1), new_skeleton(a).contour(:,2),'k');hold on;
%             plot(new_skeleton_point(1), new_skeleton_point(2),'k>'); hold on;
%             
%             plot(new_shape(1:count-1,1),new_shape(1:count-1,2),'bx'); hold on;
%             plot(new_shape_point(1),new_shape_point(2),'ro');
%             pdist2(new_shape_point, new_skeleton_point)
%            pdist2(uvnskel1,uvskel)
%             pdist2(uvnskel2,uvskel)
%          %  close all;
%         end
    end
   % new_shape(j,:) = new_shape_point;

end;