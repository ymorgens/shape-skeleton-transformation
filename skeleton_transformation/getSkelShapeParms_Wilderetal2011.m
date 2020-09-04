function [parms, skelfacts] = getSkelShapeParms_Wilderetal2011(skeleton)

% get 8 shape parameters per skeleton.  The parameters are listed in
% Wilder, Feldman, and Singh (2011) pg. 327 and depicted in their Figure 3.


% parms(1) =  total number of skeletal branches
% parms(2) =  maximum depth of skeleton
% parms(3) =  mean skeletal depth
% parms(4) mean branch angle, i.e. the angle at which each (non-root) axis
% brances from its parent (See Fig. 3);
% parms(5) mean distance along each parent axis at which each child axis stems
% (see Fig. 3);
% parms(6) mean length of each axis relative to the root
% parms(7) total absolute (unsigned) turning angle integrated along the curve
% (see Fig. 3)
% parms(8) total (absolute value of) signed turning angle of each axis,
% integrated along the curve (See Fig. 3)

% parms(9) std of branch angles in parms(4) - not from wilder et al.
% parms(10) std of of distances from parms(5) - not from wilder et al.
% parms(11) std of of distances from parms(6) - not from wilder et al.




% close all;
% clear all;
% clc;
% 
% % set paths
% p = setpaths('skeletal2DShapeTransformation');
% 
% % intialize shape toolbox
% initializeShapeToolbox;
% 
% 
% 
% 
% load('gaparms_crit_mean','new_obj','parms');
% 
% % get and show skeleton for some obj
% i = 1;
% c = 375*(new_obj{i}.shape+0.5);
% skeleton = mapskeleton(c);%300*(c+0.5));
% coribs = compute_coribs(skeleton);
% 
% for j = 1:length(skeleton)
%     transformed_skeleton(j).contour = (skeleton(j).contour/375)-0.5;    
%     skeleton_parent(j) = skeleton(j).parent;
% end

% (1) total number of skeletal branches
parms(1) = length(skeleton); % total number of skeletal branches

% get depth of skeleton
if parms(1)>1
% count = 0;
% parent(1) = count;


% for i = 2:length(skeleton)
%     v = skeleton(i).parent;
%     u = skeleton(i-1).parent;
%     
%    if u<v %& ~ismember(parent,u) % &~ismember(u,v)
%        count = count +1;
%        parent(i) = count;
%    elseif u==v
%        parent(i) = count;
%    elseif u>v
%        if v == 1
%            parent(i) = 1;
%        else
%        count = count - 1;
%        parent(i) = count;
%        end
%    end
% end
skeleton_parent_mat = SkeletonParent2mat(skeleton);
parent = zeros(1, length(skeleton));
parent(1) = 0;

id = find(skeleton_parent_mat==1);
if any(id)
parent(id) = 1;
end

id = find(skeleton_parent_mat>1);

if any(id)
    for i = 1:length(id)
     p_id = find(skeleton_parent_mat(id(i)) == skeleton_parent_mat); 
    parent(id(i)) = parent(p_id(1)-1)+1;
    
    end
end

skelfacts.parent = parent; 

% (2) maximum depth of skeleton
parms(2) = length(unique(parent)); % maximum depth of skeleton
else
    parent = 1;
    parms(2) = 1;
end
% (3) % mean skeletal depth
parms(3) = mean(parent); % mean skeletal depth


% (4) mean branch angle, i.e. the angle at which each (non-root) axis
% brances from its parent (See Fig. 3);


%clear theta;
%clear thetaa;
%clear thetab;
if parms(1)>1
for i = 2:length(skeleton)
    
    % get branch contour
    bc = skeleton(i).contour;
    
    % get branch connecting point to stem
    bcp1 = bc(1,:);
    
    % get stem contour
    h = parent(i);
    ph = h-1;
    ids = find(parent==ph);
    id = find(ids<i);
    id = ids(id(end));
    pc = skeleton(id).contour;
    
    % find where branch intersects stem contour
    id = find(sum(ismember(pc,bcp1),2)==2);
    id = id(1);
    
    % compute unit vector on branch
    bcp2 = bc(2,:);
    bv = bcp2 - bcp1;
    ubv = bv./sqrt(sum(bv.^2));
    
    % compute unit vector on stem
    
    % pt on first side of contour
    if id>=length(pc)
        id = length(pc)-2;
    end
    
    if id==1
        id = 2;
    end
    
    pcp2a = pc(id+1,:);
    pva = pcp2a - bcp1;
    upca = pva./sqrt(sum(pva.^2));
    
    % pt on second side of contour
    pcp2b = pc(id-1,:);
    pvb = pcp2b - bcp1;
    upcb = pvb./sqrt(sum(pvb.^2));
    
    % compute angles between vectors
    thetaa(i-1) = acosd(dot(upca,ubv));
    thetab(i-1) = acosd(dot(upcb,ubv));
    
    theta(i-1) = min(thetaa(i-1), thetab(i-1));
    
%     figure;
%     subplot(1,2,1);
%     % plot branch and stem contours
%     plot(bc(:,1),bc(:,2),'r') ; hold on;
%     plot(pc(:,1),pc(:,2),'g')
%     plot(pc(1,1),pc(1,2),'go')
%     
%     % plot stem/branch intersection
%     plot(pc(id,1),pc(id,2),'ko')
%     plot(bcp1(1),bcp1(2),'ko')
%     
%     % plot next points
%     plot(bcp2(1),bcp2(2),'ko')
%     plot(pcp2a(1),pcp2a(2),'ko')
%     plot(pcp2b(1),pcp2b(2),'ko')
%     
%     % plot lines showing vectors
%     plot([bcp1(1) bcp2(1)], [bcp1(2) bcp2(2)],'k');
%     plot([bcp1(1) pcp2a(1)], [bcp1(2) pcp2a(2)],'k');
%     plot([bcp1(1) pcp2b(1)], [bcp1(2) pcp2b(2)],'k');
%     
%     title(theta(i-1))
%     axis equal tight;
%     % axis square;
%     grid on;
%     subplot(1,2,2);
%     plot(bc(:,1),bc(:,2),'r') ; hold on;
%     plot(pc(:,1),pc(:,2),'g');
%     axis equal tight;
%     grid on;
%     pause
%     close all;
    
end


parms(4) = rad2deg(circ_mean(deg2rad(theta')));
parms(9) = rad2deg(circ_std(deg2rad(theta'))); % std of branch angles in parms(4)
skelfacts.mintheta = [NaN theta];
skelfacts.thetaa = [NaN thetaa];
skelfacts.thetab = [NaN thetab];
else
    parms(4) = NaN;
    parms(9) = NaN;
end

% (5) mean distance along each parent axis at whcih each child axis stems
% (see Fig. 3);

distfun = @(x,y) sum(sqrt(sum((x-y).^2,2)));
 if parms(1)>1 
for i = 2:length(skeleton)
    
     % get branch contour
    bc = skeleton(i).contour;
    
    % get branch connecting point to stem
    bcp1 = bc(1,:);
    
    % get stem contour
    h = parent(i);
    ph = h-1;
    ids = find(parent==ph);
    id = find(ids<i);
    id = ids(id(end));
    pc = skeleton(id).contour;
    pc_ulength = pc./curvedist(pc); % length of contour normalized to 1
    
     % find where branch intersects stem contour
    id = find(sum(ismember(pc,bcp1),2)==2);
    
    % calulate length of branch by summing the distances between the points
    clear len;
    if id>1
    for j = 1:id-1
        pt1 = pc(j,:);
        pt2 = pc(j+1,:);
        len(j) = distfun(pt1, pt2);
    end 
    p_len(i-1) = sum(len(:));
     p_nlen(i-1) = curvedist(pc_ulength(1:id,:)); % find the relative distance on parent
    else
        p_len(i-1) = 0;
         p_nlen(i-1) = 0; % find the relative distance on parent
    end
    
   
    
%     % plot branch and stem contours
%     subplot(1,2,1);
%     draw_shape(c);
%   hold on;
%     plot(bc(:,1),bc(:,2),'r') ; hold on;     
%     plot(pc(:,1),pc(:,2),'g')
%     plot(pc(1:id,1),pc(1:id,2),'go')
%     axis equal tight;
%     title(p_len(i-1))
%     
%     subplot(1,2,2);
%     draw_shape(c);
%     draw_skeleton(skeleton);
%     axis equal tight;
%  pause;
% close all;
end
skelfacts.p_len = [NaN p_len];
skelfacts.p_nlen = [NaN p_nlen];

% (5) mean distance along each parent axis at whcih each child axis stems
% (see Fig. 3);
parms(5) = mean(p_len);
parms(10) = std(p_len);

% find mean length for stems from each level of the hierarchy.
h_l = unique(parent);
clear hl_len;
for i = 2:length(h_l)
    id = find(parent==h_l(i));
    hl_len(i-1) = mean(p_len(id-1));
end

 else
     parms(5) = NaN;
     parms(10) = NaN;
 end


% (6) mean length of each axis relative to the root
if parms(1)>1
rt_c = skeleton(1).contour;
clear len;
for j = 1:(size(rt_c,1)-1);
    pt1 = rt_c(j,:);
    pt2 = rt_c(j+1,:);
    len(j) = distfun(pt1, pt2);
end
rt_dist = sum(len);

for i = 2:length(skeleton)
    
     % get branch contour
    bc = skeleton(i).contour;
    
    % calulate length of branch by summing the distances between the points
    clear len;
    for j = 1:(size(bc,1)-1)
        pt1 = bc(j,:);
        pt2 = bc(j+1,:);
        len(j) = distfun(pt1, pt2);
    end 
    axis_dist = sum(len(:));
    dist_rel_rt(i-1) = axis_dist./rt_dist;

end
skelfacts.curvedist_rel_rt = [NaN dist_rel_rt];
% (6) mean length of each axis relative to the root
parms(6) = mean(dist_rel_rt);
parms(11) = std(dist_rel_rt);
else
    parms(6) = NaN;
    parms(11) = NaN;
end
% (7) total absolute (unsigned) turning angle integrated along the curve
% (see Fig. 3)

% (8) total (absolute value of) signed turning angle of each axis,
% integrated along the curve (See Fig. 3)
angle = @(x,y) mod(atan2(x(1)*y(2)-x(2)*y(1),x(1)*x(2)+y(1)*y(2)), 2*pi);
for i = 1:length(skeleton)
    
    % get branch contour
    con = skeleton(i).contour;
    
    for j = 1:(size(con,1)-2)
        pt1 = con(j,:);
        pt2 = con(j+1,:);
        pt3 = con(j+2,:);
        
        % compute unit vector from first point to second point
        v1 = pt2 - pt1;
        uv1 = v1./sqrt(sum(v1.^2));
        
        % compute unit vector froms second point to third point
        v2 = pt3 - pt2;
        uv2 = v2./sqrt(sum(v2.^2));
        
       % theta(j) = acosd(dot(uv1,uv2));
        t = rad2deg(angle([uv1(1) uv2(1)],[uv1(2) uv2(2)]));
        theta(j) = mod(t+180,360) - 180;
        
%         %figure;
%         subplot(1,2,1);
%         % plot contours and points
%         plot(con(:,1),con(:,2),'k') ; hold on;
%         plot(pt1(1),pt1(2),'go')
%         plot(pt2(1),pt2(2),'ro')
%         plot(pt3(1),pt3(2),'bo')
%         
%         % plot lines showing vectors
%         plot([pt1(1) pt2(1)], [pt1(2) pt2(2)],'c');
%         plot([pt2(1) pt3(1)], [pt2(2) pt3(2)],'y');
%         
%         hold on;
%        % draw_shape(c);
%         title(theta(j))
%         
%         axis equal tight;
%         
%         
%         subplot(1,2,2);
%                 % plot points
%      hold off;
%         plot(pt1(1),pt1(2),'go'); hold on;
%         plot(pt2(1),pt2(2),'ro')
%         plot(pt3(1),pt3(2),'bo')
%         
%         % plot lines showing vectors
%         plot([pt1(1) pt2(1)], [pt1(2) pt2(2)],'c');
%         plot([pt2(1) pt3(1)], [pt2(2) pt3(2)],'y');
%         axis equal tight;
%        % title(theta2(j));
%         pause
%        % close all;
        
    end
    if (size(con,1)-2)>0
        p7(i) =  sum(abs(theta));
        p8(i) = sum(theta);
        
    else
        p7(i) = NaN;
        p8(i) = NaN;
    end
end

% (7) total absolute (unsigned) turning angle integrated along the curve
% (see Fig. 3)
parms(7) = nansum(p7);

% (8) total (absolute value of) signed turning angle of each axis,
% integrated along the curve (See Fig. 3)
parms(8) = nansum(p8);



