function [new_skeleton, parms, new_shape] = gen_new_skeleton(winning_skeleton, parms)

if nargin<2
%     parms.curve = 1; % magnitude of transformation, controls curvature
%     parms.curve_val = 0.1;
    parms.scale = 1; % scale of transformation, controls length
    parms.scale_val = [0.9 1];
end


parent = parms.parent;
skelpar_index = SkeletonParent2mat(winning_skeleton);
nCounters = length(winning_skeleton);
new_shape = parms.original_shape;
missing_skel = parms.missing_skel_part;
coribs = parms.winning_coribs(1:length(parms.original_shape),:);
distfun = @(x,y) sqrt(sum((x-y).^2,2));
%[p skelfacts] = getSkelShapeParms_Wilderetal2011(winning_skeleton);

clear new_skeleton;
for i = 1:nCounters
    % if ~ismember(i,missing_skel) % remove skeletal branch 'a'' from new skeleton
    xt = winning_skeleton(i).contour(:,1);
    yt = winning_skeleton(i).contour(:,2);
    new_skeleton(i).contour = [xt yt];
    new_skeleton(i).index = winning_skeleton(i).index;
    new_skeleton(i).parent = winning_skeleton(i).parent;
    
%     if parms.curve==1
%         
%         yt = gweightedTransformContour(xt, yt, parms.curve_val(i));
%         xt = gweightedTransformContour(yt, xt, parms.curve_val(i));
%         new_skeleton(i).contour = [xt yt];
%         new_skeleton(i).index = winning_skeleton(i).index;
%         new_skeleton(i).parent = winning_skeleton(i).parent;
%     end
%     
    if parms.scale==1
        
        mult = parms.scale_val(i);% abs(parms.scale_val(1)*rand+parms.scale_val(2)); % scale factor for contour i
        xt = mult*xt;
        yt = mult*yt;
        
        % transform shape
        new_shape(coribs(coribs(:,2)==i,1),:) = [new_shape(coribs(coribs(:,2)==i,1),1)*mult ...
            new_shape(coribs(coribs(:,2)==i,1),2)*mult];
        
        
        
        %  parms.mult(i) = mult;
        if skelpar_index(i) == -1 % oldest relative, longest branch
            
            new_skeleton(i).contour = [xt yt];
            new_skeleton(i).index = winning_skeleton(i).index;
            new_skeleton(i).parent = winning_skeleton(i).parent;
            
        elseif skelpar_index(i) ~= -1; % translate new_skelton branches to approptriate location on medial axis
            attached_pts = winning_skeleton(i).contour(1,:);
            contour_id = skelpar_index(i);
            parent_id = find(skelpar_index(1:i)<contour_id,1,'last');
            parent_pts = winning_skeleton(parent_id).contour;
            intersect_id = find(sum(ismember(parent_pts, attached_pts),2)==2);
            intersect_loc = new_skeleton(parent_id).contour(intersect_id,:);
            delta = intersect_loc(1,:) - [xt(1) yt(1)];
            new_skeleton(i).contour = [xt+delta(1) yt+delta(2)];
            new_skeleton(i).index = winning_skeleton(i).index;
            new_skeleton(i).parent = winning_skeleton(i).parent;
            
            % transform shape
            new_shape(coribs(coribs(:,2)==i,1),:) = [new_shape(coribs(coribs(:,2)==i,1),1) + delta(1) ...
                new_shape(coribs(coribs(:,2)==i,1),2) + delta(2)];
            % parms.delta(i,:) = delta;
            
            
            
        end
    end
    
    
    % adjust spatial position for non-root skeletal branch
    if parms.spatialpos == 1
        if i==1
            new_root_dist = curvedist([xt yt]);
            
            %         plot(winning_skeleton(i).contour(:,1), winning_skeleton(i).contour(:,2))
            % hold on; plot(new_skeleton(i).contour(:,1), new_skeleton(i).contour(:,2),'r')
            % pause;
        end
        
        if skelpar_index(i)~=-1
            
            new_sp_nlen = parms.sp_npos(i);
            new_sp_len = new_sp_nlen*new_root_dist;
            
            contour_id = skelpar_index(i);
            parent_id = find(skelpar_index(1:i)<contour_id,1,'last');
            
            
            parent_pts = new_skeleton(parent_id).contour;
            if new_sp_len>0
                [d, new_intersect_id] = curvedist(parent_pts, new_sp_len);
            else
                new_intersect_id = 1;
            end
            intersect_loc = new_skeleton(parent_id).contour(new_intersect_id,:);
            delta = intersect_loc(1,:) - [new_skeleton(i).contour(1,1) new_skeleton(i).contour(1,2)];
            new_skeleton(i).contour = [new_skeleton(i).contour(:,1)+delta(1) new_skeleton(i).contour(:,2)+delta(2)];
            
            % transform shape
            new_shape(coribs(coribs(:,2)==i,1),:) = [new_shape(coribs(coribs(:,2)==i,1),1)+delta(1) ...
                new_shape(coribs(coribs(:,2)==i,1),2) + delta(2)];
            %         winning_skeleton(i).contour(1,:)
            %         new_skeleton(i).contour(1,:)
            %         plot(winning_skeleton(i).contour(:,1), winning_skeleton(i).contour(:,2))
            % hold on; plot(new_skeleton(i).contour(:,1), new_skeleton(i).contour(:,2),'Color', rand(1,3))
            %
            %          pause;
            %
        end
    end
    
    % adjust orientation between child and parent
    if parms.orient == 1  & parms.scale_val(i)~=0
        
        if skelpar_index(i)~=-1
            
            % get branch contour
            bc = new_skeleton(i).contour;
            
            % get branch connecting point to stem
            bcp1 = bc(1,:);
            
            % get stem contour
            h = parent(i);
            ph = h-1;
            ids = find(parent==ph);
            id = find(ids<i);
            id = ids(id(end));
            pc = new_skeleton(id).contour;
            
            % find where branch intersects stem contour
            id = find(sum(ismember(pc,bcp1),2)==2,1,'first');
            
            
            if ~any(id)
                id = bsxfun(@minus,bcp1,pc);
                id = sum(abs(id),2);
                id = find(id == min(id),1,'first');
            end
            
            % compute unit vector on branch
            bcp2 = bc(2,:);
            bv = bcp2 - bcp1;
            ubv = bv./sqrt(sum(bv.^2));
         %   if ~any(isnan(ubv)) % no vector on child branch when it is not there, so don't compute new_shape, as orientation won't matter
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
                
                %             % pt on second side of contour
                %             pcp2b = pc(id-1,:);
                %             pvb = pcp2b - bcp1;
                %             upcb = pvb./sqrt(sum(pvb.^2));
                
                % compute angles between vectors
                thetaa = acosd(dot(upca,ubv));
                % thetab = acosd(dot(upcb,ubv));
                
                % center of rotation
                center = repmat(bcp1,length(bc),1);
                
                % find angle
                thetad = thetaa - parms.sthetas(i);
                R = [cosd(thetad) -sind(thetad); sind(thetad) cosd(thetad)];
                
                % do the rotation...
                s = bc - center;     % shift points in the plane so that the center of rotation is at the origin
                so = R*s';           % apply the rotation about the origin
                vo = so' + center;   % shift again so the origin goes back to the desired center of rotation
                new_skeleton(i).contour = vo;
                
                % transform shape
                ns = new_shape(coribs(coribs(:,2)==i,1),:);
                center = repmat(bcp1,length(ns),1);
                s = new_shape(coribs(coribs(:,2)==i,1),:) - center;
                so = R*s';
                vo = so' + center;
                
                ribs = coribs(coribs(:,2)==i,:);
                %              if abs(parms.sthetas(i) - skelfacts.thetaa(i))>90
                % %                 vo = flipud(vo);
                %                  ribs = flipud(ribs);
                %              end
                
                new_shape(ribs(:,1),:) = vo;
          %  end
        end
    end
    if parms.rib_scale(i)~=1 & parms.scale_val(i)~=0
        ribs = coribs(coribs(:,2)==i,:);
        init_dist = distfun(new_shape(ribs(:,1),:), new_skeleton(i).contour(ribs(:,3),:));
    
       % new_dist = parms.rib_scale(i)*parms.scale_val(i)*init_dist;
        new_dist = parms.rib_scale(i)*parms.scale_val(i)*init_dist;
        new_shape(ribs(:,1),:) = tr_pnts_w_simTriROP(new_skeleton(i).contour(ribs(:,3),:), new_shape(ribs(:,1),:), init_dist, new_dist);
      
    end
    

end

