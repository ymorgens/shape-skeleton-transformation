function newPnts_tr = tr_pnts_w_simTriROP(basePnts, oPnts, oDist, newDist);

% similar triangles, ratio of parts.
% 
% basePnts = skeleton_point;
% 
% oPnts = shape_point;
% 
% oDist = init_dist;
% 
% newDist = new_dist;

deltaPts = oPnts - basePnts;

ratio = newDist./oDist;
newPnts = bsxfun(@times,ratio,deltaPts);
%newPnts = ratio.*deltaPts;

newPnts_tr = newPnts+basePnts;

% can do this with solving quadratic too (See below)

%             X = [ shape_point(1) skeleton_point(1)]';
%             Y = [ shape_point(2) skeleton_point(2)]';
%             B = regress(Y,[repmat(1,[length(X) 1]) X]);


%  x2 = skeleton_point(1);
%             y2 = skeleton_point(2);
%             
%             % finding new point by solving quadratic equation (can change
%             % to use triangle similarity properties);
%             p1 = 1+ B(2).^2;
%             p2 = -2*x2 + 2*B(2)*[B(1)-y2];
%             p3 = -(new_dist.^2 - x2.^2 - [B(1)-y2].^2);
%             xnews = roots([p1 p2 p3]);
%             ynews = B(1) + B(2).*xnews;
%             
%             % quadratic gives two points on either side of skeleton
%             % pick new shape point that is closer to shape_point
%             d = pdist2(shape_point, [xnews ynews]);
%             new_id = find(d==min(d));
%            trans_shape_point = [xnews(new_id) ynews(new_id)];