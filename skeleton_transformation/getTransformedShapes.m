function [tshapes, winning_skeleton, winning_coribs] = getTransformedShapes(shape, trans, sdparms, nsamps,shape_blur_factor)
% getTransformedShapes Transform a shape's skeletal representation
%   tshapes = getTransformedShapes(shape, trans, sdparms, nsamps,shape_blur_factor)
%   tshapes is a structure with transfromed shapes and their transfromed skeletal structures
%   trans is how the subpart is transformed (1 is length, 2 is position, 3
%   is angle, and 4 is width)
%   sdparms is deviation parameter that controls how much to tranform the
%   subpart. sdparms is 1x4, where sdparms(1) is the length deviation, sdparms(2) is position deviation,
%   sdparms(3) is angle deviation, and sdparms(4) is width deviation
%   nsamps is the number of samples to produce
%   shape_blur_factor (1 = no blur; higher values show greater blur)

if nargin<3
    nsamps = 1;
    shape_blur_factor = 2;
elseif nargin<4
    shape_blur_factor = 2;
end

length_sig = sdparms(1);
pos_sig = sdparms(2);
theta_sig = sdparms(3);
width_sig = sdparms(4);

% fit shape skeleton to shape
winning_skeleton = mapskeleton(shape);
winning_coribs = compute_coribs(winning_skeleton);

% get skeletal summary.
% 1. branch lengths
for j = 1:length(winning_skeleton)
    % calculate length of skeletal branch
    c = winning_skeleton(j).contour;
    c_dist(j) = curvedist(c);
end

% 2. get other stats
[p skelfacts] = getSkelShapeParms_Wilderetal2011(winning_skeleton);

% structure of skelfacts
% skelfacts.parent = hierarchal position of stem (0 is grandparent; 1 =
% parent; 2 = grandparent)

% skelfacts.theta = branch angle; the angle at which each (non-root) axis
% branches from its parent)

% skelfacts.p_len = find the distance along each parent axis at which each
% child axis stems

% skelfacts.p_nlen = find the relative distance along each parent axis at
% which each child stems

% skelfacts.curvedist_rel_rt = length of axis relative to root (i.e.,
% grandfather)



%% synthesise

% set up initial synthesis parameters

% which factors to transform?  (0 is OFF, 1 is ON)
parms.scale = 1; % width of ribs
parms.spatialpos = 1; % position along root contour
parms.missing_skel_part = 0; % missing skeleton part
parms.orient = 1; % orientation

%  shape and matching skeleton
parms.parent = skelfacts.parent;
parms.winning_skeleton = winning_skeleton;
parms.winning_coribs = winning_coribs;
parms.original_shape = shape;

%% synthesize shape

j = 2:length(winning_skeleton); % j are the subparts of the skeleton that are changed, where 1 is the main branch

% get sythesis parameters by adding a bit of noise to the original
% parameters.  Noise levels determined by variables above (length_sig
% pos_sig theta_sig width_sig)
for n = 1:nsamps
    int = 1; % shape has self intersection
    while any(int)
        % for length
        sc_dist = c_dist;% skeletal length of each branch
        if any(ismember(trans,[1]))
            sc_dist(j) = sc_dist(j).*abs(normrnd(ones(1,length(j)), length_sig));
        end
        scale_val = sc_dist./c_dist;
        
        % for spatial position (in this case only change rib width of branch j)
        sp_npos = skelfacts.p_nlen; % spatial position on main branch
        if any(ismember(trans,[2]))
            sp_npos(j) = normrnd(sp_npos(j), pos_sig);
        end
        
        % for angle
        sthetas = skelfacts.thetaa; % angle relative to main branch
        if any(ismember(trans,[3]))
            sthetas(j) = normrnd(skelfacts.thetaa(j), theta_sig);
        end
        
        % for width
        rib_scale = ones(length(scale_val));%  width of branch
        if any(ismember(trans,[4]))
            rib_scale(j) = rib_scale(j).*abs(normrnd(ones(1,length(j)), width_sig));
            
        end
        
        % put the new parameters in parms structure
        parms.scale_val = scale_val; % scale of skeletal branch size relative to original
        parms.rib_scale = rib_scale;
        parms.sp_npos = sp_npos; % new child branch location on parent
        parms.sthetas = sthetas; % new angle relative to parent
        parms.missing_skel_part = 0; % missing skeleton part
        
        % get new skeleton
        [new_skeleton, parms] = gen_new_skeleton(winning_skeleton, parms);
        
        
        
        % pass new parms to get new_shape
        [new_shape, new_coribs] = get_new_shape_coordinates(new_skeleton, parms);
        smooth_new_shape = smooth_shape(new_shape,shape_blur_factor,'gaussian');
        
        
        % check for self-intersections
        % better self intersection script available in Matlab file exchange from geom2d toolbox: polygonSelfIntersections
        [X0,Y0,SEGMENTS] = selfintersect(smooth_new_shape(:,1),smooth_new_shape(:,2));
        int = any(X0);
    end
    tshapes(n).smooth_new_shape = smooth_new_shape;
    tshapes(n).new_shape = new_shape;
    tshapes(n).new_skeleton = new_skeleton;
    tshapes(n).new_coribs = new_coribs;
end