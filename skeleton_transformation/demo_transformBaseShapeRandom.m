%% Skeleton transformation demo

% this demo transforms the limbs of an object's shape skeleton
% here, we transform the skeleton by changing the width, angle, length, 
% and spatial position of a skeltal subpart by a small random factor, and 
% find shape with no self intersection

% Please cite the following:

% (1) Morgenstern, Y., Schmidt, F., & Fleming, R. W. (2019). One-shot categorization of novel object classes in humans. Vision research, 165, 98-108.

% This demo relies on the shape skeleton implementation by
% Feldman and Singh (2006).  The ShapeToolbox implementation
% was downloaded at http://ruccs.rutgers.edu/images/ShapeToolbox1.0.zip

% (2) Feldman, J., & Singh, M. (2006). Bayesian estimation of the shape skeleton. Proceedings of the National Academy of Sciences, 103(47), 18014-18019.

close all;
clear all;
clc;

%% initial transformation parameters
% standard deviation of values on changing the subparts, larger values = more change.
length_sig = 0.2; % length
pos_sig = 0.05; % spatial position
theta_sig = 10; % angle
width_sig = 0.2; % width

% vary the parameters above to create classes of shapes with varying
% statistics (e.g.,more variance between parts vs. less)

% choose kind of transformation for demo, in the following...
% trans = 1; % length transfromation
% trans = 2; % spatial position transform
% trans = 3; % angle transfrom
% trans = 4; % width transform

% for example,
trans = [1 2 3 4 ]; % which transformations to include. use this to include all
%trans = [ 2 3]; % use this to inlude spatial position (2) and angle (3) only

% blur shape to smooth over discontinuities from transfromation
shape_blur_factor = 2; %

%% set up paths, load example shapes, pick object
% add paths
addpath('dependencies/circStat');
addpath(genpath('dependencies/ShapeToolbox1.0'));

% intialize shape toolbox (Feldman and Singh, 2006)
initializeShapeToolbox;

% load skeletons for shapes
load('shapes/novobjs.mat','obj');

%% extract skeletal summary for object that will be transformed
id = 350; % object id (1 to 1000)
shape = obj{id}.shape;

%test = interparc(500,shape(:,1),shape(:,2)); % change shape resolution

% prefit object skeleton
%winning_skeleton = obj{id}.skeleton;
%winning_coribs = obj{id}.coribs;

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

% show shape, skeleton, ribs
h = figure('units','normalized','outerposition',[0 0 1 1]); set(h,'Color','w');

% transformed skeleton and contour
subplot(2,2,1)
draw_shape(new_shape);draw_skeleton(new_skeleton);...
    global current_shape;
current_shape = new_shape;
draw_coribs_y(new_coribs,'rainbow',new_skeleton,new_shape)
axis square;
axis tight;
hold on;
h = title('transformed skeleton'); set(h,'fontsize',30)


% original skeleton and contour
subplot(2,2,2)
draw_shape(shape);draw_skeleton(winning_skeleton);...
    global current_shape;
current_shape = shape;
draw_coribs_y(winning_coribs,'rainbow',winning_skeleton,shape)
axis square;
axis tight;
hold on;
h = title('original skeleton'); set(h,'fontsize',30)

% draw smoothed shape

% transformed shape blurred
subplot(2,2,3);
patch(smooth_new_shape(:,1),smooth_new_shape(:,2),[int*1 1 int*1]);
axis tight;
axis square;
axis off;
h = title('transformed shape - smoothed'); set(h,'fontsize',30)

% original shape blurred
subplot(2,2,4);
smooth_shape = smooth_shape(shape,shape_blur_factor,'gaussian');
patch(smooth_shape(:,1),smooth_shape(:,2),[0 1 0]);
axis tight;
axis square;
axis off;
h = title('original shape - smoothed'); set(h,'fontsize',30)
