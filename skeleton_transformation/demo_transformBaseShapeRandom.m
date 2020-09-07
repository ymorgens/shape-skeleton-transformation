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

%% set up paths, load example shapes, pick object, and other factors
% add paths
addpath('dependencies/circStat');
addpath(genpath('dependencies/ShapeToolbox1.0'));

% intialize shape toolbox (Feldman and Singh, 2006)
initializeShapeToolbox;

% set up shapes and other non transform related parameters
% load shapes and select one
load('shapes/novobjs.mat','obj');
id = 8; % object id (1 to 1000)
shape = obj{id}.shape;

shape_blur_factor = 2; % blur shape to smooth over discontinuities from transfromation
nsamps = 1; % number of samples to produce

%% initial transformation parameters

% choose kind of transformation for demo, in the following...
% trans = 1; % length transfromation
% trans = 2; % spatial position transform
% trans = 3; % angle transfrom
% trans = 4; % width transform

% for example,
trans = [1 2 3 4 ]; % which transformations to include. use this to include all
%trans = [ 2 3]; % use this to inlude spatial position (2) and angle (3) only

%% standard deviations for environment 1
% select standard deviation of values on changing the subparts, larger values = more change.
length_sig = 0.2; % length, some good starting values: 0.25:0.5:2
pos_sig = 0.05; % spatial position, some good starting values: 0:0.01:0.05
theta_sig = 10; % angle, some good starting values: 0:5:25;
width_sig = 0.2; % width, some good starting values: 0.25:0.5:2
sdparms = [length_sig, pos_sig, theta_sig, width_sig];

% get Transformed shape
[newshapes, winning_skeleton, winning_coribs] = getTransformedShapes(shape, trans, sdparms, nsamps, shape_blur_factor); 
new_shape = newshapes(1).new_shape;
smooth_new_shape = newshapes(1).smooth_new_shape;
new_skeleton = newshapes(1).new_skeleton;
new_coribs = newshapes(1).new_coribs;

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
patch(smooth_new_shape(:,1),smooth_new_shape(:,2),[0 1 0]);
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
