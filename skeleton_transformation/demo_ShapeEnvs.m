%% Skeleton transformation demo (demo_ShapeEnvs.m)

% this demo tranforms the subparts of a shape's skeleton representation 
% to produces classes of objects with different statistics

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
nsamps = 9; % number of samples to produce

%% initial transformation parameters

% choose kind of transformation for demo, in the following...
% trans = 1; % length transfromation
% trans = 2; % spatial position transform
% trans = 3; % angle transfrom
% trans = 4; % width transform

% for example,
%trans = [1 2 3 4 ]; % which transformations to include. use this to include all
%trans = [ 2 3]; % use this to inlude spatial position (2) and angle (3) only

%% standard deviations for environment 1
% select standard deviation of values on changing the subparts, larger values = more change.
length_sig = 1.2; % length , some good starting values: 0.25:0.5:2
pos_sig = 0;%0.05; % spatial position, some good starting values: 0:0.01:0.05
theta_sig = 0;%10; % angle, some good starting values: 0:5:25;
width_sig = 0;%0.2; % width, some good starting values: 0.25:0.5:2
trans = 1; % only length is transformed; ignores other parameters
sdparms = [length_sig, pos_sig, theta_sig, width_sig];

env1shapes = getTransformedShapes(shape, trans, sdparms, nsamps, shape_blur_factor); 

%% environment 2
% vary the parameters above to create classes of shapes with varying
% statistics (e.g.,more variance between parts vs. less)
length_sig = 0.2; % length
pos_sig = 0.05; % spatial position ;
theta_sig = 10; % angle
width_sig = .4; % width
sdparms = [length_sig, pos_sig, theta_sig, width_sig];
trans = [3 4];  % use this to inlude spatial position (2) and angle (3) only
env2shapes = getTransformedShapes(shape, trans, sdparms, nsamps, shape_blur_factor); 

%% show 9 samples from env1 and env2
h = figure('units','normalized','outerposition',[0 0 1 1]); set(h,'Color','w');

% subplot index
spindex = reshape(1:18, 6,3).';

% plot environment 1
spind1 = spindex(:,1:3);
for n = 1:9
    subplot(3,6,spind1(n))
    shape = env1shapes(n).smooth_new_shape;
    
    h = patch(shape(:,1), shape(:,2),'k'); h.FaceColor = [127 201 127]/255; 
    axis image; axis off;
    
end

% plot environment 2
spind2= spindex(:,4:6);
hold on;
for n = 1:9
    subplot(3,6,spind2(n))
    shape = env2shapes(n).smooth_new_shape;
    
    h = patch(shape(:,1), shape(:,2),'k'); h.FaceColor = [253 192 134]/255; 
    axis image; axis off;
    
end

