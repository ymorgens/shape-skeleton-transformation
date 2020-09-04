function articulate_shape(shape,sign,looping,filming)
% articulate_shape(shape,sign,looping,filming)
%
% "Wiggles" the shape, either at concavities (sign = -1) or 
% convexities (+1). This induces a tendency to invert figure/ground
% (border ownership), because articulation at negative curvature minima
% (concavities) moves the positive "parts", consistent with perceiving
% interior as figure (based on Hoffman & Richards' minima rule); while 
% articulation at positive curvature extrema (convexities) moves the 
% negative parts, which then tend to become quasi-figural. 
% See Barenholtz & Feldman, Cognition 2006. 
%
% If looping is 1, goes back and forth; otherwise one direction only. You usually want this on (looping = 1);

%............................................
if nargin == 2
    looping = 1;
end;
%............................................
n_frames = 7; % number of DISTINCT frames; should be an odd number
n_cycles = 3;
duration = n_cycles * n_frames;
angle_increment = pi/48;  % incremental angle the articulartion changes in each frame;

frames = -(n_frames-1)/2:(n_frames-1)/2; % sequence of factors

counter = 1;
extrema_indices = find_extrema(shape,sign);
if length(extrema_indices) > 1

if looping
    reversed_frames = frames(end-1:-1:2);
    frames_to_use = [frames reversed_frames];
else
    frames_to_use = frames; % original
end;
% Execute loop
while counter < duration;
   draw_one_frame(counter,frames_to_use,shape,extrema_indices,angle_increment,filming)
   counter = counter + 1;
end;

else
    fprintf('Only one extremum of desired sign; articulation undefined.\n');
end;

%..................
function draw_one_frame(counter,frames,shape,extrema_indices,angle_increment,filming)
% Draw one frame from "frames", subject to the modulus of "counter"
% frames is the set of multiplicands to use in the succession of frames
global current_movie;
n_frames = length(frames);
this_frame_index = mod(counter,n_frames) + 1;
step_angle = angle_increment * frames(this_frame_index);
shape_articulated = shape; % initialize;
wiggle_direction_counter = 1;
n_extrema = length(extrema_indices);
for i=1:n_extrema;
    this_extremum = extrema_indices(i);
    if i < n_extrema
        next_extremum = extrema_indices(i+1);
    else
        next_extremum = extrema_indices(1);
    end;
    wiggle_direction_counter = wiggle_direction_counter + 1;
    wiggle_direction = mod(wiggle_direction_counter,2);
    if wiggle_direction ==0   % even extremum
        this_angle = step_angle;
    else                    %odd extremum
        this_angle = -1 * step_angle;
    end;
    v = this_extremum;  % starting point
    reached_next_extremum = 0;
    while ~reached_next_extremum
        if v < size(shape,1)  
            v = v+1; % step counter
        else
            v = 1;  % wrap around
        end;
        fulcrum1 = shape(this_extremum,:);
        fulcrum2 = shape(next_extremum,:);
        p1 = rotate_point_about_fulcrum(shape(v,:),fulcrum1,this_angle);
        p2 = rotate_point_about_fulcrum(shape(v,:),fulcrum2,this_angle);
        d1 = norm(fulcrum1 - shape(v,:));
        d2 = norm(fulcrum2 - shape(v,:));
        w1 = d2/(d1+d2);  % mixing proportion
        w2 = 1- w1;
        new_point = w1 * p1 + w2 * p2;
        shape_articulated(v,:) = new_point;
        reached_next_extremum = v == next_extremum;
    end;
end;
cla; 
draw_shape(shape_articulated); %,color);
soa = .03;  %  Pause between frames. Very system dependent. 
pause(soa);
drawnow;
if filming
	new_frame = getframe;  % This line is very slow, so when filming is ON the articulation won't look right (but should store right).
	current_movie = [current_movie new_frame];
end;
    
%.........
function extrema = find_extrema(shape,desired_sign);
normals = jfnormals(shape(:,1),shape(:,2));
% First, make magnitude of normal change
n = size(shape,1);
curvature = zeros(1,n-1);
sign_of_curvature = zeros(1,n-1);
for i=2:n-2
    curvature(i) = ...
        0.5 + norm(normals(i,:) - normals(i+1,:)) +  ...   % difference between successive normals
        0.5 + norm(normals(i-1,:) - normals(i,:));
    cross_product = cross([normals(i,:),0],[normals(i+1,:),0]);
    sign_of_curvature(i) = sign(cross_product(3));;

end;
curvature(n-1) = norm(normals(n-1,:) - normals(1,:)); % wrap around
cross_product_last = cross([normals(n-1,:),0],[normals(1,:),0]);
sign_of_curvature(n-1) = cross_product_last(3)/abs(cross_product_last(3));
extrema = [];
for i=2:n-2
    switch desired_sign
        case 1
            if curvature(i) > curvature(i-1) & curvature(i) > curvature(i+1) & sign_of_curvature(i) == -1;% 
                extrema = [extrema i];
            end;
        case -1
            if curvature(i) > curvature(i-1) & curvature(i) > curvature(i+1) & sign_of_curvature(i) == 1;% 
                   extrema = [extrema i];
            end;
    end;
end;

%.........
function p_rotated = rotate_point_about_fulcrum(p,fulcrum,angle);
% rotates point 
vector = p - fulcrum;
rotated_vector = vector * [cos(angle) sin(angle); -sin(angle) cos(angle)];
p_rotated = fulcrum + rotated_vector;