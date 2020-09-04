function shape = random_attneave_shape(n_spokes,varargin)
% shape = random_attneave_shape(n_spokes,varargin)
%
% Makes a random blob shape using a procedure that Attneave did not actually
% suggest. 
% 
% Options are "spokestyle" = {fixed, random}, 'smoothing' = integer
default_smoothing = 1;
default_spokestyle = 'fixed';
max_radius = 75;
center = [50 50];
default_position = center;
n_interpolants = 10; % number of points interpolated between each spoke

%..........................
% Get values for options
spokestyle = option_value(varargin,'spokestyle');
if isequal(spokestyle,'null') spokestyle = default_spokestyle; end;
%..
smoothing = option_value(varargin,'smoothing');
if isequal(smoothing,'null') smoothing = default_smoothing; end;
%..
position = option_value(varargin,'position');
if isequal(position,'null') position = default_position; end;
%---------------------------------------------------------------------
shape_notsmooth = [];
protoshape = [];
switch spokestyle
    case 'fixed'  % spokes at even increments around center
        directions = [1:n_spokes]*2*pi/n_spokes;
    case 'random'
        directions = sort(rand(1,n_spokes)) * 2*pi;
    end;
    for i=n_spokes:-1:1 % goes down to make shape come out with correct figural parity given conventions in skeleton routines
        direction = directions(i);
        radius = rand * max_radius;
        protoshape = [protoshape; radius * cos(direction) radius * sin(direction)];                
    end;
    % next, interpolate some points in between the spokes
    for i=1:n_spokes
        if i==n_spokes 
                next_index = 1;
        else
                next_index = i+1;
        end;
        span = protoshape(next_index,:) - protoshape(i,:);
        unit_vector = span/norm(span);
        for v = 0:n_interpolants-1
            new_point = protoshape(i,:) + v * (norm(span)/n_interpolants) * unit_vector;
            shape_notsmooth = [shape_notsmooth; new_point];
        end;
    end;
max_height = 80;
max_width = 80;
shape_notsmooth =  [shape_notsmooth; shape_notsmooth(1,:)]; % close off 
shape_notsmooth = normalize_shape(shape_notsmooth,center,[max_height,max_width]);
shape = smooth_shape(shape_notsmooth,smoothing);
shape =  [shape; shape(1,:)]; % close off again

%-----------------------------------------------------
function value = option_value(options,name)
% Given an arg cell array options, returns the value corresponding to option called 'name'
value = 'null'; % returns if 'name' is not set in options
for i=1:length(options)
    if isequal(options{i},name)
        value = options{i+1};
    end;
end;