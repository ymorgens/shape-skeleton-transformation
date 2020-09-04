function draw_skeleton(skeleton,style)
% draw_skeleton(skeleton,style)
%
% Draws the skeleton (only); options are 'rainbow' and 'plain'.
% Use draw_ribs to draw ribs as well, and draw_shape to draw shape.

global current_augmented_skeleton; % Only needs to use this in 'parts' style
global current_shape; % Only needs to use this in 'parts' style
if nargin == 1
	style = 'rainbow';
end;
if ~exist('line_width')
    line_width = 1;
else
    global line_width;
end;
hold on;
switch style
    case 'rainbow' 
    for index = 1:length(skeleton)
            color_this = pick_color(index);
            axis = skeleton(index);
            plot(axis.contour(:,1),axis.contour(:,2),color_this,'LineWidth',line_width);  % plots axis
    end;
    case 'plain'
        for index = 1:length(skeleton)
            color_this = pick_color(0); % blue
            axis = skeleton(index);
            plot(axis.contour(:,1),axis.contour(:,2),color_this,'LineWidth',line_width);  % plots axis
        end;
    case 'parts'
        if ~isempty(current_augmented_skeleton)
            decompose_shape(current_shape,current_augmented_skeleton,1);
        else
            [dl_total,current_augmented_skeleton] = description_length(skeleton); % Create augmented_skeleton -- XX really should check if it exists first
            decompose_shape(current_shape,current_augmented_skeleton,1);
        end;
end;            