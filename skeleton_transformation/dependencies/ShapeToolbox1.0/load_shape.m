function load_shape(argument)
% function load_shape(argument)
% 
% Loads a shape into the shape_tool GUI.
%
% If argument is a filename, loads the shape from the file. The file should
% be a text file containing and nx2 array of numbers. 
%
% Alternatively, argument can be the shape itself, i.e. an nx2 array of numbers. 
%
% Either way, the shape will be drawn into the GUI, and necessary
% additional data structures created so that shape can be analyzed. 

global current_shape;
global shape_spacing;
global current_shape_normals;
global shape_axis_handle;

hold on; % Holds the properties of the current shape axis
if ischar(argument)    % If it's a filename
    current_shape = load(argument);
else  % presmed to be a shape
    current_shape = normalize_shape(argument,[80 80],[120 120]);
end;
finize_current_shape(shape_spacing);
current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));
axes(shape_axis_handle);
cla;
draw_shape(current_shape);