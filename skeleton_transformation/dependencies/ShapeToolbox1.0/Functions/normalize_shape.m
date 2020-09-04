function shape = normalize_shape(shape_in,center,box_size)
% shape = normalize_shape(shape_in,center,box_size)
%
% Normalize location and size of shape, placing it at center fitting in box
% of size box_size.

max_height = box_size(1);
max_width = box_size(2);
scale_factor_x = max_height/(max(shape_in(:,1)) - min(shape_in(:,1)));
scale_factor_y = max_width/(max(shape_in(:,2)) - min(shape_in(:,2)));
scale_factor = min([scale_factor_x scale_factor_y]);
shape = [scale_factor * shape_in(:,1) scale_factor * shape_in(:,2)];  % scale so it fits in desired rectangle
shape = [center(1) + shape(:,1) - mean(shape(:,1)) center(2) + shape(:,2) - mean(shape(:,2))]; % move to desired location
shape = insideright_shape(shape);