function finize_current_contour(spacing);
% Resamples the current shape at "spacing" (and recomputes normals, etc).

global current_shape;
global current_shape_normals;
global shape_spacing;
% remove duplicates:
new_shape(1,:) = current_shape(1,:);
for i=2:size(current_shape,1)
	this_point = current_shape(i,:);
	previous_point = current_shape(i-1,:);
	if ~isequal(this_point,previous_point)
		new_shape = [new_shape; this_point];
	end;
end;
current_shape = new_shape;
[xs,ys] = finespacedcontour(current_shape(:,1),current_shape(:,2),spacing,0);
current_shape = [xs ys];
current_shape_normals = jfnormals(xs,ys);
shape_spacing = spacing; % for use by other functions later