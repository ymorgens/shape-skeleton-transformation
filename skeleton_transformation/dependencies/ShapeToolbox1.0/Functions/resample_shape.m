function shape_out = resample_shape(shape,n_points)
% shape = resample_shape(shape,n_points)
%
% Resample the shape with the requested number of points. 
% If n_points is larger than the current number of points, the resulting
% shape should have the same trace. If it much smaller, the returned shape 
% may look different. 
arclength = 0;
for i=2:size(shape,1)
    arclength = arclength + norm(shape(i,:)-shape(i-1,:));
end;
spacing = arclength/(n_points-1);
[xs,ys] = finespacedcontour(shape(:,1),shape(:,2),spacing,0);
shape_out = [xs ys];