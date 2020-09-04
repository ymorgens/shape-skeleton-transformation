function shape_out = fractalize_shape(shape,noise_sd)
% shape_out = fractalize_shape(shape,noise_sd)
%
% Add noise to shape contour, by adding a little Gaussian deviate to each
% contour point. 
if nargin==1
    noise_sd = 0.7;  % noise magnitude parameter
end;
for i=1:size(shape,1) 
    shape_out(i,:) = shape(i,:) + normrnd(0,noise_sd,1,2); 
end;
