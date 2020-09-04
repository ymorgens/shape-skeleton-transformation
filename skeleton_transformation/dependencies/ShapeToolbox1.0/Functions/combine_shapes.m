function shape = combine_shapes(shape1,shape2,combination_type)
% shape = combine_shapes(shape1,shape2,combination_type)
% [combination_type = 'union' or 'intersection']
% Creates union or intersection as requested of shape1 and shape2.
% Input and output shapes are lists of points (n x 2 matrices) in drawable order.

% Basically a wrapper for bwboundaries and inpolygon, which do the hard work.

%...................................................
% First, set size of image-- large enough to encompass both shape1 and shape2
left = min([min(shape1(:,1)) min(shape2(:,1))]);
right = max([max(shape1(:,1)) max(shape2(:,1))]);
top = min([min(shape1(:,2)) min(shape2(:,2))]);
bottom = max([max(shape1(:,2)) max(shape2(:,2))]);
%...................................................
% Next, make the two images.
[x,y] = meshgrid(left:right,top:bottom);
image1 = inpolygon(x,y,shape1(:,1),shape1(:,2));  
image2 = inpolygon(x,y,shape2(:,1),shape2(:,2));
switch combination_type
	case 'intersection'
		combined_image = 255 * (image1 & image2);    % 255 is the value of white in a bw image
	case 'union'
		combined_image = 255 * (image1 | image2);
	otherwise 
		error('Unknown combination type [combine_shapes]');
end;
% Now, combined_image has 255s in the desired region, 0s elsewhere.
%...................................................
% Turn combined image back into a shape
b = bwboundaries(combined_image);
s1 = b{1};
shape = [left - 1 + s1(:,2) top - 1 + s1(:,1)];  
