function shape = image2shape(image,desired_box_size)
% shape = image2shape(image,desired_box_size)
% 
% Given image, convert it to an outline, normalized in various ways.
% Scale so fits into desired_box_size, if given; otherwise no scaling.
%
% image is in matlab format, eg the result of image = imread(

% if not given
% Basic work done by built-in bwboundaries

if nargin == 1
	desired_box_size = [300 120];
end;
center = [50 50];
b = bwboundaries(imread(image));
s1 = b{2};  % b{1} is the frame itself usually, so b{2} is the actual shape
s2 = [s1(:,2) -s1(:,1)];
shape_dimension = [abs(max(s2(:,1)) - min(s2(:,1))) abs(max(s2(:,2)) - min(s2(:,2)))]; % maxX maxY
if shape_dimension(1) > shape_dimension(2) % wider than tall
	scale_factor = desired_box_size(1)/shape_dimension(1);
else
	scale_factor = desired_box_size(2)/shape_dimension(2);
end;
%
shape = [scale_factor * s2(:,1) scale_factor * s2(:,2)];  % scale so it fits in desired rectangle
shape = [center(1) + shape(:,1) - mean(shape(:,1)) center(2) + shape(:,2) - mean(shape(:,2))]; % move to desired location
%shape = [center(1) + shape(:,1) center(2) + shape(:,2) - mean(shape(:,2))]; % this version doesn't normalize x, so skel can move
shape(size(shape,1)+1,:) = shape(1,:);  %close shape
% remove duplicate points ..............................
shape_without_duplicates = [];
last_point = [];
for i=1:size(shape,1)
    next_point = shape(i,:);
    if ~isequal(last_point,next_point)
        shape_without_duplicates = [shape_without_duplicates; next_point];
    end;
    last_point = next_point;
end;
%shape = shape_without_duplicates;
shape = insideright_shape(shape_without_duplicates);  % checks that is has right figural polarity