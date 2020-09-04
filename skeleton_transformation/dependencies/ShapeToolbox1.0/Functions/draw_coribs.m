function draw_coribs(coribs,style,skeleton)
% draw_coribs(coribs,style,skeleton)
%
% Draw the ribs (coribs), e.g. as produced by compute_coribs.
%
% coribs is an n-by-4 matrix, with one row for each correspondence. 
% Columns are 
%
%   shape_point axis_number vertex_number sign 
%
% Shape_point indexes current_shape;
% axis_number and vertex_number index current_skeleton.
% Sign is not used when drawing. 
% skeleton used is usually current_skeleton, but uses given one if given one.
global current_shape;
global current_skeleton;
%-----------------------------------------------------------
if nargin == 1
    style = 'rainbow';
    skeleton = current_skeleton;
elseif nargin == 2
    skeleton = current_skeleton;
end;
%-----------------------------------------------------------
switch style
    case 'none' 
        % do nothing -- no ribs (presumably just skeleton)
    case 'plain'
        for i=1:size(coribs,1)
            c = coribs(i,1);
            a = coribs(i,2);
            v = coribs(i,3);
            skeleton_point = skeleton(a).contour(v,:);
            shape_point = current_shape(c,:);
            private_drawline(skeleton_point,shape_point,'k',1);
        end;
    case 'rainbow' % color coded by axis
        for i=1:size(coribs,1)
            c = coribs(i,1);
            a = coribs(i,2);
            v = coribs(i,3);
            %sign = coribs(i,4);
            color = pick_color(a);  % + sign to make the left and right different colors
            skeleton_point = skeleton(a).contour(v,:);
            shape_point = current_shape(c,:);
            private_drawline(skeleton_point,shape_point,color,1);
        end;
     case 'parts' % entire 'part' color coded by axis  -- a bit klugy. Leave blank regions that are bordered by different-colored ribs ("ambiguous").
        for i=1:size(coribs,1)
			c_this = coribs(i,1);
            a_this = coribs(i,2);
            v_this = coribs(i,3);
			%s_this = coribs(i,4);
			if c_this < size(current_shape,1)-1  % ignoring very last point because it duplicates first
					c_next = c_this + 1;
			else
					c_next = 1;
			end;
			adjacent_ribs = coribs(find(coribs(:,1) == c_next),:);  % ribs that connect to the next shape point
			for i=1:size(adjacent_ribs,1)
				this_adjacent_corib = adjacent_ribs(i,:);
				a_adjacent = this_adjacent_corib(2);
				v_adjacent = this_adjacent_corib(3);
				fill_this_region = true;
				if a_adjacent == a_this  % if adjacent shape points connect to SAME axis
					part_color = pick_color(a_this);  % pick the (common) color of the axes
				else
					if skeleton(a_adjacent).parent == a_this % if a_this "dominates"
						part_color = pick_color(a_this);
					elseif skeleton(a_this).parent == a_adjacent % if a_adjacent dominates
						part_color = pick_color(a_adjacent);
					else
						fill_this_region = false; % in this case, turn filling "off"
					end;
				end;
				%if (a_adjacent == a_this) & (v_adjacent == v_this + 1) % adjacent shape point connects to THIS axis, NEXT point
				if fill_this_region
					A = current_shape(c_this,:);
					B = current_shape(c_next,:);
					C = current_skeleton(a_adjacent).contour(v_adjacent,:);
					D = current_skeleton(a_this).contour(v_this,:);
					color = pick_color(a_this);
					fill([A(1),B(1),C(1),D(1)],[A(2),B(2),C(2),D(2)],part_color,'edgecolor','none');
				end;
			end;
		end;
end;

