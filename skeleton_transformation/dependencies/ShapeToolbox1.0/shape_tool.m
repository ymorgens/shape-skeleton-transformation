function shape_tool(arg)
% shape_tool(arg)
%
% Shape toolbox GUI. Inside the GUI you can:
% - draw a shape and modify it in various ways; 
% - randomly generate a shape in various ways;
% - estimate the MAP skeleton;
% - draw its medial axis transform (MAT; Blum, 1973)
% - draw a skeleton manusally and analyze it in various ways
% - Make the shape "wiggle" in order to illustrate the principle described
%    in Barenholtz & Feldman, Cognition 2006
%
% Note that each button has a "tooltip" (mouseover to see it) which
% explains a bit more of what it does. 
%
% Click 'reset consts' to reload the constants set in Functions/set_skeleton_constants.
%
% shape_tool uses many global variables, which can be accessed from the
% command line by declaring them as global. For example, if you create a
% shape and then estimate its MAPskeleton in the GUI, and then want to
% access the skeleton that was estimated, type
% 
% >> global current_skeleton;
% 
% ..and then current_skeleton will contain the skeleton that was just
% estimated. Similarly for current_shape which always contains the shape 
% currently displayed in the GUI. 

%---------------------------------------------------------------------------
% NOTE ON DATASTRUCTURES
% A skeleton is a vector of axes, where an axis
% is a struct with fields
%   index (an integer)
%   parent (index of parent axis)
%   contour (vector of points)
%---------------------------------------------------------------------------
% To save whatever is in fig window (shape, skeleton, ribs, whatever):
% global shape_axis_handle;
% f = getframe(shape_axis_handle); imwrite(f.cdata, 'fig1.png');
%---------------------------------------------------------------------------
% Declarations:
global shape_axis_handle;
global animals;
global window_width;
global current_shape;
global default_skeleton_drawing_style; 
global part_model;
global current_skeleton;
global currently_collecting; % 1 or 0; determines whether mousemotions have any effect
global which_collecting;  % 'shape','skeleton'
                %if we are collecting, determines whether we're adding to
global last_point; % last point collected on whatever we're collecting
global min_distance_new_point;  % how far mouse must travel before new point is added
global current_skeleton_index;  % which axis are we currently adding
global current_axis_contour;  % buffer for the contour of the current axis we're working on
global current_axis;            % current axis nwe're working on
global current_parentindex;          % parent of axis contoru we're working on
global current_parentvertex;    %   ... and vertx on it.
global axis_handle_dls;         % handle for plot of dl's. 
global label_column;
global current_preskeleton_axis;  % does this need to be global?
%global n_initial_parameters;    % num of parameters preceding the first per-axis parameter set (now 1)
%global n_parameters; % number of parameters. Currently 4 but can be changed in future
global n_knotpoints_per_axis;       % n KNOT POINTS per axis -- default 3current_skeleton_parameters
global current_skeleton_parameters;  % contains parameters of current skeleton, n_parameters per axis
global n_points_per_axis;  % number of points (vertices, NOT knotpoints) on the splines that are fit to axis segments.
global current_axis_index;
global line_width;
global current_axis_hierarchy; % vector of indices; the i-th one is the parent of the i-th axis
global shape_spacing;
global distance_decay;  % ... These two constants control both "responsibility" computation and likelihood computation. 
global vonMises_b; 
global branching_penalty;  % cost of one new branch in skeleton
global show_candidates;
global rerib_frequency;  % how often (num opt steps) we recompute the ribs
global optimization_counter; % how many steps we've taken
global comet_trail;  % when optimizing, saves intermediate candidate skeleton parameters
global h_comet_trail;
global button_left;
global button_width;
global button_height;
global button_bottom;
global current_coribs;
global filming;
global n_attneave_axes;  % for random attneave shapes
global attneave_smoothing;
global current_movie;
global current_shape_normals;  % store these globally so compute_coribs can use them
global current_parent;     % used in parametric interface; the index of the axis that contains the
global current_image;
global optimization_tolerance; 
global max_optimization_iterations;
global h_figure;
global light_gray;
global parameter_label_altitude;
%global h_branching_penalty;  
global h_graphical_output;
global show_progressive_computation;
global show_ribs; % Determines whether display of ribs is included when skeletons are drawn.
global h_show_ribs; 
global h_usegradientdescent;
global usegradientdescent;
global optimization_scale_factor; % Controls step size in optimization proess
global current_augmented_skeleton; % Separate current skeleton, so we can keep the regular one "clean"

%---------------------------------------------------------------------------
if nargin==0
    arg='init';
end
switch arg
case 'init' %---------------------------------------------------------------
addpath('Functions');
set_skeleton_constants;  % sets globals 
rand('state',sum(100*clock));   % initialize the random number generator
h_figure = figure('Position',[80 80 1500 600],...
    'Name','Shape tool',...
    'Menubar','figure',...
    'windowbuttonmotionfcn','shape_tool handle_mousemotion',...
    'windowbuttondownfcn','shape_tool handle_mousedown',...
    'windowbuttonupfcn','shape_tool handle_mouseup',...
    'DoubleBuffer','on'); % sets double-buffering on for faster animation?
shape_axis_handle =  axes('Parent',h_figure,...
    'units','pixels',...
    'Position',[label_column + 10*button_left -1*button_height 900 650],... 
    'XTick',[],'Ytick',[],...
    'ButtonDownFcn','shape_tool handle_mousedown');


load Functions/animals.mat animals;
comet_trail = [];
%... Create a shape: ........................................................................
uicontrol('style','text',...
          'position',[label_column-button_width button_bottom + 13 * button_height 2*button_width button_height],...
          'string','Create a shape','BackgroundColor',[0.8 0.8 0.8]); 
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left button_bottom + 13 * button_height 2*button_width 2*button_height],...
    'string','Draw a new shape',...
	'ToolTipString','Draw a new shape by hand',...
    'Callback','shape_tool new_shape');
uicontrol('style','pushbutton',...          % Create a random axial shape (with one axis; press add axis to add more)
    'position',[label_column + button_left + 2*button_width button_bottom + 14 * button_height 1.5*button_width button_height],...
    'string','Random axial shape',...
	'ToolTipString','Make a random shape about a random axial segment',...
    'Callback','shape_tool random_axial_shape');
uicontrol('style','pushbutton',...          % Adds one random axial part to the current shape.
    'position',[label_column + button_left + 3.5*button_width button_bottom + 14 * button_height 0.5*button_width button_height],...
    'string','+axis',...
	'ToolTipString','Add one random axial part to current shape (skeleton generates a new branch)',...
    'Callback','shape_tool add_random_axis');
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left + 2*button_width button_bottom + 13 * button_height 2*button_width button_height],...
    'string','Random blob',...
	'ToolTipString','Make a random blob using a procedure not really due to Attneave',...
    'Callback','shape_tool attneave');
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left + 2*button_width button_bottom + 12 * button_height 2*button_width button_height],...
    'string','Random animal',...
	'ToolTipString','Animalus randomus',...
    'Callback','shape_tool animal');
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left + button_width button_bottom + 12 * button_height button_width button_height],...
    'string','Fractalize',...
	'ToolTipString','Add noise to current shape (adds circular gaussian deviate to each vertex)',...
    'Callback','shape_tool fractalize');
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left button_bottom + 12 * button_height button_width button_height],...
    'string','Smooth',...
	'ToolTipString','Smooth the current shape. Press repeatedly to smooth more.',...
    'Callback','shape_tool smooth');
%... Create a skeleton manually .........................................................................
manual_skeleton_altitude = 4;
uicontrol('style','text',...
          'position',[label_column-button_width-3 button_bottom+ manual_skeleton_altitude * button_height 2*button_width button_height],...
          'string','Create a skeleton','BackgroundColor',[0.8 0.8 0.8]); 
uicontrol('style','pushbutton',...          % Clears current par skeleton, allowing start over with SAME shape    
    'position',[label_column + button_left button_bottom + manual_skeleton_altitude * button_height 2*button_width 2*button_height],...
    'string',' Draw a skeleton manually',...
	'ToolTipString','Draw a new (hierarchical) skeleton by hand. New axes are attached to previous axes at closest point.',...
    'Callback','shape_tool new_skeleton');
uicontrol('style','pushbutton',...          % Attempts to optimize current skeleton
    'position',[label_column + button_left + 2 * button_width button_bottom + manual_skeleton_altitude*button_height 2*button_width button_height],...
    'string','Optimize current skeleton',...
    'Callback','shape_tool optimize',...
    'ToolTipString','Start an optimization process using current skeleton as initialization');

uicontrol('style','pushbutton',...          % Compute complexity and likelihood of current skeleton
    'position',[label_column + button_left + 2 * button_width button_bottom + (manual_skeleton_altitude+1) * button_height 2*button_width button_height],...
    'string','Description length',...
	'ToolTipString','Compute the description length of the current shape + skeleton (also shows ribs)',...
    'Callback','shape_tool analyze'); 
%............................................................................
% MAP etc 
map_altitude = 8;
uicontrol('style','text',...
          'position',[label_column-button_width-8 button_bottom+map_altitude*button_height 2*button_width button_height],...
          'string','Estimate skeleton','BackgroundColor',[0.8 0.8 0.8]); 
uicontrol('style','pushbutton',...          % Estimate skeleton, starting with the (pruned) medial axis
    'position',[label_column + button_left button_bottom+map_altitude*button_height 2 * button_width 2*button_height],...
    'string','MAP Skeleton',...
    'Callback','shape_tool estimate',...
    'ToolTipString','Estimate the Maximum A Posteriori (MAP) skeleton of the current shape.');

uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + button_left+2.5*button_width button_bottom + map_altitude * button_height button_width button_height],...
    'string','Medial axis',...
	'ToolTipString','Show conventional voronoi-based medial axis skeleton',...
    'Callback','shape_tool medial_axis',...
    'ToolTipString','Estimate the Maximum A Posteriori (MAP) skeleton');

%...................................................................................................
uicontrol('style','text',...
          'position',[label_column-button_width+20  button_bottom  button_width button_height],...
          'string','Animate shape','BackgroundColor',[0.8 0.8 0.8]); 
uicontrol('style','pushbutton',...         % Articulate the shape around connvex maxima
    'position',[label_column + button_left button_bottom button_width button_height],...
    'string','Articulate +',...
    'Callback','shape_tool articulate_plus',...
	'ToolTipString','Set the shape in motion, articulating at convexities');
uicontrol('style','pushbutton',...          % Articulate the shape around concave minima
    'position',[label_column + button_left+ button_width button_bottom button_width button_height],...
    'string','Articulate -',...
	'ToolTipString','Set the shape in motion, articulating at concavities',...
    'Callback','shape_tool articulate_minus');

%...................................................................................................
uicontrol('style','pushbutton',...          % Clears the current_shape and current_skeleton, starts over
    'position',[label_column + -2.5*button_left 1 * button_height button_width button_height],...
    'string','Reset consts',...
    'Callback','shape_tool reinit',...
    'ToolTipString','Click this button after changing any of the global parameters in shape_tool.m');

%...................................................................................................
parameter_button_width = 120;
h_usegradientdescent = ...
    uicontrol('style','popupmenu',...
          'position',[label_column + 0*parameter_button_width 1 * button_height parameter_button_width button_height],...
                'string','Include optimization step (slower) | Don"t include optimization step (faster)',...
          'ToolTipString','Turn on or off time-consuming gradient descent procedure.',...
          'Callback','shape_tool change_usegradientdescent');
    set(h_usegradientdescent,'Value',2);  
    uicontrol('style','text',...
          'position',[label_column + 0*parameter_button_width 2*button_height parameter_button_width button_height],...
          'string','Optimization','BackgroundColor',light_gray);
      
h_graphical_output = ...
    uicontrol('style','popupmenu',...
          'position',[label_column + 1*parameter_button_width 1 * button_height parameter_button_width button_height],...
                'string','Do show progressive graphics | Don"t show progressive graphics',...
          'ToolTipString','Turn on or off fancy graphical updates on computation.',...
          'Callback','shape_tool change_graphical_output');
    set(h_graphical_output,'Value',1);  
    uicontrol('style','text',...
          'position',[label_column + 1*parameter_button_width 2*button_height parameter_button_width button_height],...
          'string','Graphics','BackgroundColor',light_gray);
      
h_show_ribs = ...
    uicontrol('style','popupmenu',...
          'position',[label_column + 2*parameter_button_width 1 * button_height parameter_button_width button_height],...
                'string','Show ribs | Don"t show ribs',...
          'ToolTipString','Display ribs with skeletal estimates',...
          'Callback','shape_tool change_show_ribs');
    set(h_show_ribs,'Value',1);  
    uicontrol('style','text',...
          'position',[label_column + 2*parameter_button_width 2*button_height parameter_button_width button_height],...
          'string','Show ribs?','BackgroundColor',light_gray);
 
axes(shape_axis_handle);        % makes shape axis current
cla;                            % clear axes
current_shape = [];             % clear shape
current_skeleton = [];          % clear skeleton
current_skeleton_index = 1; % start collecting first axis again
which_collecting = 'shape';

axis equal;
axis([0 window_width 0 window_width/2]);
axes(shape_axis_handle);
cla;
case 'reinit' %---------------------------------------------------------------
 shape = current_shape;  % Preserve the shape, otherwise erase all constants. 
 set_skeleton_constants;  % Resest constants to those contained in the main constants file. 
 current_shape = shape;%
 
case 'new_shape' %---------------------------------------------------------------
	    axes(shape_axis_handle);        % makes shape axis current
        cla;
		current_shape = [];             % clear shape
        current_skeleton = [];          % clear skeleton
        current_skeleton_index = 1; % start collecting first axis again
        which_collecting = 'shape';     % switch to collecting shape mode
    case 'new_skeleton' %---------------------------------------------------------------
        axes(shape_axis_handle);        % makes shape axis current
        if length(current_skeleton) > 0 % there is a previous skeleton
            current_skeleton = [];      % erase the previous skeleton
            current_skeleton_index = 1; % start collecting first axis again
            cla;                        % clear axes
            current_shape_handle = draw_shape(current_shape);  % redraw just the shape
        end;
        which_collecting = 'skeleton'; % switch to collecting skeleton mode
    case 'new_parametric_skeleton' %---------------------------------------------------------------
        axes(shape_axis_handle);        % makes shape axis current
        current_preskeleton_axis = [];
        current_skeleton_parameters = [];
        current_skeleton = [];      % erase the previous skeleton
        current_skeleton_index = 1; % start collecting first axis again
        cla;                        % clear axes
        draw_shape(current_shape);  % redraw just the shape
        which_collecting = 'preskeleton';
        current_axis_index = 0;
        current_axis_hierarchy = [];
        current_parentvertex = -1;
    case 'new_parametric_axis' %---------------------------------------------------------------
        % Collect one axis parametrically and add it to current par axis
        current_axis_index = current_axis_index + 1;
        [knots_x knots_y] = ginput(n_knotpoints_per_axis); % collect knots of spline to be axis
        % now possibly replace first point with closest on a parent axis, if there is one:
        if current_axis_index == 1
            angle_to_first_point = atan2(knots_y(1),knots_x(1)); % this becomes the first parameter:
            current_skeleton_parameters = [angle_to_first_point]; % here start the param list with "fixed" parameter
            current_axis_hierarchy = [1 -1];  % means 1st axis has parent -1
            %distance<_to_first_point = norm([axis_xs(1) axis_ys(1]);
        else
            [parent_index,vertex_on_parent] = local_find_closest_axis(knots_x(1),knots_y(1));
            current_parentindex = parent_index;
            current_parentvertex  = vertex_on_parent;
            current_axis_hierarchy = [current_axis_hierarchy; current_axis_index current_parentindex]; % add new axis to current_axis hierarchy
 %           current_axis_hierarchy = [current_axis_hierarchy parent_index];
            knots_x(1) = current_skeleton(parent_index).contour(vertex_on_parent,1);
            knots_y(1) = current_skeleton(parent_index).contour(vertex_on_parent,2);  % this replaces first collected point with closest on existing axis
        end;
        parameters = compute_parameters_of_axis([knots_x knots_y]);   % problem here
        current_skeleton_parameters = [current_skeleton_parameters parameters];  % add to existing params
        new_axis = parameters2axis3([knots_x(1) knots_y(1)],parameters);
        for vertex=1:size(new_axis,1)-1
            private_drawline(new_axis(vertex,:),new_axis(vertex+1,:),'b');
        end;
        current_skeleton(current_axis_index).contour = new_axis;
        current_skeleton(current_axis_index).index = current_axis_index;
        if current_axis_index == 1
            current_skeleton(current_axis_index).parent = -1;
        else
            current_skeleton(current_axis_index).parent = current_parentindex;
        end;

    case 'handle_mousedown' %---------------------------------------------------------------
        %  Start accepting points in current_shape or current_skeleton
        %  or, if parametric skeleton, accept 2nd or 3rd point.
        cursor = get(shape_axis_handle,'currentpoint');  % gets current mouse location
        currently_collecting = 1; % turn collecting on
        % this applies to shape and skeleton only, collecting differently in each case
        x = cursor(1);
        y = cursor(3);  % first point
        %axes(shape_axis_handle);
        xl = xlim;
        yl = ylim; 
        inside_axis = x > xl(1) & x < xl(2) & y > yl(1) & y < yl(2);
        if inside_axis
            last_point = [x y]; % for reference later
        end;
        switch which_collecting
            case 'shape'
                current_shape = [x y];  % adds current point to beginning of list
            case 'skeleton'
                if inside_axis
                    if current_skeleton_index==1; % root axis
                        current_parent = -1;
                        %current_parentvertex = -1;
                        current_axis_contour = [x y];  % add first poitn to
                    else  %  non-root axis
                        [parent_index,vertex_on_parent] = local_find_closest_axis(x,y); % index of closest existing axis, and vertex along it
                        current_parent = parent_index;
                        closest_point = current_skeleton(parent_index).contour(vertex_on_parent,:);
                        private_drawline([x y],closest_point,'b');
                        current_axis_contour = [closest_point; x y];  % closest_point (on parent) and current point are first on new contour
                    end;
                end;
            case 'none'; % do nothing
        end; % switch which_collecting
    case 'handle_mousemotion' %---------------------------------------------------------------
        if currently_collecting,
            cursor = get(shape_axis_handle,'currentpoint');  %? gets current mouse location
            x = cursor(1);
            y = cursor(3);  % these are the x and y coords
            %axes(shape_axis_handle);
            xl = xlim;
            yl = ylim;            
            inside_axis = x > xl(1) & x < xl(2) & y > yl(1) & y < yl(2);
            if inside_axis
                candidate = [x y];
                switch which_collecting
                    case 'shape'
                        distance = norm(last_point - candidate);
                        if distance > min_distance_new_point % min distance before accepting a new point
                            private_drawline(last_point,candidate,'k');
                            current_shape = [current_shape; candidate];
                            last_point = candidate;
                        end;
                    case 'skeleton'
                        distance = norm(last_point - candidate);
                        if distance > min_distance_new_point % min distance before accepting a new point
                            private_drawline(last_point,candidate,'b');
                            current_axis_contour = [current_axis_contour; candidate];
                            last_point = candidate;
                        end; % if
                        %case 'preskeleton' % do nothing.
                    case 'none'; % do nothing
                end; % switch;
            end; % if inside axis
        end; % if

    case 'handle_mouseup' %---------------------------------------------------------------
        currently_collecting = 0;  % turn collecting off
        xl = xlim;
        yl = ylim;
        cursor = get(shape_axis_handle,'currentpoint');  %? gets current mouse location

        x = cursor(1);
        y = cursor(3);  % these are the x and y coords
        inside_axis = x > xl(1) & x < xl(2) & y > yl(1) & y < yl(2);

        switch which_collecting
            case 'shape'
                which_collecting = 'skeleton'; % automatically start getting a skeleton
                Xs = [current_shape(:,1); current_shape(1,1)];
                Ys = [current_shape(:,2); current_shape(1,2)];  % close the shape
                [smoothX smoothY] = finespacedcontour(Xs,Ys,shape_spacing,0);
                current_shape = [smoothX smoothY];
                axes(shape_axis_handle);        % makes shape axis current
                cla;                           % clear axes
                current_shape_handle = draw_shape(current_shape);  % redraw now smoothed shaep
                current_shape_normals = jfnormals(smoothX,smoothY);  % install this shape's normals

            case 'skeleton'
                cursor = get(shape_axis_handle,'currentpoint');
                x = cursor(1);
                y = cursor(3);  % current cursor
                inside_axis = x > xl(1) & x < xl(2) & y > yl(1) & y < yl(2);
                if inside_axis
                    current_axis.index = current_skeleton_index;
                    current_axis.parent = current_parent;
                    skeleton_spacing = 3;
                    %current_axis.parentvertex = current_parentvertex;
                    current_axis_contour = [current_axis_contour; x y];
                    [smoothX smoothY] = finespacedcontour(current_axis_contour(:,1),current_axis_contour(:,2),skeleton_spacing,0);
                    current_axis.contour = [smoothX smoothY];
                    current_skeleton = [current_skeleton; current_axis];
                    current_skeleton_index = current_skeleton_index + 1; % advance to be ready for next skelaxis
                    %current_skeleton_parameters = skeleton2parameters(current_skeleton);
                end;
            case 'preskeleton'
                % do nothing
            case 'none'; % do nothing
             
end;
case 'random_axial_shape' %-------------------------------------------------------------------------------------------------------------------
[current_shape,current_skeleton] = random_axial_shape(1);  % generates random shape and resets globals current_shape and current_skeleton
hold on;
cla;
current_shape = normalize_shape(current_shape,[50 50],[100 100]);
draw_shape(current_shape);
current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));
draw_shape(current_shape);


case 'articulate_plus' %-------------------------------------------------------------------------------------------------------------------
current_movie = [];  %
current_skeleton = [];  % The skeleton gets erased, so for consistency we delete it. 
articulate_shape(current_shape,1,1,filming);  
	
case 'articulate_minus' %-------------------------------------------------------------------------------------------------------------------
current_movie = [];
% The skeleton gets erased, so for consistency we delete it. 
articulate_shape(current_shape,-1,1,filming);
    
case 'add_random_axis' %-------------------------------------------------------------------------------------------------------------------
% Add one random axis, but keep the rest of the shape the same
% First, change current_skeleton if it exists
root_width = 10; % ??? same as random axial shape uses as base width
if ~isempty(current_skeleton) % if there's alread a shape and a skeleton defined
    n = length(current_skeleton);
    perm = randperm(n);
    theparent = perm(1); % random number between 1 and n-1 
    current_skeleton(n+1).parent = theparent;
    current_skeleton(n+1).index = n+1;
    parentaxis = current_skeleton(theparent).contour;
    parentnormals = jfnormals(parentaxis(:,1),parentaxis(:,2));
    perm = randperm(n_points_per_axis);
    start_index = perm(1);
    start_location = parentaxis(start_index,:);
    d = rand; sign = (d-0.5)/abs(d-0.5); % set random sign (kluge)
    start_vector = sign* parentnormals(start_index,:);
    %....
    start_direction = atan2(start_vector(2),start_vector(1));
    [Xs Ys] = random_curve(start_direction,start_location,n_points_per_axis);
    current_skeleton(n+1).contour = [Xs Ys];
    % Next, create a new standalone "part" from the new branch
    mean_rib_to_use = compute_mean_rib_length(root_width,n+1,current_skeleton); % based on generation of n+1-th axis in current_skeleton.
    new_part = skeleton2shape(current_skeleton(n+1),'root_width',mean_rib_to_use); % feeds the desired with to skeleton2shape as root width, because it's just making one blob
    % Now, union the new part and the old shape, then draw.
    current_shape = combine_shapes(current_shape,new_part,'union');
	finize_current_shape(shape_spacing);
	current_shape = smooth_shape(current_shape,1,'gaussian');
    current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));
    cla;
    draw_shape(current_shape);
    hold on
  else % if no shape, just make a random axial shape
    [current_shape,current_skeleton] = random_axial_shape(1);  % generates random shape and resets globals current_shape and current_skeleton
    hold on;
	cla;
    current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));
    draw_shape(current_shape); 
  end;

case 'attneave' %-------------------------------------------------------------------------------------------------------------------
cla;
hold on;
load_shape(random_attneave_shape(n_attneave_axes,'smoothing',attneave_smoothing,'spokestyle','random'));  % generates random shape and resets globals current_shape and current_skeleton

case 'animal' %-------------------------------------------------------------------------------------------------------------------
% Assumes animals.mat has been loaded (see init).
perm = randperm(14);
animal_index = perm(1);
cla;
%hold on;
load_shape(animals{animal_index})

case 'ribs' %-------------------------------------------------------------------------------------------------------------------
% compute ribs of current_skeleton wrt current_shape, and put them in current_ribs
% this draws just the shape-interior part of the voronoi points
cla;
draw_shape(current_shape);
current_shape_normals = jfnormals(current_shape(:,1),current_shape(:,2));
draw_skeleton(current_skeleton);
coribs = compute_coribs(current_skeleton);
draw_coribs(coribs);
case 'medial_axis' %----------------------------------------------------------------------------------------------------------------
% compute ribs of current_skeleton wrt current_shape, and put them in current_ribs

medial_skeleton = make_medial_axis(current_shape);
current_skeleton = medial_skeleton;
cla;
draw_shape(current_shape);
draw_skeleton(current_skeleton,'plain'); 

case 'analyze' %-------------------------------------------------------------------------------------------------------------------
% Compute description length for current skeleton/current_shape,
% and draws ribs if draw_ribs is on. 

%current_skeleton = compute_coribs(current_skeleton); % adds ribs to current_skeleton
if isequal(current_shape(1,:),current_shape(end,:))
    current_shape = current_shape(1:end-1,:);
end;
current_coribs = compute_coribs(current_skeleton); % set initial coribs
if show_ribs
    cla;
    draw_shape(current_shape);
    draw_coribs(current_coribs,'rainbow');
    draw_skeleton(current_skeleton,'rainbow');
end;
[dl,current_augmented_skeleton] = description_length(current_skeleton);
fprintf('Description length (-log posterior):   %f5\n',dl);
% Augmented_skeleton is original plus new field smoothedaxisriblengths,
% which is used by decompose_shape


%%%
    %XX debugging -- putting this here to test it, later will move it
    
    decomposition = decompose_shape(current_shape,current_augmented_skeleton,1);

case 'optimize' %---------------------------------------------------------------
% Start a Nelder-Mead process from current skeleton (parameterized).

comet_trail = [];
optimization_counter = 0;
%current_movie = [];  
if ~isempty(current_skeleton)
    optimize_skeleton(current_skeleton,current_shape)
else
    error('No skeleton currently defined.');
end;

case 'estimate' %---------------------------------------------------------------
% Take current_shape, compute medial axis, prune, convert to parameters, and optimize (ie, the whole shabang)
cla;
draw_shape(current_shape);
if usegradientdescent
    winning_skeleton = mapskeleton(current_shape,'gradientdescent','using_gui'); 
    %draw_skeleton(winning_skeleton,'rainbow'); 
else
    if show_progressive_computation
        winning_skeleton = mapskeleton(current_shape,'show_progressive_computation','using_gui');  
    else
        winning_skeleton = mapskeleton(current_shape,'using_gui');  % all the fastest options
    end;
    draw_skeleton(winning_skeleton,'rainbow'); 
end;
if show_ribs
    draw_coribs(compute_coribs(winning_skeleton),'rainbow');
    
end;

case 'fractalize' %---------------------------------------------------------------
% Add noise to current shape and redraw
noise_sd = 0.7;  % noise magnitude parameter
for i=1:size(current_shape,1) current_shape(i,:) = current_shape(i,:) + normrnd(0,noise_sd,1,2); end;
current_shape(end,:) = current_shape(1,:); % re-close shape
axes(shape_axis_handle);        % 
cla;
draw_shape(current_shape);

case 'smooth' %---------------------------------------------------------------
% Smooth the current shape.
% Note repeatedly pressing this button makes it smoother and smoother.
mask = 2; % 
current_shape = smooth_shape(current_shape,mask,'gaussian');
current_shape(end,:) = current_shape(1,:); % re-close shape
axes(shape_axis_handle);        % 
cla;
draw_shape(current_shape);

%------------------------------------------------------------------------------    
case 'change_graphical_output'
    graphical_output_level = get(h_graphical_output,'Value');
    if graphical_output_level == 1
       show_progressive_computation = 1;
    else
       show_progressive_computation = 0;
    end;
   
%------------------------------------------------------------------------------    
case 'change_show_ribs'   
    show_ribs_level = get(h_show_ribs,'Value');
    if show_ribs_level == 1
       show_ribs = 1;
    else
       show_ribs = 0;
    end;
%------------------------------------------------------------------------------    
case 'change_usegradientdescent'
    usegradientdescent_level = get(h_usegradientdescent,'Value');
    if usegradientdescent_level == 1
       usegradientdescent = 1;
    else
       usegradientdescent = 0;
    end;    
%------------------------------------------------------------------------------    
otherwise error('Bad argument to shape tool!');   

        
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [winning_index,winning_vertex] = local_find_closest_axis(x,y)
% Finds index of axis and vertex along axis that is closest to [x,y]
global current_skeleton; % vector of axes
best_index_so_far = 1;
best_vertex_so_far = 1;
best_distance_so_far = norm([x y] - current_skeleton(1).contour(1));
for axis_number=1:length(current_skeleton)
    this_axis_contour = current_skeleton(axis_number).contour;
    for vertex=1:size(this_axis_contour,1)
        distance = norm([x y] - this_axis_contour(vertex,:));
        if distance < best_distance_so_far
            best_distance_so_far = distance;
            best_index_so_far = axis_number;
            best_vertex_so_far = vertex;
        end;
    end;
end;
winning_index = best_index_so_far;
winning_vertex = best_vertex_so_far;


%------------------------------------------------------------------------
function this_mean_rib = compute_mean_rib_length(root_width,axis_index,skeleton)
% determines mean rib based on generation of a within skeleton
% TAKEN FROM skeleton2shape--- but we need it in shape_tool to generate random axial shapes part by part.
reduction_fraction = 0.5; % reduction in width for each generation
generation = determine_generation(axis_index,skeleton);
this_mean_rib  = root_width * reduction_fraction ^ generation;

%.........
function generation = determine_generation(axis_index,skeleton);
if axis_index == 1
    generation = 0;
else
    parents_generation = determine_generation(skeleton(axis_index).parent,skeleton);
    generation = parents_generation + 1;
end;


