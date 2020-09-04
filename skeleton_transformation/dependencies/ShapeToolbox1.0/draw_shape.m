function handle = draw_shape(shape,varargin)
% handle = draw_shape(shape,varargin)
%
% Draws a shape. Optional set of arguments, varargin =
% {closed,color,line_width}
%   closed: 1 = closed (default), 0 = open
%   color: color of shape boundary, e.g. 'k','r', etc.
%   line_width: thickness of boundary

% Basically a wrapper for plot.

if isempty(varargin)   % Default is to assume the shape is closed. If closed = 0, draws an open curve. 
    closed = true;
    color = 'k';  % black
    line_width = 1;
else
    closed = varargin{1};
    color = varargin{2};
    line_width = varargin{3};
end;
if closed 
    if ~isequal(shape(1,:),shape(end,:))  % if shape is not already closed
        shape = [shape; shape(1,:)];      % add first point to the end.
    end;
end;
if ~isempty(shape)
    handle = plot(shape(:,1),shape(:,2),color,'LineWidth',line_width);
else
    fprintf('Warning: empty shape.\n');
end;
