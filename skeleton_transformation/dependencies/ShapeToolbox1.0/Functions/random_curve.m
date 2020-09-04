function [Xs Ys] = random_curve(starting_orientation,starting_location,n_points)
% [Xs Ys] = random_curve(starting_orientation,starting_location,n_points)
step = 2; % length of each leg
sd = pi/50;
Xs = zeros(n_points,1);
Ys = zeros(n_points,1);
Xs(1) = starting_location(1);
Ys(1) = starting_location(2);
Xs(2) = Xs(1) + step * cos(starting_orientation);
Ys(2) = Ys(1) + step * sin(starting_orientation);
for i=3:n_points
    direction_change = pi * normrnd(0,sd); % norm with sd 60deg; kluge for von Mises
    previous_direction = atan2(Ys(i-1) - Ys(i-2),Xs(i-1) - Xs(i-2));
    Xs(i) = Xs(i-1) + step * cos(previous_direction + direction_change);
    Ys(i) = Ys(i-1) + step * sin(previous_direction + direction_change);
end;