function color = pick_color(index_raw)
% color = pick_color(index_raw);
%
% Picks a color based on an integer; used in color-coding shape regions
index = mod(index_raw,7); % there are 7 colors
switch index
    case 0
        color = 'b';
    case 1
        color = 'g';   
    case 2
        color = 'r';
    case 3
        color = 'c';  
    case 4
        color = 'm';
    case 5
        color = 'y';  
    case 6
        color = 'k';
end;