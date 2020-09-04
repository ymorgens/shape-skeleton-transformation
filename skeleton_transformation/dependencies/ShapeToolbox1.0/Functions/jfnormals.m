function normals = jfnormals(x,y)
% normals = jfnormals(x,y)
%
% calculate some normal vector for each point by its adjacent points'
% coordinates. 
%
% A variant of MS's normals(). 

n = length(x);
switch n
    case 1 
        normals = [1 0]; % arbitrary
    case 2  % for curve with two points, simply give the same normal twice
        p = [x(1) y(1)];
        q = [x(2) y(2)];
        tan = (q-p)/norm(q-p); % [x-normal y-normal]
        normal = [0  -1 ; 1  0] * [tan(1); tan(2)]; %[y-normal;-(x-normal)]
        normals = [normal' ; normal']; %
    otherwise
        % Typical cases follow, now assuming n>2
        if x(end) == x(1) && y(end) == y(1)  % closed
            for i=1:n-1
                if i==1
                    p = [x(1)   y(1)];   % special case for first point
                else
                    p = [x(i-1) y(i-1)];
                end;

                if i==n-1
                    q = [x(n-1) y(n-1)];  % special case for last point
                else
                    q = [x(i+1) y(i+1)];
                end;

                tan = (q-p)/norm(q-p);
                thisnormal= [0  -1 ; 1  0] * [tan(1); tan(2)];
                normals(i,:) = [thisnormal(1) thisnormal(2)];
            end;

        else                                % open
            for i=1:n-1
                if i==1
                    p = [x(1)   y(1)];   % special case for first point
                else
                    p = [x(i-1) y(i-1)];
                end;
                if i==n-1
                    q = [x(n-1) y(n-1)];  % special case for last point
                else
                    q = [x(i+1) y(i+1)];
                end;
                if norm(q-p) == 0
                    thisnormal = [0 1];  % dummy
                else
                    tan = (q-p)/norm(q-p);
                    thisnormal= [0  -1 ; 1  0] * [tan(1); tan(2)];
                end;

                normals(i,:) = [thisnormal(1) thisnormal(2)];

            end;
            normals(n,:) = normals(n-1,:);  % give last value same as second to last
        end;
    
end;