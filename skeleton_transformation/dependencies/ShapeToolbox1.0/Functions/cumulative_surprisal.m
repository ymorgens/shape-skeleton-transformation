function info = cumulative_surprisal(x, y, signed, resolution)
  % Integrates the surprisal over the length of the curve x,y
  %-----------------
 x = x(:); y = y(:);

 if length(x) ~= length(y)
 error('the x and y vectors defining the contour must have equal lengths!');
 end;

 N = length(x);
 b = 1; % spread term for the Von Mises distribution

 if (x(1) == x(N) & y(1) == y(N))
 N = N-1;
 open = 0;
 else
 open = 1;
 end

 %normalmap = [];
 surprisalmap = [];

 for j = 1:N

 % define the previous and next points relative to the current one
 xprev = 0;
 yprev = 0;
 xnext = 0;
 ynext = 0;

 for k = 1:resolution
 prev(k) = j-k;
 next(k) = j+k;
 if(prev(k) < 1) prev(k) = prev(k) + N; end
 if(next(k) > N) next(k) = next(k) - N; end

 xprev = xprev + x(prev(k));
 yprev = yprev + y(prev(k));
 xnext = xnext + x(next(k));
 ynext = ynext + y(next(k));
 end

 xprev = xprev/k;
 yprev = yprev/k;
 xnext = xnext/k;
 ynext = ynext/k;

 % vector _from_ the previous point; vector _to_ the next point
 vecprev = [x(j), y(j)] - [xprev, yprev];
 vecnext = [xnext, ynext] - [x(j), y(j)];

 % compute the magnitude of the turning angle using dot product
  
 ratio = dot(vecprev, vecnext) / (norm(vecprev)*norm(vecnext));
 if ratio > 1
     ratio_ = 1;
 elseif ratio < -1
     ratio_ = -1;
 else
     ratio_ = ratio;
 end;
 alpha = acos(ratio_); 

    
% compute the sign of turning using cross product
% (assumes a counterclockwise sampling of the contour)
 cp = cross( [vecprev, 0], [vecnext,0] );
 alpha = sign(cp(3))*alpha;

 % compute the surprisal using -log von mises
 if(signed)
 surprisal = -log( exp(b*cos(alpha - (2*pi/(N/resolution))) )/( 2*pi*besseli(0,b) ) );
 else
 surprisal = -log( exp(b*cos(alpha) )/( 2*pi*besseli(0,b) ) );
 end;

 % turning angles not defined near the end points of an open curve
 if(open & (j - resolution < 1 | j + resolution > N)) surprisal = 0; end
 
 surprisalmap = [surprisalmap; surprisal];  % add this one to cum list
 
end; % j loop
  
info = sum(surprisalmap);  % integrate over the curve
 